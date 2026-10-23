#include "candidate_spill.h"

#include <unistd.h>

#include <algorithm>
#include <cerrno>
#include <cstring>
#include <stdexcept>

namespace raven {

namespace {

// bytes buffered per group before a write
constexpr std::size_t kWriteBuffer = 8U << 20;
// bytes read at a time when loading
constexpr std::size_t kReadBuffer = 64U << 20;

std::runtime_error SpillError(const std::string& what, const std::string& path) {
  return std::runtime_error("[raven::CandidateSpill] error: " + what + " " +
                            path + ": " + std::strerror(errno));
}

}  // namespace

CandidateSpill::CandidateSpill(const std::string& dir, std::uint32_t num_groups,
                               std::uint32_t num_piles)
    : dir_(dir), num_piles_(num_piles), groups_(std::max<std::uint32_t>(1, num_groups)) {
  const std::string prefix =
      dir + "/draven_candidates." + std::to_string(::getpid()) + ".";
  for (std::uint32_t g = 0; g < groups_.size(); ++g) {
    Group& group = groups_[g];
    group.path = prefix + std::to_string(g) + ".bin";
    group.file = std::fopen(group.path.c_str(), "wb");
    if (group.file == nullptr) {
      throw SpillError("cannot create", group.path);
    }
    group.buffer.reserve(kWriteBuffer);
  }
}

CandidateSpill::~CandidateSpill() {
  for (auto& group : groups_) {
    if (group.file != nullptr) {
      std::fclose(group.file);
    }
    if (!group.path.empty()) {
      std::remove(group.path.c_str());
    }
  }
}

// Group g starts at ceil(g * num_piles / num_groups), so pile p is in
// floor(p * num_groups / num_piles).
std::uint32_t CandidateSpill::First(std::uint32_t group) const {
  const std::uint64_t n = groups_.size();
  return static_cast<std::uint32_t>((group * static_cast<std::uint64_t>(num_piles_) + n - 1) / n);
}

std::uint32_t CandidateSpill::GroupOf(std::uint32_t pile) const {
  return static_cast<std::uint32_t>(
      pile * static_cast<std::uint64_t>(groups_.size()) / num_piles_);
}

void CandidateSpill::Add(const correction_candidate& c) {
  Group& group = groups_[GroupOf(c.lhs_id)];
  const correction_candidate_fields& fields = c;
  const std::uint32_t num_anchors = c.anchors.size();
  const std::size_t anchor_bytes = num_anchors * sizeof(c.anchors[0]);
  const std::size_t size = sizeof(fields) + sizeof(num_anchors) + anchor_bytes;
  if (group.buffer.size() + size > kWriteBuffer) {
    Flush(group);
  }
  const char* parts[] = {reinterpret_cast<const char*>(&fields),
                         reinterpret_cast<const char*>(&num_anchors),
                         reinterpret_cast<const char*>(c.anchors.data())};
  const std::size_t sizes[] = {sizeof(fields), sizeof(num_anchors), anchor_bytes};
  for (int i = 0; i < 3; ++i) {
    group.buffer.insert(group.buffer.end(), parts[i], parts[i] + sizes[i]);
  }
}

void CandidateSpill::Flush(Group& group) {
  if (group.buffer.empty()) {
    return;
  }
  if (std::fwrite(group.buffer.data(), 1, group.buffer.size(), group.file) !=
      group.buffer.size()) {
    throw SpillError("cannot write", group.path);
  }
  bytes_written_ += group.buffer.size();
  group.buffer.clear();
}

void CandidateSpill::Finish() {
  for (auto& group : groups_) {
    Flush(group);
    std::vector<char>().swap(group.buffer);
    if (std::fclose(group.file) != 0) {
      group.file = nullptr;
      throw SpillError("cannot close", group.path);
    }
    group.file = nullptr;
  }
}

std::uint64_t CandidateSpill::Load(
    std::uint32_t g,
    std::vector<std::vector<correction_candidate>>& piles,
    const std::function<bool(std::uint32_t)>& keep) {
  Group& group = groups_[g];
  std::FILE* file = std::fopen(group.path.c_str(), "rb");
  if (file == nullptr) {
    throw SpillError("cannot open", group.path);
  }
  std::vector<char> buffer(kReadBuffer);
  std::size_t begin = 0;
  std::size_t end = 0;
  std::uint64_t bytes_read = 0;
  // the next n bytes (to dst, or skipped if dst is null); false at the end
  // of the file before any byte of them
  auto read = [&](void* dst, std::size_t n) -> bool {
    char* out = static_cast<char*>(dst);
    bool any = false;
    while (n > 0) {
      if (begin == end) {
        end = std::fread(buffer.data(), 1, buffer.size(), file);
        begin = 0;
        bytes_read += end;
        if (end == 0) {
          if (std::ferror(file)) {
            throw SpillError("cannot read", group.path);
          }
          if (any) {
            throw std::runtime_error("[raven::CandidateSpill] error: truncated " + group.path);
          }
          return false;
        }
      }
      const std::size_t k = std::min(n, end - begin);
      if (out != nullptr) {
        std::memcpy(out, buffer.data() + begin, k);
        out += k;
      }
      begin += k;
      n -= k;
      any = true;
    }
    return true;
  };

  correction_candidate c;
  correction_candidate_fields& fields = c;
  std::uint32_t num_anchors = 0;
  while (read(&fields, sizeof(fields))) {
    if (!read(&num_anchors, sizeof(num_anchors))) {
      throw std::runtime_error("[raven::CandidateSpill] error: truncated " + group.path);
    }
    const std::size_t anchor_bytes = num_anchors * sizeof(c.anchors[0]);
    if (c.lhs_id >= piles.size() || !keep(c.lhs_id)) {
      if (anchor_bytes > 0 && !read(nullptr, anchor_bytes)) {
        throw std::runtime_error("[raven::CandidateSpill] error: truncated " + group.path);
      }
      continue;
    }
    c.anchors.resize(num_anchors);
    if (anchor_bytes > 0 && !read(c.anchors.data(), anchor_bytes)) {
      throw std::runtime_error("[raven::CandidateSpill] error: truncated " + group.path);
    }
    piles[c.lhs_id].emplace_back(std::move(c));
    c.anchors = {};
  }
  std::fclose(file);
  std::remove(group.path.c_str());
  group.path.clear();
  return bytes_read;
}

}  // namespace raven
