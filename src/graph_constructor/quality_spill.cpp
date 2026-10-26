#include "quality_spill.h"

#include <fcntl.h>
#include <unistd.h>

#include <cerrno>
#include <cstring>
#include <stdexcept>

namespace raven {

namespace {

constexpr std::size_t kWriteBuffer = 64U << 20;

std::runtime_error SpillError(const std::string& what, const std::string& path) {
  return std::runtime_error("[raven::QualitySpill] error: " + what + " " +
                            path + ": " + std::strerror(errno));
}

}  // namespace

QualitySpill::QualitySpill(const std::string& dir)
    : path_(dir + "/draven_qualities." + std::to_string(::getpid()) + ".bin"),
      buffer_(kWriteBuffer),
      offsets_(1, 0) {
  file_ = std::fopen(path_.c_str(), "wb");
  if (file_ == nullptr) {
    throw SpillError("cannot create", path_);
  }
  std::setvbuf(file_, buffer_.data(), _IOFBF, buffer_.size());
}

QualitySpill::~QualitySpill() {
  if (file_ != nullptr) {
    std::fclose(file_);
  }
  if (fd_ >= 0) {
    ::close(fd_);
  }
  std::remove(path_.c_str());
}

void QualitySpill::Add(biosoup::NucleicAcid& read) {
  const std::size_t size = read.quality.size();
  if (size > 0 && std::fwrite(read.quality.data(), 1, size, file_) != size) {
    throw SpillError("cannot write", path_);
  }
  offsets_.push_back(offsets_.back() + size);
  std::vector<std::uint8_t>().swap(read.quality);
}

void QualitySpill::Finish() {
  if (std::fclose(file_) != 0) {
    file_ = nullptr;
    throw SpillError("cannot close", path_);
  }
  file_ = nullptr;
  std::vector<char>().swap(buffer_);
  fd_ = ::open(path_.c_str(), O_RDONLY);
  if (fd_ < 0) {
    throw SpillError("cannot open", path_);
  }
}

void QualitySpill::Load(std::uint32_t index, biosoup::NucleicAcid& read) const {
  const std::uint64_t begin = offsets_[index];
  const std::size_t size = offsets_[index + 1] - begin;
  read.quality.resize(size);
  std::size_t done = 0;
  while (done < size) {
    const ssize_t n = ::pread(fd_, read.quality.data() + done, size - done,
                              static_cast<off_t>(begin + done));
    if (n < 0 && errno == EINTR) {
      continue;
    }
    if (n <= 0) {
      throw SpillError(n < 0 ? "cannot read" : "truncated", path_);
    }
    done += n;
  }
}

}  // namespace raven
