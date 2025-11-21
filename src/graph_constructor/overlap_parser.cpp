
#include <stdexcept>
#include "overlap.h"
#include "overlap_parser.h"

namespace raven {
  void OverlapParser::AddPafOverlapToList(std::vector<std::unique_ptr<Overlap>>& dst, bool shorten_names) {
    const char* q_name = nullptr;
    std::uint32_t q_name_len = 0;
    std::uint32_t q_len = 0;
    std::uint32_t q_begin = 0;
    std::uint32_t q_end = 0;
    const char* t_name = nullptr;
    std::uint32_t t_name_len = 0;
    std::uint32_t t_len = 0;
    std::uint32_t t_begin = 0;
    std::uint32_t t_end = 0;
    std::uint32_t num_matches = 0;
    std::uint32_t overlap_len = 0;
    std::uint32_t quality = 0;
    char orientation = '\0';

    auto storage_ptr = RightStrip(
      line_.data(),
      line_ptr_);
    TerminateLine(storage_ptr);

    std::uint32_t num_values = 0;
    std::uint32_t begin_ptr = 0;
    std::vector<OverlapDescriptor> descriptors;

    while (true) {
      auto end_ptr = begin_ptr;
      while (end_ptr < storage_ptr && line_[end_ptr] != '\t') {
        ++end_ptr;
      }
      TerminateLine(end_ptr);

      switch (num_values) {
        case 0:
          q_name = line_.data() + begin_ptr;
          q_name_len = end_ptr - begin_ptr;
          break;
        case 1: q_len = std::atoi(line_.data() + begin_ptr); break;
        case 2: q_begin = std::atoi(line_.data() + begin_ptr); break;  // NOLINT
        case 3: q_end = std::atoi(line_.data() + begin_ptr); break;
        case 4: orientation = line_[begin_ptr]; break;
        case 5:
          t_name = line_.data() + begin_ptr;
          t_name_len = end_ptr - begin_ptr;
          break;
        case 6: t_len = std::atoi(line_.data() + begin_ptr); break;
        case 7: t_begin = std::atoi(line_.data() + begin_ptr); break;  // NOLINT
        case 8: t_end = std::atoi(line_.data() + begin_ptr); break;
        case 9: num_matches = std::atoi(line_.data() + begin_ptr); break;  // NOLINT
        case 10: overlap_len = std::atoi(line_.data() + begin_ptr); break;  // NOLINT
        case 11: quality = std::atoi(line_.data() + begin_ptr); break;  // NOLINT
      }

      if (num_values == 12) {
        //parse descriptors
        std::string desc_str(line_.data() + begin_ptr);
        int pos1 = desc_str.find(':', 0);
        int pos2 = desc_str.find(':', pos1 + 1);
        OverlapDescriptor desc {
          desc_str.substr(0, pos1),
          desc_str.substr(pos1 + 1, pos2 - pos1 - 1),
          desc_str.substr(pos2 + 1)
        };
        descriptors.emplace_back(desc);
      } else {
        ++num_values;
      }

      if (end_ptr == storage_ptr) {
        break;
      }
      begin_ptr = end_ptr + 1;
    }

    if (num_values != 12) {
      throw std::invalid_argument(
        "[bioparser::PafParser] error: invalid file format");
    }

    q_name_len = shorten_names ?
                 Shorten(q_name, q_name_len) :
                 RightStrip(q_name, q_name_len);

    t_name_len = shorten_names ?
                 Shorten(t_name, t_name_len) :
                 RightStrip(t_name, t_name_len);

    if (q_name_len == 0 || t_name_len == 0) {
      throw std::invalid_argument(
        "[bioparser::PafParser] error: invalid file format");
    }

    dst.emplace_back(std::unique_ptr<Overlap>(new Overlap(
      q_name, q_name_len, q_len, q_begin, q_end,
      orientation,
      t_name, t_name_len, t_len, t_begin, t_end,
      num_matches,
      overlap_len,
      quality,
      descriptors)));
  }

  void OverlapParser::AddSamOverlapToList(std::vector<std::unique_ptr<Overlap>>& dst, bool shorten_names)
  {
    const char* q_name = nullptr;
    std::uint32_t q_name_len = 0;
    std::uint32_t flag = 0;
    const char* t_name = nullptr;
    std::uint32_t t_name_len = 0;
    std::uint32_t t_begin = 0;
    std::uint32_t map_quality = 0;
    const char* cigar = nullptr;
    std::uint32_t cigar_len = 0;
    const char* t_next_name = nullptr;
    std::uint32_t t_next_name_len = 0;
    std::uint32_t t_next_begin = 0;
    std::uint32_t template_len = 0;
    const char* data = nullptr;
    std::uint32_t data_len = 0;
    const char* quality = nullptr;
    std::uint32_t quality_len = 0;

    if (line_[0] == '@') {  // file header
      ClearLine();
      return;
    }

    auto storage_ptr = RightStrip(
      line_.data(),
      line_ptr_);
    TerminateLine(storage_ptr);

    std::uint32_t num_values = 0;
    std::uint32_t begin_ptr = 0;
    std::vector<OverlapDescriptor> descriptors;

    while (true) {
      auto end_ptr = begin_ptr;
      while (end_ptr < storage_ptr && line_[end_ptr] != '\t') {
        ++end_ptr;
      }
      TerminateLine(end_ptr);

      switch (num_values) {
        case 0:
          q_name = line_.data() + begin_ptr;
          q_name_len = end_ptr - begin_ptr;
          break;
        case 1: flag = std::atoi(line_.data() + begin_ptr); break;
        case 2:
          t_name = line_.data() + begin_ptr;
          t_name_len = end_ptr - begin_ptr;
          break;
        case 3: t_begin = std::atoi(line_.data() + begin_ptr); break;  // NOLINT
        case 4: map_quality = std::atoi(line_.data() + begin_ptr); break;  // NOLINT
        case 5:
          cigar = line_.data() + begin_ptr;
          cigar_len = end_ptr - begin_ptr;
          break;
        case 6:
          t_next_name = line_.data() + begin_ptr;
          t_next_name_len = end_ptr - begin_ptr;
          break;
        case 7: t_next_begin = std::atoi(line_.data() + begin_ptr); break;  // NOLINT
        case 8: template_len = std::atoi(line_.data() + begin_ptr); break;  // NOLINT
        case 9:
          data = line_.data() + begin_ptr;
          data_len = end_ptr - begin_ptr;
          break;
        case 10:
          quality = line_.data() + begin_ptr;
          quality_len = end_ptr - begin_ptr;
          break;
      }

      if (num_values == 11) {
        //parse descriptors
        std::string desc_str(line_.data() + begin_ptr);
        int pos1 = desc_str.find(':', 0);
        int pos2 = desc_str.find(':', pos1 + 1);
        OverlapDescriptor desc {
          desc_str.substr(0, pos1),
          desc_str.substr(pos1 + 1, pos2 - pos1 - 1),
          desc_str.substr(pos2 + 1)
        };
        descriptors.emplace_back(desc);
      } else {
        ++num_values;
      }

      if (end_ptr == storage_ptr) {
        break;
      }
      begin_ptr = end_ptr + 1;
    }

    if (num_values != 11) {
      throw std::invalid_argument(
        "[bioparser::SamParser] error: invalid file format");
    }

    q_name_len = shorten_names ?
                 Shorten(q_name, q_name_len) :
                 RightStrip(q_name, q_name_len);

    t_name_len = shorten_names ?
                 Shorten(t_name, t_name_len) :
                 RightStrip(t_name, t_name_len);

    cigar_len = RightStrip(cigar, cigar_len);

    t_next_name_len = shorten_names ?
                      Shorten(t_next_name, t_next_name_len) :
                      RightStrip(t_next_name, t_next_name_len);

    data_len = RightStrip(data, data_len);
    quality_len = RightStrip(quality, quality_len);

    if (q_name_len == 0 || t_name_len == 0 || cigar_len == 0 ||
        t_next_name_len == 0 || data_len == 0 || quality_len == 0 ||
        (data_len > 1 && quality_len > 1 && data_len != quality_len)) {
      throw std::invalid_argument(
        "[bioparser::SamParser] error: invalid file format");
    }

    dst.emplace_back(std::unique_ptr<Overlap>(new Overlap(
      q_name, q_name_len,
      flag,
      t_name, t_name_len, t_begin,
      map_quality,
      cigar, cigar_len,
      t_next_name, t_next_name_len, t_next_begin,
      template_len,
      data, data_len,
      quality, quality_len,
      descriptors)));
  }

  std::vector<std::unique_ptr<Overlap>> OverlapParser::ParsePaf(std::uint64_t bytes, bool shorten_names)
  {
    return Parse(bytes, shorten_names, &OverlapParser::AddPafOverlapToList);
  }

  std::vector<std::unique_ptr<Overlap>> OverlapParser::ParseSam(std::uint64_t bytes, bool shorten_names)
  {
    return Parse(bytes, shorten_names, &OverlapParser::AddSamOverlapToList);
  }

  std::vector<std::unique_ptr<Overlap>> OverlapParser::Parse(std::uint64_t bytes, bool shorten_names,
                                                       void (OverlapParser::*add_overlap_to_list)(std::vector<std::unique_ptr<Overlap>>& dst, bool shorten_names))
  {
    std::vector<std::unique_ptr<Overlap>> dst;
    std::uint64_t parsed_bytes = 0;
    bool is_eof = false;

    while (true) {
      auto buffer_ptr = buffer_ptr_;
      for (; buffer_ptr < buffer_bytes_; ++buffer_ptr) {
        auto c = buffer_[buffer_ptr];
        if (c == '\n') {
          ReadLineFromBuffer(buffer_ptr - buffer_ptr_);
          (this->*add_overlap_to_list)(dst, shorten_names);
          parsed_bytes += line_ptr_;
          ClearLine();

          if (parsed_bytes >= bytes) {
            return dst;
          }
        }
      }
      if (buffer_ptr_ < buffer_ptr) {
        ReadLineFromBuffer(buffer_ptr - buffer_ptr_);
      }

      if (is_eof) {
        break;
      }
      is_eof = ReadBuffer();
    }

    if (line_ptr_ != 0) {
      (this->*add_overlap_to_list)(dst, shorten_names);
      parsed_bytes += line_ptr_;
      ClearLine();
    }

    return dst;
  }

  std::unique_ptr<OverlapParser> OverlapParser::Create(const std::string& path) {
    auto file = gzopen(path.c_str(), "r");
    if (file == nullptr) {
      throw std::invalid_argument("[bioparser::Parser::Create] error: unable to open file " + path);
    }
    return std::unique_ptr<OverlapParser>(new OverlapParser(file));
  }

  void OverlapParser::ResetFilePointer() {
    gzseek(file_.get(), 0, SEEK_SET);
    buffer_ptr_ = 0;
    buffer_bytes_ = 0;
  }

  bool OverlapParser::ReadBuffer() {
    buffer_ptr_ = 0;
    buffer_bytes_ = gzread(file_.get(), buffer_.data(), buffer_.size());
    return buffer_bytes_ < buffer_.size();
  }

  void OverlapParser::ReadLineFromBuffer(std::uint32_t count, bool strip) {
    if (buffer_ptr_ + count > buffer_.size()) {
      throw std::invalid_argument(
          "[bioparser::Parser::ReadLineFromBuffer] error: buffer overflow");
    }
    if (line_ptr_ + count > line_.size()) {
      line_.resize(2 * line_.size());
    }
    std::memcpy(&line_[line_ptr_], &buffer_[buffer_ptr_], count);
    line_ptr_ += strip ? RightStrip(&line_[line_ptr_], count) : count;
    buffer_ptr_ += count + 1;  // ignore sought character
  }

  void OverlapParser::TerminateLine(std::uint32_t i) {
    line_[i] = '\0';
  }

  void OverlapParser::ClearLine() {
    line_ptr_ = 0;
  }

  std::uint32_t OverlapParser::RightStrip(const char* str, std::uint32_t str_len) {
    while (str_len > 0 && std::isspace(str[str_len - 1])) {
      --str_len;
    }
    return str_len;
  }

  std::uint32_t OverlapParser::Shorten(const char* str, std::uint32_t str_len) {
    for (std::uint32_t i = 0; i < str_len; ++i) {
      if (std::isspace(str[i])) {
        return i;
      }
    }
    return str_len;
  }
} // raven
