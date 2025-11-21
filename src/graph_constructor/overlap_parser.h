
#ifndef RAVEN_OVERLAP_PARSER_H
#define RAVEN_OVERLAP_PARSER_H

#include <memory>
#include <vector>
#include <cstring>
#include "zlib.h"

#define STORAGE_SIZE 65536

namespace raven {

  struct Overlap;

  class OverlapParser {
  public:
    OverlapParser(const OverlapParser &) = delete;
    OverlapParser &operator=(const OverlapParser &) = delete;
    OverlapParser(OverlapParser &&) = delete;
    OverlapParser &operator=(OverlapParser &&) = delete;

    OverlapParser(gzFile file)
      : file_(file, gzclose),
        buffer_(STORAGE_SIZE, 0),  // 64 kB
        buffer_ptr_(0),
        buffer_bytes_(0),
        line_(STORAGE_SIZE, 0),
        line_ptr_(0) {}

    virtual ~OverlapParser() {}

    static std::unique_ptr<OverlapParser> Create(const std::string &path);

    void ResetFilePointer();

    std::vector<std::unique_ptr<Overlap>> ParsePaf(std::uint64_t bytes, bool shorten_names = true);
    std::vector<std::unique_ptr<Overlap>> ParseSam(std::uint64_t bytes, bool shorten_names = true);

  private:
    std::vector<std::unique_ptr<Overlap>> Parse(std::uint64_t bytes, bool shorten_names,
                                          void (OverlapParser::*add_overlap_to_list)(std::vector<std::unique_ptr<Overlap>> &dst, bool shorten_names));

    void AddPafOverlapToList(std::vector<std::unique_ptr<Overlap>> &dst, bool shorten_names);
    void AddSamOverlapToList(std::vector<std::unique_ptr<Overlap>> &dst, bool shorten_names);

    bool ReadBuffer();

    void ReadLineFromBuffer(std::uint32_t count, bool strip = false);

    void TerminateLine(std::uint32_t i);

    void ClearLine();

    static std::uint32_t RightStrip(const char *str, std::uint32_t str_len);

    static std::uint32_t Shorten(const char *str, std::uint32_t str_len);

    std::unique_ptr<gzFile_s, int (*)(gzFile)> file_;
    std::vector<char> buffer_;
    std::uint32_t buffer_ptr_;
    std::uint32_t buffer_bytes_;
    std::vector<char> line_;
    std::uint32_t line_ptr_;
  };

} // raven

#endif //RAVEN_OVERLAP_PARSER_H
