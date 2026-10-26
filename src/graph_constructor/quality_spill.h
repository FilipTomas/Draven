#ifndef RAVEN_QUALITY_SPILL_H
#define RAVEN_QUALITY_SPILL_H

#include <cstdint>
#include <cstdio>
#include <string>
#include <vector>

#include "biosoup/nucleic_acid.hpp"

namespace raven {

// Per-base qualities written to a file as reads are parsed and read back one
// read at a time when it is sketched (the only step that reads them on the
// correction path), so they are not all in memory while k-mers are counted.
class QualitySpill {
 public:
  // The file goes to dir (which must exist), named after the process id.
  explicit QualitySpill(const std::string& dir);
  ~QualitySpill();  // closes and removes the file

  QualitySpill(const QualitySpill&) = delete;
  QualitySpill& operator=(const QualitySpill&) = delete;

  // Write the qualities of the next read (the first Add is read 0) and free
  // them. Not thread-safe.
  void Add(biosoup::NucleicAcid& read);

  // Flush the file; Load may then be called from any thread.
  void Finish();
  std::uint64_t bytes_written() const { return offsets_.back(); }
  const std::string& path() const { return path_; }

  // Put the qualities of read `index` back into read.quality.
  void Load(std::uint32_t index, biosoup::NucleicAcid& read) const;

 private:
  std::string path_;
  std::FILE* file_ = nullptr;
  int fd_ = -1;
  std::vector<char> buffer_;  // stdio buffer of file_
  // where read i's qualities start; the last entry is the end of the file
  std::vector<std::uint64_t> offsets_;
};

}  // namespace raven

#endif  // RAVEN_QUALITY_SPILL_H
