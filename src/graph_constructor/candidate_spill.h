#ifndef RAVEN_CANDIDATE_SPILL_H
#define RAVEN_CANDIDATE_SPILL_H

#include <cstdint>
#include <cstdio>
#include <functional>
#include <string>
#include <vector>

#include "correction_candidate.h"

namespace raven {

// Correction candidates written to files during mapping and read back one
// group of piles at a time, so selection and alignment hold one group in
// memory instead of every candidate. Group g holds the piles
// [First(g), First(g + 1)); its file keeps each pile's candidates in the order
// they were added. A record is the candidate's fields, its anchor count
// (32 bits) and its anchors.
class CandidateSpill {
 public:
  // Files go to dir (which must exist), named after the process id.
  CandidateSpill(const std::string& dir, std::uint32_t num_groups,
                 std::uint32_t num_piles);
  ~CandidateSpill();  // removes the files not read back

  CandidateSpill(const CandidateSpill&) = delete;
  CandidateSpill& operator=(const CandidateSpill&) = delete;

  std::uint32_t num_groups() const { return groups_.size(); }
  // the first pile of `group`; First(num_groups()) is the number of piles
  std::uint32_t First(std::uint32_t group) const;

  // Append a candidate to the group of its pile (c.lhs_id). Not thread-safe.
  void Add(const correction_candidate& c);

  // Flush and close every file; call once after the last Add.
  void Finish();
  std::uint64_t bytes_written() const { return bytes_written_; }
  const std::string& dir() const { return dir_; }

  // Append the candidates of `group` whose pile keep(pile) accepts to
  // piles[pile], in the order they were added, then delete its file.
  // Returns the number of bytes read.
  std::uint64_t Load(std::uint32_t group,
                     std::vector<std::vector<correction_candidate>>& piles,
                     const std::function<bool(std::uint32_t)>& keep);

 private:
  struct Group {
    std::string path;
    std::FILE* file = nullptr;
    std::vector<char> buffer;
  };

  std::uint32_t GroupOf(std::uint32_t pile) const;
  void Flush(Group& group);

  std::string dir_;
  std::uint32_t num_piles_;
  std::vector<Group> groups_;
  std::uint64_t bytes_written_ = 0;
};

}  // namespace raven

#endif  // RAVEN_CANDIDATE_SPILL_H
