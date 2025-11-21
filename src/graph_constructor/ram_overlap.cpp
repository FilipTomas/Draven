#include "ram_overlap.h"
#include <algorithm>



extended_overlap ram_overlap::resolve() {
    // -------- helpers --------
    auto normalize = [](biosoup::Overlap& o) {
        if (o.lhs_end < o.lhs_begin) std::swap(o.lhs_begin, o.lhs_end);
        if (o.rhs_end < o.rhs_begin) std::swap(o.rhs_begin, o.rhs_end);
    };

    auto qlen = [](const extended_overlap& e) -> std::uint32_t {
        auto q0 = std::min(e.overlap.lhs_begin, e.overlap.lhs_end);
        auto q1 = std::max(e.overlap.lhs_begin, e.overlap.lhs_end);
        return (q1 >= q0) ? (q1 - q0) : 0u;
    };

    auto interval_gap = [](std::uint32_t a0, std::uint32_t a1,
                           std::uint32_t b0, std::uint32_t b1) -> std::uint32_t {
        if (a0 > a1) std::swap(a0, a1);
        if (b0 > b1) std::swap(b0, b1);
        if (a1 >= b0 && b1 >= a0) return 0;           // overlap or touch
        if (a1 < b0) return b0 - a1;                  // a before b
        return a0 - b1;                                // b before a
    };

    auto rect_gap_2d = [&](const biosoup::Overlap& A, const biosoup::Overlap& B) -> std::uint32_t {
        std::uint32_t qgap = interval_gap(A.lhs_begin, A.lhs_end, B.lhs_begin, B.lhs_end);
        std::uint32_t tgap = interval_gap(A.rhs_begin, A.rhs_end, B.rhs_begin, B.rhs_end);
        return std::max(qgap, tgap); // Chebyshev
    };

    auto overlap_both_axes = [&](const biosoup::Overlap& a, const biosoup::Overlap& b,
                                 std::int64_t tol = 0) -> bool {
        auto q_ok = (std::max<std::int64_t>(a.lhs_begin, b.lhs_begin) <=
                     std::min<std::int64_t>(a.lhs_end,   b.lhs_end)   + tol);
        auto t_ok = (std::max<std::int64_t>(a.rhs_begin, b.rhs_begin) <=
                     std::min<std::int64_t>(a.rhs_end,   b.rhs_end)   + tol);
        return q_ok && t_ok;
    };

    auto merge_into = [](biosoup::Overlap& anchor, const biosoup::Overlap& other) {
        anchor.lhs_begin = std::min(anchor.lhs_begin, other.lhs_begin);
        anchor.lhs_end   = std::max(anchor.lhs_end,   other.lhs_end);
        anchor.rhs_begin = std::min(anchor.rhs_begin, other.rhs_begin);
        anchor.rhs_end   = std::max(anchor.rhs_end,   other.rhs_end);
        anchor.score    += other.score; // adjust if score not additive
    };

    // -------- input guard / normalize --------
    if (overlap_fragments.empty()) {
        n_fragments = 0;
        fragment_gaps.clear();
        return extended_overlap{};
    }
    for (auto& e : overlap_fragments) normalize(e.overlap);

    if (overlap_fragments.size() == 1) {
        const auto& o = overlap_fragments.front().overlap;
        lhs_begin = o.lhs_begin; lhs_end = o.lhs_end;
        rhs_begin = o.rhs_begin; rhs_end = o.rhs_end;
        n_fragments = 1;
        fragment_gaps.clear();
        return overlap_fragments.front();
    }

    // -------- sort by query start --------
    std::sort(overlap_fragments.begin(), overlap_fragments.end(),
              [](const extended_overlap& a, const extended_overlap& b) {
                  return a.overlap.lhs_begin < b.overlap.lhs_begin;
              });

    // -------- REMOVE CONTAINED (kept) --------
    auto is_contained_q = [](const biosoup::Overlap& a, const biosoup::Overlap& b) {
        return a.lhs_begin >= b.lhs_begin && a.lhs_end <= b.lhs_end;
    };
    {
        std::vector<extended_overlap> kept;
        kept.reserve(overlap_fragments.size());
        for (const auto& ov : overlap_fragments) {
            bool contained = false;
            for (const auto& kv : kept) {
                if (is_contained_q(ov.overlap, kv.overlap)) { contained = true; break; }
            }
            if (!contained) kept.emplace_back(ov);
        }
        overlap_fragments.swap(kept);
    }
    if (overlap_fragments.empty()) {
        n_fragments = 0;
        fragment_gaps.clear();
        return extended_overlap{};
    }

    // -------- MERGE (overlap OR within proximity) --------
    constexpr std::uint32_t LINK_TOL = 5000;  // set as you like

    // Parallel cluster membership to compute gaps later
    std::vector<std::vector<extended_overlap>> clusters(overlap_fragments.size());
    for (std::size_t i = 0; i < overlap_fragments.size(); ++i)
        clusters[i].push_back(overlap_fragments[i]);

    for (std::size_t i = 0; i < overlap_fragments.size(); ) {
        auto& anchor = overlap_fragments[i].overlap;

        for (std::size_t j = i + 1; j < overlap_fragments.size(); ) {
            const auto& cand = overlap_fragments[j].overlap;

            const bool ovlp = overlap_both_axes(anchor, cand, /*tol=*/0);
            const std::uint32_t d = rect_gap_2d(anchor, cand);

            if (ovlp || d <= LINK_TOL) {
                merge_into(anchor, cand);
                clusters[i].insert(clusters[i].end(), clusters[j].begin(), clusters[j].end());
                overlap_fragments.erase(overlap_fragments.begin() + j);
                clusters.erase(clusters.begin() + j);
            } else {
                ++j;
            }
        }
        ++i;
    }

    // -------- pick longest merged union --------
    auto it_longest = std::max_element(
        overlap_fragments.begin(), overlap_fragments.end(),
        [&](const extended_overlap& a, const extended_overlap& b){ return qlen(a) < qlen(b); }
    );
    if (it_longest == overlap_fragments.end()) {
        fragment_gaps.clear();
        return extended_overlap{};
    }
    std::size_t idx_longest = static_cast<std::size_t>(
        std::distance(overlap_fragments.begin(), it_longest)
    );

    // -------- compute fragment gaps for the chosen cluster --------
    auto& members = clusters[idx_longest];

    // Sort members by query start, then end
    std::sort(members.begin(), members.end(),
        [](const extended_overlap& a, const extended_overlap& b){
            auto a0 = std::min(a.overlap.lhs_begin, a.overlap.lhs_end);
            auto b0 = std::min(b.overlap.lhs_begin, b.overlap.lhs_end);
            if (a0 != b0) return a0 < b0;
            auto a1 = std::max(a.overlap.lhs_begin, a.overlap.lhs_end);
            auto b1 = std::max(b.overlap.lhs_begin, b.overlap.lhs_end);
            return a1 < b1;
        });

    fragment_gaps.clear();
    if (!members.empty()) {
        auto prev_q0 = std::min(members[0].overlap.lhs_begin, members[0].overlap.lhs_end);
        auto prev_q1 = std::max(members[0].overlap.lhs_begin, members[0].overlap.lhs_end);
        auto prev_t0 = std::min(members[0].overlap.rhs_begin, members[0].overlap.rhs_end);
        auto prev_t1 = std::max(members[0].overlap.rhs_begin, members[0].overlap.rhs_end);

        for (std::size_t k = 1; k < members.size(); ++k) {
            auto cur_q0 = std::min(members[k].overlap.lhs_begin, members[k].overlap.lhs_end);
            auto cur_q1 = std::max(members[k].overlap.lhs_begin, members[k].overlap.lhs_end);
            auto cur_t0 = std::min(members[k].overlap.rhs_begin, members[k].overlap.rhs_end);
            auto cur_t1 = std::max(members[k].overlap.rhs_begin, members[k].overlap.rhs_end);

            // Record a gap only if BOTH axes have a real gap
            if (cur_q0 > prev_q1 && cur_t0 > prev_t1) {
                fragment_gaps.emplace_back(
                    static_cast<std::uint16_t>(prev_q1),  // q gap start
                    static_cast<std::uint16_t>(prev_t1)   // t gap start
                );
                fragment_gaps.emplace_back(
                    static_cast<std::uint16_t>(cur_q0),   // q gap end
                    static_cast<std::uint16_t>(cur_t0)    // t gap end
                );
            }
            prev_q1 = std::max(prev_q1, cur_q1);
            prev_t1 = std::max(prev_t1, cur_t1);
        }

        // Sort gaps by (q_start, then t_start)
        std::stable_sort(fragment_gaps.begin(), fragment_gaps.end(),
            [](const auto& a, const auto& b){
                if (a.first != b.first) return a.first < b.first;
                return a.second < b.second;
            });
    }

    // -------- update class fields & return --------
    normalize(it_longest->overlap);

    lhs_begin = it_longest->overlap.lhs_begin;
    lhs_end   = it_longest->overlap.lhs_end;
    rhs_begin = it_longest->overlap.rhs_begin;
    rhs_end   = it_longest->overlap.rhs_end;

    n_fragments = static_cast<std::uint16_t>(members.size());

    return *it_longest;
}


// extended_overlap ram_overlap::resolve() {
//     // is a contained in b
//     auto is_contained = [](const biosoup::Overlap& a, const biosoup::Overlap& b) -> bool {
//         return (a.lhs_begin >= b.lhs_begin && a.lhs_end <= b.lhs_end); //&&
//                //(a.rhs_begin >= b.rhs_begin && a.rhs_end <= b.rhs_end);
//     };

//     std::sort(overlap_fragments.begin(), overlap_fragments.end(),
//         [](const extended_overlap& lhs, const extended_overlap& rhs) -> bool {
//             return lhs.overlap.lhs_begin < rhs.overlap.lhs_begin;
//         }
//     );

//     // check if two overlaps are overlapping, also check for overlaps on the target side
//     auto is_overlapping = [](const biosoup::Overlap& a, const biosoup::Overlap& b) -> bool {
//         if((std::max(a.lhs_begin, b.lhs_begin) <= std::min(a.lhs_end, b.lhs_end)) &&
//                (std::max(a.rhs_begin, b.rhs_begin) <= std::min(a.rhs_end, b.rhs_end))){
//                 return 1;
//             } else if((std::max(a.lhs_begin, b.lhs_begin) <= std::min(a.lhs_end, b.lhs_end)) &&
//                !(std::max(a.rhs_begin, b.rhs_begin) <= std::min(a.rhs_end, b.rhs_end))){
//                 return 2;
//             } else {
//                 return 0;
//                }
//     };

//     auto normalize = [](biosoup::Overlap& o) -> void {
//         if (o.lhs_end < o.lhs_begin) std::swap(o.lhs_begin, o.lhs_end);
//         if (o.rhs_end < o.rhs_begin) std::swap(o.rhs_begin, o.rhs_end);
//     };

//     auto axis_overlap = [](std::int64_t a0, std::int64_t a1,
//                         std::int64_t b0, std::int64_t b1,
//                         std::int64_t tol = 0) -> bool{
//     // overlap if max(start) <= min(end) + tol
//         return std::max(a0, b0) <= std::min(a1, b1) + tol;
//     };

//     auto merge_into = [](biosoup::Overlap& anchor,
//                        const biosoup::Overlap& other) -> void {
//     // Prefix/suffix (or general) merge: expand to cover both on both axes
//         anchor.lhs_begin = std::min(anchor.lhs_begin, other.lhs_begin);
//         anchor.lhs_end   = std::max(anchor.lhs_end,   other.lhs_end);
//         anchor.rhs_begin = std::min(anchor.rhs_begin, other.rhs_begin);
//         anchor.rhs_end   = std::max(anchor.rhs_end,   other.rhs_end);
//         anchor.score += other.score;
//     // If you keep additional fields (score, matches, length, etc.), update them here.
//     // e.g., anchor.score += other.score; or anchor.len = (anchor.lhs_end - anchor.lhs_begin) ...
//     };

//     enum class OverlapType : std::uint8_t {
//         None = 0,
//         QueryOnly = 2,
//         Both = 1
//     };

//     auto classify_overlap = [&axis_overlap](const biosoup::Overlap& a,
//                                     const biosoup::Overlap& b,
//                                     std::int64_t tol = 0) -> OverlapType {
//         bool q = axis_overlap(a.lhs_begin, a.lhs_end, b.lhs_begin, b.lhs_end, tol);
//         bool t = axis_overlap(a.rhs_begin, a.rhs_end, b.rhs_begin, b.rhs_end, tol);
//         if (q && t) return OverlapType::Both;
//         if (q && !t) return OverlapType::QueryOnly;
//         return OverlapType::None;
//     };

//     auto merge_prefix_suffix_overlaps = [&normalize, &merge_into, &classify_overlap](std::vector<extended_overlap>& overlap_fragments,
//                                   std::int64_t tol = 0) -> void {
//     // Optional: normalize all intervals once
//     for (auto& e : overlap_fragments) normalize(e.overlap);

//     // We’ll treat each item as a potential anchor and greedily merge any
//     // follower that overlaps on both axes.
//     for (std::size_t i = 0; i < overlap_fragments.size(); /* no ++ here */) {
//         auto& anchor = overlap_fragments[i].overlap;
//         bool merged_any = false;

//         for (std::size_t j = i + 1; j < overlap_fragments.size(); /* no ++ here */) {
//             auto type = classify_overlap(anchor, overlap_fragments[j].overlap, tol);

//             if (type == OverlapType::Both) {
//                 // Merge j into i, erase j, and keep scanning j at same index
//                 merge_into(overlap_fragments[i].overlap, overlap_fragments[j].overlap);
//                 overlap_fragments.erase(overlap_fragments.begin() + j);
//                 merged_any = true;
//                 // Do not increment j; the next element just slid into position j
//             } else {
//                 // Skip non-mergeable (QueryOnly or None) for now
//                 ++j;
//             }
//         }

//         // if (merged_any) {
//         //     // We expanded anchor; optionally, you could re-run normalization (not needed here).
//         //     // Keep i at same position only if you need to re-scan prior elements;
//         //     // since we only ever look forward, just advance i.
//         // }
//         ++i;
//     }
//     };


//     // resolve contained overlaps
//     std::vector<extended_overlap> resolved_fragments;
//     resolved_fragments.reserve(overlap_fragments.size());

//     resolved_fragments.emplace_back(overlap_fragments.front());
//     if (n_fragments == 1){
//         extended_overlap eo = resolved_fragments.front();
//         return eo;
//     } else {
//         for(auto& ov: overlap_fragments){
//             bool is_contained_flag = false;
//             for(auto& res_ov: resolved_fragments){
//                 if(is_contained(ov.overlap, res_ov.overlap)){
//                     is_contained_flag = true;
//                     break;
//                 }
//             }
//             if(!is_contained_flag){
//                 resolved_fragments.emplace_back(ov);
//             }
//         }
//     }

//     overlap_fragments.swap(resolved_fragments);
//     n_fragments = overlap_fragments.size();

    
//     extended_overlap eo = overlap_fragments.front();
//     return eo;

//     if(n_fragments == 1){
//         return overlap_fragments.front();
//     }
//     // // merge fragments
//     // // Resolve overlapping
//     merge_prefix_suffix_overlaps(overlap_fragments, 100); // tolerance of 100 bp


// };


// class ram_overlap {
// public:
//     std::vector<biosoup::Overlap> overlaps;
//     std::uint32_t lhs_id;
//     std::uint32_t rhs_id;

//     std::uint32_t lhs_begin;
//     std::uint32_t lhs_end;
//     std::uint32_t rhs_begin;
//     std::uint32_t rhs_end;

//     std::uint16_t n_fragments = 0; // number of fragments in the overlap
//     std::vector<std::uint16_t> fragment_gaps; //gaps between fragments

// };