#pragma once

#include <string>
#include "interval.hpp"
#include <map>
#include "variant.hpp"
#include "cigar.hpp"
#include <cmath>

namespace hc
{

struct Haplotype
{
    std::string bases;
    std::map<std::size_t, Variant> event_map;
    Cigar cigar;
    std::size_t alignment_begin_wrt_ref = 0;
    double score = std::numeric_limits<double>::lowest();
    std::size_t rank = 0;

    Haplotype() = default;
    Haplotype(std::string bases, double score)
        : bases(std::move(bases)), score(score) {}

    auto size() const noexcept { return bases.size(); }

    auto get_overlapping_events(std::size_t begin) const
    {
        std::vector<Variant> events;
        auto upper_bound = event_map.upper_bound(begin);
        for (auto it = event_map.begin(); it != upper_bound; ++it)
            if (it->second.location.end > begin)
                events.push_back(it->second);
        return events;
    }
};

} // hc
