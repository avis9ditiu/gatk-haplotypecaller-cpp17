#pragma once

#include <string>
#include "../utils/interval.hpp"
#include <map>
#include "../variant/variant.hpp"
#include "../sam/cigar.hpp"
#include <cmath>
#include <boost/serialization/string.hpp>
#include <boost/serialization/map.hpp>

namespace hc
{

struct Haplotype
{
    std::string bases;
    Interval location;
    std::map<std::size_t, Variant> event_map;
    Cigar cigar;
    std::size_t alignment_begin_wrt_ref = 0;
    double score = std::numeric_limits<double>::lowest();
    std::size_t rank;

    Haplotype() = default;
    Haplotype(std::string bases, double score)
        : bases(std::move(bases)), score(score) {}

    std::size_t size() const noexcept { return bases.size(); }

    auto get_overlapping_events(std::size_t begin) const
    {
        std::vector<Variant> events;
        auto upper_bound = event_map.upper_bound(begin);
        for (auto it = event_map.begin(); it != upper_bound; ++it)
            if (it->second.location.end > begin)
                events.push_back(it->second);
        return events;
    }

    template <typename Archive>
    void serialize(Archive& ar, const unsigned int version)
    {
        ar & bases;
        ar & location;
        ar & event_map;
        ar & cigar;
        ar & alignment_begin_wrt_ref;
        ar & score;
        ar & rank;
    }
};

} // hc