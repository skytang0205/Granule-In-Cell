#pragma once

#include "Common.h"

namespace Pivot {
    class StopWatch {
    private:
        using Clock = std::chrono::steady_clock;

    public:
        explicit StopWatch(std::string_view name): m_Name { name }, m_BeginTime { Clock::now() } {}

        double Stop() {
            auto const sec = std::chrono::duration<double>(Clock::now() - m_BeginTime).count();
            if (auto it = s_Stats.find(m_Name); it != s_Stats.end()) {
                auto const i = it->second;
                s_AvgTime[i] += (sec - s_AvgTime[i]) / ++s_Counts[i];
                s_MaxTime[i] = std::max(s_MaxTime[i], sec);
            } else {
                s_Stats[m_Name] = s_Names.size();
                s_Names.push_back(m_Name);
                s_Counts.push_back(1);
                s_AvgTime.push_back(sec);
                s_MaxTime.push_back(sec);
            }
            return sec;
        }

        static void PrintStats();
        static void PrintStats(std::filesystem::path const & filename);

        static void Save(std::ostream & out);
        static void Load(std::istream & in);

    private:
        std::string       m_Name;
        Clock::time_point m_BeginTime;

        static inline std::unordered_map<std::string, std::size_t> s_Stats;
        static inline std::vector<std::string>                     s_Names;
        static inline std::vector<std::uint32_t>                   s_Counts;
        static inline std::vector<double>                          s_AvgTime;
        static inline std::vector<double>                          s_MaxTime;
    };
} // namespace Pivot