#include "StopWatch.h"

namespace Pivot {
    void StopWatch::PrintStats() {
        fmt::print(fmt::fg(fmt::color::yellow_green), "[Statistics]\n");
        for (std::size_t i = 0; i < s_Counts.size(); i++) {
            fmt::print(
                "{:>3}: {:>15}, {:>12} times,     average = {:.3f}s,     maximum = {:.3f}s\n",
                i + 1,
                s_Names[i],
                s_Counts[i],
                s_AvgTime[i],
                s_MaxTime[i]);
        }
    }

    void StopWatch::PrintStats(std::filesystem::path const & filename) {
        std::ofstream fout(filename);
        for (std::size_t i = 0; i < s_Counts.size(); i++) {
            fout << fmt::format(
                "{:>3}: {:>15}, {:>12} times,     average = {:.3f}s,     maximum = {:.3f}s",
                i + 1,
                s_Names[i],
                s_Counts[i],
                s_AvgTime[i],
                s_MaxTime[i])
                 << std::endl;
        }
    }

    void StopWatch::Save(std::ostream & out) {
        IO::Write(out, s_Stats.size());
        for (auto const & name : s_Names) { out << name << std::ends; }
        IO::Write(out, s_Counts);
        IO::Write(out, s_AvgTime);
        IO::Write(out, s_MaxTime);
    }

    void StopWatch::Load(std::istream & in) {
        std::size_t size;
        IO::Read(in, size);

        s_Names.resize(size);
        s_Stats.clear();
        for (std::size_t i = 0; i < size; i++) {
            std::getline(in, s_Names[i], '\0');
            s_Stats[s_Names[i]] = i;
        }

        s_Counts.resize(size);
        IO::Read(in, s_Counts);

        s_AvgTime.resize(size);
        IO::Read(in, s_AvgTime);

        s_MaxTime.resize(size);
        IO::Read(in, s_MaxTime);
    }
} // namespace Pivot