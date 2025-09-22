#include "Driver.h"
#include "SimBuilder.h"

#include <cxxopts.hpp>

Pivot::Simulation::Scene ParseSceneName(std::string_view name)
{
    using namespace Pivot;
    static std::unordered_map<std::string, Simulation::Scene> const s_SceneFromName = {
        {"falling", Simulation::Scene::Falling},
        {"bigball", Simulation::Scene::BigBall},
    };
    if (auto iter = s_SceneFromName.find(name.data()); iter != s_SceneFromName.end())
    {
        return iter->second;
    }
    else
    {
        spdlog::critical("Failed to parse test case name");
        std::exit(EXIT_FAILURE);
    }
}

auto ParseArgs(int argc, char **argv)
{
    try
    {
        cxxopts::Options argParser("demo", "The demo of Particle-In-Cell liquid simulation");
        argParser.add_options()("n,dirname", "Directory name", cxxopts::value<std::string>()->default_value("output"))(
            "b,begin", "Begin frame", cxxopts::value<std::uint32_t>()->default_value("0"))(
            "e,end", "End frame", cxxopts::value<std::uint32_t>()->default_value("200"))(
            "t,test", "Test case", cxxopts::value<std::string>()->default_value("falling"))(
            "s,scale", "Size scale", cxxopts::value<int>()->default_value("128"))(
            "R,radius", "Particle radius", cxxopts::value<double>()->default_value("-1"))(
            "r,rate", "Frame rate", cxxopts::value<double>()->default_value("50"))(
            "c,cfl", "Courant number", cxxopts::value<double>()->default_value(".05"))("h,help", "Print usage");
        auto result = argParser.parse(argc, argv);
        if (result.count("help"))
        {
            std::cout << argParser.help() << std::endl;
            std::exit(EXIT_SUCCESS);
        }
        Pivot::DriverCreateOptions driverOpt = {
            .Dirname = result["dirname"].as<std::string>(),
            .BeginFrame = result["begin"].as<std::uint32_t>(),
            .EndFrame = result["end"].as<std::uint32_t>(),
            .FrameRate = result["rate"].as<double>(),
            .CourantNumber = result["cfl"].as<double>(),
        };
        Pivot::SimBuildOptions simOpt = {
            .Scene = ParseSceneName(result["test"].as<std::string>()),
            .Scale = result["scale"].as<int>(),
            .ParticleRadius = result["radius"].as<double>(),
        };
        return std::pair(driverOpt, simOpt);
    }
    catch (cxxopts::exceptions::exception const &e)
    {
        spdlog::critical("Failed to parse command line: {}", e.what());
        std::exit(EXIT_FAILURE);
    }
}

int main(int argc, char **argv)
{
    // Initialize logger
    spdlog::set_pattern("[%T] %^[%l]%$ %v");
    spdlog::set_level(spdlog::level::trace);
    spdlog::flush_on(spdlog::level::trace);
    // Main process
    auto [driverOpt, simOpt] = ParseArgs(argc, argv);
    auto driver = std::make_unique<Pivot::Driver>(driverOpt);
    auto simulation = Pivot::SimBuilder::Build(simOpt);
    driver->Run(simulation.get());

    return EXIT_SUCCESS;
}
