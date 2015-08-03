#ifndef PANORAMIX_TESTS_CONFIG_HPP
#define PANORAMIX_TESTS_CONFIG_HPP
 
#include <iostream>
#include <chrono>
#include <random>

#include "gtest/gtest.h"


namespace pano {
    namespace test {

        struct ProjectDataDirStrings {
            static std::string Base;
            static std::string Normal;
            static std::string PanoramaIndoor;
            static std::string PanoramaOutdoor;
            static std::string Serialization;
            static std::string LocalManhattan;
            static std::string MeshSMF;
            static std::string Scripts;
            static std::string BPTests;
        };

    }
}

#endif