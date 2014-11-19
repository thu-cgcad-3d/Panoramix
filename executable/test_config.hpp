#ifndef TEST_CONFIG_HPP
#define TEST_CONFIG_HPP
 
#include <iostream>
#include <chrono>
#include <random>

#include "gtest/gtest.h"

#define DEBUG_TEST(...) void run()
#define DEBUG_RUN_ALL_TESTS() (run(),0)

#define DEBUG_ASSERT_EQ(a, b) assert((a) == (b))

namespace test {

    namespace ProjectDataDirStrings {
        static const std::string Base = PROJECT_TEST_DATA_DIR_STR;
        static const std::string Normal = PROJECT_TEST_DATA_DIR_STR"/normal";
        static const std::string PanoramaIndoor = PROJECT_TEST_DATA_DIR_STR"/panorama/indoor";
        static const std::string PanoramaOutdoor = PROJECT_TEST_DATA_DIR_STR"/panorama/outdoor";
        static const std::string Serialization = PROJECT_TEST_DATA_DIR_STR"/serialization";
        static const std::string LocalManhattan = PROJECT_TEST_DATA_DIR_STR"/localmanh";
        static const std::string MeshSMF = PROJECT_TEST_DATA_DIR_STR"/mesh/smf";
        static const std::string Scripts = PROJECT_TEST_DATA_DIR_STR"/scripts";
    }
 
}

#endif