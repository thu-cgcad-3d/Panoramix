#include "../src/core/regions_net.hpp"
#include "../src/vis/feature_visualize.hpp"
#include "../src/vis/regions_net_visualize.hpp"
#include "gtest/gtest.h"

#include <iostream>
#include <string>
#include <random>

using namespace panoramix;

// PROJECT_TEST_DATA_DIR_STR is predefined using CMake
static const std::string ProjectTestDataDirStr = PROJECT_TEST_DATA_DIR_STR;

TEST(RegionsNet, RegionsNet) {

    for (int i = 0; i < 4; i++) {
        std::string name = ProjectTestDataDirStr + "/" + "sampled_" + std::to_string(i) + ".png";
        cv::Mat im = cv::imread(name);
        core::RegionsNet regNet(im);
        regNet.buildNetAndComputeGeometricFeatures();
        regNet.computeImageFeatures();
        vis::ImageFeatureVisualizer()
            << regNet
            << vis::manip::Show();
    }
    
}



int main(int argc, char * argv[], char * envp[])
{
    srand(clock());
    testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}
