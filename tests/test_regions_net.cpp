#include "../src/rec/regions_net.hpp"
#include "../src/vis/visualize2d.hpp"
#include "../src/rec/regions_net_visualize.hpp"
#include "gtest/gtest.h"

#include <iostream>
#include <string>
#include <random>

using namespace panoramix;

// PROJECT_TEST_DATA_DIR_STR is predefined using CMake
static const std::string ProjectTestDataDirStr = PROJECT_TEST_DATA_DIR_STR;
static const std::string ProjectTestDataDirStr_Normal = ProjectTestDataDirStr + "/normal";
static const std::string ProjectTestDataDirStr_PanoramaIndoor = ProjectTestDataDirStr + "/panorama/indoor";
static const std::string ProjectTestDataDirStr_PanoramaOutdoor = ProjectTestDataDirStr + "/panorama/outdoor";

TEST(RegionsNet, RegionsNet) {

    for (int i = 0; i < 4; i++) {
        std::string name = ProjectTestDataDirStr_Normal + "/" + "sampled_" + std::to_string(i) + ".png";
        cv::Mat im = cv::imread(name);
        rec::RegionsNet regNet(im);
        regNet.buildNetAndComputeGeometricFeatures();
        regNet.computeImageFeatures();
        vis::Visualizer2D()
            << regNet
            << vis::manip2d::Show();
    }
    
}



int main(int argc, char * argv[], char * envp[])
{
    srand(clock());
    testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}
