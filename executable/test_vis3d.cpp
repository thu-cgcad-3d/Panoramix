#include "../src/vis/renderable_object_tree.hpp"
//#include "../src/vis/visualize3d.hpp"
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

TEST(Visualizer3D, RenderableObject) { 

    std::vector<core::Line3> lines = {
        { {-1, 1}, {5, 5} },
        { { -2, -2 }, {8, 6} }
    };
    auto lines1 = vis::MakeRenderable(lines);
    vis::RenderableObjectTree renderTree(lines1);


}

int main(int argc, char * argv[], char * envp[])
{
    srand(clock());
    testing::InitGoogleTest(&argc, argv);
    //testing::FLAGS_gtest_filter = "Texture";
    return RUN_ALL_TESTS();
    //run();
    //return 0;
}
