#include "../src/rec/regions_net.hpp"
#include "../src/vis/visualize2d.hpp"
#include "../src/rec/regions_net_visualize.hpp"
#include "../src/rec/lines_net.hpp"

#include "test_config.hpp"

using namespace panoramix;
using namespace test;

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

TEST(LinesNet, LinesNet){

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
