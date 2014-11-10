#include "../src/vis/visualize2d.hpp"
#include "gtest/gtest.h"

#include <iostream>
#include <string>
#include <random>

#include "test_config.hpp"

using namespace panoramix;
using namespace test;

TEST(Visualizer2D, Visualizer2D) {

    using namespace core;
    using namespace vis;

    std::vector<Visualizer2D> vis(1);

    {
        Image im = cv::imread(ProjectDataDirStrings::PanoramaOutdoor + "/panohk.png");
        vis[0].setImage(im);
        vis[0] << Line2({ 0, 10 }, { 100, 500 });
    }

    auto v = vis.front();
    vis.clear();
    Image lim = v.image();
    cv::imshow("Vis", lim);
    cv::waitKey();

    auto red = vis::ColorFromTag(vis::ColorTag::Red);
    ImageWithType<Vec3b> im = ImageWithType<Vec3b>::zeros(100, 100);
    for (auto & p : im) {
        p = Vec3b(red[0], red[1], red[2]);
    }
    cv::line(im, PixelLoc(0, 0), PixelLoc(100, 100), vis::ColorFromTag(vis::ColorTag::Blue), 2);
    cv::imshow("Red", im);
    cv::waitKey();
}

int main(int argc, char * argv[], char * envp[])
{
    srand(clock());
    testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}
