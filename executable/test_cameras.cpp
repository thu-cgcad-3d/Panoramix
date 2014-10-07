#include "../src/core/cameras.hpp"
#include "../src/vis/visualize2d.hpp"

#include "test_config.hpp"

using namespace panoramix;
using namespace test;

TEST(Feature, PerspectiveCamera){
    core::PerspectiveCamera cam(1000, 1000, 500, 
        core::Vec3(0, 0, 0), 
        core::Vec3(1, 0, 0));

    for (int i = 0; i < 100; i++){
        core::Vec2 v(abs(rand()), abs(rand()));
        auto p = cam.spatialDirection(v);
        auto v2 = cam.screenProjection(p);
        double dist = core::norm(v - v2);
        ASSERT_LT(dist, 0.01);
    }
    auto c = cam.screenProjection(cam.center());
    double dist = core::norm(c - core::Vec2(cam.screenSize().width / 2, cam.screenSize().height / 2));
    ASSERT_LT(dist, 2);
}

TEST(Feature, PerspectiveCameraRandom){
    for (int k = 0; k < 100; k++) {
        core::PerspectiveCamera cam(abs(rand()) % 500, abs(rand()) % 400, abs(rand()) % 600,
            core::Vec3(rand(), rand(), rand()),
            core::Vec3(rand(), rand(), rand()));
        for (int i = 0; i < 100; i++){
            core::Vec2 v(rand(), rand());
            auto p = cam.spatialDirection(v);
            auto v2 = cam.screenProjection(p);
            double dist = core::norm(v - v2);
            if (!std::isnan(dist) && !std::isinf(dist))
                ASSERT_LT(dist, 0.01);
        }
        auto c = cam.screenProjection(cam.center());
        double dist = core::norm(c - core::Vec2(cam.screenSize().width / 2, cam.screenSize().height / 2));
        if (!std::isnan(dist) && !std::isinf(dist))
            ASSERT_LT(dist, 2);
    }
}


TEST(Feature, CameraSampler) {
    cv::Mat im = cv::imread(ProjectDataDirStrings::PanoramaIndoor + "/13.jpg");

    EXPECT_EQ(2000, im.cols);
    EXPECT_EQ(1000, im.rows);
    cv::resize(im, im, cv::Size(1000, 500));
    vis::Visualizer2D viz(im);
    viz.params.alphaForNewImage = 0.5;

    core::PanoramicCamera originCam(im.cols / M_PI / 2.0);
    core::PanoramicCamera newCam(im.cols / M_PI / 2.0, 
        core::Vec3(0, 0, 0), 
        core::Vec3(0, 0, 1),
        core::Vec3(0, 1, 0));
    viz << core::CameraSampler<core::PanoramicCamera, core::PanoramicCamera>(newCam, originCam)(im)
        << vis::manip2d::Show();

    float camPositions[4][3] = {
        {1, 0, 0},
        {0, 1, 0},
        {-1, 0, 0},
        {0, -1, 0}
    };

    for (int i = 0; i < 4; i++){
        core::PerspectiveCamera cam(500, 600, 150,
            core::Vec3(0, 0, 0),
            core::Vec3(camPositions[i][0], camPositions[i][1], camPositions[i][2]),
            core::Vec3(0, 0, -1));
        core::CameraSampler<core::PerspectiveCamera, core::PanoramicCamera> sampler(cam,
            core::PanoramicCamera(im.cols / M_PI / 2.0));
        cv::Mat sampledIm = sampler(im);
        cv::imwrite(ProjectDataDirStrings::Normal + "/13-" + std::to_string(i) + ".png", sampledIm);
        vis::Visualizer2D(sampledIm) << vis::manip2d::Show();
    }
}


int main(int argc, char * argv[], char * envp[])
{
    testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}

