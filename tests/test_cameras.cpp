#include "../src/core/cameras.hpp"
#include "../src/gui/visualize2d.hpp"

#include "config.hpp"

using namespace panoramix;
using namespace test;

static_assert(core::IsCamera<core::PerspectiveCamera>::value, "");
static_assert(core::IsCamera<core::PanoramicCamera>::value, "");
static_assert(!core::IsCamera<core::Line3>::value, "");

TEST(Camera, PerspectiveCamera){
    core::PerspectiveCamera cam(1000, 1000, core::Point2(500, 500), 500, 
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

TEST(Camera, PerspectiveCameraRandom){
    for (int k = 0; k < 100; k++) {
        int w = abs(rand()) % 500;
        int h = abs(rand()) % 400;
        core::PerspectiveCamera cam(w, h, core::Point2(w, h)/2.0, abs(rand()) % 600,
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


TEST(Camera, CameraSampler) {
    cv::Mat im = cv::imread(ProjectDataDirStrings::PanoramaIndoor + "/13.jpg");

    EXPECT_EQ(2000, im.cols);
    EXPECT_EQ(1000, im.rows);
    cv::resize(im, im, cv::Size(1000, 500));
    gui::Visualizer2D viz(im);
    viz.params.alphaForNewImage = 0.5;

    core::PanoramicCamera originCam(im.cols / M_PI / 2.0);
    core::PanoramicCamera newCam(im.cols / M_PI / 2.0, 
        core::Vec3(0, 0, 0), 
        core::Vec3(0, 0, 1),
        core::Vec3(0, 1, 0));
    viz << core::MakeCameraSampler(newCam, originCam)(im)
        << gui::manip2d::Show(false);

    core::PartialPanoramicCamera newCam2(newCam);
    gui::Visualizer2D(core::MakeCameraSampler(newCam2, originCam)(im))
        << gui::manip2d::Show();


    float camPositions[4][3] = {
        {1, 0, 0},
        {0, 1, 0},
        {-1, 0, 0},
        {0, -1, 0}
    };

    for (int i = 0; i < 4; i++){
        core::PerspectiveCamera cam(500, 600, core::Point2(250, 300), 150,
            core::Vec3(0, 0, 0),
            core::Vec3(camPositions[i][0], camPositions[i][1], camPositions[i][2]),
            core::Vec3(0, 0, -1));
        core::CameraSampler<core::PerspectiveCamera, core::PanoramicCamera> sampler(cam,
            core::PanoramicCamera(im.cols / M_PI / 2.0));
        cv::Mat sampledIm = sampler(im);
        cv::imwrite(ProjectDataDirStrings::Normal + "/13-" + std::to_string(i) + ".png", sampledIm);
        gui::Visualizer2D(sampledIm) << gui::manip2d::Show();
    }
}