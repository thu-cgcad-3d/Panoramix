#include "../src/core/views_net.hpp"
#include "../src/vis/views_net_visualize.hpp"
#include "../src/vis/regions_net_visualize.hpp"
#include "../src/vis/qt_resources.hpp"
#include "gtest/gtest.h"

#include <iostream>
#include <string>
#include <random>

#include <QApplication>

using namespace panoramix;

// PROJECT_TEST_DATA_DIR_STR is predefined using CMake
static const std::string ProjectTestDataDirStr = PROJECT_TEST_DATA_DIR_STR;

TEST(ViewsNet, ViewsNet) {
    /// get perspective pieces

    cv::Mat panorama = cv::imread(ProjectTestDataDirStr + "/panofactory.jpg");

    cv::resize(panorama, panorama, cv::Size(1000, 500));
    vis::Visualizer2D(panorama) << vis::manip2d::Show();
    core::PanoramicCamera originCam(panorama.cols / M_PI / 2.0);

    std::vector<core::PerspectiveCamera> cams(10);
    std::generate(cams.begin(), cams.end(), [](){
        return core::PerspectiveCamera(500, 500, 200, 
            core::Vec3(0, 0, 0), 
            core::Vec3((rand() % 10000)-5000, (rand() % 10000)-5000, (rand() % 10000-5000)),
            core::Vec3(0, 0, -1));
    });
    
    std::vector<core::Image> ims(cams.size());
    std::transform(cams.begin(), cams.end(), ims.begin(),
        [&panorama, &originCam](const core::PerspectiveCamera & pcam){
        core::Image im;
        return core::CameraSampler<core::PerspectiveCamera, core::PanoramicCamera>(pcam, originCam)(panorama);
    });

    for (auto & im : ims) {
        vis::Visualizer2D(im) << vis::manip2d::Show();
    }

    /// insert into views net

    core::ViewsNet net;
    for (int i = 0; i < cams.size(); i++){
        qDebug() << "photo: " << i;

        auto & camera = cams[i];
        const auto & im = ims[i];
        auto viewHandle = net.insertPhoto(im, camera);
       
        qDebug() << "extracting features ...";

        net.computeFeatures(viewHandle);

        vis::Visualizer2D(im)
            << net.views().data(viewHandle).SIFTs
            << net.views().data(viewHandle).SURFs
            << vis::manip2d::Show();

        vis::Visualizer2D() 
            << *(net.views().data(viewHandle).regionNet) 
            << vis::manip2d::Show();

        vis::Visualizer2D (im)
            << vis::manip2d::SetColor(core::Color(0, 0, 255))
            << vis::manip2d::SetThickness(2)
            << net.views().data(viewHandle).lineSegments
            << vis::manip2d::SetColor(core::Color(255, 0, 0))
            << vis::manip2d::SetThickness(1)
            << net.views().data(viewHandle).lineSegmentIntersections
            << vis::manip2d::Show();

        net.updateConnections(viewHandle);
        net.computeTransformationOnConnections(viewHandle);
        net.calibrateCamera(viewHandle);
        net.calibrateAllCameras();

        if (net.isTooCloseToAnyExistingView(viewHandle).isValid()){
            continue;
        }

        qDebug() << "calibrating camera and classifying lines ...";

        // estimate vanishing points and classify lines
        net.estimateVanishingPointsAndClassifyLines();
        auto vps = net.globalData().vanishingPoints;
        for (auto & vp : vps)
            vp /= cv::norm(vp);
        double ortho = cv::norm(core::Vec3(vps[0].dot(vps[1]), vps[1].dot(vps[2]), vps[2].dot(vps[0])));
        EXPECT_LT(ortho, 1e-1);
        
        auto antivps = vps;
        for (auto & p : antivps)
            p = -p;

        std::vector<core::Vec3> allvps(vps.begin(), vps.end());
        allvps.insert(allvps.end(), antivps.begin(), antivps.end());

        std::vector<core::Point2> vp2s(allvps.size());
        std::transform(allvps.begin(), allvps.end(), vp2s.begin(),
            [&originCam](const core::Vec3 & p3){
            return originCam.screenProjection(p3);
        });

        vis::Visualizer2D(panorama)
            << vis::manip2d::SetWindowName("Vanishing points")
            << vis::manip2d::SetColor(core::Color(0, 0, 255))
            << vis::manip2d::SetThickness(3)
            << vp2s
            << vis::manip2d::Show();

        vis::Visualizer2D()
            << net.views().data(viewHandle)
            << vis::manip2d::Show();

        net.rectifySpatialLines();
        vis::Visualizer3D()
            << vis::manip3d::SetCamera(core::PerspectiveCamera(700, 700, 200, core::Vec3(1, 1, 1), core::Vec3(0, 0, 0), core::Vec3(0, 0, -1)))
            << vis::manip3d::SetBackgroundColor(core::Color(200, 200, 200))
            << vis::manip3d::SetLineWidth(2.0)
            << vis::manip3d::SetColorTableDescriptor(core::ColorTableDescriptor::RGB)
            << vis::manip3d::SetRanderMode(vis::RenderModeFlag::Lines)
            << net.globalData()
            << vis::manip3d::AutoSetCamera
            << vis::manip3d::Show();

    }
    
}

int main(int argc, char * argv[], char * envp[])
{
    srand(clock());
    testing::InitGoogleTest(&argc, argv);
    vis::InitGui(argc, argv);
    return RUN_ALL_TESTS();
}
