#include "../src/core/mesh_maker.hpp"
#include "../src/core/utilities.hpp"
#include "../src/core/views_net.hpp"
#include "../src/vis/views_net_visualize.hpp"
#include "../src/vis/regions_net_visualize.hpp"
#include "../src/vis/qt_resources.hpp"
#include "gtest/gtest.h"

#include <iostream>
#include <string>
#include <random>
#include <thread>

#include <QApplication>

using namespace panoramix;

// PROJECT_TEST_DATA_DIR_STR is predefined using CMake
static const std::string ProjectTestDataDirStr = PROJECT_TEST_DATA_DIR_STR;

TEST(ViewsNet, ViewsNet) {
    cv::Mat panorama = cv::imread(ProjectTestDataDirStr + "/13.jpg");
    cv::resize(panorama, panorama, cv::Size(2000, 1000));

    //vis::Visualizer2D(panorama) << vis::manip2d::Show();
    core::PanoramicCamera originCam(panorama.cols / M_PI / 2.0);

    std::vector<core::PerspectiveCamera> cams;
    core::Mesh<core::Vec3> cameraStand;
    core::MakeQuadFacedSphere(cameraStand, 10, 20);
    for (auto & v : cameraStand.vertices()){
        core::Vec3 direction = v.data;
        if (core::AngleBetweenDirections(direction, core::Vec3(0, 0, 1)) <= 0.1 ||
            core::AngleBetweenDirections(direction, core::Vec3(0, 0, -1)) <= 0.1){
            //cams.emplace_back(700, 700, originCam.focal(), core::Vec3(0, 0, 0), direction, core::Vec3(0, 1, 0));
            continue;
        }
        else{
            cams.emplace_back(700, 700, originCam.focal(), core::Vec3(0, 0, 0), direction, core::Vec3(0, 0, -1));
        }
    }

    // sample photos
    std::vector<core::Image> ims(cams.size());
    std::transform(cams.begin(), cams.end(), ims.begin(),
        [&panorama, &originCam](const core::PerspectiveCamera & pcam){
        core::Image im;
        qDebug() << "sampling photo ...";
        return core::CameraSampler<core::PerspectiveCamera, core::PanoramicCamera>(pcam, originCam)(panorama);
    });

    /// insert all into views net
    core::ViewsNet net;
    
    std::vector<core::ViewsNet::VertHandle> viewHandles;
    for (int i = 0; i < cams.size(); i++){
        auto & camera = cams[i];
        const auto & im = ims[i];
        viewHandles.push_back(net.insertPhoto(im, camera));
    }

    auto computeFea = [](core::ViewsNet* netptr, core::ViewsNet::VertHandle vh){
        qDebug() << "photo " << vh.id;
        qDebug() << "computing features ...";
        netptr->computeFeatures(vh);
        qDebug() << "done " << vh.id;
    };

    int vid = 0;
    while(vid < net.views().internalVertices().size()){
        std::vector<std::thread> t4(std::min(net.views().internalVertices().size() - vid, 4ull));
        for (auto & t : t4)
            t = std::thread(computeFea, &net, core::ViewsNet::VertHandle(vid++));
        for (auto & t : t4)
            t.join();
    }

    for (auto & vh : viewHandles){
        qDebug() << "photo " << vh.id;
        //net.computeFeatures(vh);
        net.updateConnections(vh);
        net.computeTransformationOnConnections(vh);
        net.calibrateCamera(vh);
    }

    {
        qDebug() << "estimating vanishing points ...";
        // estimate vanishing points and classify lines
        net.estimateVanishingPointsAndClassifyLines();
        auto vps = net.globalData().vanishingPoints;
        for (auto & vp : vps)
            vp /= core::norm(vp);
        double ortho = core::norm(core::Vec3(vps[0].dot(vps[1]), vps[1].dot(vps[2]), vps[2].dot(vps[0])));
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
    }

    {
        vis::Visualizer3D viz;
        viz << vis::manip3d::SetCamera(core::PerspectiveCamera(700, 700, 200, core::Vec3(1, 1, 1) / 4, core::Vec3(0, 0, 0), core::Vec3(0, 0, -1)))
            << vis::manip3d::SetBackgroundColor(core::Black)
            << vis::manip3d::SetColorTableDescriptor(core::ColorTableDescriptor::RGB)
            << net.globalData().spatialLineSegments
            << vis::manip3d::AutoSetCamera
            << vis::manip3d::SetRenderMode(vis::RenderModeFlag::All)
            << vis::manip3d::Show();

        net.rectifySpatialLines();

        vis::Visualizer3D viz2;
        viz2 << vis::manip3d::SetCamera(core::PerspectiveCamera(700, 700, 200, core::Vec3(1, 1, 1) / 4, core::Vec3(0, 0, 0), core::Vec3(0, 0, -1)))
            << vis::manip3d::SetBackgroundColor(core::Black)
            << vis::manip3d::SetColorTableDescriptor(core::ColorTableDescriptor::RGB)
            << net.globalData().mergedSpatialLineSegments
            << vis::manip3d::AutoSetCamera
            << vis::manip3d::SetRenderMode(vis::RenderModeFlag::All)
            << vis::manip3d::Show();
    }

    
}

int main(int argc, char * argv[], char * envp[])
{
    srand(clock());
    testing::InitGoogleTest(&argc, argv);
    vis::InitGui(argc, argv);
    return RUN_ALL_TESTS();
    //return app->exec();
}
