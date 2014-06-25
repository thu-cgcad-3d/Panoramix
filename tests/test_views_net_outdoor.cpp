#include "../src/core/mesh_maker.hpp"
#include "../src/core/utilities.hpp"
#include "../src/rec/reconstruction_engine.hpp"
#include "../src/rec/reconstruction_engine_visualize.hpp"
#include "../src/rec/regions_net_visualize.hpp"
#include "gtest/gtest.h"

#include <iostream>
#include <string>
#include <random>
#include <thread>

#include <QApplication>

using namespace panoramix;

// PROJECT_TEST_DATA_DIR_STR is predefined using CMake
static const std::string ProjectTestDataDirStr = PROJECT_TEST_DATA_DIR_STR;
static const std::string ProjectTestDataDirStr_Normal = ProjectTestDataDirStr + "/normal";
static const std::string ProjectTestDataDirStr_PanoramaIndoor = ProjectTestDataDirStr + "/panorama/indoor";
static const std::string ProjectTestDataDirStr_PanoramaOutdoor = ProjectTestDataDirStr + "/panorama/outdoor";

TEST(ViewsNet, ViewsNet) {
    cv::Mat panorama = cv::imread(ProjectTestDataDirStr_PanoramaOutdoor + "/panohk.png");
    cv::resize(panorama, panorama, cv::Size(2000, 1000));

    //vis::Visualizer2D(panorama) << vis::manip2d::Show();
    core::PanoramicCamera originCam(panorama.cols / M_PI / 2.0);

    std::vector<core::PerspectiveCamera> cams;
    core::Mesh<core::Vec3> cameraStand;
    core::MakeQuadFacedSphere(cameraStand, 6, 12);
    for (auto & v : cameraStand.vertices()){
        core::Vec3 direction = v.data;
        if (core::AngleBetweenDirections(direction, core::Vec3(0, 0, 1)) <= 0.1 ||
            core::AngleBetweenDirections(direction, core::Vec3(0, 0, -1)) <= 0.1){
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
        std::cout << "sampling photo ..." << std::endl;
        return core::CameraSampler<core::PerspectiveCamera, core::PanoramicCamera>(pcam, originCam)(panorama);
    });

    /// insert all into views net
    rec::ReconstructionEngine::Params params;
    params.mjWeightT = 2.0;
    params.intersectionConstraintLineDistanceAngleThreshold = 0.05;
    params.incidenceConstraintLineDistanceAngleThreshold = 0.2;
    params.mergeLineDistanceAngleThreshold = 0.05;
    rec::ReconstructionEngine net(params);

    std::vector<rec::ReconstructionEngine::ViewHandle> viewHandles;
    for (int i = 0; i < cams.size(); i++){
        auto & camera = cams[i];
        const auto & im = ims[i];
        viewHandles.push_back(net.insertPhoto(im, camera));
    }

    auto computeFea = [](rec::ReconstructionEngine* netptr, rec::ReconstructionEngine::ViewHandle vh){
        std::cout << "photo " << vh.id << std::endl;
        std::cout << "computing features ..." << std::endl;
        netptr->computeFeatures(vh);
        std::cout << "done " << vh.id << std::endl;
    };

    int vid = 0;
    while (vid < net.views().internalElements<0>().size()){
        std::vector<std::thread> t4(std::min(net.views().internalElements<0>().size() - vid, 4ull));
        for (auto & t : t4)
            t = std::thread(computeFea, &net, rec::ReconstructionEngine::ViewHandle(vid++));
        for (auto & t : t4)
            t.join();
    }

    for (auto & vh : viewHandles){
        std::cout << "photo " << vh.id << std::endl;
        net.updateConnections(vh);
        //net.calibrateCamera(vh);
    }

    {
        std::cout << "estimating vanishing points ..." << std::endl;
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
            << vis::manip3d::SetBackgroundColor(vis::ColorTag::Black)
            << vis::manip3d::SetColorTableDescriptor(vis::ColorTableDescriptor::RGB)
            << net.globalData().spatialLineSegments
            << vis::manip3d::AutoSetCamera
            << vis::manip3d::SetRenderMode(vis::RenderModeFlag::All)
            << vis::manip3d::Show();

        net.rectifySpatialLines();

        vis::Visualizer3D viz2;
        viz2 << vis::manip3d::SetCamera(core::PerspectiveCamera(700, 700, 200, core::Vec3(1, 1, 1) / 4, core::Vec3(0, 0, 0), core::Vec3(0, 0, -1)))
            << vis::manip3d::SetBackgroundColor(vis::ColorTag::Black)
            << vis::manip3d::SetColorTableDescriptor(vis::ColorTableDescriptor::RGB)
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
    return RUN_ALL_TESTS();
}
