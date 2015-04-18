#include "../src/core/feature.hpp"
#include "../src/core/utilities.hpp"
#include "../src/core/algorithms.hpp"
#include "../src/ml/data_set.hpp"
#include "../src/misc/matlab.hpp"
#include "../src/experimental/engine.hpp"
#include "../src/gui/visualizers.hpp"

#include "config.hpp"

using namespace panoramix;
using namespace test;
using namespace experimental;

#define LOG(...) std::cout << "[Engine Optimization] #######  " << __VA_ARGS__ << " #######"<< std::endl


static_assert(IsCamera<PerspectiveCamera>::value, "");


std::string nyu2dir = "F:\\DataSets\\NYU2\\";

struct CameraParams {
    Vec2 f_rgb, f_d;
    Vec2 c_rgb, c_d;
    Vec3 k_rgb, k_d;
    Vec2 p_rgb, p_d;
    Vec3 t;
    Mat3 R;
    double depthParam1, depthParam2;
    double maxDepth;
};



TEST(Engine, Basic){

    // load test data
    misc::Matlab::CDAndAddAllSubfolders(nyu2dir + "tinyworkspace\\");
    misc::Matlab matlab;
    matlab << "load occ edges occscore image gc depth;";
    
    // get edges
    double edgenum = 0;
    matlab << "edgenum = length(edges.indices);";
    matlab.GetVariable("edgenum", edgenum);
    assert(edgenum > 0);

    std::vector<std::vector<int32_t>> indices(edgenum);
    for (int i = 0; i < edgenum; i++){
        matlab << ("indices = edges.indices{" + std::to_string(i + 1) + "}';");
        matlab.GetVariable("indices", indices[i]);
    }

    // get occ scores
    std::vector<double> occscore;
    matlab << "occscore = occscore';";
    matlab.GetVariable("occscore", occscore);
    assert(occscore.size() == (int)edgenum);

    // get image
    Image image;
    matlab.GetVariable("image", image, true);

    // get gc
    Image7d gc;
    matlab.GetVariable("gc", gc, true);

    // get depth
    Imagef depth;
    matlab.GetVariable("depth", depth, false);

    auto c = ml::annotations::nyu2::c_rgb();
    auto f = ml::annotations::nyu2::f_rgb();
    
    Engine mg(image, c, (f(0) + f(1)) / 2.0);
    mg.installOcclusionResponce(indices, occscore);
    mg.installGCResponse(gc);

    //mg.showBoundaryJunctions();
    //mg.showDetachableRegionLineConnections();
    //mg.showOcclusionResponse();

    mg.solve();
    
}


TEST(Engine, Batch){

    ForEachCase([](int id, const std::vector<std::vector<int32_t>> & indices, const std::vector<double> & occscore,
        const Image & image, const Image7d & gc, const Imagef & depth, const CameraParams & cameraParams){
        
        //if (id < 4)
        //    return;
       
        core::GeneralPerspectiveCamera cam_rgb(640, 480, cameraParams.c_rgb, cameraParams.f_rgb);
        core::GeneralPerspectiveCamera cam_d(640, 480, cameraParams.c_d, cameraParams.f_d);


        Engine mg(image, cameraParams.c_rgb, Mean(cameraParams.f_rgb(0), cameraParams.f_rgb(1)));

        std::vector<Point3> gtPoints;
        gtPoints.reserve(depth.cols * depth.rows);
        for (auto it = depth.begin(); it != depth.end(); ++it){
            auto d = mg.view.camera.spatialDirection(it.pos()) * depth(it.pos()) / 10000.0; 
            gtPoints.emplace_back(d);
        }
        gui::Visualizer vis;
        vis.installingOptions.pointSize = 1.0;
        vis.renderOptions.renderMode |= gui::RenderModeFlag::Points;
        vis.renderOptions.backgroundColor = gui::White;
        vis.camera(core::PerspectiveCamera(1000, 800, core::Point2(500, 400), 800, Point3(-1, 1, 1), Point3(0, 0, 0), Point3(0, 0, -1)));
        vis.add(gtPoints).show(true, true);
        
        mg.installOcclusionResponce(indices, occscore);
        mg.installGCResponse(gc);
        mg.depthCandidates.push_back(depth);

        mg.showSegmentations();

        mg.solve();


    });

}