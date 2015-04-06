#include "../src/core/feature.hpp"
#include "../src/core/utilities.hpp"
#include "../src/core/algorithms.hpp"
#include "../src/ml/data_set.hpp"
#include "../src/misc/matlab.hpp"
#include "../src/experimental/mixed_graph.hpp"
#include "../src/gui/visualizers.hpp"

#include "config.hpp"

using namespace panoramix;
using namespace test;
using namespace experimental;

#define LOG(...) std::cout << "[MixedGraph Optimization] #######  " << __VA_ARGS__ << " #######"<< std::endl


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

template <class FunT>
void ForEachCase(FunT && fun){

    // load test data
    misc::Matlab::CDAndAddAllSubfolders(nyu2dir);
    misc::Matlab matlab;

    // get camera params
    matlab << "camera_params;";
    matlab
        << "f_rgb = [fx_rgb;fy_rgb];"
        << "f_d = [fx_d;fy_d];"
        << "c_rgb = [cx_rgb;cy_rgb];"
        << "c_d = [cx_d;cy_d];"
        << "k_rgb = [k1_rgb;k2_rgb;k3_rgb];"
        << "k_d = [k1_d;k2_d;k3_d];"
        << "p_rgb = [p1_rgb;p2_rgb];"
        << "p_d = [p1_d;p2_d];"
        << "t = [t_x;t_y;t_z];";
    
    CameraParams cameraParams;
    matlab.GetVariable("f_rgb", cameraParams.f_rgb);
    matlab.GetVariable("f_d", cameraParams.f_d);
    matlab.GetVariable("c_rgb", cameraParams.c_rgb);
    matlab.GetVariable("c_d", cameraParams.c_d);
    matlab.GetVariable("k_rgb", cameraParams.k_rgb);
    matlab.GetVariable("k_d", cameraParams.k_d);
    matlab.GetVariable("p_rgb", cameraParams.p_rgb);
    matlab.GetVariable("p_d", cameraParams.p_d);
    matlab.GetVariable("t", cameraParams.t);
    matlab.GetVariable("R", cameraParams.R);
    matlab.GetVariable("depthParam1", cameraParams.depthParam1);
    matlab.GetVariable("depthParam2", cameraParams.depthParam2);
    matlab.GetVariable("maxDepth", cameraParams.maxDepth);

    int offsetper30 = 0;
    int offsetper100 = 0;

    for (int id = 1; id <= 1449; id++){

        std::cout << "ID: " << id << std::endl;

        if (id % 30 == 1){
            int first = id;
            int last = id + 30 - 1;
            std::cout << ("load occscores_[" + std::to_string(first) + "-" + std::to_string(last) + "] bndinfos occscores;") << std::endl;
            matlab << ("load occscores_[" + std::to_string(first) + "-" + std::to_string(last) + "] bndinfos occscores;");
            offsetper30 = id - 1;
        }
        if (id % 100 == 1){
            int first = id;
            int last = id + 100 - 1;
            std::cout << ("load geometric_contexts_[" + std::to_string(first) + "-" + std::to_string(last) + "] gc;") << std::endl;
            std::cout << ("load imdps_[" + std::to_string(first) + "-" + std::to_string(last) + "] ims dps;") << std::endl;
            matlab << ("load geometric_contexts_[" + std::to_string(first) + "-" + std::to_string(last) + "] gc;");
            matlab << ("load imdps_[" + std::to_string(first) + "-" + std::to_string(last) + "] ims dps;");
            offsetper100 = id - 1;
        }

        int trueIdPer30 = id - offsetper30;
        int trueIdPer100 = id - offsetper100;

        // get edges
        double edgenum = 0;
        matlab << ("inds = bndinfos(" + std::to_string(trueIdPer30) + ").edges.indices;");
        matlab << "edgenum = length(inds);";
        matlab.GetVariable("edgenum", edgenum);
        assert(edgenum > 0);

        std::vector<std::vector<int32_t>> indices(edgenum);
        for (int i = 0; i < edgenum; i++){
            matlab << ("indices = inds{" + std::to_string(i + 1) + "}';");
            matlab.GetVariable("indices", indices[i]);
        }

        // get occ scores
        std::vector<double> occscore;
        matlab << ("occscore = occscores{" + std::to_string(trueIdPer30) + "}';");
        matlab.GetVariable("occscore", occscore);
        assert(occscore.size() == (int)edgenum);



        // get image
        Image image;
        matlab << ("image = crop_image(ims(:,:,:, " + std::to_string(trueIdPer100) + "));");
        matlab.GetVariable("image", image, true);

        // get gc
        Image7d gc;
        matlab << ("g = crop_image(gc(:,:,:," + std::to_string(trueIdPer100) + "));");
        matlab.GetVariable("g", gc, true);

        // get depth
        Imagef depth;
        matlab << ("depth = crop_image(dps(:,:," + std::to_string(trueIdPer100) + "));");
        matlab.GetVariable("depth", depth, false);

        assert(image.size() == gc.size());
        assert(image.size() == depth.size());

        fun(id, indices, occscore, image, gc, depth, cameraParams);

    }
}


TEST(MixedGraph, Basic){

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
    
    MixedGraph mg(image, c, (f(0) + f(1)) / 2.0);
    mg.installOcclusionResponce(indices, occscore);
    mg.installGCResponse(gc);

    //mg.showBoundaryJunctions();
    //mg.showDetachableRegionLineConnections();
    //mg.showOcclusionResponse();

    mg.solve();
    
}


TEST(MixedGraph, Batch){

    ForEachCase([](int id, const std::vector<std::vector<int32_t>> & indices, const std::vector<double> & occscore,
        const Image & image, const Image7d & gc, const Imagef & depth, const CameraParams & cameraParams){
        
        if (id < 5)
            return;
       
        core::GeneralPerspectiveCamera cam_rgb(640, 480, cameraParams.c_rgb, cameraParams.f_rgb);
        core::GeneralPerspectiveCamera cam_d(640, 480, cameraParams.c_d, cameraParams.f_d);


        MixedGraph mg(image, cameraParams.c_rgb, Mean(cameraParams.f_rgb(0), cameraParams.f_rgb(1)));

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

        mg.solve();


    });

}