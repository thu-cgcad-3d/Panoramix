#include "../src/core/feature.hpp"
#include "../src/core/utilities.hpp"
#include "../src/core/algorithms.hpp"
#include "../src/ml/data_set.hpp"
#include "../src/misc/matlab.hpp"
#include "../src/experimental/mixed_graph.hpp"

#include "config.hpp"

using namespace panoramix;
using namespace test;
using namespace experimental;


std::string nyu2dir = "F:\\DataSets\\NYU2\\tinyworkspace\\";

TEST(MixedGraph, Basic){

    static_assert(core::IsIteratorOfType<int*, int>::value, "");

    // load test data
    misc::Matlab::CDAndAddAllSubfolders(nyu2dir);
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
