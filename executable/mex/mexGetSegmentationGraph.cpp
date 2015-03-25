#include "class_handle.hpp"
#include "common.hpp"

void mexFunction(int nlhs, mxArray* plhs[],
    const int nrhs, const mxArray* prhs[]){

    using namespace mex;
   
    if (nrhs != 1 || nlhs != 1){
        mexErrMsgTxt("Usage: \n"
            "   seggraphhandle = mexGetSegmentationGraph(image)\n");
        return;
    }

    Image image;
    misc::Matlab::GetVariable(prhs[0], image);
    if (image.channels() >= 3)
        cv::mixChannels(image, image, { 0, 2, 1, 1, 2, 0 });
   
    int segsNum = 0;
    Imagei segs;
    SegmentationExtractor segmenter;
    segmenter.params().algorithm = SegmentationExtractor::GraphCut;
    std::tie(segs, segsNum) = segmenter(image);

    std::map<std::pair<int, int>, std::vector<std::vector<PixelLoc>>> boundaryEdges;
    FindContoursOfRegionsAndBoundaries(segs, segsNum, boundaryEdges, 3);

    Graph * g = new Graph;
    auto & graph = *g;
    graph.internalElements<0>().reserve(segsNum);
    graph.internalElements<1>().reserve(boundaryEdges.size());

    for (int i = 0; i < segsNum; i++){
        Region rd;
        // todo
        graph.add(std::move(rd));
    }

    for (auto & bep : boundaryEdges){
        auto & rids = bep.first;
        Boundary bd;
        // todo
        bd.pixels = std::move(bep.second);
        graph.add<1>({ RegionHandle(rids.first), RegionHandle(rids.second) }, std::move(bd));
    }

    plhs[0] = convertPtr2Mat(g);

}