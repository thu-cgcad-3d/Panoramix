#include "../class_handle.hpp"

#include "../../src/misc/mxarray.hpp"

#include "../../src/core/algorithms.hpp"
#include "../../src/core/utility.hpp"
#include "../../src/core/cons_graph.hpp"
#include "../../src/ml/factor_graph.hpp"

#include "../../src/gui/visualize2d.hpp"
#include "../../src/gui/visualizers.hpp"

using namespace panoramix;
using namespace panoramix::core;


template <class T>
std::valarray<T> GenRange(T first, T step, T last){
    std::valarray<T> r(std::ceil(last - first / step));
    for (int i = 0; i < r.size(); i++){
        r[i] = first + i * step;
    }
    return r;
}

template <class T1, class T2>
std::pair<std::valarray<T1>, std::valarray<T2>> MeshGrid(const std::valarray<T1> & xs, const std::valarray<T2> & ys){
    std::pair<std::valarray<T1>, std::valarray<T2>> results;
    results.first = std::valarray<T1>(xs.size() * ys.size());
    results.second = std::valarray<T2>(xs.size() * ys.size());
    int o = 0;
    for (int i = 0; i < xs.size(); i++){
        for (int j = 0; j < ys.size(); j++){
            results.first[o] = xs[o];
            results.second[o] = ys[o];
            o++;
        }
    }
    return results;
}


struct Patch {
    Point2 center;
    double scale;
    void feature(const Image3f & im, float * fea) const {

    }
};


struct PatchTable {

};


/// perspective pattern finder
void FindPerspectivePattern(const Image3f & oim){

    Image3f im = oim.clone();
    ResizeToMakeWidthUnder(im, 500);
    PerspectiveCamera cam(im.cols, im.rows, Point2(im.cols / 2, im.rows / 2), std::min(im.cols, im.rows));

    struct Patch {
        Point2 center;
        double radius; // radius in image
        Plane3 plane; // to be optimized
    };

    struct PatchCon {
        double similarity;
        Mat3 transform;
    };

    using PatchGraph = HomogeneousGraph02<Patch, PatchCon>;

    // for each patch, find neighbor patches that are similiar to this after certain transform
    // construct a patch graph
    auto radiusTable = GenRange(2.0, 0.5, 10.0);
    auto radiusRatioTable = std::log(GenRange(-1.0, 0.1 , 1.0));
    auto posTable = MeshGrid(GenRange(-8.0, 0.1, 8.0), GenRange(-8.0, 0.1, 8.0));
    auto angleTable = GenRange(-M_PI_2, 0.1, +M_PI_2);

    for (auto it = im.begin(); it != im.end(); ++it){
        Point2 pos = ecast<double>(it.pos());
        for (double r : radiusTable){
            Patch p1;
            p1.center = pos;
            p1.radius = r;

            

            for (double rr : radiusRatioTable){
                double r2 = rr * r;
                for (int i = 0; i < posTable.first.size(); i++){
                    Patch p2;
                    p2.center = p1.center + Point2(posTable.first[i], posTable.second[i]);
                    p2.radius = r2;

                    // compare p1 and p2;
                    
                }
            }
        }
    }

}




void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]) {

    auto ins = misc::ToMXArray(nrhs, prhs);

    // Get the command string
    if (ins.empty() || ins.front().toString().empty())
		mexErrMsgTxt("First input should be a command string.");
    std::string cmd = ins.front().toString();

    int argc = ins.size() - 1;
    int outc = nlhs;
    auto argv = &(ins[1]);
    auto outv = plhs;

    if (cmd == "saveMatToPanoramix"){
        if (argc != 3){
            mexErrMsgTxt("Wrong inputs num!");
            return;
        }  

        std::string filename = argv[0].toString();
        double lastDimAsChannel = argv[2].scalar();
        core::Image mat;
        misc::FromMXArray(argv[1], mat, lastDimAsChannel != 0);
        core::SaveToDisk(filename, mat);
        return;
    }

    if (cmd == "loadMatFromPanoramix"){
        if (argc != 2 || outc != 1){
            mexErrMsgTxt("Wrong inputs/outputs num!");
            return;
        }
        std::string filename = argv[0].toString();
        core::Image mat;
        core::LoadFromDisk(filename, mat);
        outv[0] = misc::ToMXArray(mat);
        return;
    }

    if (cmd == "segmentGraphCut"){
        if (argc == 0){
            mexErrMsgTxt("Wrong inputs/outputs num!");
            return;
        }
        core::SegmentationExtractor segmenter;
        segmenter.params().algorithm = core::SegmentationExtractor::GraphCut;
        if (argc > 1)
            segmenter.params().sigma = argv[1];
        if (argc > 2)
            segmenter.params().c = argv[2];
        core::Image image;
        misc::FromMXArray(argv[0], image);
        core::Imagei segs;
        int segnum;
        std::tie(segs, segnum) = segmenter(image, false);
        if (outc > 0){
            outv[0] = misc::ToMXArray(segs);
        }
        if (outc > 1){
            outv[1] = mxCreateDoubleScalar(segnum);
        }
        return;
    }

    if (cmd == "segmentGraphCutPano"){
        if (argc == 0){
            mexErrMsgTxt("Wrong inputs/outputs num!");
            return;
        }
        core::SegmentationExtractor segmenter;
        segmenter.params().algorithm = core::SegmentationExtractor::GraphCut;
        if (argc > 1)
            segmenter.params().sigma = mxGetScalar(argv[1]);
        if (argc > 2)
            segmenter.params().c = mxGetScalar(argv[2]);
        core::Image image;
        misc::FromMXArray(argv[0], image);
        core::Imagei segs;
        int segnum;
        std::tie(segs, segnum) = segmenter(image, true);
        if (outc > 0){
            outv[0] = misc::ToMXArray(segs);
        }
        if (outc > 1){
            outv[1] = mxCreateDoubleScalar(segnum);
        }
        return;
    }

    if (cmd == "segmentSLIC"){
        if (argc == 0){
            mexErrMsgTxt("Wrong inputs/outputs num!");
            return;
        }
        core::SegmentationExtractor segmenter;
        segmenter.params().algorithm = core::SegmentationExtractor::SLIC;
        if (argc > 1)
            segmenter.params().superpixelSizeSuggestion = mxGetScalar(argv[1]);
        if (argc > 2)
            segmenter.params().superpixelNumberSuggestion = mxGetScalar(argv[2]);
        core::Image image;
        misc::FromMXArray(argv[0], image);
        core::Imagei segs;
        int segnum;
        std::tie(segs, segnum) = segmenter(image);
        if (outc > 0){
            outv[0] = misc::ToMXArray(segs);
        }
        if (outc > 1){
            outv[1] = mxCreateDoubleScalar(segnum);
        }
        return;
    }

    if (cmd == "SIFT"){
        if (argc == 0){
            mexErrMsgTxt("Wrong inputs/outputs num!");
            return;
        }
        core::Image image;
        misc::FromMXArray(argv[0], image);
        cv::SIFT sift;
        std::vector<cv::KeyPoint> kps;
        core::DenseMat<float> descriptors;
        sift(image, cv::noArray(), kps, descriptors);
        if (outc > 0){
            outv[0] = misc::ToMXArray(kps);
        }
        if (outc > 1){
            outv[1] = misc::ToMXArray(descriptors);
        }
        return;
    }

    if (cmd == "SURF"){
        if (argc == 0){
            mexErrMsgTxt("Wrong inputs/outputs num!");
            return;
        }
        core::Image image;
        misc::FromMXArray(argv[0], image);
        cv::SURF surf;
        std::vector<cv::KeyPoint> kps;
        core::DenseMat<float> descriptors;
        surf(image, cv::noArray(), kps, descriptors);
        if (outc > 0){
            outv[0] = misc::ToMXArray(kps);
        }
        if (outc > 1){
            outv[1] = misc::ToMXArray(descriptors);
        }
        return;
    }
    

    if (cmd == "estimatePerspectivePattern"){
        if (argc == 0){
            mexErrMsgTxt("Wrong inputs/outputs num!");
            return;
        }

        core::Image image;
        misc::FromMXArray(argv[0], image);

        if (outc > 0){ // depth
            //outv[0]
        }
        else{ // visualize

        }

    }
    


    // Got here, so command not recognized
    mexErrMsgTxt("Command not recognized.");
}
