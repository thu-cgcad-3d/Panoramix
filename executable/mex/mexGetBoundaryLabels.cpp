#include "class_handle.hpp"
#include "common.hpp"


void mexFunction(int nlhs, mxArray* plhs[],
    const int nrhs, const mxArray* prhs[]){
   
    using namespace mex;

    if (nrhs != 2 || nlhs != 1){
        mexErrMsgTxt("Usage: \n"
            "   seggraphhandle = mexGetBoundaryLabels(seggraphhandle, depth)");
        return;
    }

    // get graph
    auto seggraphhandle = convertMat2HandlePtr<Graph>(prhs[0]);
    if (!seggraphhandle->isValid()){
        mexErrMsgTxt("Invalid seggraphhandle!");
        return;
    }

    Graph & graph = **seggraphhandle;
    printf("graph has %d regions and %d boundaries\n", 
        graph.internalElements<0>().size(), graph.internalElements<1>().size());

    Imaged depth;
    Matlab::GetVariable(prhs[1], depth);

    depth /= MinMaxValOfImage(depth).second;
    cv::imshow("depth", depth);
    cv::waitKey();



}