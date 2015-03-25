#include "class_handle.hpp"
#include "common.hpp"

void mexFunction(int nlhs, mxArray* plhs[],
    const int nrhs, const mxArray* prhs[]){

    using namespace mex;

    if (nrhs != 3 || nlhs != 1){
        mexErrMsgTxt("Usage: \n"
            "   seggraphhandle = mexGetBoundaryFeas(seggraphhandle, image, gc)");
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

    // get image
    // get gc
    Image image;
    ImageOfType<Vec<double, 7>> gc;
    Matlab::GetVariable(prhs[1], image);
    if (image.channels() >= 3)
        cv::mixChannels(image, image, { 0, 2, 1, 1, 2, 0 });
    Matlab::GetVariable(prhs[2], gc);

    
    plhs[0] = convertPtr2Mat(seggraphhandle->ptr());

}