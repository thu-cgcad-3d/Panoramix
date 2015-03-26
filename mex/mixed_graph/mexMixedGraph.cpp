#include "../class_handle.hpp"

#include "../../src/misc/matlab.hpp"
#include "../../src/core/mixed_graph.hpp"

using namespace panoramix;
using namespace panoramix::core;

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{	
    // Get the command string
    char cmdbuffer[64];
    if (nrhs < 1 || mxGetString(prhs[0], cmdbuffer, sizeof(cmdbuffer)))
		mexErrMsgTxt("First input should be a command string less than 64 characters long.");
    std::string cmd = cmdbuffer;

    const mxArray ** argv = prhs + 1;
    int argc = nrhs - 1;

    // New
    if (cmd == "new") {
        if (argc != 1){
            mexErrMsgTxt("Wrong Arguments Num.");
            return;
        }
        auto imageMXA = argv[0];
        Image image;
        misc::Matlab::GetVariable(imageMXA, image);
        if (image.channels() >= 3)
            cv::mixChannels(image, image, { 0, 2, 1, 1, 2, 0 });
        
        cv::imshow("image", image);
        cv::waitKey();

        plhs[0] = convertPtr2Mat<MixedGraph>(new MixedGraph);        
        return;
    }
    
    // Delete
    if (cmd == "delete") {
        destroyObject<MixedGraph>(prhs[1]);
        printf("BoundaryGraph destroyed\n");
        return;
    }
    
    // Get the class instance pointer from the second input
    MixedGraph *instance = convertMat2Ptr<MixedGraph>(prhs[1]);



    // Got here, so command not recognized
    mexErrMsgTxt("Command not recognized.");
}
