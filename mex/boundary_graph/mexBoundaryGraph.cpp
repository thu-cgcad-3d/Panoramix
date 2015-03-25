#include "../class_handle.hpp"

#include "boundary_graph.hpp"

using namespace mex;

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{	
    // Get the command string
    char cmd[64];
	if (nrhs < 1 || mxGetString(prhs[0], cmd, sizeof(cmd)))
		mexErrMsgTxt("First input should be a command string less than 64 characters long.");
        
    // New
    if (!strcmp("new", cmd)) {
        // Check parameters
        if (nlhs != 1)
            mexErrMsgTxt("New: One output expected.");
        // Return a handle to a new C++ instance
        if (nrhs == 1){
            plhs[0] = convertPtr2Mat<mex::BoundaryGraph>(new mex::BoundaryGraph);
        }
        else if (nrhs >= 2){
            Image image;
            misc::Matlab::GetVariable(prhs[1], image);
            if (image.channels() >= 3)
                cv::mixChannels(image, image, { 0, 2, 1, 1, 2, 0 });
            cv::imshow("image", image);
            cv::waitKey();
            Imaged7 gc;
            if (nrhs >= 3){
                misc::Matlab::GetVariable(prhs[2], gc);
            }
            Imaged depth;
            if (nrhs >= 4){
                misc::Matlab::GetVariable(prhs[3], depth);
            }
            plhs[0] = convertPtr2Mat<mex::BoundaryGraph>(new mex::BoundaryGraph(image, gc, depth));
        }
        else{
            mexErrMsgTxt("New: Invalid Arguments Num");
            return;
        }
        printf("BoundaryGraph created");
        return;
    }
    
    // Check there is a second input, which should be the class instance handle
    if (nrhs < 2)
		mexErrMsgTxt("Second input should be a class instance handle.");
    
    // Delete
    if (!strcmp("delete", cmd)) {
        // Destroy the C++ object
        destroyObject<mex::BoundaryGraph>(prhs[1]);
        // Warn if other commands were ignored
        if (nlhs != 0 || nrhs != 2)
            mexWarnMsgTxt("Delete: Unexpected arguments ignored.");
        printf("BoundaryGraph destroyed\n");
        return;
    }
    
    // Get the class instance pointer from the second input
    mex::BoundaryGraph *instance = convertMat2Ptr<mex::BoundaryGraph>(prhs[1]);

    // Call the various class methods
    // Train    
    if (!strcmp("print", cmd)) {
        // Check parameters
        if (nlhs < 0 || nrhs < 2)
            mexErrMsgTxt("Print: Unexpected arguments.");
        // Call the method
        instance->print();
        return;
    }




    // Got here, so command not recognized
    mexErrMsgTxt("Command not recognized.");
}
