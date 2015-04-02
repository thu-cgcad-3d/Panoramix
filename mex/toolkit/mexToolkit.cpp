#include "../class_handle.hpp"

#include "../../src/misc/matlab.hpp"
#include "../../src/core/algorithms.hpp"
#include "../../src/core/utilities.hpp"
#include "../../src/gui/visualize2d.hpp"
#include "../../src/gui/visualizers.hpp"

using namespace panoramix;
using namespace panoramix::core;


std::string GetString(const mxArray * mxa){
    size_t len = mxGetN(mxa);
    char * buffer = new char[len + 1];
    std::fill(buffer, buffer + len + 1, '\0');
    mxGetString(mxa, buffer, len + 1);
    std::string str(buffer);
    delete[] buffer;
    return str;
}


void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]) {

    // Get the command string
    char cmdbuffer[64];
    if (nrhs < 1 || mxGetString(prhs[0], cmdbuffer, sizeof(cmdbuffer)))
		mexErrMsgTxt("First input should be a command string less than 64 characters long.");
    std::string cmd = cmdbuffer;

    const mxArray ** argv = prhs + 1;
    int argc = nrhs - 1;
    int outc = nlhs;
    mxArray ** outv = plhs;

    if (cmd == "saveMatToPanoramix"){
        if (argc != 3){
            mexErrMsgTxt("Wrong inputs num!");
            return;
        }  
        std::string filename = GetString(argv[0]);
        double lastDimAsChannel = mxGetScalar(argv[2]);
        core::Image mat;
        misc::Matlab::GetVariable(argv[1], mat, lastDimAsChannel != 0);
        core::SaveToDisk(filename, mat);
        return;
    }

    if (cmd == "loadMatFromPanoramix"){
        if (argc != 2 || outc != 1){
            mexErrMsgTxt("Wrong inputs/outputs num!");
            return;
        }
        std::string filename = GetString(argv[0]);
        core::Image mat;
        core::LoadFromDisk(filename, mat);
        auto mxa = static_cast<mxArray*>(misc::Matlab::PutVariable(mat));
        outv[0] = mxa;
        return;
    }

    // Got here, so command not recognized
    mexErrMsgTxt("Command not recognized.");
}
