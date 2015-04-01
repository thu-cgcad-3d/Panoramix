#include "../class_handle.hpp"

#include "../../src/misc/matlab.hpp"
#include "../../src/core/algorithms.hpp"
#include "../../src/core/utilities.hpp"
#include "../../src/core/mixed_graph.hpp"
#include "../../src/gui/visualize2d.hpp"
#include "../../src/gui/visualizers.hpp"

using namespace panoramix;
using namespace panoramix::core;


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

    if (cmd == "showEdges"){
        
        return;
    }

    // Got here, so command not recognized
    mexErrMsgTxt("Command not recognized.");
}
