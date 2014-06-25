#ifndef PANORAMIX_REC_RECONSTRUCTION_ENGINE_VISUALIZE_HPP
#define PANORAMIX_REC_RECONSTRUCTION_ENGINE_VISUALIZE_HPP

#include "../vis/visualize2d.hpp"
#include "../vis/visualize3d.hpp"
#include "reconstruction_engine.hpp"
 
namespace panoramix {
    namespace rec {

        using vis::Visualizer2D;
        using vis::Visualizer3D;

        // 2d vis for vertdata
        Visualizer2D operator << (Visualizer2D viz, const ReconstructionEngine::ViewData & vd);


        // 3d vis for globaldata
        Visualizer3D operator << (Visualizer3D viz, const ReconstructionEngine::GlobalData & netgb);

        // 2d for whole net
        Visualizer2D operator << (Visualizer2D viz, const ReconstructionEngine & net);

    }
}
 
#endif