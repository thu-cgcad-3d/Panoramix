#include "pi_graph.hpp"

namespace pano {
    namespace experimental {

        void DetectOcclusions(PIGraph & mg);

        // just set all bndpiece and linepice connected/attached
        void AssumeThereAreNoOcclusions(PIGraph & mg);

    }
}