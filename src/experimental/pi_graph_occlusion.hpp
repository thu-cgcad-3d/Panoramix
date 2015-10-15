#pragma once

#include "pi_graph.hpp"

namespace pano {
    namespace experimental {

        void DetectOcclusions(PIGraph & mg);



        struct LineLabel {
            bool connectLeft, connectRight;
            bool operator == (LineLabel ll) const { return connectLeft == ll.connectLeft && connectRight == ll.connectRight; }
            bool operator != (LineLabel ll) const { return !(*this == ll); }
            LineLabel operator !() const { return LineLabel{ !connectLeft, !connectRight }; }
            LineLabel leftRightSwapped() const { return LineLabel{ connectRight, connectLeft }; }
        };
        struct LineLabelCost {
            double connectLeftConnectRight;
            double connectLeftDisconnectRight;
            double disconnectLeftConnectRight;
            double disconnectAll;
        };

        void DetectOcclusions2(PIGraph & mg, 
            double minAngleSizeOfLineInTJunction = DegreesToRadians(3),
            double lambdaShrinkForHLineDetectionInTJunction = 0.2,
            double lambdaShrinkForVLineDetectionInTJunction = 0.1,
            double angleSizeForPixelsNearLines = DegreesToRadians(2));

    }
}