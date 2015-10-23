#pragma once

#include "pi_graph.hpp"

namespace pano {
    namespace experimental {



        void DetectOcclusions(PIGraph & mg,
            double minAngleSizeOfLineInTJunction = DegreesToRadians(3),
            double lambdaShrinkForHLineDetectionInTJunction = 0.2,
            double lambdaShrinkForVLineDetectionInTJunction = 0.1,
            double angleSizeForPixelsNearLines = DegreesToRadians(5));


        struct LineSidingWeight {
            double leftWeightRatio;
            double rightWeightRatio;
            bool isOcclusion() const { return leftWeightRatio == 0 || rightWeightRatio == 0; }
            bool isLineDetached() const { return leftWeightRatio == 0 && rightWeightRatio == 0; }
            bool onlyConnectLeft() const { return leftWeightRatio > 0 && rightWeightRatio == 0; }
            bool onlyConnectRight() const { return leftWeightRatio == 0 && rightWeightRatio > 0; }
            template <class Archiver>
            void serialize(Archiver & ar) {
                ar(leftWeightRatio, rightWeightRatio);
            }
        };

        std::vector<LineSidingWeight> ComputeLinesSidingWeights(const PIGraph & mg,
            double minAngleSizeOfLineInTJunction = DegreesToRadians(3),
            double lambdaShrinkForHLineDetectionInTJunction = 0.2,
            double lambdaShrinkForVLineDetectionInTJunction = 0.1,
            double angleSizeForPixelsNearLines = DegreesToRadians(2),
            std::vector<std::array<std::set<int>, 2>> * line2leftRightSegsPtr = nullptr);

        void ApplyLinesSidingWeights(PIGraph & mg, const std::vector<LineSidingWeight> & lsw, 
            const std::vector<std::array<std::set<int>, 2>> & line2leftRightSegs);



    }
}