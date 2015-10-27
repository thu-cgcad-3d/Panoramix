#pragma once

#include "pi_graph.hpp"
#include "pi_graph_occlusion.hpp"

namespace pano {
    namespace experimental {

        class PILayoutAnnotation;

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
            bool connectLeft() const { return leftWeightRatio > 0; }
            bool connectRight() const { return rightWeightRatio > 0; }
            bool onlyConnectLeft() const { return leftWeightRatio > 0 && rightWeightRatio == 0; }
            bool onlyConnectRight() const { return leftWeightRatio == 0 && rightWeightRatio > 0; }
            double minWeightRatio() const { return std::min(leftWeightRatio, rightWeightRatio); }
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
            std::vector<std::map<int, double>> * line2leftSegsWithWeightPtr = nullptr,
            std::vector<std::map<int, double>> * line2rightSegsWithWeightPtr = nullptr);

        std::vector<LineSidingWeight> ComputeLinesSidingWeightsFromAnnotation(const PIGraph & mg,
            const PILayoutAnnotation & anno, 
            double sampleAngleStep = DegreesToRadians(0.5),
            double angleThres = DegreesToRadians(2),
            double ratioThres = 0.6);

        std::vector<std::array<std::set<int>, 2>> CollectSegsNearLines(const PIGraph & mg, 
            double angleSizeForPixelsNearLines = DegreesToRadians(2));

        void ApplyLinesSidingWeights(PIGraph & mg, const std::vector<LineSidingWeight> & lsw, 
            const std::vector<std::array<std::set<int>, 2>> & line2leftRightSegs,
            bool connectSegsOnDanglingLine);



    }
}