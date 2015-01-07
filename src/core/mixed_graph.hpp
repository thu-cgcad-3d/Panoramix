#ifndef PANORAMIX_CORE_MIXED_GRAPH_HPP
#define PANORAMIX_CORE_MIXED_GRAPH_HPP

#include "view.hpp"

namespace panoramix {
    namespace core {
     
        // the mixed graph
        // unary variable
        struct MGUnaryVariable {
            int claz;
            double depthOfCenter;
            template <class Archive> void serialize(Archive & ar) {
                ar(claz, depthOfCenter);
            }
        };

        // unary entity
        struct MGUnary {
            enum Type { Region, Line } type;
            std::vector<Vec3> normalizedCorners;
            Vec3 normalizedCenter;
            template <class Archive> void serialize(Archive & ar) {
                ar(type, normalizedCorners, normalizedCenter);
            }
        };

        Plane3 PlaneOfMGUnary(const MGUnary & unary, const std::vector<Vec3> & vps,
            const MGUnaryVariable & var);
        Line3 LineOfMGUnary(const MGUnary & unary, const std::vector<Vec3> & vps,
            const MGUnaryVariable & var);

        double DepthRatioOnMGUnary(const Vec3 & direction, const MGUnary & unary,
            const std::vector<Vec3> & vps, int claz);
        Point3 LocationOnMGUnary(const Vec3 & direction, const MGUnary & unary,
            const std::vector<Vec3> & vps, const MGUnaryVariable & var);

        // binary variable
        struct MGBinaryVariable{
            std::array<std::vector<double>, 2> sampleDepthsOnRelatedUnaries;
            template <class Archive> void serialize(Archive & ar) {
                ar(sampleDepthsOnRelatedUnaries);
            }
        };

        // binary entity
        struct MGBinary {
            enum Type {
                RegionRegionConnection,
                RegionRegionOverlapping,
                RegionLineConnection,
                LineLineIntersection,
                LineLineIncidence
            } type;
            double weight;
            std::vector<Vec3> normalizedAnchors;
            std::array<double, 2> importanceRatioInRelatedUnaries;
            template <class Archive> void serialize(Archive & ar) {
                ar(type, weight, normalizedAnchors, importanceRatioInRelatedUnaries);
            }
        };

        using MixedGraph = HomogeneousGraph02<MGUnary, MGBinary>;
        using MGUnaryHandle = HandleAtLevel<0>;
        using MGBinaryHandle = HandleAtLevel<1>;

        using MGUnaryVarTable = std::unordered_map<MGUnaryHandle, MGUnaryVariable>;
        using MGBinaryVarTable = std::unordered_map<MGBinaryHandle, MGBinaryVariable>;

        void InitializeUnaryVarDepths(MGUnaryVarTable & unaryVars, double depth);
        void UpdateBinaryVars(const MixedGraph & mg, const std::vector<Vec3> & vps,
            const MGUnaryVarTable & unaryVars, MGBinaryVarTable & binaryVars);

        // 0.999995
        double FeasibilityOfBinary(const MGBinary & b, int uv1Claz, int uv2Claz,
            const std::vector<Vec3> & vps);
        inline double FeasibilityOfBinary(const MixedGraph & mg, const MGBinaryHandle & bh,
            const MGUnaryVarTable & unaryVars, const std::vector<Vec3> & vps){
            return FeasibilityOfBinary(mg.data(bh),
                unaryVars.at(mg.topo(bh).lowers[0]).claz,
                unaryVars.at(mg.topo(bh).lowers[1]).claz, vps);
        }
        inline bool IsBadBinary(const MixedGraph & mg, const MGBinaryHandle & bh,
            const MGUnaryVarTable & unaryVars, const std::vector<Vec3> & vps){
            return FeasibilityOfBinary(mg, bh, unaryVars, vps) < 0.999;
        }
        inline bool IsGoodBinary(const MixedGraph & mg, const MGBinaryHandle & bh,
            const MGUnaryVarTable & unaryVars, const std::vector<Vec3> & vps){
            return FeasibilityOfBinary(mg, bh, unaryVars, vps) >= 0.999;
        }



        // build mixed graph
        MixedGraph BuildMixedGraph(const std::vector<View<PerspectiveCamera>> & views,
            const std::vector<RegionsGraph> & regionsGraphs, const std::vector<LinesGraph> & linesGraphs,
            const ComponentIndexHashMap<std::pair<RegionIndex, RegionIndex>, double> & regionOverlappingsAcrossViews,
            const ComponentIndexHashMap<std::pair<LineIndex, LineIndex>, Vec3> & lineIncidencesAcrossViews,
            const std::vector<std::map<std::pair<RegionHandle, LineHandle>, std::vector<Point2>>> & regionLineConnections,
            const std::vector<Vec3> & vps,
            MGUnaryVarTable & unaryVars, MGBinaryVarTable & binaryVars,
            double initialDepth = 1.0);



        //// patch on mixed graph
        struct MGPatch {
            MGUnaryVarTable uhs;
            MGBinaryVarTable bhs;

            inline void updateBinaryVars(const MixedGraph & mg, const std::vector<Vec3> & vps) {
                UpdateBinaryVars(mg, vps, uhs, bhs);
            }
            MGPatch firstUnary() const { return MGPatch{ { *uhs.begin() }, {} }; }
            template <class Archive> void serialize(Archive & ar) {
                ar(uhs, bhs);
            }
        };

        inline bool Contains(const MGPatch & p, const MGUnaryHandle & uh) { return Contains(p.uhs, uh); }
        inline bool Contains(const MGPatch & p, const MGBinaryHandle & bh) { return Contains(p.bhs, bh); }

        bool BinaryHandlesAreValidInPatch(const MixedGraph & mg, const MGPatch & patch);
        bool UnariesAreConnectedInPatch(const MixedGraph & mg, const MGPatch & patch);
        double FeasibilityOfPatch(const MixedGraph & mg, const MGPatch & patch, const std::vector<Vec3> & vps);
        inline bool IsGoodPatch(const MixedGraph & mg, const MGPatch & patch, const std::vector<Vec3> & vps){
            return FeasibilityOfPatch(mg, patch, vps) < 0.999;
        }

        double AnchorDistanceSumOnBinaryOfPatch(const MGBinaryHandle & bh, const MGPatch & patch);
        double AnchorDistanceSumOnBinary(const MixedGraph & mg, const MGBinaryHandle & bh, const MGPatch & patch, const std::vector<Vec3> & vps);
        double BinaryDistanceOfPatch(const MGBinaryHandle & bh, const MGPatch & patch);
        double AverageBinaryDistanceOfPatch(const MGPatch & patch, int power = 1);
        double AverageDepthOfPatch(const MGPatch & patch);


        std::vector<MGPatch> SplitMixedGraphIntoPatches(const MixedGraph & mg,
            const MGUnaryVarTable & unaryVars, const MGBinaryVarTable & binaryVars);

        MGPatch MakePatchOnBinary(const MixedGraph & mg, const MGBinaryHandle & bh,
            const MGUnaryVarTable & unaryVars, const MGBinaryVarTable & binaryVars);
        MGPatch MakeStarPatchAroundUnary(const MixedGraph & mg, const MGUnaryHandle & uh,
            const MGUnaryVarTable & unaryVars, const MGBinaryVarTable & binaryVars);

        void CommitPatchToVariableTable(const MGPatch & patch,
            MGUnaryVarTable & unaryVars, MGBinaryVarTable & binaryVars);


        std::vector<MGPatch> SplitPatch(const MixedGraph & mg, const MGPatch & patch,
            std::function<bool(MGBinaryHandle)> useBh);
        inline std::vector<MGPatch> SplitIntoGoodPatches(const MixedGraph & mg, const MGPatch & patch,
            const std::vector<Vec3> & vps){
            return SplitPatch(mg, patch, [&mg, &patch, &vps](MGBinaryHandle bh){
                return IsGoodBinary(mg, bh, patch.uhs, vps);
            });
        }

        MGPatch MinimumSpanningTreePatch(const MixedGraph & mg, const MGPatch & patch,
            std::function<bool(MGBinaryHandle, MGBinaryHandle)> compareBh);
        inline MGPatch MinimumSpanningTreePatch(const MixedGraph & mg, const MGPatch & patch){
            return core::MinimumSpanningTreePatch(mg, patch,
                [&patch](core::MGBinaryHandle a, core::MGBinaryHandle b){
                return core::BinaryDistanceOfPatch(a, patch) <
                    core::BinaryDistanceOfPatch(b, patch);
            });
        }

        bool IsTreePatch(const MixedGraph & mg, const MGPatch & patch);

        void CompletePatch(const MixedGraph & mg, MGPatch & patch, const MGBinaryVarTable & binaryVars);



        void ScalePatch(MGPatch & patch, double scale);
        inline MGPatch & operator *= (MGPatch & patch, double scale){
            ScalePatch(patch, scale);
            return patch;
        }
        inline MGPatch operator * (const MGPatch & patch, double scale) {
            MGPatch result = patch;
            ScalePatch(result, scale);
            return result;
        }
        inline MGPatch & operator /= (MGPatch & patch, double scale){
            patch *= (1.0 / scale); return patch;
        }




        // patch depth optimizer
        class MGPatchDepthsOptimizer {
        public:
            enum AlgorithmType {
                MosekLinearProgramming,
                EigenSparseQR
            };

            MGPatchDepthsOptimizer(const MixedGraph & mg, MGPatch & patch,
                const std::vector<Vec3> & vanishingPoints,
                bool useWeights = true,
                AlgorithmType at = MosekLinearProgramming);
            MGPatchDepthsOptimizer(const MGPatchDepthsOptimizer &) = delete;
            ~MGPatchDepthsOptimizer();

        public:
            void setDepthBounds(double depthLb, double depthUb);
            void setDepthsAllGreaterThan(double lob);
            void setUnaryClass(const MGUnaryHandle & uh, int claz);
            bool optimize();

        private:
            const AlgorithmType _at;
            const MixedGraph & _mg;
            MGPatch & _patch;
            const std::vector<Vec3> & _vanishingPoints;
            void * _internal;
        };


        void ExtandPatch(const MixedGraph & mg, MGPatch & patch,
            const MGUnaryVarTable & unaryVars, const MGBinaryVarTable & binaryVars,
            const std::vector<Vec3> & vps,
            HandledTable<MGUnaryHandle, double> & uhScores,
            HandledTable<MGBinaryHandle, double> & bhScores,
            std::vector<MGUnaryHandle> & uhsOrder,
            double scoreThreshold,
            const std::function<bool(MGUnaryHandle)> & uhIsValid);

        // reconstruction
        /*std::vector<MGPatch> Reconstruct(const MixedGraph & mg,
        const MGUnaryVarTable & unaryVars, const MGBinaryVarTable & binaryVars,
        const std::vector<Vec3> & vps);*/

    }
}
 
#endif