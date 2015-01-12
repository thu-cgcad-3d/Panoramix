#ifndef PANORAMIX_CORE_MIXED_GRAPH_HPP
#define PANORAMIX_CORE_MIXED_GRAPH_HPP
 
#include "view.hpp"

namespace panoramix {
    namespace core {


        struct MGUnary;
        struct MGBinary;


        // the mixed graph
        // unary variable
        struct MGUnaryVariable {
            bool fixed; // fixed -> not optimizable
            std::vector<double> variables; // (a, b, c) for region plane ax+by+c=1, 1/centerDepth for line
            
            double rawDepth() const;
            Plane3 interpretAsPlane() const;
            Line3 interpretAsLine(const MGUnary & unary, const std::vector<Vec3> & vps) const;
            std::vector<double> variableCoeffsForInverseDepthAtDirection(const Vec3 & direction,
                const MGUnary & unary, const std::vector<Vec3> & vps) const;
            double inverseDepthAtDirection(const Vec3 & direction, const MGUnary & unary, const std::vector<Vec3> & vps) const;
            double depthAtCenter(const MGUnary & unary, const std::vector<Vec3> & vps) const;
            std::vector<double> variableUpperBounds(const MGUnary & unary, const std::vector<Vec3> & vps) const;
            std::vector<double> variableLowerBounds(const MGUnary & unary, const std::vector<Vec3> & vps) const;

            template <class Archive> void serialize(Archive & ar) {
                ar(fixed, variables);
            }
        };

        // unary entity
        struct MGUnary {
            enum Type { Region, Line } type;            
            std::vector<Vec3> normalizedCorners;
            Vec3 normalizedCenter;
            int lineClaz;
            template <class Archive> void serialize(Archive & ar) {
                ar(type, normalizedCorners, normalizedCenter, lineClaz);
            }
        };        
      

        struct MGBinaryVariable {
            bool enabled;
            template <class Archive> void serialize(Archive & ar) {
                ar(enabled);
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



        // build mixed graph
        MixedGraph BuildMixedGraph(const std::vector<View<PerspectiveCamera>> & views,
            const std::vector<RegionsGraph> & regionsGraphs, const std::vector<LinesGraph> & linesGraphs,
            const ComponentIndexHashMap<std::pair<RegionIndex, RegionIndex>, double> & regionOverlappingsAcrossViews,
            const ComponentIndexHashMap<std::pair<LineIndex, LineIndex>, Vec3> & lineIncidencesAcrossViews,
            const std::vector<std::map<std::pair<RegionHandle, LineHandle>, std::vector<Point2>>> & regionLineConnections,
            const std::vector<Vec3> & vps,
            MGUnaryVarTable & unaryVars, MGBinaryVarTable& binaryVars,
            double initialDepth = 1.0);

        MixedGraph BuildMixedGraph(const std::vector<View<PerspectiveCamera>> & views, 
            std::vector<Vec3> & vps,
            MGUnaryVarTable & unaryVars, MGBinaryVarTable& binaryVars,

            double initialDepth = 1.0,
            const core::LineSegmentExtractor & lineseger = core::LineSegmentExtractor(),
            double intersectionDistanceThreshold = 10,
            double incidenceDistanceAlongDirectionThreshold = 100,
            double incidenceDistanceVerticalDirectionThreshold = 5,
            
            const core::SegmentationExtractor & segmenter = core::SegmentationExtractor(),
            double samplingStepLengthOnBoundary = 10.0, 
            double samplingStepLengthOnLines = 5.0,
            int dilationSize = 3,

            double interViewIncidenceAngleAlongDirectionThreshold = M_PI_4,
            double interViewIncidenceAngleVerticalDirectionThreshold = 20.0 / 300.0);





        //// patch on mixed graph
        struct MGPatch {
            MGUnaryVarTable uhs;
            MGBinaryVarTable bhs;
            MGPatch firstUnary() const { return MGPatch{ { *uhs.begin() }, {} }; }
            template <class Archive> void serialize(Archive & ar) {
                ar(uhs, bhs);
            }
        };

        inline bool Contains(const MGPatch & p, const MGUnaryHandle & uh) { return Contains(p.uhs, uh); }
        inline bool Contains(const MGPatch & p, const MGBinaryHandle & bh) { return Contains(p.bhs, bh); }


        bool BinaryHandlesAreValidInPatch(const MixedGraph & mg, const MGPatch & patch);
        bool UnariesAreConnectedInPatch(const MixedGraph & mg, const MGPatch & patch);

        MGPatch MakePatchOnBinary(const MixedGraph & mg, const MGBinaryHandle & bh,
            const MGUnaryVarTable & unaryVars, const MGBinaryVarTable & binaryVars);
        MGPatch MakeStarPatchAroundUnary(const MixedGraph & mg, const MGUnaryHandle & uh,
            const MGUnaryVarTable & unaryVars, const MGBinaryVarTable & binaryVars);


        double AnchorDistanceSumOnBinaryOfPatch(const MixedGraph & mg, const MGBinaryHandle & bh,
            const MGPatch & patch, const std::vector<Vec3> & vps);
        double BinaryDistanceOfPatch(const MixedGraph & mg, const MGBinaryHandle & bh,
            const MGPatch & patch, const std::vector<Vec3> & vps);
        double AverageBinaryDistanceOfPatch(const MixedGraph & mg,
            const MGPatch & patch, const std::vector<Vec3> & vps);
        double AverageUnaryCenterDepthOfPatch(const MixedGraph & mg,
            const MGPatch & patch, const std::vector<Vec3> & vps);
        double AverageRawDepthOfPatch(const MGPatch & patch);

        
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




        std::vector<MGPatch> SplitMixedGraphIntoPatches(const MixedGraph & mg,
            const MGUnaryVarTable & unaryVars, const MGBinaryVarTable & binaryVars);

        std::vector<MGPatch> SplitPatch(const MixedGraph & mg, const MGPatch & patch,
            std::function<bool(MGBinaryHandle)> useBh);

        MGPatch MinimumSpanningTreePatch(const MixedGraph & mg, const MGPatch & patch,
            std::function<bool(MGBinaryHandle, MGBinaryHandle)> compareBh);
        inline MGPatch MinimumSpanningTreePatch(const MixedGraph & mg, const MGPatch & patch, 
            const std::vector<Vec3> & vps){
            std::unordered_map<MGBinaryHandle, double> binaryDistances;
            for (auto & bhv : patch.bhs){
                binaryDistances[bhv.first] = core::BinaryDistanceOfPatch(mg, bhv.first, patch, vps);
            }
            return core::MinimumSpanningTreePatch(mg, patch,
                [&binaryDistances](core::MGBinaryHandle a, core::MGBinaryHandle b){
                return binaryDistances.at(a) < binaryDistances.at(b);
            });
        }



        // patch depth optimizer
        class MGPatchDepthsOptimizer {
        public:
            enum Algorithm {
                MosekLinearProgramming,
                MosekLinearProgrammingSimplified,
                EigenSparseQR,
                EigenSparseQRSimplified,
                MATLAB_CVX
            };

            MGPatchDepthsOptimizer(const MixedGraph & mg, MGPatch & patch,
                const std::vector<Vec3> & vanishingPoints,
                bool useWeights = false,
                Algorithm at = EigenSparseQR);
            MGPatchDepthsOptimizer(const MGPatchDepthsOptimizer &) = delete;
            MGPatchDepthsOptimizer(MGPatchDepthsOptimizer && pdo);
            ~MGPatchDepthsOptimizer();

        public:
            bool optimize();

        private:
            const Algorithm _at;
            const MixedGraph & _mg;
            MGPatch & _patch;
            const std::vector<Vec3> & _vanishingPoints;
            void * _internal;
        };






    }
}


#endif