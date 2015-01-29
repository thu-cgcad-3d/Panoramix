#ifndef PANORAMIX_CORE_MIXED_GRAPH_HPP
#define PANORAMIX_CORE_MIXED_GRAPH_HPP


#include "basic_types.hpp"
#include "generic_topo.hpp"
#include "cons_graph.hpp"
#include "cameras.hpp"

//#include "surface_labels.hpp"

namespace panoramix {
    namespace core {

        // the mixe graph definition
        struct RegionData;
        struct RegionBoundaryData;
        struct LineData;
        struct LineRelationData;
        struct RegionLineConnectionData;


        using MixedGraph = ConstraintGraph<std::tuple<RegionData, LineData>,
            std::tuple<
            ConstraintConfig<RegionBoundaryData, RegionData, RegionData>,
            ConstraintConfig<LineRelationData, LineData, LineData>,
            ConstraintConfig<RegionLineConnectionData, RegionData, LineData>
            >
        >;   



        // regions
        struct RegionData {
            Vec3 normalizedCenter;
            double area; // projected onto the plane rooted at normalizedCenter
            std::vector<std::vector<Vec3>> normalizedContours;
            template <class Archiver>
            void serialize(Archiver & ar) {
                ar(normalizedCenter, area, normalizedContours);
            }
        };
        using RegionHandle = ComponentHandle<RegionData>;
        struct RegionBoundaryData {
            std::vector<std::vector<Vec3>> normalizedEdges;
            double length;
            std::vector<std::vector<Vec3>> normalizedSampledPoints;
            template <class Archiver>
            void serialize(Archiver & ar) {
                ar(edges, length, normalizedSampledPoints);
            }
        };
        using RegionBoundaryHandle = ConstraintHandle<RegionBoundaryData>;



        // lines
        struct LineData {
            Classified<Line3> line;
            template <class Archiver>
            void serialize(Archiver & ar) {
                ar(line);
            }
        };
        using LineHandle = ComponentHandle<LineData>;
        struct LineRelationData {
            Vec3 normalizedRelationCenter;
            float junctionWeight;
            enum Type {
                Incidence,
                Intersection
            } type;
            template <class Archiver>
            void serialize(Archiver & ar) {
                ar(normalizedRelationCenter, junctionWeight, type);
            }
        };
        using LineRelationHandle = ConstraintHandle<LineRelationData>;


        // region and line
        struct RegionLineConnectionData {
            std::vector<Vec3> normalizedAnchors;
            template <class Archiver>
            void serialize(Archiver & ar) {
                ar(normalizedAnchors);
            }
        };
        using RegionLineConnectionHandle = ConstraintHandle<RegionLineConnectionData>;





        // junction weight 
        float IncidenceJunctionWeight(bool acrossViews);
        float OutsiderIntersectionJunctionWeight();
        float ComputeIntersectionJunctionWeightWithLinesVotes(const Mat<float, 3, 2> & votes);


        // estimate vanishing points
        std::vector<Vec3> EstimateVanishingPointsAndClassifyLines(const std::vector<PerspectiveCamera> & cams,
            std::vector<std::vector<Classified<Line2>>> & lineSegments);



        // add lines to mixed graph from classified line segments
        void AddLines(MixedGraph & mg, const std::vector<Classified<Line2>> & lineSegments,
            const PerspectiveCamera & cam,
            const std::vector<Vec3> & vps,
            double intersectionAngleThreshold = 0.04,
            double incidenceAngleAlongDirectionThreshold = 0.1,
            double incidenceAngleVerticalDirectionThreshold = 0.02,
            double interViewIncidenceAngleAlongDirectionThreshold = 0.15, // for new line-line incidence recognition
            double interViewIncidenceAngleVerticalDirectionThreshold = 0.03,
            bool includeUnclassifiedLines = false);


        // add more regions and related constraints to mixed graph
        void AddRegions(MixedGraph & mg, const Imagei & segmentedRegions, const PerspectiveCamera & cam,
            double samplingStepAngleOnBoundary, double samplingStepAngleOnLine);
        void AddRegions(MixedGraph & mg, const Imagei & segmentedRegions, const PanoramicCamera & cam,
            double samplingStepAngleOnBoundary, double samplingStepAngleOnLine);






    }
}


#if 0

namespace panoramix {
    namespace core {

        struct MGUnary;
        struct MGBinary;


        // the mixed graph
        // unary variable
        struct MGUnaryVariable {
            bool fixed; // fixed -> not optimizable, variables stay the same during optimization

            // (a, b, c) for RegionFree ax+by+c=1, 
            // (a, b) for RegionAlongFixedAxis  ax+by+c=1,
            // 1/centerDepth for RegionWithFixedNormal
            // (1/cornerDepth1, 1/cornerDepth2) for LineFree
            // 1/centerDepth for LineOriented,
            std::vector<double> variables; 

            inline MGUnaryVariable(bool f = false) : fixed(f) {}
            explicit MGUnaryVariable(const MGUnary & u, bool f = false);
            
            double rawDepth() const;
            Plane3 interpretAsPlane(const MGUnary & unary) const;
            Line3 interpretAsLine(const MGUnary & unary) const;
            std::vector<double> variableCoeffsForInverseDepthAtDirection(const Vec3 & direction, const MGUnary & unary) const;
            double inverseDepthAtDirection(const Vec3 & direction, const MGUnary & unary) const;         
            double depthAtCenter(const MGUnary & unary) const;

            void fitToClosestOrientation(MGUnary & unary, const std::vector<Vec3> & vps, double angleThreshold);

            template <class Archive> void serialize(Archive & ar) {
                ar(fixed, variables);
            }
        };


        // unary entity
        struct MGUnary {
            enum Type { 
                Undefined = -1,
                RegionFree, RegionAlongFixedAxis, RegionWithFixedNormal, 
                LineFree, LineOriented,
            } type;            
            
            std::vector<Vec3> normalizedCorners;
            Vec3 normalizedCenter;            
            Vec3 normalizedOrientation;

            enum Material {
                Unknown = -1, Planar, NonPlanar, Void
            } material;

            bool hasOrientationConstraints() const;
            bool isRegion() const;
            bool isLine() const;

            template <class Archive> void serialize(Archive & ar) {
                ar(type, normalizedCorners, normalizedCenter, normalizedOrientation, material);
            }
        };  
      

        struct MGBinaryVariable {
            enum Type {
                Disconnected,
                Connected
            } type;

            inline MGBinaryVariable() : type(Connected) {}
            explicit MGBinaryVariable(const MGBinary & b) : type(Connected) {}
            template <class Archive> void serialize(Archive & ar) {
                ar(type);
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




        // build mixed graph
        MixedGraph BuildMixedGraph(const std::vector<View<PerspectiveCamera>> & views,
            const std::vector<RegionsGraph> & regionsGraphs, const std::vector<LinesGraph> & linesGraphs,
            const ComponentIndexHashMap<std::pair<RegionIndex, RegionIndex>, double> & regionOverlappingsAcrossViews,
            const ComponentIndexHashMap<std::pair<LineIndex, LineIndex>, Vec3> & lineIncidencesAcrossViews,
            const std::vector<std::map<std::pair<RegionHandle, LineHandle>, std::vector<Point2>>> & regionLineConnections,
            const std::vector<Vec3> & vps,
            double initialDepth = 1.0);

        MixedGraph BuildMixedGraph(const std::vector<View<PerspectiveCamera>> & views,
            std::vector<Vec3> & vps,

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
        using MGUnaryVarTable = std::unordered_map<MGUnaryHandle, MGUnaryVariable>;
        using MGBinaryVarTable = std::unordered_map<MGBinaryHandle, MGBinaryVariable>;
        struct MGPatch {
            MGUnaryVarTable uhs;
            MGBinaryVarTable bhs;
            bool empty() const { return uhs.empty(); }
            MGPatch firstUnary() const { return MGPatch{ { *uhs.begin() }, {} }; }
            template <class Archive> void serialize(Archive & ar) {
                ar(uhs, bhs);
            }
        };

        inline bool Contains(const MGPatch & p, const MGUnaryHandle & uh) { return Contains(p.uhs, uh); }
        inline bool Contains(const MGPatch & p, const MGBinaryHandle & bh) { return Contains(p.bhs, bh); }


        bool BinaryHandlesAreValidInPatch(const MixedGraph & mg, const MGPatch & patch);
        bool UnariesAreConnectedInPatch(const MixedGraph & mg, const MGPatch & patch);

        MGPatch MakePatchOnBinary(const MixedGraph & mg, const MGBinaryHandle & bh, const std::vector<Vec3> & vps);
        MGPatch MakeStarPatchAroundUnary(const MixedGraph & mg, const MGUnaryHandle & uh, const std::vector<Vec3> & vps);


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




        std::vector<MGPatch> SplitMixedGraphIntoPatches(const MixedGraph & mg, const std::vector<Vec3> & vps);

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
                EigenSparseQR,
                MATLAB_CVX,
                MATLAB_CVX_v2
            };

            MGPatchDepthsOptimizer(const MixedGraph & mg, MGPatch & patch,
                const std::vector<Vec3> & vanishingPoints,
                bool useWeights = false,
                Algorithm at = EigenSparseQR);
            MGPatchDepthsOptimizer(const MGPatchDepthsOptimizer &) = delete;
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



        void FitUnariesToClosestOrientations(const MixedGraph & mg, MGPatch & patch,
            const std::vector<Vec3> & vps, double angleThreshold);


    }
}

#endif

#endif