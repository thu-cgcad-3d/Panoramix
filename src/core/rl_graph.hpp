#ifndef PANORAMIX_CORE_RL_GRAPH_HPP
#define PANORAMIX_CORE_RL_GRAPH_HPP


#include "basic_types.hpp"
#include "generic_topo.hpp"
#include "cons_graph.hpp"
#include "cameras.hpp"

namespace panoramix {
    namespace core {

        // the mixe region-line graph definition
        struct RegionData;
        struct RegionBoundaryData;
        struct LineData;
        struct LineRelationData;
        struct RegionLineConnectionData;


        using RLGraph = ConstraintGraph<std::tuple<RegionData, LineData>,
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
                ar(normalizedEdges, length, normalizedSampledPoints);
            }
        };
        using RegionBoundaryHandle = ConstraintHandle<RegionBoundaryData>;



        // lines
        struct LineData {
            Line3 line;
            int initialClaz;
            template <class Archiver>
            void serialize(Archiver & ar) {
                ar(line, initialClaz);
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
            bool detachable; // if the line lies on the edge of a region, it may be detachable from the region
            template <class Archiver>
            void serialize(Archiver & ar) {
                ar(normalizedAnchors, detachable);
            }
        };
        using RegionLineConnectionHandle = ConstraintHandle<RegionLineConnectionData>;




        // junction weight 
        float IncidenceJunctionWeight(bool acrossViews);
        float OutsiderIntersectionJunctionWeight();
        float ComputeIntersectionJunctionWeightWithLinesVotes(const Mat<float, 3, 2> & votes);


        // estimate vanishing points
        std::vector<Vec3> EstimateVanishingPointsAndClassifyLines(const PerspectiveCamera & cams,
            std::vector<Classified<Line2>> & lineSegments);
        std::vector<Vec3> EstimateVanishingPointsAndClassifyLines(const std::vector<PerspectiveCamera> & cams,
            std::vector<std::vector<Classified<Line2>>> & lineSegments);


        // add lines to rl graph from classified line segments
        void AppendLines(RLGraph & mg, const std::vector<Classified<Line2>> & lineSegments,
            const PerspectiveCamera & cam,
            const std::vector<Vec3> & vps,
            double intersectionAngleThreshold = 0.04,
            double incidenceAngleAlongDirectionThreshold = 0.1,
            double incidenceAngleVerticalDirectionThreshold = 0.02,
            double interViewIncidenceAngleAlongDirectionThreshold = 0.15, // for new line-line incidence recognition
            double interViewIncidenceAngleVerticalDirectionThreshold = 0.03);


        // add more regions and related constraints to rl graph
        std::vector<RegionHandle> AppendRegions(RLGraph & mg, const Imagei & segmentedRegions, const PerspectiveCamera & cam,
            double samplingStepAngleOnBoundary, double samplingStepAngleOnLine,
            int samplerSizeOnBoundary = 3, int samplerSizeOnLine = 3);
        std::vector<RegionHandle> AppendRegions(RLGraph & mg, const Imagei & segmentedRegions, const PanoramicCamera & cam,
            double samplingStepAngleOnBoundary, double samplingStepAngleOnLine, 
            int samplerSizeOnBoundary = 3, int samplerSizeOnLine = 3);




        // mixed grpah property
        struct RLGraphComponentProperty {
            bool used; // not used for void/non-planar areas
            int orientationClaz;
            int orientationNotClaz; // if region is tangential with some vp ?
            std::vector<double> variables;
            template <class Archiver>
            inline void serialize(Archiver & ar) {
                ar(used, orientationClaz, orientationNotClaz, variables);
            }
        };
        struct RLGraphConstraintProperty {
            bool used;
            template <class Archiver>
            inline void serialize(Archiver & ar) { 
                ar(used);
            }
        };
        
        using RLGraphComponentPropertyTable = 
            ComponentHandledTableFromConstraintGraph<RLGraphComponentProperty, RLGraph>::type;
        using RLGraphConstraintPropertyTable =
            ConstraintHandledTableFromConstraintGraph<RLGraphConstraintProperty, RLGraph>::type;
        
        // property table
        struct RLGraphPropertyTable {
            std::vector<Vec3> vanishingPoints;
            RLGraphComponentPropertyTable componentProperties;
            RLGraphConstraintPropertyTable constraintProperties;
            template <class DataT> 
            inline RLGraphComponentProperty & operator[](ComponentHandle<DataT> h) { 
                return componentProperties[h]; 
            }
            template <class DataT>
            inline const RLGraphComponentProperty & operator[](ComponentHandle<DataT> h) const {
                return componentProperties[h]; 
            }
            template <class DataT>
            inline RLGraphConstraintProperty & operator[](ConstraintHandle<DataT> h) {
                return constraintProperties[h]; 
            }
            template <class DataT>
            inline const RLGraphConstraintProperty & operator[](ConstraintHandle<DataT> h) const {
                return constraintProperties[h]; 
            }

            template <class Archiver>
            inline void serialize(Archiver & ar) {
                ar(vanishingPoints, componentProperties, constraintProperties); 
            }
        };

        // make property table
        RLGraphPropertyTable MakeRLGraphPropertyTable(const RLGraph & mg, const std::vector<Vec3> & vps);

        // reset variables
        void ResetVariables(const RLGraph & mg, RLGraphPropertyTable & props);

        // make dangling constraints unused
        void UpdateConstraintUsabilities(const RLGraph & mg, RLGraphPropertyTable & props,
            bool enableDisabledConstraints = true);

        // get component instances
        Line3 Instance(const RLGraph & mg, const RLGraphPropertyTable & props, const LineHandle & lh);
        Plane3 Instance(const RLGraph & mg, const RLGraphPropertyTable & props, const RegionHandle & rh);



        // solve equations
        void SolveVariablesUsingInversedDepths(const RLGraph & mg, RLGraphPropertyTable & props, bool useWeights = true);
        void SolveVariablesUsingNormalDepths(const RLGraph & mg, RLGraphPropertyTable & props, bool useWeights = true);
        
        // model properties
        double ComponentMedianCenterDepth(const RLGraph & mg, const RLGraphPropertyTable & props);
        double ComputeScore(const RLGraph & mg, const RLGraphPropertyTable & props);

        void NormalizeVariables(const RLGraph & mg, RLGraphPropertyTable & props);





        // adjust constraints
        void AttachPrincipleDirectionConstraints(const RLGraph & mg, RLGraphPropertyTable & props, 
            double rangeAngle = M_PI / 100.0, bool avoidLineConflictions = true);
        void AttachWallConstriants(const RLGraph & mg, RLGraphPropertyTable & props,
            double rangeAngle = M_PI / 100.0, const Vec3 & verticalSeed = Vec3(0, 0, 1));
        void AttachFloorAndCeilingConstraints(const RLGraph & mg, RLGraphPropertyTable & props,
            double eyeHeightRatioLowerBound = 0.1, double eyeHeightRatioUpperBound = 0.6,
            double angleThreshold = M_PI / 100.0, const Vec3 & verticalSeed = Vec3(0, 0, 1));

        void AttachGeometricContextConstraints(const RLGraph & mg, RLGraphPropertyTable & props, 
            const std::vector<GeometricContextEstimator::Feature> & perspectiveGCs,
            const std::vector<PerspectiveCamera> & gcCameras,
            int shrinkRegionOrientationIteration = 1, bool considerGCVerticalConstraint = false);


        void LooseOrientationConstraintsOnComponents(const RLGraph & mg, RLGraphPropertyTable & props,
            double linesLoosableRatio = 0.2, double regionsLoosableRatio = 0.05, double distThresRatio = 0.12);


        void LooseMaybeOcclusionBoundaryConstraints(const RLGraph & mg, RLGraphPropertyTable & props);





        // visualize current mixed graph
        void Visualize(const View<PanoramicCamera> & texture, 
            const RLGraph & mg, RLGraphPropertyTable & props);
        void Visualize(const View<PerspectiveCamera> & texture,
            const RLGraph & mg, RLGraphPropertyTable & props);




        // refine the mixed graph variables to a structure
        HandledTable<RegionHandle, int> ClusterRegions(const RLGraph & mg, const RLGraphPropertyTable & props);



    }
}


#endif