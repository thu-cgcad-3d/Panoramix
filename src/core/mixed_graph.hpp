#ifndef PANORAMIX_CORE_MIXED_GRAPH_HPP
#define PANORAMIX_CORE_MIXED_GRAPH_HPP


#include "basic_types.hpp"
#include "generic_topo.hpp"
#include "cons_graph.hpp"
#include "cameras.hpp"

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
        void AppendLines(MixedGraph & mg, const std::vector<Classified<Line2>> & lineSegments,
            const PerspectiveCamera & cam,
            const std::vector<Vec3> & vps,
            double intersectionAngleThreshold = 0.04,
            double incidenceAngleAlongDirectionThreshold = 0.1,
            double incidenceAngleVerticalDirectionThreshold = 0.02,
            double interViewIncidenceAngleAlongDirectionThreshold = 0.15, // for new line-line incidence recognition
            double interViewIncidenceAngleVerticalDirectionThreshold = 0.03);


        // add more regions and related constraints to mixed graph
        void AppendRegions(MixedGraph & mg, const Imagei & segmentedRegions, const PerspectiveCamera & cam,
            double samplingStepAngleOnBoundary, double samplingStepAngleOnLine,
            int samplerSizeOnBoundary = 3, int samplerSizeOnLine = 3);
        void AppendRegions(MixedGraph & mg, const Imagei & segmentedRegions, const PanoramicCamera & cam,
            double samplingStepAngleOnBoundary, double samplingStepAngleOnLine, 
            int samplerSizeOnBoundary = 3, int samplerSizeOnLine = 3);





        struct MixedGraphComponentProperty {
            bool used; // not used for void/non-planar areas
            int orientationClaz;
            int orientationNotClaz; // if region is tangential with some vp ?
            std::vector<double> variables;
            template <class Archiver>
            inline void serialize(Archiver & ar) {
                ar(used, orientationClaz, orientationNotClaz, variables);
            }
        };
        struct MixedGraphConstraintProperty {
            bool used;
            template <class Archiver>
            inline void serialize(Archiver & ar) { 
                ar(used);
            }
        };
        
        using MixedGraphComponentPropertyTable = 
            ComponentHandledTableFromConstraintGraph<MixedGraphComponentProperty, MixedGraph>::type;
        using MixedGraphConstraintPropertyTable =
            ConstraintHandledTableFromConstraintGraph<MixedGraphConstraintProperty, MixedGraph>::type;
        
        struct MixedGraphPropertyTable {
            std::vector<Vec3> vanishingPoints;
            MixedGraphComponentPropertyTable componentProperties;
            MixedGraphConstraintPropertyTable constraintProperties;
            template <class DataT> 
            inline MixedGraphComponentProperty & operator[](ComponentHandle<DataT> h) { 
                return componentProperties[h]; 
            }
            template <class DataT>
            inline const MixedGraphComponentProperty & operator[](ComponentHandle<DataT> h) const {
                return componentProperties[h]; 
            }
            template <class DataT>
            inline MixedGraphConstraintProperty & operator[](ConstraintHandle<DataT> h) {
                return constraintProperties[h]; 
            }
            template <class DataT>
            inline const MixedGraphConstraintProperty & operator[](ConstraintHandle<DataT> h) const {
                return constraintProperties[h]; 
            }

            template <class Archiver>
            inline void serialize(Archiver & ar) {
                ar(vanishingPoints, componentProperties, constraintProperties); 
            }
        };
        MixedGraphPropertyTable MakeMixedGraphPropertyTable(const MixedGraph & mg, const std::vector<Vec3> & vps);
        MixedGraphPropertyTable MakeMixedGraphPropertyTable(const MixedGraph & mg, const std::vector<Vec3> & vps,
            const std::vector<GeometricContextEstimator::Feature> & perspectiveGCs,
            const std::vector<PerspectiveCamera> & gcCameras, int shrinkRegionOrientationIteration = 1);


        void InitializeVariables(const MixedGraph & mg, MixedGraphPropertyTable & props);

        Line3 Instance(const MixedGraph & mg, const MixedGraphPropertyTable & props, const LineHandle & lh);
        Plane3 Instance(const MixedGraph & mg, const MixedGraphPropertyTable & props, const RegionHandle & rh);     

        std::vector<double> VariableCoefficientsForInverseDepthAtDirection(const MixedGraph & mg, 
            const MixedGraphPropertyTable & props, const Vec3 & direction, const LineHandle & lh);
        std::vector<double> VariableCoefficientsForInverseDepthAtDirection(const MixedGraph & mg,
            const MixedGraphPropertyTable & props, const Vec3 & direction, const RegionHandle & rh);

        double DepthAtDirectionGivenVariables(const MixedGraph & mg, const double * variables,
            const MixedGraphPropertyTable & props, const Vec3 & direction, const LineHandle & lh);
        double DepthAtDirectionGivenVariables(const MixedGraph & mg, const double * variables,
            const MixedGraphPropertyTable & props, const Vec3 & direction, const RegionHandle & rh);

        void SolveVariablesUsingInversedDepths(const MixedGraph & mg, MixedGraphPropertyTable & props);
        void SolveVariablesUsingNormalDepths(const MixedGraph & mg, MixedGraphPropertyTable & props);
        
        double ComponentMedianCenterDepth(const MixedGraph & mg, const MixedGraphPropertyTable & props);
        void NormalizeVariables(const MixedGraph & mg, MixedGraphPropertyTable & props);
        double ComputeScore(const MixedGraph & mg, const MixedGraphPropertyTable & props);

        void LooseOrientationConstraintsOnLines(const MixedGraph & mg, MixedGraphPropertyTable & props,
            double loosableRatio = 0.2);



        void Visualize(const View<PanoramicCamera> & texture, 
            const MixedGraph & mg, MixedGraphPropertyTable & props);
        void Visualize(const View<PerspectiveCamera> & texture,
            const MixedGraph & mg, MixedGraphPropertyTable & props);


    }
}


#endif