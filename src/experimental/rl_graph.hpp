#ifndef PANORAMIX_EXPERIMENTAL_RL_GRAPH_HPP
#define PANORAMIX_EXPERIMENTAL_RL_GRAPH_HPP


#include "../core/basic_types.hpp"
#include "../core/utility.hpp"
#include "../core/generic_topo.hpp"
#include "../core/cons_graph.hpp"
#include "../core/cameras.hpp"



namespace panoramix {
    namespace experimental {

        using namespace core;

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

        template <class T>
        using RLGraphComponentTable = typename ComponentHandledTableFromConstraintGraph<T, RLGraph>::type;
        template <class T>
        using RLGraphConstraintTable = typename ConstraintHandledTableFromConstraintGraph<T, RLGraph>::type;



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
            double length;
            bool detachable; // if the line lies on the edge of a region, it may be detachable from the region
            template <class Archiver>
            void serialize(Archiver & ar) {
                ar(normalizedAnchors, length, detachable);
            }
        };
        using RegionLineConnectionHandle = ConstraintHandle<RegionLineConnectionData>;




        // junction weight 
        // 7.0, 10.0
        float IncidenceJunctionWeight(bool acrossViews);
        float OutsiderIntersectionJunctionWeight();
        // [0.0 ~ 5.0]
        float ComputeIntersectionJunctionWeightWithLinesVotes(const Mat<float, 3, 2> & votes);




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
            int samplerSizeOnBoundary = 3, int samplerSizeOnLine = 3, bool noBoundaryUnderLines = false);
        std::vector<RegionHandle> AppendRegions(RLGraph & mg, const Imagei & segmentedRegions, const PanoramicCamera & cam,
            double samplingStepAngleOnBoundary, double samplingStepAngleOnLine, 
            int samplerSizeOnBoundary = 3, int samplerSizeOnLine = 3, bool noBoundaryUnderLines = false);



        // get a perfect mask view for a region
        View<PartialPanoramicCamera, Imageub> PerfectRegionMaskView(const RLGraph & mg, RegionHandle rh, double focal = 100.0);







        namespace {
            template <class HandleT> using Decomposed = std::pair<int, HandleT>;
            template <class T> struct Original_ {};
            template <class HandleT> struct Original_<Decomposed<HandleT>> { using type = HandleT; };
            template <class T> using Original = typename Original_<T>::type;
            template <class T> using HandleMapOldToNew = std::map<T, Decomposed<T>>;
            template <class T> using HandleMapNewToOld = std::map<T, Original<T>>;
        }
        using RLGraphOldToNew = MetaBind<HandleMapOldToNew, 
            RegionHandle, 
            LineHandle,
            RegionBoundaryHandle, 
            LineRelationHandle, 
            RegionLineConnectionHandle>;
        using RLGraphNewToOld = MetaBind<HandleMapNewToOld, 
            Decomposed<RegionHandle>, 
            Decomposed<LineHandle>,
            Decomposed<RegionBoundaryHandle>, 
            Decomposed<LineRelationHandle>, 
            Decomposed<RegionLineConnectionHandle>>;

        // connected component ids
        int ConnectedComponents(const RLGraph & mg, RLGraphComponentTable<int> & ccids);
        std::vector<RLGraph> Decompose(const RLGraph & mg, const RLGraphComponentTable<int> & ccids, int ccnum,
            RLGraphOldToNew * old2new = nullptr,
            RLGraphNewToOld * new2old = nullptr);

        
        // rl graph controlers
        struct RLGraphComponentControl {
            bool used; // not used for void/non-planar areas
            int orientationClaz;
            int orientationNotClaz; // if region is tangential with some vp ?
            std::vector<Weighted<Point3>> weightedAnchors;
            std::vector<Bounded<Vec3>> boundedAnchors;
            template <class Archiver>
            inline void serialize(Archiver & ar) {
                ar(used, orientationClaz, orientationNotClaz, weightedAnchors, boundedAnchors);
            }
        };
        using RLGraphComponentControls = RLGraphComponentTable<RLGraphComponentControl>;

        struct RLGraphConstraintControl {
            bool used;
            std::vector<Weighted<Vec3>> weightedAnchors;
            std::vector<Bounded<Vec3>> boundedAnchors;
            template <class Archiver>
            inline void serialize(Archiver & ar) {
                ar(used, weightedAnchors, boundedAnchors); 
            }
        };
        using RLGraphConstraintControls = RLGraphConstraintTable<RLGraphConstraintControl>;


        struct RLGraphControls {
            std::vector<Vec3> vanishingPoints;
            RLGraphComponentControls componentControls;
            RLGraphConstraintControls constraintControls;

            RLGraphControls(){}
            RLGraphControls(const RLGraph & mg, const std::vector<Vec3> & vps);
            RLGraphControls(RLGraphControls && c)
                : vanishingPoints(std::move(c.vanishingPoints)),
                componentControls(std::move(c.componentControls)),
                constraintControls(std::move(c.constraintControls)){
            }
            RLGraphControls & operator = (RLGraphControls && c){
                std::swap(vanishingPoints, c.vanishingPoints);
                std::swap(componentControls, c.componentControls);
                std::swap(constraintControls, c.constraintControls);
                return *this;
            }

            template <class DataT>
            inline RLGraphComponentControl & operator[](ComponentHandle<DataT> h) {
                return componentControls[h];
            }
            template <class DataT>
            inline const RLGraphComponentControl & operator[](ComponentHandle<DataT> h) const {
                return componentControls[h];
            }
            template <class DataT>
            inline RLGraphConstraintControl & operator[](ConstraintHandle<DataT> h) {
                return constraintControls[h];
            }
            template <class DataT>
            inline const RLGraphConstraintControl & operator[](ConstraintHandle<DataT> h) const {
                return constraintControls[h];
            }

            void disable(RegionHandle h, const RLGraph & mg);
            void disable(LineHandle h, const RLGraph & mg);
            void enable(RegionHandle h, const RLGraph & mg);
            void enable(LineHandle h, const RLGraph & mg);
            
            template <class DataT>
            void disable(ConstraintHandle<DataT> h, const RLGraph & mg) { constraintControls[h].used = false; }
            template <class DataT>
            bool enable(ConstraintHandle<DataT> h, const RLGraph & mg) {
                constraintControls[h].used = componentControls[mg.topo(h).constraint<0>()].used &&
                    componentControls[mg.topo(h).constraint<1>()].used;
                return constraintControls[h].used;
            }

            void disableAllInvalidConstraints(const RLGraph & mg);
            void enableAll();

            template <class Archiver>
            inline void serialize(Archiver & ar) {
                ar(vanishingPoints, componentControls, constraintControls);
            }
        };

        
        void SetNecessaryConstraintWeightedAnchors(const RLGraph & mg, RLGraphControls & controls);
        void SetFullConstraintWeightedAnchors(const RLGraph & mg, RLGraphControls & controls);


        struct RLGraphVar {
            std::vector<double> variables;
            template <class Archiver>
            inline void serialize(Archiver & ar) {
                ar(variables);
            }
        };
        using RLGraphVars = RLGraphComponentTable<RLGraphVar>;



        // connected components
        int ConnectedComponents(const RLGraph & mg, const RLGraphControls & controls,
            RLGraphComponentTable<int> & ccids,
            const std::function<bool(const RLGraphConstraintControl &)> & constraintAsConnected);
        std::vector<RLGraphControls> Decompose(const RLGraph & mg, const RLGraphControls & controls,
            const RLGraphComponentTable<int> & ccids, int ccnum);




        // reset all weighted anchor's weights to 1.0
        template <class ConstraintDataT, class FunT>
        inline void SetConstraintControls(RLGraphControls & controls, FunT && fun){
            auto & data = controls.constraintControls.dataOfType<ConstraintHandle<ConstraintDataT>>();
            for (int i = 0; i < data.size(); i++){
                fun(ConstraintHandle<ConstraintDataT>(i), data[i]);
            }
        }
        template <class ComponentDataT, class FunT>
        inline void SetComponentControl(RLGraphControls & controls, FunT && fun){
            auto & data = controls.componentControls.dataOfType<ComponentHandle<ComponentDataT>>();
            for (int i = 0; i < data.size(); i++){
                fun(ComponentHandle<ComponentDataT>(i), data[i]);
            }
        }



        // component instances and related tools
        Line3 Instance(const RLGraph & mg, const RLGraphControls & controls, 
            const RLGraphVars & vars, const LineHandle & lh);
        Plane3 Instance(const RLGraph & mg, const RLGraphControls & controls, 
            const RLGraphVars & vars, const RegionHandle & rh);


        template <class ComponentT>
        using InstanceType = decltype(Instance(std::declval<RLGraph>(), std::declval<RLGraphControls>(),
            std::declval<RLGraphVars>(), std::declval<ComponentHandle<ComponentT>>()));
        template <class ComponentT>
        using InstanceTable = HandledTable<ComponentHandle<ComponentT>, InstanceType<ComponentT>>;

        template <class ComponentT>
        inline InstanceTable<ComponentT> Instances(
            const RLGraph & mg, const RLGraphControls & controls,
            const RLGraphVars & vars) {
            auto instances = mg.createComponentTable<ComponentT, InstanceType<ComponentT>>();
            for (auto & c : mg.components<ComponentT>()){
                if (controls[c.topo.hd].used)
                    instances[c.topo.hd] = Instance(mg, controls, vars, c.topo.hd);
            }
            return instances;
        }

        inline double DepthAt(const Vec3 & direction, const Plane3 & plane, const Point3 & eye = Point3(0, 0, 0)){
           return norm(IntersectionOfLineAndPlane(Ray3(eye, direction), plane).position);
        }
        inline double DepthAt(const Vec3 & direction, const Line3 & line, const Point3 & eye = Point3(0, 0, 0)){
            return norm(DistanceBetweenTwoLines(Ray3(Point3(0, 0, 0), direction), line.ray()).second.first);
        }

        template <class InstanceT>
        inline Point3 PointAt(const Vec3 & direction, const InstanceT & inst, const Point3 & eye = Point3(0, 0, 0)){
            return normalize(direction) * DepthAt(direction, inst, eye);
        }



        // region polygons
        std::vector<Polygon3> RegionPolygon(const RLGraph & mg, const RLGraphControls & controls,
            const RLGraphVars & vars, RegionHandle rh);
        HandledTable<RegionHandle, std::vector<Polygon3>> RegionPolygons(const RLGraph & mg, 
            const RLGraphControls & controls,
            const RLGraphVars & vars);
        double MedianCenterDepth(const RLGraph & mg, const RLGraphControls & controls,
            const RLGraphVars & vars);
        double Score(const RLGraph & mg, const RLGraphControls & controls,
            const RLGraphVars & vars);



        // anchors
        int NumberOfComponentWeightedAnchors(const RLGraphControls & controls);
        int NumberOfComponentBoundedAnchors(const RLGraphControls & controls);
        int NumberOfConstraintWeightedAnchors(const RLGraphControls & controls);
        int NumberOfConstraintBoundedAnchors(const RLGraphControls & controls);

        bool AttachWeightedAnchorToCenterOfLargestLineIfNoExists(const RLGraph & mg,
            RLGraphControls & controls,
            double depth = 1.0, double weight = 1.0, bool orientedOnly = true);
        bool AttachWeightedAnchorToCenterOfLargestRegionIfNoExists(const RLGraph & mg,
            RLGraphControls & controls,
            double depth = 1.0, double weight = 1.0, bool orientedOnly = true);

        void ClearAllComponentAnchors(RLGraphControls & controls);
        void ClearAllConstraintAnchors(RLGraphControls & controls);




        template <class T, int N, class CameraT>
        HandledTable<RegionHandle, Vec<T, N>> CollectFeatureMeanOnRegions(const RLGraph & mg,
            const CameraT & pcam, const ImageOf<Vec<T, N>> & feature){
            HandledTable<RegionHandle, Vec<T, N>> featureMeanTable = mg.createComponentTable<RegionData, Vec<T, N>>();
            for (auto & r : mg.components<RegionData>()){
                auto rh = r.topo.hd;
                auto regionMaskView = PerfectRegionMaskView(mg, rh);
                auto sampler = MakeCameraSampler(regionMaskView.camera, pcam);
                auto featureOnRegion = sampler(feature);
                int votes = 0;
                Vec<T, N> featureSum;
                for (auto it = regionMaskView.image.begin(); it != regionMaskView.image.end(); ++it){
                    if (!*it){
                        continue;
                    }
                    featureSum += featureOnRegion(it.pos());
                    votes += 1;
                }
                auto featureMean = featureSum / std::max(votes, 1);
                featureMeanTable[rh] = featureMean;
            }
            return featureMeanTable;
        }






        void AttachPrincipleDirectionConstraints(const RLGraph & mg, RLGraphControls & controls,
            double rangeAngle = M_PI / 30.0, bool avoidLineConflictions = true);
        void AttachWallConstriants(const RLGraph & mg, RLGraphControls & controls,
            double rangeAngle = M_PI / 60.0, const Vec3 & verticalSeed = Vec3(0, 0, 1));
        void AttachFloorAndCeilingConstraints(const RLGraph & mg,
            RLGraphControls & controls, const RLGraphVars & vars,
            double eyeHeightRatioLowerBound = 0.1, double eyeHeightRatioUpperBound = 0.6,
            double angleThreshold = M_PI / 100.0, const Vec3 & verticalSeed = Vec3(0, 0, 1));



        void LooseOrientationConstraintsOnComponents(const RLGraph & mg,
            RLGraphControls & controls, const RLGraphVars & vars,
            double linesLoosableRatio = 0.2, double regionsLoosableRatio = 0.05,
            double distThresRatio = 0.12);




        RLGraphVars MakeVariables(const RLGraph & mg, const RLGraphControls & controls, bool randomized = false);
        bool HasInfOrNaNValue(const RLGraphVars & v);




        // solve equations
        RLGraphVars SolveVariablesWithoutBoundedAnchors(const RLGraph & mg, 
            const RLGraphControls & controls,
            bool useWeights = false);

        void ResetToFullArmorAnchors(const RLGraph & mg, RLGraphControls & controls);
        void ResetToSampledArmorAnchors(const RLGraph & mg, RLGraphControls & controls, double sampleStepAngle);
        
        RLGraphVars SolveVariablesWithBoundedAnchors(const RLGraph & mg, const RLGraphControls & controls,
            bool useWeights = false, int tryNum = 10);

        void OptimizeVariablesWithBoundedAnchors(const RLGraph & mg, const RLGraphControls & controls, RLGraphVars & vars,
            bool useWeights = true, int tryNum = 100, int maxOptimizeNum = 10,
            const std::function<bool(const RLGraphVars &)> & callback = nullptr);


        void NormalizeVariables(const RLGraph & mg, const RLGraphControls & controls,
            RLGraphVars & vars);

     






       

    }
}


#endif