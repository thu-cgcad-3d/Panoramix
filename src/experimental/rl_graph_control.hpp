#pragma once

#include "rl_graph.hpp"

namespace pano {
    namespace experimental {

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
            Decomposed<RegionLineConnectionHandle >> ;

        // connected component ids
        int ConnectedComponents(const RLGraph & mg, RLGraphComponentTable<int> & ccids);
        std::vector<RLGraph> Decompose(const RLGraph & mg, const RLGraphComponentTable<int> & ccids, int ccnum,
            RLGraphOldToNew * old2new = nullptr,
            RLGraphNewToOld * new2old = nullptr);


        // rl graph controls
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

            RLGraphControls() {}
            RLGraphControls(const RLGraph & mg, const std::vector<Vec3> & vps);
            RLGraphControls(RLGraphControls && c)
                : vanishingPoints(std::move(c.vanishingPoints)),
                componentControls(std::move(c.componentControls)),
                constraintControls(std::move(c.constraintControls)) {}
            RLGraphControls & operator = (RLGraphControls && c) {
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
        void SetMoreConstraintWeightedAnchors(const RLGraph & mg, RLGraphControls & controls, double expandAngle);


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
        inline void SetConstraintControls(RLGraphControls & controls, FunT && fun) {
            auto & data = controls.constraintControls.dataOfType<ConstraintHandle<ConstraintDataT>>();
            for (int i = 0; i < data.size(); i++) {
                fun(ConstraintHandle<ConstraintDataT>(i), data[i]);
            }
        }
        template <class ComponentDataT, class FunT>
        inline void SetComponentControl(RLGraphControls & controls, FunT && fun) {
            auto & data = controls.componentControls.dataOfType<ComponentHandle<ComponentDataT>>();
            for (int i = 0; i < data.size(); i++) {
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
            for (auto & c : mg.components<ComponentT>()) {
                if (controls[c.topo.hd].used)
                    instances[c.topo.hd] = Instance(mg, controls, vars, c.topo.hd);
            }
            return instances;
        }

        inline double DepthAt(const Vec3 & direction, const Plane3 & plane, const Point3 & eye = Point3(0, 0, 0)) {
            return norm(IntersectionOfLineAndPlane(Ray3(eye, direction), plane).position);
        }
        inline double DepthAt(const Vec3 & direction, const Line3 & line, const Point3 & eye = Point3(0, 0, 0)) {
            return norm(DistanceBetweenTwoLines(Ray3(Point3(0, 0, 0), direction), line.ray()).second.first);
        }

        template <class InstanceT>
        inline Point3 PointAt(const Vec3 & direction, const InstanceT & inst, const Point3 & eye = Point3(0, 0, 0)) {
            return normalize(direction) * DepthAt(direction, inst, eye);
        }



        // region polygons
        std::vector<Polygon3> RegionPolygon(const RLGraph & mg, const RLGraphControls & controls,
            const RLGraphVars & vars, RegionHandle rh);
        HandledTable<RegionHandle, std::vector<Polygon3>> RegionPolygons(const RLGraph & mg,
            const RLGraphControls & controls,
            const RLGraphVars & vars);




        std::vector<double> VariableCoefficientsForInverseDepthAtDirection(const RLGraph & mg,
            const RLGraphControls & controls, const RLGraphVars & vars, const Vec3 & direction, const LineHandle & lh);
        std::vector<double> VariableCoefficientsForInverseDepthAtDirection(const RLGraph & mg,
            const RLGraphControls & controls, const RLGraphVars & vars, const Vec3 & direction, const RegionHandle & rh);

        double DepthAtDirectionGivenVariables(const RLGraph & mg, const double * variables,
            const RLGraphControls & controls, const Vec3 & direction, const LineHandle & lh);
        double DepthAtDirectionGivenVariables(const RLGraph & mg, const double * variables,
            const RLGraphControls & controls, const Vec3 & direction, const RegionHandle & rh);






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






        std::vector<RegionHandle> CollectRegionsIntersectingDirection(const Vec3 & dirction, bool alsoConsiderBackward,
            const RLGraph & mg, double rangeAngle = M_PI / 30.0);
        void AttachPrincipleDirectionConstraints(const RLGraph & mg, RLGraphControls & controls,
            double rangeAngle = M_PI / 30.0, bool avoidLineConflictions = true);
        void AttachPrincipleDirectionConstraints2(const RLGraph & mg, RLGraphControls & controls, const Vec3 & up,
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






        double MedianCenterDepth(const RLGraph & mg, const RLGraphControls & controls,
            const RLGraphVars & vars);
        double Score(const RLGraph & mg, const RLGraphControls & controls,
            const RLGraphVars & vars);

        void NormalizeVariables(const RLGraph & mg, const RLGraphControls & controls,
            RLGraphVars & vars);



    }
}

