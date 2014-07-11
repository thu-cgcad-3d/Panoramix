#ifndef PANORAMIX_CORE_RECONSTRUCTION_ENGINE_HPP
#define PANORAMIX_CORE_RECONSTRUCTION_ENGINE_HPP

#include "../core/basic_types.hpp"
#include "../core/feature.hpp"
#include "../core/utilities.hpp"

#include "../deriv/derivative.hpp"

#include "regions_net.hpp"

namespace panoramix {
    namespace rec { 

        using namespace core;


        template <class T>
        class DisableableExpression {
        public:
            inline DisableableExpression() : _enabled(nullptr) {}
            inline explicit DisableableExpression(const deriv::Expression<T> & rawExpr, deriv::ExpressionGraph & g)
                : _enabled(std::make_shared<bool>(true)) {
                auto enabledExpr = deriv::composeFunction(g, [this]()
                    -> double {return *_enabled ? 1.0 : -1.0; });
                _expr = deriv::cwiseSelect(enabledExpr, rawExpr, 0.0);
            }

            inline deriv::Expression<T> toExpression() const { return _expr; }
            inline void setEnabled(bool b) { *_enabled = b; }
            inline void enable() { *_enabled = true; }
            inline void disable() { *_enabled = false; }

        private:
            deriv::Expression<T> _expr;
            std::shared_ptr<bool> _enabled;
        };

        using EHandleTable = std::vector < deriv::EHandle > ;

        template <class T>
        class OptimizibleExpression {
        public:
            inline OptimizibleExpression() :_data(0), _lastChange(0) {}
            inline explicit OptimizibleExpression(const T & d) : _data(d), _lastChange(0) {}

            void registerHandleTable(EHandleTable & table) { 
                table.push_back(_expr.handle()); 
                _positionInHandleTable = table.size() - 1;
            }

            void getDerivative(const EHandleTable & derivTable) {
                _dexpr = _expr.g()->asDerived<T>(derivTable[_positionInHandleTable]);
            }

            void optimizeData(double delta, double momentum){

            }

        private:
            deriv::Expression<T> _expr;
            deriv::DerivativeExpression<T> _dexpr;
            int _positionInHandleTable;
            T _data;
            T _lastChange;
        };


        // engine
        // ReconstructionEngine
        class ReconstructionEngine {
        public:
            struct Params {
                Params();
                PanoramicCamera camera; // camera for generating the panorama
                double lineSegmentWeight;
                double siftWeight;
                double surfWeight;
                LineSegmentExtractor lineSegmentExtractor;
                CVFeatureExtractor<cv::SIFT> siftExtractor;
                CVFeatureExtractor<cv::SURF> surfExtractor;
                SegmentationExtractor segmenter;

                // angle scalar to judge whether two views may share certain common features
                double cameraAngleScaler; 
                // angle scalar to judge whether two views are too close
                double smallCameraAngleScalar; 
                // angle threshold to judge whether two lines are constrained (intersection/incidence), 
                // used for building the Constraint Graph
                double intersectionConstraintLineDistanceAngleThreshold;
                double incidenceConstraintLineDistanceAngleThreshold; // for incidence constraits
                // angle threshold to judge whether two colineared lines should be merged
                double mergeLineDistanceAngleThreshold;
                // manhattan junction weights
                // mjWeightTriplet includes the Y W and K junctions
                double mjWeightTriplet, mjWeightX, mjWeightT, mjWeightL, mjWeightI;
            };

            struct ViewData;
            struct ViewConnectionData;
            using ViewsGraph = GraphicalModel02<ViewData, ViewConnectionData>;
            using ViewHandle = HandleAtLevel<0>;
            using ViewConnectionHandle = HandleAtLevel<1>;

            struct ComponentData;
            struct ConstraintData;
            using ConstraintGraph = GraphicalModel02<ComponentData, ConstraintData>;
            using ComponentHandle = HandleAtLevel<0>;
            using ConstraintHandle = HandleAtLevel<1>;

        public:
            inline explicit ReconstructionEngine(const Params params = Params()) : _params(params) {}
            inline const Params & params() const { return _params; }
            
            inline ViewHandle insertView(const ViewData & vd) { return _views.add(vd); }

            // insert a new photo, with known parameters
            ViewHandle insertPhoto(const Image & im, const PerspectiveCamera & cam, 
                double cameraDirectionErrorScale = 0.0);

            // insert a new panoramic image and connect close views and set panorama
            void insertPanorama(const Image & panorama, const std::vector<PerspectiveCamera> & viewCams, 
                const PanoramicCamera & panCam);

            // compute features for a single view
            void computeFeatures(ViewHandle h);

            // segment view image and build net of regions for a single view
            // after computeFeatures(h)
            void buildRegionNet(ViewHandle h);

            // connect this view to neighbor views who may overlap with h
            size_t updateConnections(ViewHandle h);

            // whether this view overlaps some existing views a lot, measured by the smallCameraAngleScalar parameter
            ViewHandle isTooCloseToAnyExistingView(ViewHandle h) const;

            // find the feature matches to connected views
            // after computeFeatures(h) and buildRTrees(h)
            void findMatchesToConnectedViews(ViewHandle h); 

            // calibrate all cameras
            // after findMatchesToConnectedViews(h)
            void calibrateAllCameras();

            // stitch panorama
            void stitchPanorama();
            
            // estimate vanishing points using lines extract from all views, classify this lines and lift them all to space
            void estimateVanishingPointsAndClassifyLines();

            // build constraints on spatial lines and rectify their parameters to build a more reasonable 3D sketch
            void rectifySpatialLines();

            // reconstruct regions
            void reconstructFaces();
                
        public:

            //// views
            // view data
            struct ViewData {
                // cameras
                PerspectiveCamera originalCamera, camera;
                double cameraDirectionErrorScale;
                
                // image and 2d image features
                Image image;
                std::vector<Classified<Line2>> lineSegments;                
                std::vector<HPoint2> lineSegmentIntersections;
                std::vector<std::pair<int, int>> lineSegmentIntersectionLineIDs;                
                
                std::vector<KeyPoint> keypointsForMatching;
                cv::Mat descriptorsForMatching;

                // regions
                std::shared_ptr<RegionsNet> regionNet;
            };

            // view connection data
            struct ViewConnectionData {
                cv::detail::MatchesInfo matchInfo;
            };



            //// components and constraints
            // component data
            struct LineStructureComponentData {
                double eta;
                deriv::Expression<const double &> etaExpr;
                int mergedSpatialLineSegmentId;
            };

            struct RegionComponentData {
                Vec3 theta;
                deriv::Expression<Eigen::Vector3d> thetaExpr;
                DisableableExpression<double> manhattanEnergyExpr;
                ReconstructionEngine::ViewHandle viewHandle;
                RegionsNet::RegionHandle regionHandle;
            };

            struct ComponentData {
                enum class Type {
                    UnInitialized,
                    LineStructure,
                    Region
                };
                explicit ComponentData(Type t = Type::UnInitialized);
                Type type;
                LineStructureComponentData asLineStructure;
                RegionComponentData asRegion;
            };

            // constraint data
            struct RegionOverlapConstraintData {
                double overlapRatio;
            };

            struct RegionConnectivityConstraintData {
                ReconstructionEngine::ViewHandle viewHandle;
                RegionsNet::BoundaryHandle boundaryHandle;
            };

            struct LineStructureConnectivityConstraintData {
                LineStructureConnectivityConstraintData();

                size_t mergedSpatialLineSegmentIds[2]; // corresponded mergedSpatialLineSegments ids
                PositionOnLine3 positionOnLines[2];
                Vec3 position; // location of intersecion

                // [i][0] -> line lengths with class i lying between vp[i] and position
                // [i][1] -> line lengths with class i lying between position and anti-vp[i]
                double lineVotings[3][2];
                double weight;
                struct { double I, L, X, T, Triplet; } junctionWeights;
                enum { Intersection, Incidence } type;
                double slackValue; // retreived after optimization
            };

            struct RegionLineStructureConnectivityConstraintData {
                std::vector<Vec3> sampledPoints;
            };

            struct ConstraintData {
                enum class Type {
                    UnInitialized,
                    RegionOverlap,
                    RegionConnectivity,
                    LineStructureConnectivity,
                    RegionLineStructureConnectivity
                };
                explicit ConstraintData(Type t = Type::UnInitialized);

                Type type;
                RegionOverlapConstraintData asRegionOverlap;
                RegionConnectivityConstraintData asRegionPairConsistency;
                LineStructureConnectivityConstraintData asLineStructureConnectivity;
                RegionLineStructureConnectivityConstraintData asRegionLineStructureConnectivity;

                DisableableExpression<double> constraintEnergyExpr;
            };




            // global data
            struct GlobalData {
                Image panorama;

                std::array<Vec3, 3> vanishingPoints;
                //std::vector<Image> geometricContext;
                //std::vector<Image> manhattanJunctionDistribution;
                std::vector<Classified<Line3>> spatialLineSegments;
                std::vector<Vec3> mergedSpatialLineSegmentIntersections;
                std::vector<Classified<Line3>> mergedSpatialLineSegments;
                std::vector<int> mergedSpatialLineSegmentChainIds;

                //std::vector<std::pair<int, int>> constraintsOnMergedSpatialLineSegments;

                std::map<int, std::vector<int>> spatialStructuresOfMergedSpatialLineIds;
                std::vector<Classified<Line3>> mergedSpatialLineSegmentsClassifiedWithStructureIds;
            };

            inline const ViewsGraph & views() const { return _views; }
            inline const ConstraintGraph & constraints() const { return _constraints; }
            inline const GlobalData & globalData() const { return _globalData; }

        private:
            ViewsGraph _views;
            ConstraintGraph _constraints;
            deriv::ExpressionGraph _exprGraph;
            Params _params;
            GlobalData _globalData;
        };

 
    }
}
 
#endif