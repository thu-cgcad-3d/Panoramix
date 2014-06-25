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

                //RTreeWrapper<HPoint2> lineSegmentIntersectionsRTree;
                //RTreeWrapper<KeyPoint> keypointsForMatchingRTree;
            };

            // view connection data
            struct ViewConnectionData {
                cv::detail::MatchesInfo matchInfo;
            };



            // component data
            struct LineStructureComponentData {
                double eta;
                deriv::Expression<const double &> etaExpr;
                std::vector<int> mergedSpatialLineSegmentIds;
            };

            struct RegionComponentData {
                Vec3 theta;
                deriv::Expression<Eigen::Vector3d> thetaExpr;
                deriv::Expression<double> manhattanEnergy;
                bool isVoid;
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

            struct LineStructureConnectivityConstraintData {
                LineStructureConnectivityConstraintData();

                size_t mergedSpatialLineSegmentIds[2]; // corresponded mergedSpatialLineSegments ids
                Vec3 position; // location of intersecion

                // [i][0] -> line lengths with class i lying between vp[i] and position
                // [i][1] -> line lengths with class i lying between position and anti-vp[i]
                double lineVotings[3][2];
                double weight;
                struct { double I, L, X, T, Triplet; } junctionWeights;
                enum { Intersection, Incidence } type;
                double slackValue; // retreived after optimization
            };

            struct RegionPairConsistencyConstraintData {
                ReconstructionEngine::ViewHandle viewHandle;
                RegionsNet::BoundaryHandle boundaryHandle;
                bool isOccludingBoundary;
            };

            struct ConstraintData {
                enum class Type {
                    UnInitialized,
                    RegionOverlap,
                    LineStructureConnectivity,
                    RegionPairConsistency
                };
                explicit ConstraintData(Type t = Type::UnInitialized);

                Type type;
                RegionOverlapConstraintData asRegionOverlap;
                LineStructureConnectivityConstraintData asLineStructureConnectivity;
                RegionPairConsistencyConstraintData asRegionPairConsistency;

                deriv::Expression<double> constraintEnergy;
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
            Params _params;
            GlobalData _globalData;
        };

 
    }
}
 
#endif