#ifndef PANORAMIX_CORE_VIEWS_NET_HPP
#define PANORAMIX_CORE_VIEWS_NET_HPP

#include "../core/basic_types.hpp"
#include "../core/feature.hpp"
#include "../core/utilities.hpp"
#include "regions_net.hpp"

namespace panoramix {
    namespace rec { 

        using namespace core;

        // views net
        class ViewsNet {
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

        public:
            inline explicit ViewsNet(const Params params = Params()) : _params(params){}
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

            struct ViewConnectionData {
                cv::detail::MatchesInfo matchInfo;
            };

            struct GlobalData {
                Image panorama;

                std::array<Vec3, 3> vanishingPoints;
                std::vector<Image> geometricContext;
                std::vector<Image> manhattanJunctionDistribution;
                std::vector<Classified<Line3>> spatialLineSegments;
                std::vector<Vec3> mergedSpatialLineSegmentIntersections;
                std::vector<Classified<Line3>> mergedSpatialLineSegments;
                std::vector<int> mergedSpatialLineSegmentChainIds;

                //std::vector<std::pair<int, int>> constraintsOnMergedSpatialLineSegments;

                std::map<int, std::vector<int>> spatialStructuresOfMergedSpatialLineIds;
                std::vector<Classified<Line3>> mergedSpatialLineSegmentsClassifiedWithStructureIds;
            };

            inline const ViewsGraph & views() const { return _views; }
            inline const GlobalData & globalData() const { return _globalData; }

        private:
            ViewsGraph _views;
            Params _params;
            GlobalData _globalData;
        };

 
    }
}
 
#endif