#ifndef PANORAMIX_CORE_VIEWS_NET_HPP
#define PANORAMIX_CORE_VIEWS_NET_HPP

#include "basic_types.hpp"
#include "feature.hpp"
#include "regions_net.hpp"

namespace panoramix {
    namespace core { 

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
                double cameraAngleScaler; // angle scalar to judge whether two views may share certain common features
                double smallCameraAngleScalar; // angle scalar to judge whether two views are too close
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

            struct VertData;
            struct HalfData;
            using ViewMesh = Mesh<VertData, HalfData>;
            using VertHandle = ViewMesh::VertHandle;
            using HalfHandle = ViewMesh::HalfHandle;

        public:
            inline explicit ViewsNet(const Params params = Params()) : _params(params){}
            inline const Params & params() const { return _params; }
            
            inline VertHandle insertVertex(const VertData & vd) { return _views.addVertex(vd); }

            // insert a new photo, with known parameters
            VertHandle insertPhoto(const Image & im, const PerspectiveCamera & cam, 
                double cameraDirectionErrorScale = 0.0,
                double cameraPositionErrorScale = 0.0);

            // compute features for a single view
            void computeFeatures(VertHandle h);

            // segment view image and build net of regions for a single view
            void buildRegionNet(VertHandle h);

            // connect this view to neighbor views
            int updateConnections(VertHandle h);

            // whether this view overlaps some existing views a lot, measured by the smallCameraAngleScalar parameter
            VertHandle isTooCloseToAnyExistingView(VertHandle h) const;

            // calibrate camera parameters of this view, and if an error loop is closed, calibrate all related view cameras
            void calibrateCamera(VertHandle h); 

            // calibrate all cameras
            void calibrateAllCameras();
            
            // estimate vanishing points using lines extract from all views, classify this lines and lift them all to space
            void estimateVanishingPointsAndClassifyLines();

            // build constraints on spatial lines and rectify their parameters to build a more reasonable 3D sketch
            void rectifySpatialLines();
                
        public:
            struct VertData {
                PerspectiveCamera originalCamera, camera;
                double cameraDirectionErrorScale, cameraPositionErrorScale;
                Image image;
                double weight;
                std::vector<Classified<Line2>> lineSegments;
                std::vector<core::HPoint2> lineSegmentIntersections;
                std::vector<std::pair<int, int>> lineSegmentIntersectionLineIDs;
                CVFeatureExtractor<cv::SIFT>::Feature SIFTs;
                CVFeatureExtractor<cv::SURF>::Feature SURFs;
                std::shared_ptr<RegionsNet> regionNet;
            };

            struct HalfData {
                double cameraAngleDistance;
                double weight;
                // transform //
            };

            struct ConstraintData {
                ConstraintData();
                size_t mergedSpatialLineSegmentIds[2]; // corresponded mergedSpatialLineSegments ids
                Vec3 position; // location of intersecion
                
                // [i][0] -> line lengths with class i lying between vp[i] and position
                // [i][1] -> line lengths with class i lying between position and anti-vp[i]
                double lineVotings[3][2];
                double weight;
                struct JunctionWeights { double I, L, X, T, Triplet; } junctionWeights;
                enum { Intersection,  Incidence } type;
                double slackValue; // retreived after optimization
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
                // constraints for 3D lines reconstruction
                std::vector<ConstraintData> constraints, refinedConstraints;
            };

            inline const ViewMesh & views() const { return _views; }
            inline const GlobalData & globalData() const { return _globalData; }

        private:
            ViewMesh _views;
            Params _params;
            GlobalData _globalData;
        };

 
    }
}
 
#endif