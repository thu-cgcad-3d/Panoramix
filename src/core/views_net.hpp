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
                double connectedLinesDistanceAngleThreshold;
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
            VertHandle insertPhoto(const Image & im, const PerspectiveCamera & cam);

            void computeFeatures(VertHandle h);
            int updateConnections(VertHandle h);
            VertHandle isTooCloseToAnyExistingView(VertHandle h) const;
            void computeTransformationOnConnections(VertHandle h);
            void calibrateCamera(VertHandle h); 
            void calibrateAllCameras();
            void updateExternalRegionConnections(VertHandle h); // build region connections across views
            void estimateVanishingPointsAndClassifyLines();
            void rectifySpatialLines(); // 
            //void estimateLayoutStructure();

        public:
            struct VertData {
                PerspectiveCamera originalCamera, camera;
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
            struct GlobalData {
                Image panorama;
                std::array<Vec3, 3> vanishingPoints;
                std::vector<Image> geometricContext;
                std::vector<Image> manhattanJunctionDistribution;
                std::vector<Classified<Line3>> spatialLineSegments;
                std::vector<Vec3> mergedSpatialLineSegmentIntersections;
                std::vector<Classified<Line3>> mergedSpatialLineSegments;
                std::vector<int> mergedSpatialLineSegmentChainIds;
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