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
                Params() : camera(250.0), lineSegmentWeight(1.0), siftWeight(1.0), 
                      surfWeight(1.0), cameraAngleScaler(1.8), smallCameraAngleScalar(0.05), 
                      linePieceSpanAngle(M_PI / 32.0) {}
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
                double linePieceSpanAngle;
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
                SegmentationExtractor::Feature segmentedRegions;

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
                std::vector<Classified<Line3>> mergedSpatialLineSegments;
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