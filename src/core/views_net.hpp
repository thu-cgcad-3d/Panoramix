#ifndef PANORAMIX_CORE_VIEWS_NET_HPP
#define PANORAMIX_CORE_VIEWS_NET_HPP

#include "basic_types.hpp"
#include "feature.hpp"
#include "mesh.hpp"

namespace panoramix {
    namespace core { 

        // views net
        class ViewsNet {
        public:
            struct VertData;
            struct HalfData;
            struct Params {
                Params() : camera(250.0),
                      lineSegmentWeight(1.0), siftWeight(1.0), 
                      surfWeight(1.0), cameraAngleScaler(1.8){}
                PanoramicCamera camera; // camera for generating the panorama
                double lineSegmentWeight;
                double siftWeight;
                double surfWeight;
                LineSegmentExtractor lineSegmentExtractor;
                CVFeatureExtractor<cv::SIFT> siftExtractor;
                CVFeatureExtractor<cv::SURF> surfExtractor;
                double cameraAngleScaler;
            };

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
            void computeTransformationOnConnections(VertHandle h);
            void calibrateCamera(VertHandle h); 
            void calibrateAllCameras();
            void computeGlobalFeatures();

        public:
            struct VertData {
                PerspectiveCamera originalCamera, camera;
                Image image;
                double weight;
                LineSegmentExtractor::Feature featureLineSegment;
                std::vector<core::HPoint2> featureLineIntersections;
                std::vector<std::pair<int, int>> featureLineIntersectionLineIDs;
                CVFeatureExtractor<cv::SIFT>::Feature featureSIFT;
                CVFeatureExtractor<cv::SURF>::Feature featureSURF;
            };
            struct HalfData {
                double cameraAngleDistance;
                double weight;
                Eigen::Matrix<double, 4, 4, Eigen::DontAlign> from2ToTransformation;
            };
            struct GlobalData {
                Image panorama;
                std::array<Vec3, 3> vanishingPoints;
                std::vector<Image> geometricContext;
                std::vector<Image> manhattanJunctionDistribution;
                std::vector<Line3> spatialLineSegments;
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