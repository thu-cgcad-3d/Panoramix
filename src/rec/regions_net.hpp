#ifndef PANORAMIX_CORE_REGIONS_NET_HPP
#define PANORAMIX_CORE_REGIONS_NET_HPP

#include "../core/mesh.hpp"
#include "../core/utilities.hpp"
#include "../core/feature.hpp"

namespace panoramix {
    namespace rec {

        using namespace core;

        // net of segmented image regions
        class RegionsNet {
        public:
            struct Params {
                SegmentationExtractor segmenter;
            };
            struct VertData {
                Image regionMask; // 8UC1
                Vec2 center;
                double area;
                double borderLength;
                Box2 boundingBox;
                int lineClassScore;
            };
            struct HalfData {
                double boundaryLength;
                
            };
            using RegionMesh = Mesh<VertData, HalfData>;
            using VertHandle = RegionMesh::VertHandle;
            using HalfHandle = RegionMesh::HalfHandle;

        public:
            explicit RegionsNet(const Image & image, const Params & params = Params());
            void buildNetAndComputeGeometricFeatures(const std::vector<Classified<Line2>> & classifiedLines = std::vector<Classified<Line2>>(),
                const Size & imageSizeContainingLines = Size(0, 0));
            void computeImageFeatures();

            inline const RegionMesh & regions() const { return _regions; }
            inline const Image & image() const { return _image; }
            inline const SegmentationExtractor::Feature & segmentedRegions() const { return _segmentedRegions; }

        private:
            Image _image;
            SegmentationExtractor::Feature _segmentedRegions; // 32SC1
            RegionMesh _regions;
            Params _params;

            struct VertDataBoundingBoxFunctor {
                inline explicit VertDataBoundingBoxFunctor(const RegionMesh & r) : regions(r) {}
                inline Box2 operator()(VertHandle vh) const { return regions.data(vh).boundingBox; }
                const RegionMesh & regions;
            };
            RTreeWrapper<VertHandle, VertDataBoundingBoxFunctor> _regionsRTree;
        };

    }
}
 
#endif