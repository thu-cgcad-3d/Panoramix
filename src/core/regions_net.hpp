#ifndef PANORAMIX_CORE_REGIONS_NET_HPP
#define PANORAMIX_CORE_REGIONS_NET_HPP

#include "../deriv/mesh.hpp"
#include "feature.hpp"

namespace panoramix {
    namespace core {

        using deriv::Mesh;

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
            };
            struct HalfData {
                double boundaryLength;
                
            };
            using RegionMesh = Mesh<VertData, HalfData>;
            using VertHandle = RegionMesh::VertHandle;
            using HalfHandle = RegionMesh::HalfHandle;

        public:
            explicit RegionsNet(const Image & image, const Params & params = Params());
            void buildNetAndComputeGeometricFeatures();
            void computeImageFeatures();

            inline const RegionMesh & regions() const { return _regions; }
            inline const Image & image() const { return _image; }
            inline const SegmentationExtractor::Feature & segmentedRegions() const { return _segmentedRegions; }

        private:
            Image _image;
            SegmentationExtractor::Feature _segmentedRegions; // 32SC1
            RegionMesh _regions;
            Params _params;
        };

    }
}
 
#endif