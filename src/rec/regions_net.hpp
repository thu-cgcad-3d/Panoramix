#ifndef PANORAMIX_CORE_REGIONS_NET_HPP
#define PANORAMIX_CORE_REGIONS_NET_HPP

#include "../core/graphical_model.hpp"
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
            struct RegionData {
                Image regionMask; // 8UC1
                Vec2 center;
                double area;
                double borderLength;
                Box2 boundingBox;
                int orientationClass; // -1 0 1 2
            };
            struct BoundaryData {
                double boundaryLength;                
            };
            using RegionsGraph = GraphicalModel02<RegionData, BoundaryData>;
            using RegionHandle = HandleAtLevel<0>;
            using BoundaryHandle = HandleAtLevel<1>;

        public:
            explicit RegionsNet(const Image & image, const Params & params = Params());
            void buildNetAndComputeGeometricFeatures(const std::vector<Classified<Line2>> & classifiedLines = std::vector<Classified<Line2>>(),
                const Size & imageSizeContainingLines = Size(0, 0));
            void computeImageFeatures();

            inline const RegionsGraph & regions() const { return _regions; }
            inline const Image & image() const { return _image; }
            inline const SegmentationExtractor::Feature & segmentedRegions() const { return _segmentedRegions; }

        private:
            Image _image;
            SegmentationExtractor::Feature _segmentedRegions; // 32SC1
            RegionsGraph _regions;
            Params _params;

            struct RegionDataBoundingBoxFunctor {
                inline explicit RegionDataBoundingBoxFunctor(const RegionsGraph & r) : regions(r) {}
                inline Box2 operator()(RegionHandle vh) const { return regions.data(vh).boundingBox; }
                const RegionsGraph & regions;
            };
            RTreeWrapper<RegionHandle, RegionDataBoundingBoxFunctor> _regionsRTree;
        };

    }
}
 
#endif