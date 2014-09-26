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
                Params();
                SegmentationExtractor segmenter;
                double samplingStepLengthOnBoundary;
                template <class Archiver>
                void serialize(Archiver & ar) {
                    ar(segmenter, samplingStepLengthOnBoundary);
                }
            };
            struct RegionData {
                Image regionMask; // 8UC1
                Vec2 center;
                double area;
                std::vector<std::vector<PixelLoc>> contours;
                std::vector<std::vector<PixelLoc>> dilatedContours;
                Box2 boundingBox;
                template <class Archiver>
                void serialize(Archiver & ar) {
                    ar(regionMask, center, area, contours, dilatedContours, boundingBox);
                }
            };
            struct BoundaryData {
                std::vector<std::vector<PixelLoc>> edges;
                double length;
                InfiniteLine2 fittedLine;
                double tjunctionLikelihood;
                double straightness;
                double interleavedArea;
                double interleavedLength;
                std::vector<std::vector<Point2>> sampledPoints;
                template <class Archiver>
                void serialize(Archiver & ar) {
                    ar(edges, length, fittedLine, 
                        tjunctionLikelihood, straightness, interleavedArea, 
                        interleavedLength, sampledPoints);
                }
            };
            using RegionsGraph = GraphicalModel02<RegionData, BoundaryData>;
            using RegionHandle = HandleAtLevel<0>;
            using BoundaryHandle = HandleAtLevel<1>;

        public:
            inline RegionsNet() {}
            explicit RegionsNet(const Image & image, const Params & params = Params());
            void buildNetAndComputeGeometricFeatures();
            void computeImageFeatures();            

            inline const RegionsGraph & regions() const { return _regions; }
            inline const Image & image() const { return _image; }

            // 32SC1
            inline const SegmentationExtractor::Feature & segmentedRegions() const { return _segmentedRegions; }



        private:
            Image _image;
            SegmentationExtractor::Feature _segmentedRegions; // 32SC1
            RegionsGraph _regions;
            Params _params;

            template <class Archiver>
            void serialize(Archiver & ar) {
                ar(_image, _segmentedRegions, _regions, _params);
            }
            friend class cereal::access;
        };

    }
}
 
#endif