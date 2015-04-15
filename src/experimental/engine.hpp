#ifndef PANORAMIX_EXPERIMENTAL_ENGINE_HPP
#define PANORAMIX_EXPERIMENTAL_ENGINE_HPP 

#include "../core/homo_graph.hpp"
#include "rl_graph.hpp"

namespace panoramix {
    namespace experimental {

        using namespace core;

        using Vec7 = Vec<double, 7>;
        using Image7d = ImageOfType<Vec7>;


        struct BoundaryJunction {
            std::set<RegionBoundaryHandle> bhs;
            std::vector<Vec3> positions;
        };

        class Engine {
        public:
            Engine() {}
            Engine(const Image & image);
            Engine(const Image & image,
                const Point2 & cameraCenterPosition, 
                double cameraFocal);

            void installImageAndExtractFeatures(const Image & image);
            void installImageAndExtractFeatures(const Image & image,
                const Point2 & cameraCenterPosition,
                double cameraFocal);

            void installGCResponse(const Image7d & gc);
            void installOcclusionResponce(const std::vector<std::vector<PixelLoc>> & edges, 
                const std::vector<double> & scores);
            void installOcclusionResponce(const std::vector<std::vector<int>> & edges,
                const std::vector<double> & scores);
            
            void showGCResponse() const;
            void showOcclusionResponse() const;
            void showBoundaryJunctions() const;
            void showDetachableRegionLineConnections() const;

            void showSegmentations() const;
            void showVanishingPointsAndLines() const;
            void showRegionOrientationConstriants() const;

            HandledTable<RegionHandle, core::Plane3> solve() const;

        private:
            // build all based on view
            void initialize();

        private:
            View<PerspectiveCamera> view;

            std::vector<Classified<core::Line2>> lines;
            std::vector<Vec3> vps;
            Imagei regions;
            std::vector<RegionHandle> rhs;
            std::vector<BoundaryJunction> boundaryJunctions;

            RLGraph mg;

            HandledTable<RegionHandle, Vec7> gcResponse;
            HandledTable<RegionBoundaryHandle, double> occlusionResponse;

            std::vector<Imagef> depthCandidates;

            friend class cereal::access;
            template <class Archiver>
            inline void serialize(Archiver & ar) {
                ar(view, lines, vps, regions, rhs, boundaryJunctions, 
                    mg, gcResponse, occlusionResponse, depthCandidates);
            }
        };


        

    }
}


#endif