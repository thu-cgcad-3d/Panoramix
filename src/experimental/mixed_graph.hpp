#ifndef PANORAMIX_MEX_MIXED_GRAPH_HPP
#define PANORAMIX_MEX_MIXED_GRAPH_HPP 

#include "../core/homo_graph.hpp"
#include "../core/rl_graph.hpp"

namespace panoramix {
    namespace experimental {

        using namespace panoramix::core;

        using Vec7 = Vec<double, 7>;
        using Imaged7 = ImageOfType<Vec7>;


        struct BoundaryJunction {
            std::set<RegionBoundaryHandle> bhs;
            std::vector<Vec3> positions;
        };

        struct MixedGraph {
            View<PerspectiveCamera> view;
            std::vector<Classified<core::Line2>> lines;
            std::vector<Vec3> vps;
            Imagei regions;
            std::vector<RegionHandle> regionIds2Handles;
            std::vector<BoundaryJunction> boundaryJunctions;

            RLGraph regionLineGraph;

            HandledTable<RegionHandle, Vec7> gcResponse;
            HandledTable<RegionBoundaryHandle, double> occlusionResponse;

            MixedGraph() {}
            MixedGraph(const Image & image,
                const Point2 & cameraCenterPosition, 
                double cameraFocal);

            void installGCResponse(const Imaged7 & gc);
            void installOcclusionResponce(const std::vector<std::vector<PixelLoc>> & edges, 
                const std::vector<double> & scores);
            void installOcclusionResponce(const std::vector<std::vector<int>> & edges,
                const std::vector<double> & scores);
            
            void showOcclusionResponse() const;

            void showSegmentations() const;
            void showVanishingPointsAndLines() const;
            void showRegionOrientationConstriants() const;

            HandledTable<RegionHandle, Plane3> solve() const;
        };


        

    }
}


#endif