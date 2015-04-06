#ifndef PANORAMIX_MEX_MIXED_GRAPH_HPP
#define PANORAMIX_MEX_MIXED_GRAPH_HPP 

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

            std::vector<Imagef> depthCandidates;

            MixedGraph() {}
            MixedGraph(const Image & image,
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
        };


        

    }
}


#endif