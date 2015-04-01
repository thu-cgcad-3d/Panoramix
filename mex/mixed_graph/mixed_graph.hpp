#ifndef PANORAMIX_MEX_MIXED_GRAPH_HPP
#define PANORAMIX_MEX_MIXED_GRAPH_HPP 

#include "../../src/core/homo_graph.hpp"
#include "../../src/core/rl_graph.hpp"

namespace panoramix {
    namespace experimental {

        using namespace panoramix;
        using namespace panoramix::core;

        using Vec7 = Vec<double, 7>;

        template <class T> 
        struct Pack { 
            T h; 
            Pack() {}
            Pack(T hh) : h(hh) {}
            operator T() const { return h; }
        };

        struct JunctionData {
            Vec3 center;
        };
        using BndGraph = HomogeneousGraph0x<Pack<RegionBoundaryHandle>, JunctionData>;



        struct MixedGraph {
            View<PerspectiveCamera> view;
            std::vector<Classified<core::Line2>> lines;
            std::vector<Vec3> vps;
            Imagei regions;
            std::vector<RegionHandle> regionIds2Handles;

            RLGraph regionLineGraph;
            BndGraph boundaryGraph;

            HandledTable<RegionHandle, Vec7> gcResponse;
            HandledTable<RegionBoundaryHandle, bool> occlusionResponse;

            MixedGraph() {}
            MixedGraph(const Image & image, 
                const HandledTable<RegionHandle, Vec7> & gcResponse,
                const HandledTable<RegionBoundaryHandle, bool> & occlusionResponse,
                const Point2 & cameraCenterPosition, 
                double cameraFocal);

            void showSegmentations() const;
            void showVanishingPointsAndLines() const;
            void showRegionOrientationConstriants() const;

            HandledTable<RegionHandle, Plane3> solve() const;
        };


        

    }
}


#endif