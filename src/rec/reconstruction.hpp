#ifndef PANORAMIX_REC_RECONSTRUCTION_HPP
#define PANORAMIX_REC_RECONSTRUCTION_HPP

#include "../core/view.hpp"
 
namespace panoramix {
    namespace rec {

        using namespace core;

        struct Context {
            View<PanoramicCamera> originalView;
            
            std::vector<View<PerspectiveCamera>> views;
            std::vector<RegionsGraph> regionsGraphs;
            std::vector<LinesGraph> & linesGraphs;
            std::vector<Vec3> vanishingPoints;

            ComponentIndexHashMap<std::pair<RegionIndex, RegionIndex>, double> regionOverlappings;
            ComponentIndexHashMap<std::pair<RegionIndex, LineIndex>, std::vector<Vec3>> regionLineConnections;
            ComponentIndexHashMap<std::pair<LineIndex, LineIndex>, Vec3> lineIndencesAcrossViews;

            int regionConnectedComponentsNum;
            ComponentIndexHashMap<RegionIndex, int> regionConnectedComponentIds;
            int lineConnectedComponentsNum;
            ComponentIndexHashMap<LineIndex, int> lineConnectedComponentIds;
        };


        void EstimateSpatialLineDepths(const Context & context,
            ComponentIndexHashMap<LineIndex, Line3> & reconstructedLines,
            double constantEtaForFirstLineInEachConnectedComponent,
            bool twiceEstimation);


        void EstimateSpatialRegionPlanes(const std::vector<View<PerspectiveCamera>> & views,
            const std::vector<RegionsGraph> & regionsNets, const std::vector<LinesGraph> & linesNets,
            const std::array<Vec3, 3> & vanishingPoints,
            const ComponentIndexHashMap<std::pair<RegionIndex, RegionIndex>, double> & regionOverlappings,
            const ComponentIndexHashMap<std::pair<RegionIndex, LineIndex>, std::vector<Vec3>> & regionLineConnections,
            const ComponentIndexHashMap<std::pair<LineIndex, LineIndex>, Vec3> & interViewLineIncidences,
            int regionConnectedComponentsNum, const ComponentIndexHashMap<RegionIndex, int> & regionConnectedComponentIds,
            int lineConnectedComponentsNum, const ComponentIndexHashMap<LineIndex, int> & lineConnectedComponentIds,
            ComponentIndexHashMap<LineIndex, Line3> & reconstructedLines,
            ComponentIndexHashMap<RegionIndex, Plane3> & reconstructedPlanes,
            const Image & globalTexture);


    }
}

#endif