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


        void EstimateSpatialRegionPlanes(const Context & context,
            ComponentIndexHashMap<LineIndex, Line3> & reconstructedLines,
            ComponentIndexHashMap<RegionIndex, Plane3> & reconstructedPlanes);


    }
}

#endif