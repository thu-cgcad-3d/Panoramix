#include <iostream>
#include <string>
#include <random>

#include "../src/rec/reconstruction.hpp"
#include "test_config.hpp"

using namespace panoramix;


inline std::string ComposeCacheFileName(const std::string & inputFile, const std::string & tag) {
    std::string name = inputFile;
    for (auto & c : name){
        if (c == '/' || c == '\\')
            c = '.';
        if (c == ':')
            c = '@';
    }
    std::string d = "";
    return "./cache/" + name + '.' + tag + d + ".state";
}


struct All {

	rec::View<core::PanoramicCamera> originalView;
	std::vector<rec::View<core::PerspectiveCamera>> perspectiveViews;
	std::vector<rec::RegionsNet> regionsNets;
    std::vector<rec::LinesNet> linesNets;

    std::array<core::Vec3, 3> vanishingPoints;

    rec::ComponentIndexHashMap<std::pair<rec::RegionIndex, rec::RegionIndex>, double> regionOverlappings;
    rec::ComponentIndexHashMap<std::pair<rec::RegionIndex, rec::LineIndex>, std::vector<rec::Vec3>> regionLineConnections;
    rec::ComponentIndexHashMap<std::pair<rec::LineIndex, rec::LineIndex>, rec::Vec3> interViewLineIncidences;
    int regionConnectedComponentsNum;
    rec::ComponentIndexHashMap<rec::RegionIndex, int> regionConnectedComponentIds;
    int lineConnectedComponentsNum;
    rec::ComponentIndexHashMap<rec::LineIndex, int> lineConnectedComponentIds;

    rec::ComponentIndexHashMap<rec::LineIndex, core::Line3> reconstructedLines;
    rec::ComponentIndexHashMap<rec::RegionIndex, core::Plane3> reconstructedPlanes;

    template <class Archiver>
    inline void serialize(Archiver & ar) {
        ar(originalView, perspectiveViews, regionsNets, linesNets,
            vanishingPoints,
            regionOverlappings, regionLineConnections, interViewLineIncidences, 
            regionConnectedComponentsNum, regionConnectedComponentIds,
            lineConnectedComponentsNum, lineConnectedComponentIds,
            reconstructedLines, reconstructedPlanes);
    }

};


int main(int argc, char * argv[], char * envp[]) {

    std::string originalFile = test::ProjectDataDirStrings::PanoramaIndoor + "/13.jpg";
    std::string cacheFileBeforeComputingFeatures = ComposeCacheFileName(originalFile, "1_before_fea");
    std::string cacheFileAfterComputingFeatures = ComposeCacheFileName(originalFile, "2_after_fea");
    std::string cacheFileAfterEstimatingVPs = ComposeCacheFileName(originalFile, "3_after_vp");
    std::string cacheFileAfterEstimatingLineDepths = ComposeCacheFileName(originalFile, "4_after_linedepths");

    All all;

    core::UpdateIfFileIsTooOld(originalFile, cacheFileBeforeComputingFeatures,
        [&](const std::string & in, const std::string & out) {
        core::Image im = cv::imread(in);
        core::ResizeToMakeWidthUnder(im, 2000);
        all.originalView = rec::CreatePanoramicView(im);
        all.perspectiveViews = rec::PerspectiveSampling(all.originalView, {
            core::PerspectiveCamera(700, 700, 350, { 0, 0, 0 }, { 1, 0, 0 }, { 0, 0, -1 }),
            core::PerspectiveCamera(700, 700, 350, { 0, 0, 0 }, { 0, 1, 0 }, { 0, 0, -1 }),
            core::PerspectiveCamera(700, 700, 350, { 0, 0, 0 }, { -1, 0, 0 }, { 0, 0, -1 }),
            core::PerspectiveCamera(700, 700, 350, { 0, 0, 0 }, { 0, -1, 0 }, { 0, 0, -1 }),
            core::PerspectiveCamera(700, 700, 350, { 0, 0, 0 }, { 0, 0, 1 }, { 1, 0, 0 }),
            core::PerspectiveCamera(700, 700, 350, { 0, 0, 0 }, { 0, 0, -1 }, { 1, 0, 0 })
        });
        core::SaveToDisk(out, all);
    });

    core::UpdateIfFileIsTooOld(cacheFileBeforeComputingFeatures, cacheFileAfterComputingFeatures,
        [&](const std::string & in, const std::string & out){
        core::LoadFromDisk(in, all);
        all.regionsNets.resize(all.perspectiveViews.size());
        all.linesNets.resize(all.perspectiveViews.size());
        #pragma omp for
        for (int i = 0; i < all.perspectiveViews.size(); i++){
            std::tie(all.regionsNets[i], all.linesNets[i]) = 
                rec::InitializeFeatureNets(all.perspectiveViews[i], 5.0, 20.0, 15.0, 40.0);
        }
        core::SaveToDisk(out, all);
    });

    core::UpdateIfFileIsTooOld(cacheFileAfterComputingFeatures, cacheFileAfterEstimatingVPs,
        [&](const std::string & in, const std::string & out) {
        core::LoadFromDisk(in, all);
        all.vanishingPoints = 
            rec::EstimateVanishingPointsAndClassifyLines(all.perspectiveViews, all.linesNets);
        rec::RecognizeRegionLineConstraints(all.perspectiveViews, all.regionsNets, all.linesNets,
            all.regionOverlappings, all.regionLineConnections, all.interViewLineIncidences,
            M_PI_4 / 8.0, 10.0);
        rec::ComputeConnectedComponentsUsingRegionLineConstraints(all.perspectiveViews, all.regionsNets, all.linesNets,
            all.regionOverlappings, all.regionLineConnections, all.interViewLineIncidences,
            all.regionConnectedComponentsNum, all.regionConnectedComponentIds,
            all.lineConnectedComponentsNum, all.lineConnectedComponentIds);
        core::SaveToDisk(out, all);
    });

    core::UpdateIfFileIsTooOld(cacheFileAfterEstimatingVPs, cacheFileAfterEstimatingLineDepths,
        [&](const std::string & in, const std::string & out) {
        core::LoadFromDisk(in, all);
        rec::EstimateSpatialLineDepths(all.perspectiveViews, all.linesNets, all.vanishingPoints, all.interViewLineIncidences, 
            all.lineConnectedComponentsNum, all.lineConnectedComponentIds, 
            all.reconstructedLines, 10.0, false);
        core::SaveToDisk(out, all);
    });

    core::UpdateIfFileIsTooOld(cacheFileAfterEstimatingLineDepths, "xxx",
        [&](const std::string & in, const std::string & out) {
        core::LoadFromDisk(in, all);
        rec::EstimateSpatialRegionPlanes(all.perspectiveViews, all.regionsNets, all.linesNets, all.vanishingPoints,
            all.regionOverlappings, all.regionLineConnections, all.interViewLineIncidences,
            all.regionConnectedComponentsNum, all.regionConnectedComponentIds,
            all.lineConnectedComponentsNum, all.lineConnectedComponentIds,
            all.reconstructedLines, all.reconstructedPlanes, all.originalView.image);
        //core::SaveToDisk(out, all);
    });

   
    return 0;
}