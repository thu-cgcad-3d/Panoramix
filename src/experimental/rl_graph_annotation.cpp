#include <QFileInfo>

#include "../core/single_view.hpp"

#include "rl_graph_annotation.hpp"

namespace pano {

    namespace experimental {


        std::string AnnotationFilePath(const std::string & imageFilePath) {
            QFileInfo finfo(QString::fromStdString(imageFilePath));
            if (!finfo.exists())
                return "";
            auto annoFileName = finfo.absoluteFilePath() + ".anno.cereal";
            return annoFileName.toStdString();
        }

        void ProjectOn(const Polygon3 & polygon, const PanoramicCamera & cam, Imagei & canvas, int val) {
            Vec3 center;
            for (auto & p : polygon.corners) {
                center += p;
            }
            center /= norm(center);
            double radiusAngle = 0.0;
            for (auto & c : polygon.corners) {
                radiusAngle = std::max(radiusAngle, AngleBetweenDirections(c, center));
            }
            auto view = core::PerfectRegionMaskView(polygon.corners, center, 50);
            auto panoView = view.sampled(cam);
            if (canvas.empty()) {
                canvas = Imagei(cam.screenSize(), -1);
            }
            for (auto it = canvas.begin(); it != canvas.end(); ++it) {
                if (panoView.image(it.pos())) {
                    *it = val;
                }
            }
        }



        void AttachAnnotationConstraints(const RLGraph & mg, 
            const PanoramicCamera & cam, 
            const SegmentationTopo & segtopo,
            const std::vector<RegionHandle> & rhs,
            const std::vector<RegionBoundaryHandle> & bhs,
            RLGraphControls & controls, 
            const std::vector<OrientedPolygon> & polygons, 
            const std::vector<Chain3> & occlusions) {
            
            // set region controls
            Imagei polyIds;
            for (int i = 0; i < polygons.size(); i++) {
                ProjectOn(polygons[i].polygon, cam, polyIds, i);
            }

            auto labelRatiosTable = mg.createComponentTable<RegionData, std::map<std::tuple<int, int, bool>, double>>();
            for (auto & r : mg.components<RegionData>()) {
                auto rh = r.topo.hd;
                auto regionMaskView = PerfectRegionMaskView(mg, rh);
                auto sampler = MakeCameraSampler(regionMaskView.camera, cam);
                auto polyIdsOnRegion = sampler(polyIds);
                int regionSize = 0;
                for (auto it = regionMaskView.image.begin(); it != regionMaskView.image.end(); ++it) {
                    if (!*it) {
                        continue;
                    }
                    int pid = polyIdsOnRegion(it.pos());
                    if (pid >= 0) {
                        auto & poly = polygons[pid];
                        auto label = std::make_tuple(poly.towardVPId, poly.alongVPId, poly.isClutter);
                        labelRatiosTable[rh][label] ++;
                    }
                }
                regionSize++;
                for (auto & r : labelRatiosTable[rh]) {
                    r.second /= regionSize;
                }

                auto & labelRatios = labelRatiosTable[rh];
                if (labelRatios.empty()) {
                    continue;
                }
                
                double maxRatio = 0;
                int toward = -1, along = -1, clutter = false;
                for (auto & r : labelRatios) {
                    if (r.second > maxRatio) {
                        maxRatio = r.second;
                        toward = std::get<0>(r.first);
                        along = std::get<1>(r.first);
                        clutter = std::get<2>(r.first);
                    }
                }
                controls[rh].orientationClaz = toward;
                controls[rh].orientationNotClaz = along;
                controls[rh].used = !clutter;
            }


            // set boundary controls


        }



        QuantEvaluation Evaluate(const LayeredShape3 & shape, const PanoIndoorAnnotation & anno) {
            NOT_IMPLEMENTED_YET();
        }

    }

}