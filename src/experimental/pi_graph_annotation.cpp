#include "../gui/qttools.hpp"
#include "../gui/utility.hpp"
#include "../gui/singleton.hpp"

#include "pi_graph_annotation_widgets.hpp"
#include "pi_graph_occlusion.hpp"
#include "pi_graph_annotation.hpp"

namespace pano {
    namespace experimental {

        std::string AnnotationFilePath(const std::string & imagePath) {
            QFileInfo finfo(QString::fromStdString(imagePath));
            if (!finfo.exists())
                return "";
            auto annoFileName = finfo.absoluteFilePath() + ".anno.cereal";
            return annoFileName.toStdString();
        }


        PIAnnotation LoadOrInitializeNewAnnotation(const std::string & imagePath) {
            
            assert(QFileInfo(QString::fromStdString(imagePath)).exists());

            auto annoPath = AnnotationFilePath(imagePath);
            QFileInfo annofinfo(QString::fromStdString(annoPath));

            PIAnnotation anno;

            // if exist, load it
            if (!annofinfo.exists() || !LoadFromDisk(annofinfo.absoluteFilePath().toStdString(), anno)) {
                
                // initialize new annotation
                anno.originalImage = cv::imread(imagePath);
                
                // rectify the image
                std::cout << "rectifying image" << std::endl;
                anno.rectifiedImage = anno.originalImage.clone();
                gui::MakePanoramaByHand(anno.rectifiedImage, &anno.extendedOnTop, &anno.extendedOnBottom);

                // create view
                std::cout << "creating views" << std::endl;
                Image image = anno.rectifiedImage.clone();
                ResizeToHeight(image, 700);
                anno.view = CreatePanoramicView(image);

                // collect lines
                std::cout << "collecting lines" << std::endl;
                auto cams = CreateCubicFacedCameras(anno.view.camera, image.rows, image.rows, image.rows * 0.4);
                std::vector<Line3> rawLine3s;
                for (int i = 0; i < cams.size(); i++) {
                    auto pim = anno.view.sampled(cams[i]).image;
                    LineSegmentExtractor lineExtractor;
                    lineExtractor.params().algorithm = LineSegmentExtractor::LSD;
                    auto ls = lineExtractor(pim, 3, 300); // use pyramid
                    for (auto & l : ls) {
                        rawLine3s.emplace_back(normalize(cams[i].toSpace(l.first)),
                            normalize(cams[i].toSpace(l.second)));
                    }
                }
                rawLine3s = MergeLines(rawLine3s, DegreesToRadians(1));

                // estimate vp and initially classify lines
                std::cout << "estimating vps and clasifying lines" << std::endl;
                auto line3s = ClassifyEachAs(rawLine3s, -1);
                auto vps = EstimateVanishingPointsAndClassifyLines(line3s);
                OrderVanishingPoints(vps);
                int vertVPId = NearestDirectionId(vps, Vec3(0, 0, 1));

                anno.lines = std::move(line3s);
                anno.vps = std::move(vps);
                anno.vertVPId = vertVPId;
            }

            return anno;
        }


        void EditAnnotation(PIAnnotation & anno) {
            PIAnnotation annoEdited = anno;
            gui::Singleton::InitGui();
            PIAnnotationWidget w;
            w.setCurAnnotation(&annoEdited);
            w.resize(900, 900);
            w.show();
            gui::Singleton::ContinueGui();
            int selected = gui::SelectFrom({ "Accept", "Cancel" }, "Accept the edit?", "Accept the edit?", 0, 1);
            if (selected == 0) {
                anno = std::move(annoEdited);
            }
        }


        void SaveAnnotation(const std::string & imagePath, const PIAnnotation & anno) {
            auto annoPath = AnnotationFilePath(imagePath);
            QFileInfo annofinfo(QString::fromStdString(annoPath));

            SaveToDisk(annofinfo.absoluteFilePath().toStdString(), anno);
        }



        void AttachAnnotatedPolygonsAndOcclusions(PIGraph & mg,
            const std::vector<AnnotedPolygon> & polygons, const std::vector<AnnotedOcclusion> & occs) {
            std::cout << "attaching annotated polygons and occlusions" << std::endl;
            AssumeThereAreNoOcclusions(mg);

            // attach polygons
            View<PanoramicCamera, Imagei> segView(mg.segs, mg.view.camera);
            std::vector<std::map<int, double>> seg2polygonOccupationArea(mg.nsegs);
            for (int i = 0; i < polygons.size(); i++) {
                auto & polygon = polygons[i].polygon;
                Vec3 center = Origin();
                for (auto & p : polygon.corners) {
                    center += normalize(p);
                }
                center /= norm(center);
                double maxAngleRadius = 0.0;
                auto chain = polygon.boundary();
                static const double sampleStepAngle = DegreesToRadians(1);
                std::vector<Vec3> contour;
                for (int k = 0; k < chain.size(); k++) {
                    auto & p1 = chain.at(k);
                    auto & p2 = chain.next(k);
                    for (double a = 0.0; a < AngleBetweenDirections(p1, p2); a += sampleStepAngle) {
                        auto dir = RotateDirection(p1, p2, a);
                        double angle = AngleBetweenDirections(dir, center);
                        if (maxAngleRadius < angle) {
                            maxAngleRadius = angle;
                        }
                        contour.push_back(dir);
                    }
                }

                Vec3 up;
                std::tie(std::ignore, up) = ProposeXYDirectionsFromZDirection(center);
                UniformSphericalCamera cam(mg.view.camera.focal(), maxAngleRadius + 0.01, mg.view.camera.eye(), center, up);
                Imageub mask = Imageub::zeros(cam.screenSize());
                
                // project polygon using this camera
                std::vector<Point2i> proj;
                proj.reserve(contour.size());
             
                for (auto & dir : contour) {
                    proj.push_back(cam.toScreen(dir));
                }
                cv::fillPoly(mask, std::vector<std::vector<Point2i>>{proj}, (uint8_t)1);

                auto segProj = segView.sampled(cam).image;
                for (auto it = mask.begin(); it != mask.end(); ++it) {
                    if (*it) {
                        int seg = segProj(it.pos());
                        seg2polygonOccupationArea[seg][i] += 1.0 / mg.seg2area[seg] / mg.fullArea; // TODO weight error
                    }
                }
            }

            for (int seg = 0; seg < mg.nsegs; seg++) {
                auto & polyAreas = seg2polygonOccupationArea[seg];
                double maxAreaRatio = 0.2;
                int bestPoly = -1;
                for (auto & p : polyAreas) {
                    if (p.second > maxAreaRatio) {
                        bestPoly = p.first;
                        maxAreaRatio = p.second;
                    }
                }
                if (bestPoly != -1) {
                    mg.seg2control[seg] = polygons[bestPoly].control;
                }
            }


            // attach occlusions

        }


    }
}