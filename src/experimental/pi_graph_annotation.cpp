#include "../core/containers.hpp"
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
            const std::vector<AnnotedPolygon> & polygons, const std::vector<AnnotedOcclusion> & occs,
            double polygonBoundarySampleStepAngle, double occChainSampleStepAngle, double occChainToBndPieceAngleThres) {
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
                std::vector<Vec3> contour;
                for (int k = 0; k < chain.size(); k++) {
                    auto & p1 = chain.at(k);
                    auto & p2 = chain.next(k);
                    double span = AngleBetweenDirections(p1, p2);

                    for (double a = 0.0; a < span; a += polygonBoundarySampleStepAngle) {
                        auto dir = RotateDirection(p1, p2, a);
                        assert(IsFuzzyZero(AngleBetweenDirections(p1, dir) - a, 1e-2));
                        double angle = AngleBetweenDirections(dir, center);
                        if (maxAngleRadius < angle) {
                            maxAngleRadius = angle;
                        }
                        int ss = contour.size();
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
                
                for (int k = 0; k < contour.size(); k++) {
                    auto & dir = contour[k];
                    auto p = cam.toScreen(dir);
                    proj.push_back(p);
                }
                cv::fillPoly(mask, std::vector<std::vector<Point2i>>{proj}, (uint8_t)1);

                auto segProj = segView.sampled(cam).image;
                auto imProj = mg.view.sampled(cam).image;
                for (auto it = mask.begin(); it != mask.end(); ++it) {
                    if (*it) {
                        int seg = segProj(it.pos());
                        seg2polygonOccupationArea[seg][i] += 1.0 / mg.seg2area[seg] / mg.fullArea; // TODO there may be bugs here
                    }
                }
            }

            for (int seg = 0; seg < mg.nsegs; seg++) {
                auto & polyAreas = seg2polygonOccupationArea[seg];
                double maxAreaRatio = 0.4;
                // check whether their is any clutter mask polygons
                bool hasClutterMask = false;
                for (auto & p : polyAreas) {
                    if (p.second > maxAreaRatio && !polygons[p.first].control.used) {
                        hasClutterMask = true;
                        break;
                    }
                }
                if (hasClutterMask) {
                    mg.seg2control[seg].orientationClaz = mg.seg2control[seg].orientationNotClaz = -1;
                    mg.seg2control[seg].used = false;
                    continue;
                }
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
            RTreeMap<Vec3, int> bndPiecesRTree;
            for (int bp = 0; bp < mg.nbndPieces(); bp++) {
                for (auto & d : mg.bndPiece2dirs[bp]) {
                    bndPiecesRTree.emplace(normalize(d), bp);
                }
            }
            for (int i = 0; i < occs.size(); i++) {
                auto & chain = occs[i].chain;

                // sample on each chain               
                for (int j = 1; j < chain.size(); j++) {
                    auto & p1 = chain[j - 1];
                    auto & p2 = chain[j];
                    std::map<int, int> bndPiecesMet;
                    for (double a = 0.0; a < AngleBetweenDirections(p1, p2); a += occChainSampleStepAngle) {
                        auto dir = normalize(RotateDirection(p1, p2, a));
                        // find nearest bnd piece
                        int nearestbp = -1;
                        double minAngleDist = occChainToBndPieceAngleThres;
                        bndPiecesRTree.search(BoundingBox(dir).expand(occChainToBndPieceAngleThres * 3),
                            [&dir, &minAngleDist, &nearestbp](const std::pair<Vec3, int> & bndPieceD) {
                            double angleDist = AngleBetweenDirections(bndPieceD.first, dir);
                            if (angleDist < minAngleDist) {
                                minAngleDist = angleDist;
                                nearestbp = bndPieceD.second;
                            }
                            return true;
                        });
                        if (nearestbp != -1) {
                            bndPiecesMet[nearestbp] ++;
                        }
                    }
                    for (auto & bpp : bndPiecesMet) {                        
                        int bndPiece = bpp.first;
                        auto & bndPieceDirs = mg.bndPiece2dirs[bndPiece];
                        assert(bndPieceDirs.size() > 1);
                        Vec3 bndRotNormal = bndPieceDirs.front().cross(bndPieceDirs.back());
                        Vec3 edgeRotNormal = normalize(p1.cross(p2));
                        bool sameDirection = bndRotNormal.dot(edgeRotNormal) > 0;
                        if (sameDirection) {
                            mg.bndPiece2occlusion[bndPiece] = OcclusionRelation::LeftIsFront;
                        } else {
                            mg.bndPiece2occlusion[bndPiece] = OcclusionRelation::RightIsFront;
                        }
                    }
                }

            }


        }


    }
}