
#include "eval.hpp"

#include "../../src/experimental/rl_graph_solver.hpp"
#include "../../src/experimental/rl_graph_occlusion.hpp"
#include "../../src/experimental/rl_graph_modeling.hpp"
#include "../../src/experimental/tools.hpp"

#include "../../src/experimental/pi_graph_control.hpp"
#include "../../src/experimental/pi_graph_occlusion.hpp"
//#include "../../src/experimental/pi_graph_solve.hpp"
#include "../../src/experimental/pi_graph_vis.hpp"

namespace panolyz {

    std::string PanoramixResultFilePath(const std::string & impath, int version) {
        return impath + ".panoramix." + std::to_string(version) + ".cereal";
    }

    // pi graph

    struct PIGraphModel : ReconstructedModel {
        PIGraph mg;
        explicit PIGraphModel(PIGraph && g) : mg(std::move(g)) {}

        virtual double depthAt(const Vec3 & direction) const {
            auto p = ToPixel(mg.view.camera.toScreen(direction));
            p.x = WrapBetween(p.x, 0, mg.segs.cols);
            p.y = BoundBetween(p.y, 0, mg.segs.rows - 1);
            int seg = mg.segs(p);
            return norm(IntersectionOfLineAndPlane(Ray3(Origin(), direction), mg.seg2recPlanes[seg]).position);
        }
        virtual void visualize(const std::vector<Vec3> & directions) const {
            //VisualizeReconstruction({ mg.ccidsBigToSmall.front() }, mg);
        }
    };


    std::unique_ptr<ReconstructedModel> PredictionOfPanoramix(const std::string & impath,
        const PredictOptions & options, misc::Matlab & matlab, const PILayoutAnnotation & anno) {

        auto resultPath = PanoramixResultFilePath(impath, 0);

        PIGraph mg;
        if (1 || !LoadFromDisk(resultPath, mg)) {
            // rebuild pigraph
            misc::Matlab matlab;

            Image3ub originalImage = anno.rectifiedImage;

            auto image = originalImage.clone();
            ResizeToHeight(image, 700);

            /// prepare things!
            View<PanoramicCamera, Image3ub> view;
            std::vector<PerspectiveCamera> cams;
            std::vector<Classified<Line3>> line3s;
            std::vector<Vec3> vps;
            int vertVPId;
            Imagei segs;
            int nsegs;


            bool refresh_preparation = false;
            if (refresh_preparation || !misc::LoadCache(impath, "preparation", view, cams, line3s, vps, vertVPId, segs, nsegs)) {
                view = CreatePanoramicView(image);

                // collect lines in each view
                cams = CreateCubicFacedCameras(view.camera, image.rows, image.rows, image.rows * 0.4);
                std::vector<Line3> rawLine3s;
                std::vector<Line2> rawLine2s;
                for (int i = 0; i < cams.size(); i++) {
                    auto pim = view.sampled(cams[i]).image;
                    LineSegmentExtractor lineExtractor;
                    lineExtractor.params().algorithm = LineSegmentExtractor::LSD;
                    auto ls = lineExtractor(pim, 2, 300); // use pyramid
                    for (auto & l : ls) {
                        rawLine2s.push_back(l);
                        rawLine3s.emplace_back(normalize(cams[i].toSpace(l.first)),
                            normalize(cams[i].toSpace(l.second)));
                    }
                }
                rawLine3s = MergeLines(rawLine3s, DegreesToRadians(3), DegreesToRadians(5));

                // estimate vp
                line3s = ClassifyEachAs(rawLine3s, -1);
                vps = EstimateVanishingPointsAndClassifyLines(line3s);
                vertVPId = NearestDirectionId(vps, Vec3(0, 0, 1));

                if (0) {
                    gui::ColorTable ctable = gui::RGBGreys;
                    for (int i = 0; i < cams.size(); i++) {
                        auto & cam = cams[i];
                        std::vector<Classified<Line2>> lines;
                        for (auto & l3 : line3s) {
                            if (!cam.isVisibleOnScreen(l3.component.first) || !cam.isVisibleOnScreen(l3.component.second)) {
                                continue;
                            }
                            auto p1 = cam.toScreen(l3.component.first);
                            auto p2 = cam.toScreen(l3.component.second);
                            lines.push_back(ClassifyAs(Line2(p1, p2), l3.claz));
                        }
                        auto pim = view.sampled(cams[i]).image;
                        gui::AsCanvas(pim).thickness(3).colorTable(ctable).add(lines).show();
                    }
                }


                // estimate segs
                nsegs = SegmentationForPIGraph(view, line3s, segs, DegreesToRadians(1));
                RemoveThinRegionInSegmentation(segs, 1, true);
                nsegs = DensifySegmentation(segs, true);
                assert(IsDenseSegmentation(segs));

                if (1) {
                    auto ctable = gui::CreateGreyColorTableWithSize(nsegs);
                    ctable.randomize();
                    gui::ColorTable rgb = gui::RGBGreys;
                    auto canvas = gui::MakeCanvas(view.image).alpha(0.9).add(ctable(segs));
                    for (auto & l : line3s) {
                        static const double sampleAngle = M_PI / 100.0;
                        auto & line = l.component;
                        double spanAngle = AngleBetweenDirections(line.first, line.second);
                        std::vector<Point2> ps; ps.reserve(spanAngle / sampleAngle);
                        for (double angle = 0.0; angle <= spanAngle; angle += sampleAngle) {
                            Vec3 dir = RotateDirection(line.first, line.second, angle);
                            ps.push_back(view.camera.toScreen(dir));
                        }
                        for (int i = 1; i < ps.size(); i++) {
                            auto & p1 = ps[i - 1];
                            auto & p2 = ps[i];
                            if (Distance(p1, p2) >= view.image.cols / 2) {
                                continue;
                            }
                            canvas.thickness(2);
                            canvas.colorTable(rgb).add(gui::ClassifyAs(Line2(p1, p2), l.claz));
                        }
                    }
                    canvas.show();
                }


                // save
                misc::SaveCache(impath, "preparation", view, cams, line3s, vps, vertVPId, segs, nsegs);
            }



            // gc !!!!
            std::vector<PerspectiveCamera> hcams;
            std::vector<Weighted<View<PerspectiveCamera, Image5d>>> gcs;
            Image5d gc;
            static const int hcamNum = 16;
            static const Sizei hcamScreenSize(500, 500);
            //static const Sizei hcamScreenSize(500, 700);
            static const int hcamFocal = 200;
            std::string hcamsgcsFileName;
            {
                std::stringstream ss;
                ss << "hcamsgcs_" << hcamNum << "_" << hcamScreenSize.width << "_" << hcamScreenSize.height << "_" << hcamFocal;
                hcamsgcsFileName = ss.str();
            }
            if (0 || !misc::LoadCache(impath, hcamsgcsFileName, hcams, gcs)) {
                // extract gcs
                hcams = CreateHorizontalPerspectiveCameras(view.camera, hcamNum, hcamScreenSize.width, hcamScreenSize.height, hcamFocal);
                gcs.resize(hcams.size());
                for (int i = 0; i < hcams.size(); i++) {
                    auto pim = view.sampled(hcams[i]);
                    auto pgc = ComputeGeometricContext(matlab, pim.image, false, true);
                    gcs[i].component.camera = hcams[i];
                    gcs[i].component.image = pgc;
                    gcs[i].score = abs(1.0 - normalize(hcams[i].forward()).dot(normalize(view.camera.up())));
                }
                misc::SaveCache(impath, hcamsgcsFileName, hcams, gcs);
            }
            std::string gcmergedFileName;
            {
                std::stringstream ss;
                ss << "gc_" << hcamNum << "_" << hcamScreenSize.width << "_" << hcamScreenSize.height << "_" << hcamFocal;
                gcmergedFileName = ss.str();
            }
            if (0 || !misc::LoadCache(impath, gcmergedFileName, gc)) {
                gc = Combine(view.camera, gcs).image;
                if (1) {
                    std::vector<Imaged> gcChannels;
                    cv::split(gc, gcChannels);
                    gui::AsCanvas(ConvertToImage3d(gc)).show();
                }
                misc::SaveCache(impath, gcmergedFileName, gc);
            }



            // build pigraph!
            bool refresh_mg_init = refresh_preparation || false;
            if (refresh_mg_init || !misc::LoadCache(impath, "mg_init", mg)) {
                mg = BuildPIGraph(view, vps, vertVPId, segs, line3s,
                    DegreesToRadians(1), DegreesToRadians(1), DegreesToRadians(2),
                    DegreesToRadians(5), DegreesToRadians(60), DegreesToRadians(5));

                const auto printPIGraph = [&mg](int delay) {
                    static const gui::ColorTable randColors = gui::CreateRandomColorTableWithSize(mg.nsegs);
                    static const gui::ColorTable rgbGrayColors = gui::RGBGreys;
                    auto pim = Print(mg,
                        [&mg](int seg) -> gui::Color {
                        return gui::White;
                    },
                        [&mg](int lp) {
                        return mg.linePiece2bndPiece[lp] == -1 ? gui::Red : gui::Black;
                    },
                        [&mg](int bp) -> gui::Color {
                        return gui::LightGray;
                    }, 1, 3);
                    gui::AsCanvas(pim).show(delay, "pi graph");
                };

                printPIGraph(0);
                misc::SaveCache(impath, "mg_init", mg);
            }


            const auto printPIGraphControls = [&mg](int delay) {
                const auto bpColor = [&mg](int bp) -> gui::Color {
                    auto occ = mg.bndPiece2segRelation[bp];
                    if (occ == SegRelation::Connected) {
                        return gui::LightGray;
                    } else if (occ == SegRelation::LeftIsFront) {
                        return gui::Color(gui::Green);
                    } else {
                        return gui::Color(gui::Blue);
                    }
                };
                const auto lpColor = [&mg, bpColor](int lp) -> gui::Color {
                    if (mg.linePiece2bndPiece[lp] != -1) {
                        return bpColor(mg.linePiece2bndPiece[lp]) * 0.9;
                    } else {
                        return mg.linePiece2segLineRelation[lp] == SegLineRelation::Detached ? gui::Black : gui::LightGray;
                    }
                };
                const auto segColor = [&mg](int seg) -> gui::Color {
                    static const gui::ColorTable ctable = gui::ColorTableDescriptor::RGBGreys;
                    auto & c = mg.seg2control[seg];
                    if (!c.used) {
                        return gui::Black;
                    }
                    if (c.orientationClaz != -1) {
                        return ctable[c.orientationClaz].blendWith(gui::White, 0.5);
                    }
                    if (c.orientationNotClaz != -1) {
                        return ctable[c.orientationNotClaz].blendWith(gui::White, 0.3);
                    }
                    return gui::White;
                };
                auto pim = Print(mg, segColor, ConstantFunctor<gui::Color>(gui::Transparent), bpColor, 1, 3);
                gui::AsCanvas(pim).show(delay, "pi graph controled");
            };



            // attach orientation constraints
            bool refresh_mg_oriented = refresh_mg_init || false;
            if (refresh_mg_oriented || !misc::LoadCache(impath, "mg_oriented", mg)) {
                AttachPrincipleDirectionConstraints(mg);
                AttachWallConstraints(mg, M_PI / 60.0);
                AttachGCConstraints(mg, gc);
                printPIGraphControls(0);
                misc::SaveCache(impath, "mg_oriented", mg);
            }



            // detect occlusions
            if (1 || !misc::LoadCache(impath, "mg_occdetected", mg)) {
                DetectOcclusions2(mg);
                printPIGraphControls(0);
                misc::SaveCache(impath, "mg_occdetected", mg);
            }


            PIConstraintGraph cg(mg);
            cg.reconstructLargestCC();


            // save to disk
            SaveToDisk(resultPath, mg);
        }

        return std::make_unique<PIGraphModel>(std::move(mg));
    }

}
