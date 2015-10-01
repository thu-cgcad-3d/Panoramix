
#include "eval.hpp"

#include "../../src/experimental/rl_graph_solver.hpp"
#include "../../src/experimental/rl_graph_occlusion.hpp"
#include "../../src/experimental/rl_graph_modeling.hpp"
#include "../../src/experimental/tools.hpp"

#include "../../src/experimental/pi_graph_control.hpp"
#include "../../src/experimental/pi_graph_occlusion.hpp"
#include "../../src/experimental/pi_graph_solve.hpp"
#include "../../src/experimental/pi_graph_vis.hpp"

namespace panolyz {

    std::string PanoramixResultFilePath(const std::string & impath, int version) {
        return impath + ".panoramix." + std::to_string(version) + ".cereal";
    }


    // old style rlgraph
    struct RLGraphModel : ReconstructedModel {
        View<PanoramicCamera, Image3ub> view;
        std::vector<Vec3> vps;
        int vertVPId;
        RLGraph mg;
        std::vector<RegionHandle> rhs;
        std::vector<RegionBoundaryHandle> bhs;
        RLGraphControls controls;
        RLGraphVars vars;
    };


    // pi graph

    struct PIGraphModel : ReconstructedModel {
        PIGraph mg;
        explicit PIGraphModel(PIGraph && g) : mg(std::move(g)) {}

        virtual double depthAt(const Vec3 & direction) const {
            int seg = mg.segs(ToPixel(mg.view.camera.toScreen(direction)));
            return norm(IntersectionOfLineAndPlane(Ray3(Origin(), direction), mg.seg2recPlanes[seg]).position);
        }
        virtual void visualize(const std::vector<Vec3> & directions) const {
            VisualizeReconstruction({ mg.ccidsBigToSmall.front() }, mg);
        }
    };


    std::unique_ptr<ReconstructedModel> PredictionOfPanoramix(const std::string & impath,
        const PredictOptions & options, misc::Matlab & matlab) {

        auto resultPath = PanoramixResultFilePath(impath, 0);

        PIGraph mg;
        if (1 || !LoadFromDisk(resultPath, mg)) {
            // rebuild pigraph
            misc::Matlab matlab;

            Image3ub originalImage = cv::imread(impath);

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


            if (0 || !misc::LoadCache(impath, "preparation", view, cams, line3s, vps, vertVPId, segs, nsegs)) {
                view = CreatePanoramicView(image);

                // collect lines in each view
                cams = CreateCubicFacedCameras(view.camera, image.rows, image.rows, image.rows * 0.4);
                std::vector<Line3> rawLine3s;
                for (int i = 0; i < cams.size(); i++) {
                    auto pim = view.sampled(cams[i]).image;
                    LineSegmentExtractor lineExtractor;
                    lineExtractor.params().algorithm = LineSegmentExtractor::LSD;
                    auto ls = lineExtractor(pim, 2, 300); // use pyramid
                    for (auto & l : ls) {
                        rawLine3s.emplace_back(normalize(cams[i].toSpace(l.first)),
                            normalize(cams[i].toSpace(l.second)));
                    }
                }
                rawLine3s = MergeLines(rawLine3s, DegreesToRadians(1));

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



                // segmentation
                SegmentationExtractor segmenter;
                segmenter.params().algorithm = SegmentationExtractor::GraphCut;
                segmenter.params().sigma = 10.0;
                segmenter.params().c = 5.0;
                segmenter.params().superpixelSizeSuggestion = 2000;
                std::tie(segs, nsegs) = segmenter(view.image, rawLine3s, view.camera, DegreesToRadians(1));

                RemoveThinRegionInSegmentation(segs, 1, true);
                RemoveSmallRegionInSegmentation(segs, 100, true);
                RemoveDanglingPixelsInSegmentation(segs, true);
                nsegs = DensifySegmentation(segs, true);
                assert(IsDenseSegmentation(segs));
                if (1) {
                    auto ctable = gui::CreateRandomColorTableWithSize(nsegs);
                    gui::AsCanvas(ctable(segs)).add(view.image).show();
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
            mg = BuildPIGraph(view, vps, vertVPId, segs, line3s,
                DegreesToRadians(1), DegreesToRadians(1), DegreesToRadians(1),
                0.04, M_PI_2, 0.02);

            {
                auto pim = Print(mg,
                    [&mg](int seg) -> gui::Color {
                    return gui::White;
                },
                    [&mg](int lp) {
                    return gui::Transparent;
                },
                    [&mg](int bp) -> gui::Color {
                    static gui::ColorTable ctable(std::vector<gui::Color>{ gui::Red, gui::Green, gui::Blue }, gui::Gray);
                    return mg.bndPiece2linePieces[bp].empty() ? gui::LightGray : ctable[mg.bndPiece2classes[bp]];
                }, 1, 3);
                cv::imshow("pigraph", Image3f(pim * 0.99f) + Image3f(mg.view.image / 255.0f * 0.01f));
                cv::waitKey();
            }

            // attach orientation constraints
            AttachPrincipleDirectionConstraints(mg);
            AttachWallConstraints(mg, M_PI / 100.0);
            AttachGCConstraints(mg, gc);

            // detect occlusions
            if (options.useGroundTruthOcclusions) {
                auto anno = LoadOrInitializeNewLayoutAnnotation(impath);
                // todo!!!
            } else {
                DetectOcclusions(mg);
            }
            //AssumeThereAreNoOcclusions(mg);

            {
                auto pim = Print(mg,
                    [&mg](int seg) -> gui::Color {
                    static const gui::ColorTable ctable = gui::ColorTableDescriptor::RGBGreys;
                    auto & c = mg.seg2control[seg];
                    if (!c.used) {
                        return gui::Black;
                    }
                    if (c.orientationClaz != -1) {
                        return ctable[c.orientationClaz];
                    }
                    if (c.orientationNotClaz != -1) {
                        return ctable[c.orientationNotClaz] * 0.5;
                    }
                    return gui::White;
                },
                    [&mg](int lp) {
                    static const gui::ColorTable ctable = gui::ColorTableDescriptor::RGBGreys;
                    return ctable[mg.lines[mg.linePiece2line[lp]].claz].blendWith(gui::White, 0.5) *
                        ((mg.linePiece2bndPiece[lp] == -1 && mg.linePiece2attachment[lp] == AttachmentRelation::Detached) ? 0.1 : 1.0);
                },
                    [&mg](int bp) -> gui::Color {
                    auto occ = mg.bndPiece2occlusion[bp];
                    if (occ == OcclusionRelation::Connected) {
                        return gui::DarkGray;
                    } else if (occ == OcclusionRelation::LeftIsFront) {
                        return gui::Color(gui::Black) * 0.6;
                    } else {
                        return gui::Color(gui::Black) * 0.6;
                    }
                }, 2, 3);
                cv::imshow("pigraph_controled", Image3f(pim * 0.99f) + Image3f(mg.view.image / 255.0f * 0.01f));
                cv::waitKey();
            }

            // build constraint graph
            BuildConstraintGraph(mg);
            {
                auto cctable = gui::CreateRandomColorTableWithSize(mg.nccs);
                auto pim = Print(mg,
                    [&mg, &cctable](int seg) -> gui::Color {
                    int vert = mg.seg2vert[seg];
                    if (vert == -1) return gui::Transparent;
                    return cctable[mg.vert2cc[vert]];
                },
                    [&mg, &cctable](int lp) -> gui::Color {
                    int vert = mg.line2vert[mg.linePiece2line[lp]];
                    if (vert == -1) return gui::Transparent;
                    return cctable[mg.vert2cc[vert]];
                },
                    [](int bp) -> gui::Color {
                    return gui::LightGray;
                }, 1, 2);
                cv::imshow("pigraph_ccs", Image3f(pim * 0.99f) + Image3f(mg.view.image / 255.0f * 0.01f));
                cv::waitKey();
            }

            // solve
            SolvePIGraph(mg.ccidsBigToSmall.front(), mg, matlab, 100);


            // save to disk
            SaveToDisk(resultPath, mg);
        }

        return std::make_unique<PIGraphModel>(std::move(mg));
    }

}