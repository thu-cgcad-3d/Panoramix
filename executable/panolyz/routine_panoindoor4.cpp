#include "../../src/core/basic_types.hpp"
#include "../../src/core/containers.hpp"
#include "../../src/experimental/pi_graph.hpp"
#include "../../src/experimental/pi_graph_control.hpp"
#include "../../src/experimental/pi_graph_occlusion.hpp"
#include "../../src/experimental/pi_graph_solve.hpp"
#include "../../src/experimental/pi_graph_vis.hpp"
#include "../../src/experimental/tools.hpp"
#include "../../src/misc/matlab_api.hpp"
#include "../../src/gui/scene.hpp"
#include "../../src/gui/canvas.hpp"
#include "../../src/gui/utility.hpp"

#include "routines.hpp"


namespace panolyz {

    namespace PanoramaIndoor4 {

        static const bool refresh = true;

        void Run() {

            using namespace pano;
            using namespace core;
            using namespace experimental;

            misc::Matlab matlab;

            std::vector<std::string> paths;
            std::vector<Image> images = gui::PickImages(PROJECT_TEST_DATA_DIR_STR"/panorama/indoor", &paths);
            assert(paths.size() == images.size());

            for (int k = 0; k < images.size(); k++) {
                const std::string & path = paths[k];
                const Image3ub & original = images[k];
                if (original.empty())
                    continue;


                Image3ub image;
                bool extendedOnTop = false, extendedOnBottom = false;

                static const bool automaticallyRectifyIfNeeded = true;
                if (0 || !Load(path, "rectified", image, extendedOnTop, extendedOnBottom)) {
                    image = original.clone();
                    if (!automaticallyRectifyIfNeeded) {
                        if (!gui::MakePanoramaByHand(image, &extendedOnTop, &extendedOnBottom)) {
                            WARNNING("failed making panorama");
                        }
                    } else {
                        if (!MakePanorama(image, -1, &extendedOnTop, &extendedOnBottom)) {
                            WARNNING("failed making panorama");
                        }
                    }
                    Save(path, "rectified", image, extendedOnTop, extendedOnBottom);
                }

                ResizeToHeight(image, 700);
                if (1) {
                    cv::imshow("rectified", image);
                    cv::waitKey();
                }


                View<PanoramicCamera, Image3ub> view;

                std::vector<PerspectiveCamera> cams;
                std::vector<Classified<Line3>> line3s;
                std::vector<Vec3> vps;
                int vertVPId;

                Imagei segmentedImage;

                if (1 || !Load(path, "pre2", view, cams, line3s, vps, segmentedImage, vertVPId)) {
                    view = CreatePanoramicView(image);

                    // collect lines in each view
                    cams = CreateCubicFacedCameras(view.camera, image.rows, image.rows, image.rows * 0.4);
                    std::vector<Line3> rawLine3s;
                    for (int i = 0; i < cams.size(); i++) {
                        auto pim = view.sampled(cams[i]).image;
                        LineSegmentExtractor lineExtractor;
                        lineExtractor.params().algorithm = LineSegmentExtractor::LSD;
                        auto ls = lineExtractor(pim, 3, 300); // use pyramid
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


                    // segmentation
                    SegmentationExtractor segmenter;
                    segmenter.params().algorithm = SegmentationExtractor::GraphCut;
                    segmenter.params().sigma = 10.0;
                    segmenter.params().c = 5.0;
                    segmenter.params().superpixelSizeSuggestion = 2000;
                    int segmentsNum = 0;
                    std::tie(segmentedImage, segmentsNum) = segmenter(view.image, rawLine3s, view.camera, DegreesToRadians(1));
                    //RemoveThinRegionInSegmentation(segmentedImage, true);
                    //segmentsNum = DensifySegmentation(segmentedImage, true);
                    segmentsNum = RemoveSmallRegionInSegmentation(segmentedImage, 100, true);
                    assert(IsDenseSegmentation(segmentedImage));

                    Save(path, "pre2", view, cams, line3s, vps, segmentedImage, vertVPId);
                }

                //if (refresh) {
                //    auto ctable = gui::CreateRandomColorTableWithSize(vps.size());
                //    for (int i = 0; i < cams.size(); i++) {
                //        auto pim = view.sampled(cams[i]).image;
                //        for (auto & l : line3s) {
                //            auto ll = ClassifyAs(Line2(cams[i].toScreen(l.component.first, cams[i].toScreen(l.component.second))
                //            gui::AsCanvas(pim).thickness(3).colorTable(ctable).add(lines[i]).show();
                //        }
                //    }
                //}

                if (0) {
                    auto ctable = gui::CreateRandomColorTableWithSize(MinMaxValOfImage(segmentedImage).second + 1);
                    //gui::AsCanvas(ctable(segmentedImage)).show();
                    gui::AsCanvas(ctable(segmentedImage)).add(view.image).show();
                }


                std::vector<PerspectiveCamera> hcams;
                std::vector<Weighted<View<PerspectiveCamera, Image5d>>> gcs;
                Image5d gc;
                {
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
                    if (0 || !Load(path, hcamsgcsFileName, hcams, gcs)) {
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
                        Save(path, hcamsgcsFileName, hcams, gcs);
                    }
                    std::string gcmergedFileName;
                    {
                        std::stringstream ss;
                        ss << "gc_" << hcamNum << "_" << hcamScreenSize.width << "_" << hcamScreenSize.height << "_" << hcamFocal;
                        gcmergedFileName = ss.str();
                    }
                    if (0 || !Load(path, gcmergedFileName, gc)) {
                        gc = Combine(view.camera, gcs).image;
                        Save(path, gcmergedFileName, gc);
                    }
                }


                if (1) {
                    std::vector<Imaged> gcChannels;
                    cv::split(gc, gcChannels);
                    gui::AsCanvas(ConvertToImage3d(gc)).show();
                    misc::MAT gcMat("./cache/gcMat.mat", misc::MAT::Write);
                    gcMat.setVar("gc", gc);
                }


                // pi graph
                static const bool rebuild = true;

                PIGraph mg;
                if (rebuild || !Load(path, "pigraph", mg)) {
                    mg = BuildPIGraph(view, vps, vertVPId, segmentedImage, line3s,
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

                    Save(path, "pigraph", mg);
                }

                // set pi graph constraints
                if (rebuild || !Load(path, "pigraph_controled", mg)) {
                    auto up = normalize(vps[vertVPId]);
                    if (up.dot(-view.camera.up()) < 0) {
                        up = -up;
                    }
                    AttachPrincipleDirectionConstraints(mg);
                    AttachWallConstraints(mg, M_PI / 100.0);
                    if (extendedOnTop) {
                        DisableTopSeg(mg);
                    }

                    AttachGCConstraints(mg, gc);
                    //DetectOcclusions(mg);
                    AssumeThereAreNoOcclusions(mg);

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
                        }, 1, 3);
                        cv::imshow("pigraph_controled", Image3f(pim * 0.99f) + Image3f(mg.view.image / 255.0f * 0.01f));
                        cv::waitKey();
                    }

                    Save(path, "pigraph_controled", mg);
                }

                // solve pi graph
                if (rebuild || !Load(path, "pigraph_solved", mg)) {
                    
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

                    SolvePIGraph(mg.ccidsBigToSmall.front(), mg, matlab, 100);
                    
                    Save(path, "pigraph_solved", mg);
                }


                VisualizeReconstruction({ mg.ccidsBigToSmall.front() }, mg);

            }

        }
    }
}