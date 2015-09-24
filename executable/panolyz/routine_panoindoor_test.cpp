#include "../../src/core/basic_types.hpp"
#include "../../src/core/containers.hpp"
#include "../../src/experimental/pi_graph.hpp"
#include "../../src/experimental/pi_graph_control.hpp"
#include "../../src/experimental/pi_graph_occlusion.hpp"
#include "../../src/experimental/pi_graph_solve.hpp"
#include "../../src/experimental/pi_graph_vis.hpp"
#include "../../src/experimental/pi_graph_annotation.hpp"
#include "../../src/experimental/tools.hpp"
#include "../../src/misc/matlab_api.hpp"
#include "../../src/gui/scene.hpp"
#include "../../src/gui/canvas.hpp"
#include "../../src/gui/utility.hpp"

#include "routines.hpp"


namespace panolyz {

    namespace PanoramaIndoorTest {

        static const bool refresh = true;

        void Run() {

            using namespace pano;
            using namespace core;
            using namespace experimental;

            misc::Matlab matlab;

            std::vector<std::string> paths;
            std::vector<Image> images = gui::PickImages("H:\\DataSet\\pi\\", &paths);
            assert(paths.size() == images.size());

            for (int k = 0; k < images.size(); k++) {
                // get annotation
                const std::string & path = paths[k];                
                PIAnnotation anno = LoadOrInitializeNewAnnotation(path);

                View<PanoramicCamera, Image3ub> view = anno.view;
                auto image = view.image;
                std::vector<Classified<Line3>> line3s;
                { // remove small lines and unclassified lines
                    for (auto & line : anno.lines) {
                        if (AngleBetweenDirections(line.component.first, line.component.second) < DegreesToRadians(3) || line.claz == -1) {
                            continue;
                        }
                        line3s.push_back(line);
                    }
                }
                std::vector<Vec3> vps = anno.vps;
                int vertVPId = anno.vertVPId;               
                
                Imagei segmentedImage;
                if (1 || !Load(path, "segs", segmentedImage)) {              
                    std::vector<Line3> rawLine3s;
                    for (auto & line : anno.lines) {
                        rawLine3s.push_back(line.component);
                    }
                    for (auto & occ : anno.occlusions) {
                        for (int i = 1; i < occ.chain.size(); i++) {
                            rawLine3s.emplace_back(occ.chain[i - 1], occ.chain[i]);
                        }
                    }

                    // segmentation
                    SegmentationExtractor segmenter;
                    segmenter.params().algorithm = SegmentationExtractor::GraphCut;
                    segmenter.params().sigma = 10.0;
                    segmenter.params().c = 5.0;
                    segmenter.params().superpixelSizeSuggestion = 2000;
                    int segmentsNum = 0;
                    std::tie(segmentedImage, segmentsNum) = segmenter(view.image, rawLine3s, view.camera, DegreesToRadians(1));
                    std::cout << "before removing small regions, nsegs = " << segmentsNum << std::endl;
                    segmentsNum = RemoveSmallRegionInSegmentation(segmentedImage, 500.0, true);
                    std::cout << "after removing small regions, nsegs = " << segmentsNum << std::endl;
                    //RemoveThinRegionInSegmentation(segmentedImage, true);
                    //segmentsNum = DensifySegmentation(segmentedImage, true);
                    assert(IsDenseSegmentation(segmentedImage));

                    Save(path, "segs", segmentedImage);
                }

                if (1) {
                    auto ctable = gui::CreateRandomColorTableWithSize(MinMaxValOfImage(segmentedImage).second + 1);
                    gui::AsCanvas(ctable(segmentedImage)).add(view.image).show();
                }


                // pi graph
                static const bool rebuild = true;

                PIGraph mg;
                if (1 || !Load(path, "pigraph", mg)) {
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
                if (1 || !Load(path, "pigraph_controled", mg)) {
                    auto up = normalize(vps[vertVPId]);
                    if (up.dot(-view.camera.up()) < 0) {
                        up = -up;
                    }
                   
                    AttachAnnotatedPolygonsAndOcclusions(mg, anno.polygons, anno.occlusions);

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
                if (1 || !Load(path, "pigraph_solved", mg)) {

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
                            [&mg](int bp) -> gui::Color {
                            return gui::Gray;
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