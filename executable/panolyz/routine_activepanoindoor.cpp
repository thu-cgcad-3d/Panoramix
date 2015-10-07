#include "../../src/core/basic_types.hpp"
#include "../../src/experimental/rl_graph_solver.hpp"
#include "../../src/experimental/tools.hpp"
#include "../../src/misc/matlab_api.hpp"
#include "../../src/gui/scene.hpp"
#include "../../src/gui/canvas.hpp"
#include "../../src/gui/utility.hpp"

#include "routines.hpp"


namespace panolyz {

    namespace ActivePanoramaIndoor{

        void Run(){       

            std::string path;

            using namespace pano;
            using namespace core;
            using namespace experimental;

            misc::Matlab matlab;

            /*while (1){
                int r = gui::SelectFrom({ "Yes", "No", "Cancel" }, "Confirm?");
                std::cout << "selected: " << r << std::endl;
            }*/

            core::Image3ub image = gui::PickAnImage(PROJECT_TEST_DATA_DIR_STR"/panorama/indoor", &path);

            std::vector<gui::PenConfig> penConfigs = {
                { "A", "a", 2.0, gui::Red, gui::PenStyle::SolidLine },
                { "B", "b", 2.0, gui::Blue, gui::PenStyle::DashDotDotLine }
            };


            if (image.empty())
                return;

            MakePanorama(image);
            ResizeToHeight(image, 700);

            View<PanoramicCamera, Image3ub> view;

            std::vector<PerspectiveCamera> cams;
            std::vector<std::vector<Classified<Line2>>> lines;
            std::vector<Vec3> vps;
            int vertVPId;

            Imagei segmentedImage;

#define REFRESH 1

            if (REFRESH){
                view = CreatePanoramicView(image);

                //gui::PaintWithPanorama(view, penConfigs,
                //    [](const std::vector<Point2> & polyline, int penId){
                //    return false;
                //});

                // collect lines in each view
                cams = CreateCubicFacedCameras(view.camera, image.rows, image.rows, image.rows * 0.4);
                lines.resize(cams.size());
                for (int i = 0; i < cams.size(); i++){
                    auto pim = view.sampled(cams[i]).image;
                    LineSegmentExtractor lineExtractor;
                    lineExtractor.params().algorithm = LineSegmentExtractor::LSD;
                    auto ls = lineExtractor(pim, 2, 300); // use pyramid
                    lines[i].reserve(ls.size());
                    for (auto & l : ls){
                        lines[i].push_back(ClassifyAs(l, -1));
                    }
                }

                // estimate vp
                vps = EstimateVanishingPointsAndClassifyLines(cams, lines);
                if (0){
                    auto ctable = gui::CreateRandomColorTableWithSize(vps.size());
                    for (int i = 0; i < cams.size(); i++){
                        auto pim = view.sampled(cams[i]).image;
                        gui::AsCanvas(pim).thickness(3).colorTable(ctable).add(lines[i]).show();
                    }
                }
                vertVPId = NearestDirectionId(vps, Vec3(0, 0, 1));

                // extract lines from segmentated region boundaries and classify them using estimated vps
                // make 3d lines
                std::vector<Line3> line3ds;
                for (int i = 0; i < cams.size(); i++){
                    for (auto & l : lines[i]){
                        line3ds.emplace_back(normalize(cams[i].toSpace(l.component.first)),
                            normalize(cams[i].toSpace(l.component.second)));
                    }
                }
                SegmentationExtractor segmenter;
                segmenter.params().algorithm = SegmentationExtractor::GraphCut;
                segmenter.params().sigma = 10.0;
                segmenter.params().c = 1.0;
                segmenter.params().superpixelSizeSuggestion = 2000;
                int segmentsNum = 0;
                std::tie(segmentedImage, segmentsNum) = segmenter(view.image, line3ds, view.camera, M_PI / 36.0);

                if (0){
                    auto ctable = gui::CreateRandomColorTableWithSize(segmentsNum);
                    gui::AsCanvas(ctable(segmentedImage)).show();
                    gui::AsCanvas(ctable(segmentedImage)).add(view.image).show();
                }

                misc::SaveCache(path, "pre", view, cams, lines, vps, segmentedImage, vertVPId);
            }
            else{
                misc::LoadCache(path, "pre", view, cams, lines, vps, segmentedImage, vertVPId);
            }


            std::vector<PerspectiveCamera> hcams;
            std::vector<Weighted<View<PerspectiveCamera, Image5d>>> gcs;
            if (REFRESH){
                // extract gcs
                hcams = CreateHorizontalPerspectiveCameras(view.camera, 16, 500, 400, 300);
                gcs.resize(hcams.size());
                for (int i = 0; i < hcams.size(); i++){
                    auto pim = view.sampled(hcams[i]);
                    auto pgc = ComputeGeometricContext(matlab, pim.image, false, true);
                    gcs[i].component.camera = hcams[i];
                    gcs[i].component.image = pgc;
                    gcs[i].score = sin(AngleBetweenUndirectedVectors(hcams[i].forward(), view.camera.up()));
                }
                misc::SaveCache(path, "hcamsgcs", hcams, gcs);
            }
            else{
                misc::LoadCache(path, "hcamsgcs", hcams, gcs);
            }

            Image5d gc;
            gc = Combine(view.camera, gcs).image;
            if (1){
                gui::AsCanvas(gc).show();
            }






        }

    }

}