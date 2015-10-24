#include "../../src/core/basic_types.hpp"
#include "../../src/core/containers.hpp"
#include "../../src/experimental/rl_graph_solver.hpp"
#include "../../src/experimental/rl_graph_occlusion.hpp"
#include "../../src/experimental/rl_graph_modeling.hpp"
#include "../../src/experimental/tools.hpp"
#include "../../src/misc/matlab_api.hpp"
#include "../../src/gui/scene.hpp"
#include "../../src/gui/canvas.hpp"
#include "../../src/gui/utility.hpp"

#include "routines.hpp"


namespace panolyz {

    namespace GC {

        void Run() {

            using namespace pano;
            using namespace core;
            using namespace experimental;

            std::vector<std::string> errorPaths;

            misc::Matlab matlab("", true);
            gui::ForEachImageFromAFolder("H:\\DataSet\\", [&matlab, &errorPaths](const std::string & path) -> bool {
                try {
                    std::cout << "processing image " << path << std::endl;
                    Image3ub original = cv::imread(path);
                    if (original.empty())
                        return true;

                    Image3ub image;
                    bool extendedOnTop = false, extendedOnBottom = false;

                    static const bool automaticallyRectifyIfNeeded = false;
                    if (0 || !misc::LoadCache(path, "rectified", image, extendedOnTop, extendedOnBottom)) {
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
                        misc::SaveCache(path, "rectified", image, extendedOnTop, extendedOnBottom);
                    }

                    ResizeToHeight(image, 700);
                    if (0) {
                        cv::imshow("rectified", image);
                        cv::waitKey();
                    }


                    View<PanoramicCamera, Image3ub> view;
                    view = CreatePanoramicView(image);

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
                        if (0 || !misc::LoadCache(path, hcamsgcsFileName, hcams, gcs)) {
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
                            misc::SaveCache(path, hcamsgcsFileName, hcams, gcs);
                        }
                        std::string gcmergedFileName;
                        {
                            std::stringstream ss;
                            ss << "gc_" << hcamNum << "_" << hcamScreenSize.width << "_" << hcamScreenSize.height << "_" << hcamFocal;
                            gcmergedFileName = ss.str();
                        }
                        if (0 || !misc::LoadCache(path, gcmergedFileName, gc)) {
                            gc = Combine(view.camera, gcs).image;
                            misc::SaveCache(path, gcmergedFileName, gc);
                        }
                    }

                    if (1) {
                        std::vector<Imaged> gcChannels;
                        cv::split(gc, gcChannels);
                        gui::AsCanvas(ConvertToImage3d(gc)).show(5, "gc merged");
                    }
                } catch (...) {
                    errorPaths.push_back(path);
                }
                return true;
            });

            misc::SaveCache("routine_gc", "error_paths", errorPaths);
        }
    }
}