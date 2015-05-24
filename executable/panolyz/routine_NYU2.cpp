#include <thread>

#include "../../src/core/basic_types.hpp"
#include "../../src/experimental/rl_graph.hpp"
#include "../../src/experimental/rl_graph_occlusion.hpp"
#include "../../src/experimental/tools.hpp"
#include "../../src/gui/scene.hpp"
#include "../../src/gui/canvas.hpp"
#include "../../src/misc/matlab_engine.hpp"
#include "../../src/gui/singleton.hpp"
#include "../../src/ml/factor_graph.hpp"

#include "tools.hpp"
#include "routines.hpp"


namespace panolyz {

    namespace NYU2 {


        using namespace panoramix;
        using namespace core;
        using namespace experimental;


        using Image7d = ImageOf<Vec<double, 7>>;


        void Show(const PerspectiveView & v, const RLGraph & mg, const RLGraphControls & controls, const RLGraphVars & vars){
            gui::SceneBuilder vis;
            Visualize(vis, v, mg, controls, vars);
            vis.show(true, true, gui::RenderOptions().camera(v.camera).cullBackFace(false).cullFrontFace(false));
        }

        struct CameraParams {
            Vec2 f_rgb, f_d;
            Vec2 c_rgb, c_d;
            Vec3 k_rgb, k_d;
            Vec2 p_rgb, p_d;
            Vec3 t;
            Mat3 R;
            double depthParam1, depthParam2;
            double maxDepth;
        };


        struct Data {
            CameraParams cp;
            Image image;
            Imagef depth, rawDepth;
            std::vector<std::vector<int32_t>> indices;
            std::vector<double> occscore;
            Image7d gc;
        };


        std::string nyu2dir = "F:\\DataSets\\NYU2\\";


        template <class FunT>
        void ForEachCase(FunT && fun){

            // load test data
            misc::MatlabEngine::CDAndAddAllSubfolders(nyu2dir);
            misc::MatlabEngine matlab;

            // get camera params
            matlab << "camera_params;";
            matlab
                << "f_rgb = [fx_rgb;fy_rgb];"
                << "f_d = [fx_d;fy_d];"
                << "c_rgb = [cx_rgb;cy_rgb];"
                << "c_d = [cx_d;cy_d];"
                << "k_rgb = [k1_rgb;k2_rgb;k3_rgb];"
                << "k_d = [k1_d;k2_d;k3_d];"
                << "p_rgb = [p1_rgb;p2_rgb];"
                << "p_d = [p1_d;p2_d];"
                << "t = [t_x;t_y;t_z];";

            CameraParams cameraParams;
            matlab.GetVariable("f_rgb", cameraParams.f_rgb);
            matlab.GetVariable("f_d", cameraParams.f_d);
            matlab.GetVariable("c_rgb", cameraParams.c_rgb);
            matlab.GetVariable("c_d", cameraParams.c_d);
            matlab.GetVariable("k_rgb", cameraParams.k_rgb);
            matlab.GetVariable("k_d", cameraParams.k_d);
            matlab.GetVariable("p_rgb", cameraParams.p_rgb);
            matlab.GetVariable("p_d", cameraParams.p_d);
            matlab.GetVariable("t", cameraParams.t);
            matlab.GetVariable("R", cameraParams.R);
            matlab.GetVariable("depthParam1", cameraParams.depthParam1);
            matlab.GetVariable("depthParam2", cameraParams.depthParam2);
            matlab.GetVariable("maxDepth", cameraParams.maxDepth);

            int offsetper30 = 0;
            int offsetper100 = 0;

            for (int id = 1; id <= 1449; id++){

                std::cout << "ID: " << id << std::endl;

                if (id % 30 == 1){
                    int first = id;
                    int last = id + 30 - 1;
                    std::cout << ("load occscores_[" + std::to_string(first) + "-" + std::to_string(last) + "] bndinfos occscores;") << std::endl;
                    matlab << ("load occscores_[" + std::to_string(first) + "-" + std::to_string(last) + "] bndinfos occscores;");
                    offsetper30 = id - 1;
                }
                if (id % 100 == 1){
                    int first = id;
                    int last = id + 100 - 1;
                    std::cout << ("load geometric_contexts_[" + std::to_string(first) + "-" + std::to_string(last) + "] gc;") << std::endl;
                    std::cout << ("load imdps_[" + std::to_string(first) + "-" + std::to_string(last) + "] ims dps;") << std::endl;
                    matlab << ("load geometric_contexts_[" + std::to_string(first) + "-" + std::to_string(last) + "] gc;");
                    matlab << ("load imdps_[" + std::to_string(first) + "-" + std::to_string(last) + "] ims dps rdps;");
                    offsetper100 = id - 1;
                }

                int trueIdPer30 = id - offsetper30;
                int trueIdPer100 = id - offsetper100;

                // get edges
                double edgenum = 0;
                matlab << ("inds = bndinfos(" + std::to_string(trueIdPer30) + ").edges.indices;");
                matlab << "edgenum = length(inds);";
                matlab.GetVariable("edgenum", edgenum);
                assert(edgenum > 0);

                std::vector<std::vector<int32_t>> indices(edgenum);
                for (int i = 0; i < edgenum; i++){
                    matlab << ("indices = inds{" + std::to_string(i + 1) + "}';");
                    matlab.GetVariable("indices", indices[i]);
                }

                // get occ scores
                std::vector<double> occscore;
                matlab << ("occscore = occscores{" + std::to_string(trueIdPer30) + "}';");
                matlab.GetVariable("occscore", occscore);
                assert(occscore.size() == (int)edgenum);



                // get image
                Image image;
                matlab << ("image = crop_image(ims(:,:,:, " + std::to_string(trueIdPer100) + "));");
                matlab.GetVariable("image", image, true);

                // get gc
                Image7d gc;
                matlab << ("g = crop_image(gc(:,:,:," + std::to_string(trueIdPer100) + "));");
                matlab.GetVariable("g", gc, true);

                // get depth
                Imagef depth;
                matlab << ("depth = crop_image(dps(:,:," + std::to_string(trueIdPer100) + "));");
                matlab.GetVariable("depth", depth, false);
                assert(image.size() == gc.size());
                assert(image.size() == depth.size());


                // get rawdepth
                Imagef rawDepth;
                matlab << ("rawDepth = crop_image(rdps(:,:," + std::to_string(trueIdPer100) + "));");
                matlab.GetVariable("rawDepth", rawDepth, false);
                assert(image.size() == gc.size());
                assert(image.size() == rawDepth.size());

                Data data{ cameraParams, image, depth, rawDepth, std::move(indices), std::move(occscore), gc };

                fun(id, data);

            }
        }


        template <class T, int N>
        T Mean(const Vec<T, N> & v) {
            return std::accumulate(v.val, v.val + N, 0.0) / N;
        }


        void Run(){


            ForEachCase([](int id, const Data & data){

                //if (id < 4)
                //    return;

                core::PerspectiveCamera cam_rgb(640, 480, data.cp.c_rgb, Mean(data.cp.f_rgb));
                core::PerspectiveCamera cam_d(640, 480, data.cp.c_d, Mean(data.cp.f_d));

                std::vector<Point3> gtPoints;
                gtPoints.reserve(data.depth.cols * data.depth.rows);
                for (auto it = data.depth.begin(); it != data.depth.end(); ++it){
                    auto d = cam_d.toSpace(it.pos()) * data.depth(it.pos()) / 1e5;
                    //auto d = normalize(cam_d.toSpace(it.pos())) * data.rawDepth(it.pos());
                    gtPoints.emplace_back(d);
                }

               /* gui::Visualizer vis;
                vis.installingOptions.pointSize = 1.0;
                vis.renderOptions.renderMode |= gui::RenderModeFlag::Points;
                vis.renderOptions.backgroundColor = gui::White;*/
                gui::SceneBuilder vis;
                vis.pointSize(1.0).add(gtPoints);
                vis.show(true, true, gui::RenderOptions().renderMode(gui::RenderModeFlag::All));


                Image image = data.image;


                PerspectiveView view;
                view.camera = cam_rgb;
                view.image = image;

                LineSegmentExtractor lineExtractor;
                lineExtractor.params().algorithm = core::LineSegmentExtractor::LSD;
                std::vector<Classified<Line2>> lines;
                lines = core::ClassifyEachAs(lineExtractor(view.image, 1), -1);

                // classify lines
                std::vector<Vec3> vps = EstimateVanishingPointsAndClassifyLines(view.camera, lines);
                // remove non-manhattan vps
                vps.erase(vps.begin() + 3, vps.end());
                for (auto & l : lines){
                    if (l.claz >= 3){
                        l.claz = -1;
                    }
                }

                // append regions
                SegmentationExtractor segmenter;
                segmenter.params().algorithm = SegmentationExtractor::GraphCut;
                segmenter.params().c = 100.0;
                std::vector<Line2> pureLines(lines.size());
                for (int i = 0; i < lines.size(); i++)
                    pureLines[i] = lines[i].component;
                
                Imagei segmentedImage;
                int segmentsNum;
                std::tie(segmentedImage, segmentsNum) = segmenter(view.image, pureLines, view.image.cols / 100.0);
                


                RLGraph mg;
                RLGraphControls controls;
                RLGraphVars vars;
                std::vector<RegionHandle> rhs;

                int vertVPId = NearestDirectionId(vps, view.camera.upward());
                int hVPId1 = (vertVPId + 1) % 3;
                int hVPId2 = (vertVPId + 2) % 3;

                AppendLines(mg, lines, view.camera, vps, 40.0 / view.camera.focal(), 100.0 / view.camera.focal());
                rhs = AppendRegions(mg, segmentedImage, view.camera, 0.01, 0.001, 1, 3);

                controls = RLGraphControls(mg, vps);


                if (true){
                    AttachWallConstriants(mg, controls, 5 / view.camera.focal(), view.camera.upward());
                }

                // set weights
                // set constraint weights
                SetNecessaryConstraintWeightedAnchors(mg, controls);


                // cc decompose
                auto ccids = MakeHandledTableForAllComponents(mg, -1);
                int ccnum = ConnectedComponents(mg, controls, ccids, [](const RLGraphConstraintControl & c){
                    return c.used && c.weightedAnchors.size() > 0;
                });
                RLGraphOldToNew old2new;
                auto mgs = Decompose(mg, ccids, ccnum, &old2new);
                for (auto & o2n : old2new.container<RegionHandle>()){
                    auto nrh = o2n.second.second;
                    int ccid = o2n.second.first;
                    assert(!nrh.invalid());
                    assert(nrh.id < mgs[ccid].internalComponents<RegionData>().size());
                }
                auto cs = Decompose(mg, controls, ccids, ccnum);
                assert(mgs.size() == cs.size());

                int largestCC = std::max_element(mgs.begin(), mgs.end(), [](const RLGraph & g1, const RLGraph & g2){
                    return g1.internalComponents<RegionData>().size() + g1.internalComponents<LineData>().size()
                        < g2.internalComponents<RegionData>().size() + g2.internalComponents<LineData>().size();
                }) - mgs.begin();

                mg = std::move(mgs[largestCC]);
                controls = std::move(cs[largestCC]);

                std::vector<RegionHandle> newrhs = rhs;
                for (auto & rh : newrhs){
                    if (rh.invalid()){
                        continue;
                    }
                    auto newrh = old2new.at(rh);
                    if (newrh.first == largestCC){
                        rh = newrh.second;
                    }
                    else{
                        rh.reset();
                    }
                }
                rhs = std::move(newrhs);

                // attach depths as anchors
                Imagef realDepth(image.size(), 0.0);
                for (auto it = realDepth.begin(); it != realDepth.end(); ++it){
                    *it = norm(cam_d.toSpace(it.pos()) * data.depth(it.pos())) / 1e5;
                }
                
                SetComponentControl<RegionData>(controls, [&mg, &view, &realDepth](RegionHandle rh, RLGraphComponentControl & c){
                    auto center = mg.data(rh).normalizedCenter;
                    float f = realDepth(ToPixel(view.camera.toScreen(center)));
                    c.weightedAnchors.push_back(WeightAs(center * f, 1.0));
                });

                // remove constraint anchors on occlusions
                SetConstraintControls<RegionBoundaryData>(controls, [&mg, &view, &realDepth](RegionBoundaryHandle bh, RLGraphConstraintControl & c){
                    double maxDepth = MinMaxValOfImage(realDepth).second;
                    std::vector<Weighted<Vec3>> validAnchors;
                    for (auto & wa : c.weightedAnchors){
                        auto p = ToPixel(view.camera.toScreen(wa.component));
                        std::vector<float> ds;
                        for (int xx = -2; xx <= 2; xx++){
                            for (int yy = -2; yy <= 2; yy++){
                                Pixel pp(p.x + xx, p.y + yy);
                                if (Contains(realDepth, pp)){
                                    ds.push_back(realDepth(pp));
                                }
                            }
                        }
                        std::sort(ds.begin(), ds.end());
                        if (abs(ds.front() - ds.back()) <= maxDepth / 20.0){
                            validAnchors.push_back(wa);
                        }
                    }
                    c.weightedAnchors = std::move(validAnchors);
                });

                vars = SolveVariablesWithBoundedAnchors(mg, controls, false, 100);

                Show(view, mg, controls, vars);


            });

        }

    }

}