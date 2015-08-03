#include "../../src/core/basic_types.hpp"
#include "../../src/experimental/rl_graph_solver.hpp"
#include "../../src/experimental/tools.hpp"
#include "../../src/experimental/rl_graph_occlusion.hpp"
#include "../../src/gui/scene.hpp"
#include "../../src/gui/canvas.hpp"
#include "../../src/core/cameras.hpp"
#include "../../src/core/clock.hpp"

#include "routines.hpp"



namespace panolyz {

    namespace NormalIndoor {


        using namespace pano;
        using namespace core;
        using namespace experimental;


        void Show(const PerspectiveView & v, const RLGraph & mg, const RLGraphControls & controls, const RLGraphVars & vars){
            gui::SceneBuilder vis;
            Visualize(vis, v, mg, controls, vars);
            vis.show(true, true, gui::RenderOptions().camera(v.camera).cullBackFace(false).cullFrontFace(false));
        }


        void Run(){

            misc::Matlab matlab;

            bool withGC = false;
            bool withHCons = true;
            bool useWeights = false;
            bool useAllAnchors = false;
            bool outdoor = true;

            std::string path;
            //path = PROJECT_TEST_DATA_DIR_STR"/normal/room.png"; // bingo
            //path = PROJECT_TEST_DATA_DIR_STR"/normal/room4.jpg"; // bingo
            //path = PROJECT_TEST_DATA_DIR_STR"/normal/room7.jpg"; // bingo!!!
            //path = PROJECT_TEST_DATA_DIR_STR"/normal/room10.jpg"; // bingo!
            //path = PROJECT_TEST_DATA_DIR_STR"/normal/room12.png"; // bingo!
            //path = PROJECT_TEST_DATA_DIR_STR"/normal/room13.jpg"; // bingo
            //path = PROJECT_TEST_DATA_DIR_STR"/normal/room14.jpg"; // bingo
            //path = PROJECT_TEST_DATA_DIR_STR"/normal/room21.jpg"; // bingo! ?
            //path = PROJECT_TEST_DATA_DIR_STR"/normal/room22.jpg"; // bingo
            //path = PROJECT_TEST_DATA_DIR_STR"/normal/room23.jpg"; // bingo

            //path = PROJECT_TEST_DATA_DIR_STR"/normal/ng1.png";
            //path = PROJECT_TEST_DATA_DIR_STR"/normal/room2e.jpg";
            //path = PROJECT_TEST_DATA_DIR_STR"/normal/room3.jpg";
            //path = PROJECT_TEST_DATA_DIR_STR"/normal/room5.jpg"; //!
            //path = PROJECT_TEST_DATA_DIR_STR"/normal/room6.jpg";
            //path = PROJECT_TEST_DATA_DIR_STR"/normal/room8.jpg"; //!
            //path = PROJECT_TEST_DATA_DIR_STR"/normal/room11.jpg"; // vp failed
            //path = PROJECT_TEST_DATA_DIR_STR"/normal/room15.jpg"; // vp failed
            //path = PROJECT_TEST_DATA_DIR_STR"/normal/room16.jpg";
            //path = PROJECT_TEST_DATA_DIR_STR"/normal/room17.jpg";
            //path = PROJECT_TEST_DATA_DIR_STR"/normal/room18.jpg";
            //path = PROJECT_TEST_DATA_DIR_STR"/normal/room19.jpg"; // flat...
            //path = PROJECT_TEST_DATA_DIR_STR"/normal/room20.jpg";


            //path = "F:\\DataSets\\yanghao.thumanhattan\\IMG_3490.jpg";
            //path = "F:\\DataSets\\yanghao.thumanhattan\\IMG_3489.jpg";
            //path = "F:\\DataSets\\yanghao.thumanhattan\\IMG_3567.jpg";
            //path = "F:\\DataSets\\yanghao.thumanhattan\\IMG_3512.jpg"; // bingo!
            path = "F:\\DataSets\\yanghao.thumanhattan\\IMG_3503.jpg"; // bingo!
            //path = "F:\\DataSets\\yanghao.thumanhattan\\IMG_3554.jpg";
            //path = "F:\\DataSets\\yanghao.thumanhattan\\IMG_3502.jpg";
            //path = "F:\\DataSets\\yanghao.thumanhattan\\IMG_3558.jpg";
            //path = "F:\\DataSets\\yanghao.thumanhattan\\IMG_3546.jpg"; // bin...
            //path = "F:\\DataSets\\yanghao.thumanhattan\\IMG_3487.jpg";
            //path = "F:\\DataSets\\yanghao.thumanhattan\\IMG_3571.jpg";
            //path = "F:\\DataSets\\yanghao.thumanhattan\\IMG_3499.jpg"; // bingo
            //path = "F:\\DataSets\\yanghao.thumanhattan\\IMG_3500.jpg";
            //path = "F:\\DataSets\\yanghao.thumanhattan\\IMG_3487.jpg";
            //path = "F:\\DataSets\\yanghao.thumanhattan\\IMG_3568.jpg";
            //path = "F:\\DataSets\\yanghao.thumanhattan\\IMG_3488.jpg"; // bingo!
            //path = "F:\\DataSets\\yanghao.thumanhattan\\IMG_3557.jpg";
            //path = "F:\\DataSets\\yanghao.thumanhattan\\IMG_3480.jpg";
            //path = "F:\\DataSets\\yanghao.thumanhattan\\IMG_3485.jpg";
            //path = "F:\\DataSets\\yanghao.thumanhattan\\IMG_3539.jpg";
            //path = "F:\\DataSets\\yanghao.thumanhattan\\IMG_3534.jpg";
            //path = "F:\\DataSets\\yanghao.thumanhattan\\IMG_3492.jpg";

            Image image = cv::imread(path);
            ResizeToHeight(image, 400);


            View<PerspectiveCamera> view;

            std::vector<Classified<Line2>> lines;
            std::vector<Vec3> vps;

            Imagei segmentedImage;

            bool RE_BUILD_FEATURES = true;

            if (RE_BUILD_FEATURES){

                double focal;
                std::vector<Classified<Line3>> line3s;

                VanishingPointsDetector::Params vpdParams(VanishingPointsDetector::TardifSimplified);
                view = CreatePerspectiveView(image, Point3(0, 0, 0), Point3(1, 0, 0), Point3(0, 0, -1),
                    LineSegmentExtractor(), VanishingPointsDetector(vpdParams), &line3s, &lines, &vps, &focal).unwrap();

                {
                    Clock clock("omap");
                    std::vector<HPoint2> hvps(vps.size());
                    for (int i = 0; i < vps.size(); i++){
                        hvps[i] = view.camera.toScreenInHPoint(vps[i]);
                    }
                    auto om = ComputeOrientationMaps(lines, hvps, view.image.size());
                    gui::ColorTable rgb = gui::ColorTableDescriptor::RGBGreys;
                    cv::imshow("om", rgb(om));
                    cv::waitKey();
                }

                SegmentationExtractor segmenter;
                segmenter.params().algorithm = SegmentationExtractor::GraphCut;
                segmenter.params().sigma = 10.0;
                segmenter.params().c = 1.0;
                segmenter.params().superpixelSizeSuggestion = 2000;
                int segmentsNum = 0;

                std::vector<Line2> line2s;
                for (auto & l : lines){
                    line2s.push_back(l.component);
                }

                std::tie(segmentedImage, segmentsNum) = segmenter(view.image, line2s, 5.0);

                if (0){
                    auto ctable = gui::CreateRandomColorTableWithSize(segmentsNum);
                    gui::AsCanvas(ctable(segmentedImage)).show();
                    gui::AsCanvas(ctable(segmentedImage)).show();
                }

                Save(path, "pre", view, lines, vps, segmentedImage);
            }
            else{
                Load(path, "pre", view, lines, vps, segmentedImage);
            }


            // gc
            ImageOf<Vec<double, 5>> gc;
            if (RE_BUILD_FEATURES){
                gc = ComputeGeometricContext(matlab, view.image, outdoor);
                Save(path, "gc", gc);
            }
            else{
                Load(path, "gc", gc);
            }



            int vertVPId = NearestDirectionId(vps, view.camera.upward());
            int hVPId1 = (vertVPId + 1) % 3;
            int hVPId2 = (vertVPId + 2) % 3;


            RLGraph mg;
            RLGraphControls controls;
            RLGraphVars vars;
            std::vector<RegionHandle> rhs;


            AppendLines(mg, lines, view.camera, vps, 40.0 / view.camera.focal(), 100.0 / view.camera.focal());
            rhs = AppendRegions(mg, segmentedImage, view.camera, 0.01, 0.001, 1, 3);

            controls = RLGraphControls(mg, vps);



            if (withHCons){
                AttachWallConstriants(mg, controls, 5 / view.camera.focal(), view.camera.upward());
            }
           /* AttachGeometricContextConstraints(mg, controls, view.camera, gc,
                [outdoor, vertVPId, withGC, withHCons](RLGraphComponentControl & c, const Vec<double, 5> & v, double s){
                int maxlabel = std::max_element(v.val, v.val + 5) - v.val;
                if (outdoor && s > 2){
                    switch (maxlabel){
                    case GeometricContextEstimator::OI_Ground:
                        if (withGC){
                            c.orientationClaz = vertVPId;
                            c.orientationNotClaz = -1;
                        }
                        break;
                    case GeometricContextEstimator::OI_VerticalPlanarFace:
                        if (withGC && !withHCons){
                            c.orientationClaz = -1;
                            c.orientationNotClaz = vertVPId;
                        }
                        break;
                    case GeometricContextEstimator::OI_Clutter:
                        c.used = false;
                        break;
                    case GeometricContextEstimator::OI_Porous:
                        c.used = false;
                        break;
                    case GeometricContextEstimator::OI_Sky:
                        c.used = false;
                        break;
                    default:
                        break;
                    }
                }
                else if (!outdoor && s > 3){
                    switch (maxlabel){
                    case GeometricContextEstimator::II_VerticalPlanarFace:
                        if (withGC && !withHCons){
                            c.orientationClaz = -1;
                            c.orientationNotClaz = vertVPId;
                        }
                        break;
                    case GeometricContextEstimator::II_HorizontalPlanarFace:
                        if (withGC){
                            c.orientationClaz = vertVPId;
                        }
                        break;
                    case GeometricContextEstimator::II_Clutter:
                        if (withGC && !withHCons){
                            c.orientationClaz = c.orientationNotClaz = -1;
                        }
                        break;
                    default:
                        break;
                    }
                }
            });*/


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


            if (!AttachWeightedAnchorToCenterOfLargestLineIfNoExists(mg, controls))
                return;

            ResetToSampledArmorAnchors(mg, controls, 2.0 / view.camera.focal());
            vars = MakeVariables(mg, controls);

            OptimizeVariablesWithBoundedAnchors(matlab, mg, controls, vars, useWeights, 1000, 10,
                [&view, &mg, &controls](const RLGraphVars & vs){
                Show(view, mg, controls, vs);
                return true;
            });
            assert(!HasInfOrNaNValue(vars));


            //vars = SolveVariablesWithBoundedAnchors(mg, controls, config.useWeights/*false*/,
            //    config.useAllAnchors);


            //Show(view, mg, controls, vars);

            //LooseOrientationConstraintsOnComponents(mg, controls, vars, 0.0);
            //vars = SolveVariablesBoundComponentAnchors(mg, controls, config.useWeights/*false*/,
            //    config.useAllAnchors);
            //{
            //    // save
            //    auto planes = Instances<RegionData>(mg, controls, vars);
            //    Image3f rlomap = Image3f::zeros(view.image.size());
            //    for (auto it = rlomap.begin(); it != rlomap.end(); ++it){
            //        int rid = segmentedImage(it.pos());
            //        auto rh = rhs[rid];
            //        if (rh.invalid()){
            //            continue;
            //        }
            //        auto & plane = planes[rh];
            //        for (int i = 0; i < 3; i++){
            //            (*it)[i] = abs(normalize(plane.normal).dot(normalize(vps[i])));
            //        }
            //    }
            //    matlab.PutVariable("rlomap_afterloose", rlomap);
            //}

            //matlab << ("save('" + p.fileprefix + tag + ".mat', 'rlomap');");
            Show(view, mg, controls, vars);


        }


    }


}