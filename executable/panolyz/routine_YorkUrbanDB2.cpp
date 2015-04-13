#include "../../src/core/basic_types.hpp"
#include "../../src/experimental/rl_graph.hpp"
#include "../../src/experimental/rl_graph_occlusion.hpp"
#include "../../src/experimental/tools.hpp"
#include "../../src/gui/visualize2d.hpp"
#include "../../src/gui/visualizers.hpp"
#include "../../src/misc/matlab.hpp"
#include "../../src/gui/singleton.hpp"

#include "tools.hpp"
#include "routines.hpp"

namespace panolyz {

    namespace YorkUrbanDB2 {

        using namespace panoramix;
        using namespace core;
        using namespace experimental;

        namespace gt {
            static const double focal = 6.0532;
            static const double pixelSize = 0.0090;
            static const double focalReal = focal / pixelSize;
            static const core::Point2 pp(307.5513, 251.4542);

            struct Data {
                bool isOutdoor;
                std::vector<Vec3> vps;
                std::vector<Classified<Line2>> lines;
                ImageOfType<Vec<double, 7>> gcOutdoor, gcIndoor, gcIndoor2;
            };

            Data LoadVPs(const std::string & imagepath){
                std::string base = imagepath.substr(0, imagepath.size() - 4);
                std::string linespath = base + "LinesAndVP.mat";
                std::string vporthopath = base + "GroundTruthVP_Orthogonal_CamParams.mat";
                std::string scenetype = base + "scene_type.mat";
                std::string gcpath = base + "gc.mat";

                Data gt;
                misc::Matlab matlab;
                matlab << ("load(\'" + linespath + "\')");

                DenseMatd linepts;
                matlab.GetVariable("lines", linepts);
                std::vector<double> vp_association;
                matlab << "vp_association = vp_association';";
                matlab.GetVariable("vp_association", vp_association);

                assert(linepts.rows == vp_association.size() * 2);

                gt.lines.resize(vp_association.size());
                for (int i = 0; i < vp_association.size(); i++){
                    gt.lines[i].claz = vp_association[i] - 1;
                    gt.lines[i].component.first = Point2(linepts(i * 2, 0), linepts(i * 2, 1));
                    gt.lines[i].component.second = Point2(linepts(i * 2 + 1, 0), linepts(i * 2 + 1, 1));
                }

                DenseMatd vps;
                matlab << ("load('" + vporthopath + "');");
                matlab.GetVariable("vp_orthogonal", vps);
                gt.vps = {
                    { vps(0, 0), vps(1, 0), vps(2, 0) },
                    { vps(0, 1), vps(1, 1), vps(2, 1) },
                    { vps(0, 2), vps(1, 2), vps(2, 2) }
                };

                double outdoor;
                matlab << ("load('" + scenetype + "');");
                matlab << "isOutdoor = (scene_type == 'o');";
                matlab.GetVariable("isOutdoor", outdoor);
                gt.isOutdoor = outdoor != 0;

                matlab.Load(gcpath);
                matlab.GetVariable("gc_outdoor", gt.gcOutdoor);
                matlab.GetVariable("gc_indoor", gt.gcIndoor);
                matlab.GetVariable("gc_indoor2", gt.gcIndoor2);

                assert(gt.gcIndoor.size() == gt.gcOutdoor.size());

                return gt;
            }
        };


        struct Preparation {
            View<PerspectiveCamera> view;
            std::vector<Classified<Line2>> lines;
            std::vector<Vec3> vps;
            Imagei segmentedImage;
            int segmentsNum;
            ImageOfType<Vec<double, 5>> gc;

            bool outdoor;

            void showVPLines() const {
                std::vector<core::Classified<core::Ray2>> vpRays;
                for (int i = 0; i < 3; i++){
                    std::cout << "vp[" << i << "] = " << view.camera.screenProjection(vps[i]) << std::endl;
                    for (double a = 0; a <= M_PI * 2.0; a += 0.1){
                        core::Point2 p = core::Point2(view.image.cols / 2, view.image.rows / 2) + core::Vec2(cos(a), sin(a)) * 1000.0;
                        vpRays.push_back(core::ClassifyAs(core::Ray2(p, (view.camera.screenProjectionInHPoint(vps[i]) - core::HPoint2(p, 1.0)).numerator), i));
                    }
                }
                gui::Visualizer2D(view.image)
                    << gui::manip2d::SetColorTable(gui::ColorTable(gui::ColorTableDescriptor::RGB).appendRandomizedGreyColors(vps.size() - 3))
                    << gui::manip2d::SetThickness(1)
                    << vpRays
                    << gui::manip2d::SetThickness(2)
                    << lines
                    << gui::manip2d::Show();
            }
        };

        Preparation Prepare(const std::string & name){
            auto path = "F:\\DataSets\\YorkUrbanDB\\data\\" + name + "\\" + name + ".jpg";

            View<PerspectiveCamera> view;
            std::vector<Classified<Line2>> lines;
            std::vector<Vec3> vps;
            Imagei segmentedImage;
            int segmentsNum;
            ImageOfType<Vec<double, 5>> gc;

            bool outdoor;


            auto image = cv::imread(path);
            assert(!image.empty());
            view.image = image;

            // load gt            
            auto data = gt::LoadVPs(path);
            lines = std::move(data.lines);
            vps = std::move(data.vps);
            outdoor = data.isOutdoor;
            GeometricContextEstimator gcEstimator;
            gc = gcEstimator(outdoor ? data.gcOutdoor : data.gcIndoor2,
                outdoor ? SceneClass::Outdoor : SceneClass::Indoor);

            view.camera = core::PerspectiveCamera(view.image.cols, view.image.rows, gt::pp, gt::focalReal,
                core::Point3(0, 0, 0), core::Point3(0, 0, -1), core::Point3(0, -1, 0));

            if (0) {
                std::vector<core::Classified<core::Ray2>> vpRays;
                for (int i = 0; i < 3; i++){
                    std::cout << "vp[" << i << "] = " << view.camera.screenProjection(vps[i]) << std::endl;
                    for (double a = 0; a <= M_PI * 2.0; a += 0.1){
                        core::Point2 p = core::Point2(image.cols / 2, image.rows / 2) + core::Vec2(cos(a), sin(a)) * 1000.0;
                        vpRays.push_back(core::ClassifyAs(core::Ray2(p, (view.camera.screenProjectionInHPoint(vps[i]) - core::HPoint2(p, 1.0)).numerator), i));
                    }
                }
                gui::Visualizer2D(image)
                    << gui::manip2d::SetColorTable(gui::ColorTable(gui::ColorTableDescriptor::RGB).appendRandomizedGreyColors(vps.size() - 3))
                    << gui::manip2d::SetThickness(1)
                    << vpRays
                    << gui::manip2d::SetThickness(2)
                    << lines
                    << gui::manip2d::Show();
            }

            // append regions
            SegmentationExtractor segmenter;
            segmenter.params().algorithm = SegmentationExtractor::GraphCut;
            segmenter.params().c = 100.0;
            std::vector<Line2> pureLines(lines.size());
            for (int i = 0; i < lines.size(); i++)
                pureLines[i] = lines[i].component;
            std::tie(segmentedImage, segmentsNum) = segmenter(view.image, pureLines, view.image.cols / 100.0);

            return{ view, lines, vps, segmentedImage, segmentsNum, gc, outdoor };
        }


        void DepthOrdering(const Preparation & p){
            auto & view = p.view;
            auto & lines = p.lines;
            auto & vps = p.vps;
            auto & segmentedImage = p.segmentedImage;
            int segmentsNum = p.segmentsNum;
            auto & gc = p.gc;
            bool outdoor = p.outdoor;

            //auto edges = FindContoursOfRegionsAndBoundaries()


            RLGraph mg;
            RLGraphControls controls;
            RLGraphVars vars;
            std::vector<RegionHandle> rhs;

            int vertVPId = NearestDirectionId(vps, view.camera.upward());
            AppendLines(mg, lines, view.camera, vps, 40.0 / view.camera.focal(), 100.0 / view.camera.focal());
            rhs = AppendRegions(mg, segmentedImage, view.camera, 0.01, 0.001, 1, 2, true);

            if (1){
                // show line and boundary connections
                auto ctable = gui::CreateRandomColorTableWithSize(mg.internalComponents<RegionData>().size());
                auto rlim = Print(mg, segmentedImage, view.camera, rhs, [&ctable](RegionHandle rh){
                    return ctable[rh.id];
                }, [](LineHandle lh){
                    return gui::White;
                }, [](RegionBoundaryHandle bh){
                    return gui::Black;
                }, 1, 2);
                for (auto & rl : mg.constraints<RegionLineConnectionData>()){
                    RegionHandle rh = rl.topo.component<0>();
                    LineHandle lh = rl.topo.component<1>();
                    auto rc = view.camera.screenProjection(mg.data(rh).normalizedCenter);
                    for (auto & a : { rl.data.normalizedAnchors.front(), rl.data.normalizedAnchors.back() }){
                        auto p = ToPixelLoc(view.camera.screenProjection(a));
                        cv::line(rlim, ToPixelLoc(rc), p, gui::Color(gui::Gray));
                    }
                }
                cv::imshow("", rlim);
                cv::waitKey();
            }

            if(1){
                auto ctable = gui::CreateGreyColorTableWithSize(mg.internalComponents<LineData>().size());
                Image3i im = Image3i::zeros(view.image.size());
                for (auto & l : mg.components<LineData>()){
                    auto p1 = view.camera.screenProjection(l.data.line.first);
                    auto p2 = view.camera.screenProjection(l.data.line.second);
                    cv::line(im, ToPixelLoc(p1), ToPixelLoc(p2), ctable[l.topo.hd.id]);
                }

                auto ldr = GuessLineDepthRelation(mg, 450 / view.camera.focal());
                for (auto & lr : mg.constraints<LineRelationData>()){
                    auto & lh1 = lr.topo.component<0>();
                    auto & lh2 = lr.topo.component<1>();
                    auto & ld1 = mg.data(lh1);
                    auto & ld2 = mg.data(lh2);
                    Vec3 eq1 = ld1.line.first.cross(ld1.line.second);
                    Vec3 eq2 = ld2.line.first.cross(ld2.line.second);
                    Vec3 inter = normalize(eq1.cross(eq2));
                    if (ld1.initialClaz == -1 && ld2.initialClaz == -1){
                        continue;
                    }
                    auto p = ToPixelLoc(view.camera.screenProjection(inter));
                    auto r = ldr[lr.topo.hd];
                    switch (r) {
                    case panoramix::experimental::DepthRelationGuess::FirstMaybeCloser:
                        cv::circle(im, p, 2, cv::Scalar(255, 0, 0), 3);
                        break;
                    case panoramix::experimental::DepthRelationGuess::SecondMaybeCloser:
                        cv::circle(im, p, 2, cv::Scalar(0, 0, 255), 3);
                        break;
                    default:
                        break;
                    }
                }

                cv::imshow("first closer", im * 255);
                cv::waitKey();
             }




        }
        
        void ReconstructLines(const Preparation & p){
            auto & view = p.view;
            auto & lines = p.lines;
            auto & vps = p.vps;
            auto & segmentedImage = p.segmentedImage;
            int segmentsNum = p.segmentsNum;
            auto & gc = p.gc;

            bool outdoor = p.outdoor;

            // consider only lines
            RLGraph mg;
            RLGraphControls controls;
            RLGraphVars vars;

            AppendLines(mg, lines, view.camera, vps, 40.0 / gt::focalReal, 100.0 / gt::focalReal);

            controls = RLGraphControls(mg, vps);
            SetConstraintWeights<LineRelationData>(controls, [&mg](LineRelationHandle h){
                return std::max(mg.data(h).junctionWeight, 1e-1f);
            });

            auto ccids = MakeHandledTableForAllComponents(mg, -1);
            int ccnum = ConnectedComponents(mg, controls, ccids, [](const RLGraphConstraintControl & c){
                return c.used && c.weight > 0;
            });
            auto mgs = Decompose(mg, ccids, ccnum);
            auto cs = Decompose(mg, controls, ccids, ccnum);
            assert(mgs.size() == cs.size());

            for (int i = 0; i < ccnum; i++){
                auto & mg = mgs[i];
                auto & controls = cs[i];
                AttachAnchorToCenterOfLargestLineIfNoAnchorExists(mg, controls);
                auto vars = SolveVariables(mgs[i], controls, true);
                NormalizeVariables(mg, controls, vars);

                gui::Visualizer vis;
                Visualize(vis, view, mg, controls, vars);
                vis.camera(view.camera);
                vis.show(true, true);
            }

        }

        void ReconstructModel(const Preparation & p){
            auto & view = p.view;
            auto & lines = p.lines;
            auto & vps = p.vps;
            auto & segmentedImage = p.segmentedImage;
            int segmentsNum = p.segmentsNum;
            auto & gc = p.gc;

            bool outdoor = p.outdoor;

            RLGraph mg;
            RLGraphControls controls;
            RLGraphVars vars;
            std::vector<RegionHandle> rhs;

            int vertVPId = NearestDirectionId(vps, view.camera.upward());
            AppendLines(mg, lines, view.camera, vps, 40.0 / view.camera.focal(), 100.0 / view.camera.focal());
            rhs = AppendRegions(mg, segmentedImage, view.camera, 0.01, 0.001, 1, 3, true);

            controls = RLGraphControls(mg, vps);

            AttachGeometricContextConstraints(mg, controls, view.camera, gc,
                [outdoor, vertVPId](RLGraphComponentControl & c, const Vec<double, 5> & v, double s){
                int maxlabel = std::max_element(v.val, v.val + 5) - v.val;
                if (outdoor && s > 2){
                    switch (maxlabel){
                    case GeometricContextEstimator::OI_Ground:
                        c.orientationClaz = vertVPId;
                        c.orientationNotClaz = -1;
                        break;
                    case GeometricContextEstimator::OI_VerticalPlanarFace:
                        c.orientationClaz = -1;
                        c.orientationNotClaz = vertVPId;
                        break;
                    case GeometricContextEstimator::OI_Clutter:
                        c.orientationClaz = -1;
                        c.orientationNotClaz = -1;
                        break;
                    case GeometricContextEstimator::OI_Porous:
                    case GeometricContextEstimator::OI_Sky:
                        c.used = false;
                        break;
                    default:
                        break;
                    }
                }
                else if (!outdoor && s > 4){
                    switch (maxlabel){
                    case GeometricContextEstimator::II_FrontVerticalPlanarFace:
                    case GeometricContextEstimator::II_SideVerticalPlanarFace:
                        c.orientationClaz = -1;
                        c.orientationNotClaz = vertVPId;
                        break;
                    case GeometricContextEstimator::II_HorizontalPlanarFace:
                        c.orientationClaz = vertVPId;
                        break;
                    case GeometricContextEstimator::II_Clutter:
                        c.orientationClaz = c.orientationNotClaz = -1;
                        break;
                    default:
                        break;
                    }
                }
            });

            //AttachPrincipleDirectionConstraints(mg, controls);
            //AttachWallConstriants(mg, controls);

            // set constraint weights
            SetConstraintWeights<LineRelationData>(controls, [&mg](LineRelationHandle h){
                return std::max(mg.data(h).junctionWeight * 10, 1.0f);
            });
            SetConstraintWeights<RegionBoundaryData>(controls, [&mg](RegionBoundaryHandle h){
                return std::max(mg.data(h).length / M_PI * 10.0, 1.0);
            });
            SetConstraintWeights<RegionLineConnectionData>(controls, [&mg](RegionLineConnectionHandle h){
                return std::max(mg.data(h).length / M_PI * 10.0, 1.0);
            });

            // cc decompose
            auto ccids = MakeHandledTableForAllComponents(mg, -1);
            int ccnum = ConnectedComponents(mg, controls, ccids, [](const RLGraphConstraintControl & c){
                return c.used && c.weight > 0;
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

            // for each cc
            for (int i = 0; i < ccnum; i++){
                auto & mg = mgs[i];
                auto & controls = cs[i];
                RLGraphVars vars;

                if (!AttachAnchorToCenterOfLargestLineIfNoAnchorExists(mg, controls))
                    continue;

                std::vector<RegionHandle> newrhs = rhs;
                for (auto & rh : newrhs){
                    if (rh.invalid()){
                        continue;
                    }
                    auto newrh = old2new.at(rh);
                    if (newrh.first == i){
                        rh = newrh.second;
                    }
                    else{
                        rh.reset();
                    }
                }

                // show who is anchored?
                if (0){
                    auto printed = Print(mg, segmentedImage, view.camera, newrhs, [&controls](RegionHandle rh){
                        if (!controls[rh].used){
                            return gui::Transparent;
                        }
                        if (controls[rh].weightedAnchors.size() > 0){
                            return gui::Blue;
                        }
                        return gui::Red;
                    }, [&controls](LineHandle rh){
                        if (!controls[rh].used){
                            return gui::Transparent;
                        }
                        if (controls[rh].weightedAnchors.size() > 0){
                            return gui::Black;
                        }
                        return gui::White;
                    }, [&controls](RegionBoundaryHandle bh){
                        if (!controls[bh].used){
                            return gui::Transparent;
                        }
                        return gui::Yellow;
                    });
                    cv::imshow("who is anchored?", printed);
                    cv::waitKey();
                }

                { // with inversed depth setup
                    vars = SolveVariables(mg, controls, true, true);
                    NormalizeVariables(mg, controls, vars);
                    LooseOrientationConstraintsOnComponents(mg, controls, vars, 0.0, 0.1);
                    vars = SolveVariables(mg, controls, true, true);
                    NormalizeVariables(mg, controls, vars);
                    if (1){
                        gui::Visualizer vis("inversed depth setup");
                        Visualize(vis, view, mg, controls, vars);
                        vis.camera(view.camera);
                        vis.renderOptions.cullBackFace = vis.renderOptions.cullFrontFace = false;
                        vis.show(true, true);
                    }
                    
                    return;

                    ClearAllComponentAnchors(controls);
                    OptimizeVariables(mg, controls, vars, true, false);
                    NormalizeVariables(mg, controls, vars);
                    if (1){
                        gui::Visualizer vis("inversed depth setup -> optimized");
                        Visualize(vis, view, mg, controls, vars);
                        vis.camera(view.camera);
                        vis.renderOptions.cullBackFace = vis.renderOptions.cullFrontFace = false;
                        vis.show(true, true);
                    }
                }

                {  // without inversed depth setup
                    vars = MakeVariables(mg, controls);
                    OptimizeVariables(mg, controls, vars, true, true);
                    NormalizeVariables(mg, controls, vars);
                    if (1){
                        gui::Visualizer vis("directly optimized");
                        Visualize(vis, view, mg, controls, vars);
                        vis.camera(view.camera);
                        vis.renderOptions.cullBackFace = vis.renderOptions.cullFrontFace = false;
                        vis.show(true, true);
                    }
                }

            }
        }


        

        void Show(const Preparation & p, const RLGraph & mg, const RLGraphControls & controls, const RLGraphVars & vars){
            gui::Visualizer vis("inversed depth setup");
            Visualize(vis, p.view, mg, controls, vars);
            vis.camera(p.view.camera);
            vis.renderOptions.cullBackFace = vis.renderOptions.cullFrontFace = false;
            vis.show(true, true);
        }


        void LabelGT(const Preparation & p){
            auto & view = p.view;
            auto & lines = p.lines;
            auto & vps = p.vps;
            auto & segmentedImage = p.segmentedImage;
            int segmentsNum = p.segmentsNum;
            auto & gc = p.gc;

            auto bps = FindContoursOfRegionsAndBoundaries(segmentedImage, 1, false);
            std::vector<std::string> regionLabelNames = {
                "Free Plane",
                "Toward VP1",
                "Toward VP2",
                "Toward VP3",
                "Along VP1",
                "Along VP2",
                "Along VP3",
                "Void"
            };
            std::vector<gui::Color> regionLabelColors = {
                gui::LightGray,
                gui::Red,
                gui::Yellow,
                gui::Blue,
                gui::Magenta,
                gui::Cyan,
                gui::Orange,
                gui::DarkGray
            };
            std::vector<std::string> boundaryLabelNames = {
                "Connected",
                "Disconnected"
            };
            std::vector<gui::Color> boundaryLabelColors = {
                gui::White,
                gui::Black
            };

            RLGraph mg;
            RLGraphControls controls;
            RLGraphVars vars;
            std::vector<RegionHandle> rhs;

            int vertVPId = NearestDirectionId(vps, view.camera.upward());
            AppendLines(mg, lines, view.camera, vps, 40.0 / view.camera.focal(), 100.0 / view.camera.focal());
            rhs = AppendRegions(mg, segmentedImage, view.camera, 0.01, 0.001, 1, 3, true);

            controls = RLGraphControls(mg, vps);

            p.showVPLines();

            Labels labels;

            while (true){
                LabelIt(labels, view.image, segmentedImage, bps, 
                    regionLabelNames, boundaryLabelNames, regionLabelColors, boundaryLabelColors);

                // apply labels to controls
                // regions
                for (int i = 0; i < labels.regionLabels.size(); i++){
                    RegionHandle rh = rhs[i];
                    auto & control = controls[rh];
                    switch (labels.regionLabels[i]){
                    case 0: control.used = true; control.orientationClaz = -1; control.orientationNotClaz = -1; break;
                    case 1: control.used = true; control.orientationClaz = 0; control.orientationNotClaz = -1; break;
                    case 2: control.used = true; control.orientationClaz = 1; control.orientationNotClaz = -1; break;
                    case 3: control.used = true; control.orientationClaz = 2; control.orientationNotClaz = -1; break;
                    case 4: control.used = true; control.orientationClaz = -1; control.orientationNotClaz = 0; break;
                    case 5: control.used = true; control.orientationClaz = -1; control.orientationNotClaz = 1; break;
                    case 6: control.used = true; control.orientationClaz = -1; control.orientationNotClaz = 2; break;
                    case 7: control.used = false; control.orientationClaz = -1; control.orientationNotClaz = -1; break;
                    }
                }
                // boundaries, line-region relations
                // todo

                controls.disableAllInvalidConstraints(mg);

                // cc decompose
                auto ccids = MakeHandledTableForAllComponents(mg, -1);
                int ccnum = ConnectedComponents(mg, controls, ccids, [](const RLGraphConstraintControl & c){
                    return c.used && c.weight > 0;
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

                // for each cc
                for (int i = 0; i < ccnum; i++){
                    auto & mg = mgs[i];
                    auto & controls = cs[i];
                    RLGraphVars vars;

                    if (!AttachAnchorToCenterOfLargestLineIfNoAnchorExists(mg, controls))
                        continue;

                    std::vector<RegionHandle> newrhs = rhs;
                    for (auto & rh : newrhs){
                        if (rh.invalid()){
                            continue;
                        }
                        auto newrh = old2new.at(rh);
                        if (newrh.first == i){
                            rh = newrh.second;
                        }
                        else{
                            rh.reset();
                        }
                    }


                    // initial status
                    vars = SolveVariables(mg, controls, false, true);
                    NormalizeVariables(mg, controls, vars);
                    if (1){
                        Show(p, mg, controls, vars);
                    }
                }

            }

        }


        void ExhausticReconstruct(const Preparation & p) {

            auto & view = p.view;
            auto & lines = p.lines;
            auto & vps = p.vps;
            auto & segmentedImage = p.segmentedImage;
            int segmentsNum = p.segmentsNum;
            auto & gc = p.gc;

            bool outdoor = p.outdoor;

            RLGraph mg;
            RLGraphControls controls;
            RLGraphVars vars;
            std::vector<RegionHandle> rhs;

            int vertVPId = NearestDirectionId(vps, view.camera.upward());
            AppendLines(mg, lines, view.camera, vps, 40.0 / view.camera.focal(), 100.0 / view.camera.focal());
            rhs = AppendRegions(mg, segmentedImage, view.camera, 0.01, 0.001, 1, 3, true);

            controls = RLGraphControls(mg, vps);

            AttachGeometricContextConstraints(mg, controls, view.camera, gc,
                [outdoor, vertVPId](RLGraphComponentControl & c, const Vec<double, 5> & v, double s){
                int maxlabel = std::max_element(v.val, v.val + 5) - v.val;
                if (outdoor && s > 2){
                    switch (maxlabel){
                    case GeometricContextEstimator::OI_Ground:
                        c.orientationClaz = vertVPId;
                        c.orientationNotClaz = -1;
                        break;
                    case GeometricContextEstimator::OI_VerticalPlanarFace:
                        c.orientationClaz = -1;
                        c.orientationNotClaz = vertVPId;
                        break;
                    case GeometricContextEstimator::OI_Clutter:
                        c.orientationClaz = -1;
                        c.orientationNotClaz = -1;
                        break;
                    case GeometricContextEstimator::OI_Porous:
                    case GeometricContextEstimator::OI_Sky:
                        c.used = false;
                        break;
                    default:
                        break;
                    }
                }
                else if (!outdoor && s > 4){
                    switch (maxlabel){
                    case GeometricContextEstimator::II_FrontVerticalPlanarFace:
                    case GeometricContextEstimator::II_SideVerticalPlanarFace:
                        c.orientationClaz = -1;
                        c.orientationNotClaz = vertVPId;
                        break;
                    case GeometricContextEstimator::II_HorizontalPlanarFace:
                        c.orientationClaz = vertVPId;
                        break;
                    case GeometricContextEstimator::II_Clutter:
                        c.orientationClaz = c.orientationNotClaz = -1;
                        break;
                    default:
                        break;
                    }
                }
            });


            // cc decompose
            auto ccids = MakeHandledTableForAllComponents(mg, -1);
            int ccnum = ConnectedComponents(mg, controls, ccids, [](const RLGraphConstraintControl & c){
                return c.used && c.weight > 0;
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


            // for each cc
            for (int i = 0; i < ccnum; i++){
                auto & mg = mgs[i];
                auto & controls = cs[i];
                RLGraphVars vars;

                if (!AttachAnchorToCenterOfLargestLineIfNoAnchorExists(mg, controls))
                    continue;

                std::vector<RegionHandle> newrhs = rhs;
                for (auto & rh : newrhs){
                    if (rh.invalid()){
                        continue;
                    }
                    auto newrh = old2new.at(rh);
                    if (newrh.first == i){
                        rh = newrh.second;
                    }
                    else{
                        rh.reset();
                    }
                }


                // initial status
                vars = SolveVariables(mg, controls, false, true);
                NormalizeVariables(mg, controls, vars);
                if (1){
                    Show(p, mg, controls, vars);
                }


                //
                while (true){
                    std::queue<RLGraphControls> Q;
                    Q.push(controls);

                    auto & bestControls = Q.front();
                    
                    // search in branches
                    break;

                }
                
            }
        }


        void Run(){

            std::string name;
            name = "p1020171";

            misc::Matlab matlab;
            matlab.Load("F:\\DataSets\\YorkUrbanDB\\data\\names.mat");
            double dnum = 0;
            matlab << "num = length(names);";
            matlab.GetVariable("num", dnum);
            int num = dnum;

            for (int i = 0; i < num; i++){
                matlab << ("name = names{" + std::to_string(i + 1) + "};");
                std::string name;
                matlab.GetVariable("name", name);
                std::cout << "processing " << name << std::endl;

                if (name != "p1080104"){
                    continue;
                }

                auto p = Prepare(name);
                
                LabelGT(p);

                //DepthOrdering(p);
                //ReconstructLines(p);
                //ReconstructModel(p);
                //ExhausticReconstruct(p);
            }

        }
    }

}