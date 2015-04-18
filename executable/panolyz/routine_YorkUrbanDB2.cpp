#include <thread>

#include "../../src/core/basic_types.hpp"
#include "../../src/experimental/rl_graph.hpp"
#include "../../src/experimental/rl_graph_occlusion.hpp"
#include "../../src/experimental/tools.hpp"
#include "../../src/gui/visualize2d.hpp"
#include "../../src/gui/visualizers.hpp"
#include "../../src/misc/matlab.hpp"
#include "../../src/gui/singleton.hpp"
#include "../../src/ml/factor_graph.hpp"

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
            std::string fileprefix;
            View<PerspectiveCamera> view;
            std::vector<Classified<Line2>> lines;
            std::vector<Vec3> vps;
            Imagei segmentedImage;
            int segmentsNum;
            ImageOfType<Vec<double, 5>> gc;

            std::vector<std::vector<PixelLoc>> edges;
            std::vector<double> scores;

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
            auto fileprefix = "F:\\DataSets\\YorkUrbanDB\\data\\" + name + "\\" + name;

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

            // load occlusion detection results
            //NOT_IMPLEMENTED_YET();


            return{ fileprefix, view, lines, vps, segmentedImage, segmentsNum, gc, { {} }, {}, outdoor };
        }


        void DepthOrdering(const Preparation & p){
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
            SetNecessaryConstraintWeightedAnchors(mg, controls);

            auto ccids = MakeHandledTableForAllComponents(mg, -1);
            int ccnum = ConnectedComponents(mg, controls, ccids, [](const RLGraphConstraintControl & c){
                return c.used && c.weightedAnchors.size() > 0;
            });
            auto mgs = Decompose(mg, ccids, ccnum);
            auto cs = Decompose(mg, controls, ccids, ccnum);
            assert(mgs.size() == cs.size());

            for (int i = 0; i < ccnum; i++){
                auto & mg = mgs[i];
                auto & controls = cs[i];
                AttachWeightedAnchorToCenterOfLargestLineIfNoExists(mg, controls);
                auto vars = SolveVariablesWithBoundedAnchors(mgs[i], controls, true);
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
                    case GeometricContextEstimator::II_VerticalPlanarFace:
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

            // for each cc
            for (int i = 0; i < ccnum; i++){
                auto & mg = mgs[i];
                auto & controls = cs[i];
                RLGraphVars vars;

                if (!AttachWeightedAnchorToCenterOfLargestLineIfNoExists(mg, controls))
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
                    vars = SolveVariablesWithBoundedAnchors(mg, controls, true, true);
                    NormalizeVariables(mg, controls, vars);
                    LooseOrientationConstraintsOnComponents(mg, controls, vars, 0.0, 0.1);
                    vars = SolveVariablesWithBoundedAnchors(mg, controls, true, true);
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
                    OptimizeVariablesWithBoundedAnchors(mg, controls, vars, true, false);
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
                    OptimizeVariablesWithBoundedAnchors(mg, controls, vars, true, true);
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

        void Show(const PerspectiveView & v, const RLGraph & mg, const RLGraphControls & controls, const RLGraphVars & vars){
            gui::Visualizer vis("inversed depth setup");
            Visualize(vis, v, mg, controls, vars);
            vis.camera(v.camera);
            vis.renderOptions.cullBackFace = vis.renderOptions.cullFrontFace = false;
            vis.show(true, true);
        }

        int Control2Label(const RLGraphComponentControl & control){
            if (control.used == true && control.orientationClaz == -1 && control.orientationNotClaz == -1) return 0;
            if (control.used == true && control.orientationClaz == 0 && control.orientationNotClaz == -1) return 1;
            if (control.used == true && control.orientationClaz == 1 && control.orientationNotClaz == -1) return 2;
            if (control.used == true && control.orientationClaz == 2 && control.orientationNotClaz == -1) return 3;
            if (control.used == true && control.orientationClaz == -1 && control.orientationNotClaz == 0) return 4;
            if (control.used == true && control.orientationClaz == -1 && control.orientationNotClaz == 1) return 5;
            if (control.used == true && control.orientationClaz == -1 && control.orientationNotClaz == 2) return 6;
            if (control.used == false && control.orientationClaz == -1 && control.orientationNotClaz == -1) return 7;
            return -1;
        }

        void Label2Control(int label, RLGraphComponentControl & control){
            switch (label){
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

        void LabelGT(const Preparation & p){

            auto & view = p.view;
            auto & lines = p.lines;
            auto & vps = p.vps;
            auto & segmentedImage = p.segmentedImage;
            int segmentsNum = p.segmentsNum;
            auto & gc = p.gc;
            bool outdoor = p.outdoor;

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
                    case GeometricContextEstimator::II_VerticalPlanarFace:
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


            p.showVPLines();

            Labels labels;

            {
                /// read region labels
                labels.regionLabels.resize(rhs.size());
                for (int i = 0; i < labels.regionLabels.size(); i++){
                    RegionHandle rh = rhs[i];
                    if (rh.invalid())
                        continue;
                    const auto & control = controls[rh];
                    labels.regionLabels[i] = Control2Label(control);
                }
                /// read boundary labels
                // todo
            }

            while(true){
                bool accepted = LabelIt(labels, view.image, segmentedImage, bps, 
                    regionLabelNames, boundaryLabelNames, regionLabelColors, boundaryLabelColors);
                if (accepted){
                    LOG("Accepted!");
                    SaveToDisk(p.fileprefix + "gtlabels.cereal", labels);
                    break;
                }

                // apply labels to controls
                // regions
                for (int i = 0; i < labels.regionLabels.size(); i++){
                    RegionHandle rh = rhs[i];
                    if (rh.invalid())
                        continue;
                    auto & control = controls[rh];
                    Label2Control(labels.regionLabels[i], control);
                }
                // boundaries, line-region relations
                // todo

                controls.disableAllInvalidConstraints(mg);

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

                // for each cc
                for (int i = 0; i < ccnum; i++){
                    auto & mg = mgs[i];
                    auto & controls = cs[i];
                    RLGraphVars vars;

                    if (!AttachWeightedAnchorToCenterOfLargestLineIfNoExists(mg, controls))
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
                    


                    vars = SolveVariablesWithBoundedAnchors(mg, controls, false, true);
                    NormalizeVariables(mg, controls, vars);
                    if (1){
                        Show(p, mg, controls, vars);
                    }
                }
            }

        }


        template <class T, class CostFunT, class BranchFunT>
        void BranchAndBound(T & solution, CostFunT && cost, BranchFunT && branch){
            Scored<T> bound = ScoreAs(solution, cost(solution));            
            std::priority_queue<Scored<T>> Q;
            Q.push(bound);

            while (!Q.empty()){
                Scored<T> node = std::move(Q.top());
                Q.pop();
                if (node.score < bound.score){
                    bound = std::move(node);
                    std::cout << "current bound: " << bound.score << std::endl;
                }
                else{
                    // branch
                    std::vector<Scored<T>> subNodes;
                    branch(std::move(node.component), [&subNodes](const T & subNode){
                        subNodes.push_back(ScoreAs(subNode, 0.0));
                    });

                    // calc costs
                    int pnum = std::thread::hardware_concurrency();
                    auto task = [&subNodes, &cost](int tid){
                        auto & subNode = subNodes[tid];
                        subNode.score = cost(subNode.component);
                    };
                    std::cout << "COMPUTING SCORES" << std::endl;
                    std::vector<std::thread> threads;
                    for (int i = 0; i < subNodes.size(); i++){
                        threads.push_back(std::thread(task, i));
                        if (threads.size() == pnum){
                            for (auto & t : threads){
                                t.join();
                            }
                            threads.clear();
                        }
                    }

                    for (auto & t : threads){
                        t.join();
                    }
                    threads.clear();
                    
                    std::cout << "COMPUTING SCORES DONE" << std::endl;
                    for (auto & n : subNodes){
                        if (n.score <= bound.score){
                            Q.push(std::move(n));
                        }
                    }
                    std::cout << "Q size: " << Q.size() << std::endl;
                }
            }
        }


        struct branchFun {
            template <class FunT>
            inline void operator()(RLGraphControls && curNode, FunT && processSubNode) {
                for (auto & r : mg.components<RegionData>()){
                    RLGraphComponentControl & c = curNode[r.topo.hd];
                    if (!c.used){
                        continue;
                    }
                    if (c.orientationClaz == -1 && c.orientationNotClaz == -1){
                        // free -> ground
                        c.orientationClaz = vertVPId;
                        processSubNode(curNode);
                        c.orientationClaz = -1;
                        // free -> along vertical
                        c.orientationNotClaz = vertVPId;
                        processSubNode(curNode);
                        c.orientationNotClaz = -1;
                    }
                    else if (c.orientationClaz == -1 && c.orientationNotClaz == vertVPId){
                        // along vertical -> to hvp1
                        c.orientationClaz = hVPId1;
                        c.orientationNotClaz = -1;
                        processSubNode(curNode);
                        // along vertical -> to hvp2
                        c.orientationClaz = hVPId2;
                        c.orientationNotClaz = -1;
                        processSubNode(curNode);
                        c.orientationClaz = -1;
                        c.orientationNotClaz = vertVPId;
                    }
                }

                // connected -> disconnected
                for (auto & b : mg.constraints<RegionBoundaryData>()){
                    RLGraphConstraintControl & c = curNode[b.topo.hd];
                    if (!c.used){
                        continue;
                    }
                    c.used = false;
                    processSubNode(curNode);
                    c.used = true;
                }
                for (auto & b : mg.constraints<RegionLineConnectionData>()){
                    RLGraphConstraintControl & c = curNode[b.topo.hd];
                    if (!c.used){
                        continue;
                    }
                    c.used = false;
                    processSubNode(curNode);
                    c.used = true;
                }

            }
            int vertVPId, hVPId1, hVPId2;
            const RLGraph & mg;
        };

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
            int hVPId1 = (vertVPId + 1) % 3;
            int hVPId2 = (vertVPId + 2) % 3;

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
                    case GeometricContextEstimator::II_VerticalPlanarFace:
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


            // for each cc
            for (int i = 0; i < ccnum; i++){
                auto & mg = mgs[i];
                auto & controls = cs[i];
                RLGraphVars vars;

                if (!AttachWeightedAnchorToCenterOfLargestLineIfNoExists(mg, controls))
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
                vars = SolveVariablesWithBoundedAnchors(mg, controls, false, true);
                NormalizeVariables(mg, controls, vars);
                if (1){
                    Show(p, mg, controls, vars);
                }
            }
        }
















        struct BoundaryJunction {
            std::set<RegionBoundaryHandle> bhs;
            std::vector<Vec3> positions;
        };

        std::vector<BoundaryJunction> GetJunctions(const RLGraph & mg,
            const PerspectiveCamera & cam, const Imagei & regions, 
            const std::vector<RegionHandle> & regionIds2Handles){

            std::vector<BoundaryJunction> boundaryJunctions;
            auto bjunctions = ExtractBoundaryJunctions(regions);
            for (auto & junc : bjunctions){
                auto & regionIds = junc.first;
                auto & ps = junc.second;
                std::set<RegionHandle> rhs;
                for (int regionId : regionIds){
                    RegionHandle rh = regionIds2Handles[regionId];
                    if (rh.invalid())
                        continue;
                    rhs.insert(rh);
                }
                if (rhs.size() < 3)
                    continue;
                // locate boundary handles among rhs
                std::set<RegionBoundaryHandle> bhs;
                for (RegionHandle rh : rhs){
                    auto & relatedBhs = mg.topo(rh).constraints<RegionBoundaryData>();
                    for (RegionBoundaryHandle bh : relatedBhs){
                        auto anotherRh = mg.topo(bh).component<0>();
                        if (anotherRh == rh){
                            anotherRh = mg.topo(bh).component<1>();
                        }
                        if (Contains(rhs, anotherRh)){
                            bhs.insert(bh);
                        }
                    }
                }
                if (bhs.size() < 3)
                    continue;
                std::vector<Vec3> positions(ps.size());
                for (int i = 0; i < ps.size(); i++){
                    positions[i] = normalize(cam.spatialDirection(ps[i]));
                }
                boundaryJunctions.push_back(BoundaryJunction{ std::move(bhs), std::move(positions) });
            }
            return boundaryJunctions;
        }

        double Distance(const Point2 & p, const std::vector<PixelLoc> & edge){
            double mind = std::numeric_limits<double>::max();
            for (const auto & c : edge){
                auto cc = core::vec_cast<double>(c);
                double d = core::Distance(p, cc);
                if (d < mind){
                    mind = d;
                }
            }
            return mind;
        }



        HandledTable<RegionBoundaryHandle, double> OcclusionResponce(const RLGraph & mg,
            const PerspectiveCamera & cam, const std::vector<std::vector<PixelLoc>> & edges,
            const std::vector<double> & scores){
            assert(edges.size() == scores.size());
            HandledTable<RegionBoundaryHandle, double> occlusionResponse = 
                mg.createConstraintTable<RegionBoundaryData>(0.0);
            auto getBB = [&edges](int eid)->Box2 {
                return BoundingBoxOfContainer(edges.at(eid));
            };
            std::vector<int> eids(edges.size());
            std::iota(eids.begin(), eids.end(), 0);
            RTreeWrapper<int, decltype(getBB)> rtree(eids.begin(), eids.end(), getBB);
            for (auto & bd : mg.constraints<RegionBoundaryData>()){
                double & resp = occlusionResponse[bd.topo.hd];
                double allSampleScoreSum = 0;
                int allSampleNum = 0;
                for (auto & ps : bd.data.normalizedEdges){
                    for (auto & p : ps){
                        allSampleNum++;
                        auto sample = cam.screenProjection(p);
                        double sampleScore = 0.0;
                        int nearbyEdgeNum = 0;
                        rtree.search(Box2(sample, sample).expand(5.0),
                            [&scores, &sampleScore, &nearbyEdgeNum, &edges, &sample](int eid) -> bool {
                            double dist = Distance(sample, edges[eid]);
                            if (dist <= 3.0){
                                sampleScore += scores[eid] > 0 ? 1.0 : (scores[eid] == 0.0 ? 0.0 : -1.0);
                                nearbyEdgeNum++;
                            }
                            return true;
                        });
                        sampleScore /= std::max(nearbyEdgeNum, 1);
                        allSampleScoreSum += sampleScore;
                    }
                }
                occlusionResponse[bd.topo.hd] = allSampleScoreSum / std::max(allSampleNum, 1);
            }
            return occlusionResponse;
        }










        void BPReconstruct(const Preparation & p){
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
            int hVPId1 = (vertVPId + 1) % 3;
            int hVPId2 = (vertVPId + 2) % 3;

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
                    case GeometricContextEstimator::II_VerticalPlanarFace:
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

            // initial status
            vars = SolveVariablesWithBoundedAnchors(mg, controls, false, true);
            NormalizeVariables(mg, controls, vars);
            if (1){
                Show(p, mg, controls, vars);
            }

            // 
            auto boundaryJunctions = GetJunctions(mg, view.camera, segmentedImage, rhs);
            auto gcResponse = CollectFeatureMeanOnRegions(mg, view.camera, gc);
            auto occResponse = OcclusionResponce(mg, view.camera, p.edges, p.scores);




            size_t vpnum = vps.size();

            auto lineLeftRegionConnections =
                mg.createComponentTable<LineData, std::set<RegionLineConnectionHandle>>();
            auto lineRightRegionConnections =
                mg.createComponentTable<LineData, std::set<RegionLineConnectionHandle>>();
            for (auto & rl : mg.constraints<RegionLineConnectionData>()){
                RegionHandle rh = rl.topo.component<0>();
                LineHandle lh = rl.topo.component<1>();
                const Line3 & lineProj = mg.data(lh).line;
                Vec3 rightDir = normalize(lineProj.first.cross(lineProj.second));
                const Vec3 & regionCenter = mg.data(rh).normalizedCenter;
                bool isOnRight = (regionCenter - lineProj.center()).dot(rightDir) >= 0;
                (isOnRight ? lineRightRegionConnections : lineLeftRegionConnections)[lh].insert(rl.topo.hd);
            }


            // solve!


            size_t regionNum = mg.internalComponents<RegionData>().size();
            size_t lineNum = mg.internalComponents<LineData>().size();
            size_t boundaryNum = mg.internalConstraints<RegionBoundaryData>().size();
            size_t rlConNum = mg.internalConstraints<RegionLineConnectionData>().size();
            size_t llConNum = mg.internalConstraints<LineRelationData>().size();
            size_t bjuncNum = boundaryJunctions.size();



            // data that MUST be updated after each inverse depth optimization !!!!
            struct ContinuousCaches {
                HandledTable<RegionHandle, Plane3> regionPlanes;
                HandledTable<LineHandle, Line3> spatialLines;
                HandledTable<RegionBoundaryHandle, double> rrDistances;
                HandledTable<RegionLineConnectionHandle, double> rlDistances;
                HandledTable<LineRelationHandle, double> llDistances;

                Plane3 & operator[](RegionHandle rh) { return regionPlanes[rh]; }
                Line3 & operator[](LineHandle lh) { return spatialLines[lh]; }
                double & operator[](RegionBoundaryHandle bh) { return rrDistances[bh]; }
                double & operator[](RegionLineConnectionHandle bh) { return rlDistances[bh]; }
                double & operator[](LineRelationHandle bh) { return llDistances[bh]; }

                void update(const RLGraph & mg, const RLGraphControls & controls, const RLGraphVars & vars){
                    // update continuous cache
                    for (auto & r : mg.components<RegionData>()){
                        (*this)[r.topo.hd] = Instance(mg, controls, vars, r.topo.hd);
                    }
                    for (auto & l : mg.components<LineData>()){
                        (*this)[l.topo.hd] = Instance(mg, controls, vars, l.topo.hd);
                    }

                    for (auto & c : mg.constraints<RegionBoundaryData>()){
                        double dist = 0.0;
                        int num = 0;
                        auto & samples = c.data.normalizedSampledPoints;
                        auto & inst1 = (*this)[c.topo.component<0>()];
                        auto & inst2 = (*this)[c.topo.component<1>()];
                        for (auto & ss : samples){
                            for (auto & s : ss){
                                double d1 = DepthAt(s, inst1);
                                double d2 = DepthAt(s, inst2);
                                dist += abs(d1 - d2);
                                num++;
                            }
                        }
                        dist /= num;
                        (*this)[c.topo.hd] = dist;
                    }
                    for (auto & c : mg.constraints<RegionLineConnectionData>()){
                        double dist = 0.0;
                        int num = 0;
                        auto & inst1 = (*this)[c.topo.component<0>()];
                        auto & inst2 = (*this)[c.topo.component<1>()];
                        for (auto & s : c.data.normalizedAnchors){
                            double d1 = DepthAt(s, inst1);
                            double d2 = DepthAt(s, inst2);
                            dist += abs(d1 - d2);
                            num++;
                        }
                        dist /= num;
                        (*this)[c.topo.hd] = dist;
                    }
                }
            };
            ContinuousCaches cache = {
                mg.createComponentTable<RegionData, Plane3>(),
                mg.createComponentTable<LineData, Line3>(),
                mg.createConstraintTable<RegionBoundaryData>(0.0),
                mg.createConstraintTable<RegionLineConnectionData>(0.0),
                mg.createConstraintTable<LineRelationData>(0.0)
            };


            ml::FactorGraph fg;
            fg.reserveVarCategories(5);
            fg.reserveVars(regionNum + lineNum + boundaryNum + rlConNum + llConNum);

            /// add vars

            // add region orientation constraint flags
            // 3: {vpvert, vp, free}
            auto regionOrientationVc = fg.addVarCategory(ml::FactorGraph::VarCategory{ 3, 0.5 });
            auto regionVhs = mg.createComponentTable<RegionData, ml::FactorGraph::VarHandle>();
            for (auto & r : mg.components<RegionData>()){
                regionVhs[r.topo.hd] = fg.addVar(regionOrientationVc);
            }

            // add boundary flags
            // 2: {not connected, connected}
            auto boundaryOcclusionVc = fg.addVarCategory(ml::FactorGraph::VarCategory{ 2, 0.5 });
            auto boundaryVhs = mg.createConstraintTable<RegionBoundaryData, ml::FactorGraph::VarHandle>();
            for (auto & b : mg.constraints<RegionBoundaryData>()){
                boundaryVhs[b.topo.hd] = fg.addVar(boundaryOcclusionVc);
            }

            // add line connection flags
            // 3 : {both, left, right}
            auto lineConnectionSidesVc = fg.addVarCategory(ml::FactorGraph::VarCategory{ 3, 0.5 });
            auto lineConnectionSidesVhs = mg.createComponentTable<LineData, ml::FactorGraph::VarHandle>();
            for (auto & l : mg.components<LineData>()){
                lineConnectionSidesVhs[l.topo.hd] = fg.addVar(lineConnectionSidesVc);
            }

            // add ll connection flags
            // 2: {not connected, connected}
            auto llConVc = fg.addVarCategory(ml::FactorGraph::VarCategory{ 2, 0.5 });
            auto llConVhs = mg.createConstraintTable<LineRelationData, ml::FactorGraph::VarHandle>();
            for (auto & r : mg.constraints<LineRelationData>()){
                llConVhs[r.topo.hd] = fg.addVar(llConVc);
            }


            /// add factors



        }






#define TAG(name) ("." + std::string(#name) + std::to_string(name))

        struct Config {
            bool reestVPLines;
            bool withGC;
            bool withHCons;
            bool useWeights;
            bool useAllAnchors;

            template <class A>
            void serialize(A & a) { a(reestVPLines, withGC, withHCons, useWeights, useAllAnchors); }

            std::string tag() const {
                return TAG(reestVPLines) + TAG(withGC) + TAG(withHCons) + TAG(useWeights) + TAG(useAllAnchors);
            }
        };




        void SurfaceOrientationPrediction(const Preparation & p, const Config & config){

            misc::Matlab matlab;
            std::string tag = config.tag();

            auto view = p.view;
            auto lines = p.lines;
            auto vps = p.vps;
            auto & segmentedImage = p.segmentedImage;
            int segmentsNum = p.segmentsNum;
            auto & gc = p.gc;

            if (config.reestVPLines){

                matlab << "load('" + p.fileprefix + "om_vp_lines.mat');";
                DenseMatd vpsd;
                matlab.GetVariable("vp", vpsd);
                std::vector<HPoint2> hvps(3);
                for (int i = 0; i < 3; i++){
                    hvps[i].numerator = { vpsd(i, 0), vpsd(i, 1) };
                    hvps[i].denominator = 1.0;
                }
                view = CreatePerspectiveView(view.image, hvps, view.camera.eye(), view.camera.center(), view.camera.up());
                for (int i = 0; i < 3; i++){
                    vps[i] = normalize(view.camera.spatialDirection(hvps[i].value()));
                }
                double linenum = 0.0;
                matlab << "linenum = length(lines);";
                matlab.GetVariable("linenum", linenum);
                lines.clear();
                lines.reserve(linenum);
                for (int i = 0; i < linenum; i++){
                    DenseMatd pdata1, pdata2;
                    matlab << ("p = lines(" + std::to_string(i + 1) + ").point1;");
                    matlab.GetVariable("p", pdata1);
                    matlab << ("p = lines(" + std::to_string(i + 1) + ").point2;");
                    matlab.GetVariable("p", pdata2);
                    double claz = 0;
                    matlab << ("c = lines(" + std::to_string(i + 1) + ").lineclass;");
                    matlab.GetVariable("c", claz);
                    lines.push_back(ClassifyAs(Line2(Point2(pdata1(0), pdata1(1)), Point2(pdata2(0), pdata2(1))), int(claz - 1)));
                }

                //std::vector<core::Classified<core::Ray2>> vpRays;
                //for (int i = 0; i < 3; i++){
                //    std::cout << "vp[" << i << "] = " << view.camera.screenProjection(vps[i]) << std::endl;
                //    for (double a = 0; a <= M_PI * 2.0; a += 0.1){
                //        core::Point2 p = core::Point2(view.image.cols / 2, view.image.rows / 2) + core::Vec2(cos(a), sin(a)) * 1000.0;
                //        vpRays.push_back(core::ClassifyAs(core::Ray2(p, (view.camera.screenProjectionInHPoint(vps[i]) - core::HPoint2(p, 1.0)).numerator), i));
                //    }
                //}
                //gui::Visualizer2D(view.image)
                //    << gui::manip2d::SetColorTable(gui::ColorTable(gui::ColorTableDescriptor::RGB).appendRandomizedGreyColors(vps.size() - 3))
                //    << gui::manip2d::SetThickness(1)
                //    << vpRays
                //    << gui::manip2d::SetThickness(2)
                //    << lines
                //    << gui::manip2d::Show();
            }

            bool outdoor = p.outdoor;

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

            // gc
            AttachGeometricContextConstraints(mg, controls, view.camera, gc,
                [outdoor, vertVPId, config](RLGraphComponentControl & c, const Vec<double, 5> & v, double s){
                bool withGC = config.withGC;
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
                        if (withGC && !config.withHCons){
                            c.orientationClaz = -1;
                            c.orientationNotClaz = vertVPId;
                        }
                        break;
                    case GeometricContextEstimator::OI_Clutter:
                        if (withGC && !config.withHCons){
                            c.orientationClaz = -1;
                            c.orientationNotClaz = -1;
                        }
                        break;
                    case GeometricContextEstimator::OI_Porous:
                        if (withGC && !config.withHCons){
                            c.orientationClaz = -1;
                            c.orientationNotClaz = -1;
                        }
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
                        if (withGC && !config.withHCons){
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
                        if (withGC && !config.withHCons){
                            c.orientationClaz = c.orientationNotClaz = -1;
                        }
                        break;
                    default:
                        break;
                    }
                }
            });
            if (config.withHCons){
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


            if (!AttachWeightedAnchorToCenterOfLargestLineIfNoExists(mg, controls))
                return;

            ResetToFullArmorAnchors(mg, controls);
            vars = SolveVariablesWithBoundedAnchors(mg, controls,
                    config.useWeights, 1000);
            assert(!HasInfOrNaNValue(vars));


            //vars = SolveVariablesWithBoundedAnchors(mg, controls, config.useWeights/*false*/,
            //    config.useAllAnchors);

            {
                // save
                auto planes = Instances<RegionData>(mg, controls, vars);
                Image3f rlomap = Image3f::zeros(view.image.size());
                for (auto it = rlomap.begin(); it != rlomap.end(); ++it){
                    int rid = segmentedImage(it.pos());
                    auto rh = rhs[rid];
                    if (rh.invalid()){
                        continue;
                    }
                    auto & plane = planes[rh];
                    for (int i = 0; i < 3; i++){
                        (*it)[i] = abs(normalize(plane.normal).dot(normalize(vps[i])));
                    }
                }
                matlab.PutVariable("rlomap", rlomap);
            }

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






        void Run(){

            std::string name;
            name = "p1020171";

            misc::Matlab matlab;
            matlab.Load("F:\\DataSets\\YorkUrbanDB\\data\\names.mat");
            double dnum = 0;
            matlab << "num = length(names);";
            matlab.GetVariable("num", dnum);
            int num = dnum;

            Config conf;
            conf.reestVPLines = false;
            conf.useAllAnchors = true;
            conf.useWeights = true;
            conf.withGC = true;
            conf.withHCons = true;

            std::vector<std::pair<Config, int>> errors;

            //for (bool c1 : {true}){
            //    conf.useAllAnchors = c1;
            //    for (bool c2 : {false}){
            //        conf.useWeights = c2;

            //        for (bool c3 : {false}){
            //            conf.withGC = c3;
            //            for (bool c4 : {true}){
            //                conf.withHCons = c4;
            int i = 99;
                            //for (int i = 0; i < num; i++)
                            {
                                matlab << ("name = names{" + std::to_string(i + 1) + "};");
                                std::string name;
                                matlab.GetVariable("name", name);
                                std::cout << "processing [" << i << "] " << name << std::endl;
                                std::cout << "tag: " << conf.tag() << std::endl;

                                //if (name != "p1020867"){
                                //    continue;
                                //}

                                //if (!nowstart)
                                //    continue;


                                auto p = Prepare(name);

                                //LabelGT(p);

                                //DepthOrdering(p);
                                //ReconstructLines(p);
                                //ReconstructModel(p);
                                //ExhausticReconstruct(p);
                                //BPReconstruct(p);
                                try{
                                    SurfaceOrientationPrediction(p, conf);
                                }
                                catch (std::exception e){
                                    std::cout << e.what() << std::endl;
                                    errors.emplace_back(conf, i);
                                }
                            }

            //            }
            //        }
            //    }
            //}

            //SaveToDisk("cache/errorids.cereal", errors);

        }
    }

}