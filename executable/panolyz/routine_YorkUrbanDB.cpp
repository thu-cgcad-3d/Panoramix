#include "../../src/core/basic_types.hpp"
#include "../../src/experimental/rl_graph.hpp"
#include "../../src/experimental/rl_graph_occlusion.hpp"
#include "../../src/experimental/tools.hpp"
#include "../../src/gui/visualize2d.hpp"
#include "../../src/gui/visualizers.hpp"
#include "../../src/misc/matlab.hpp"

#include "routines.hpp"

namespace panolyz {

    namespace YorkUrbanDB {

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



        void Test(const std::string & name){
            auto path = "F:\\DataSets\\YorkUrbanDB\\data\\" + name + "\\" + name + ".jpg";

            View<PerspectiveCamera> view;
            std::vector<Classified<Line2>> lines;
            std::vector<Vec3> vps;
            Imagei segmentedImage;
            int segmentsNum;
            ImageOfType<Vec<double, 5>> gc;

            bool outdoor;

            if (1){
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

                if (1) {
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
                        << gui::manip2d::Show(0);
                }

                // append regions
                SegmentationExtractor segmenter;
                segmenter.params().algorithm = SegmentationExtractor::GraphCut;
                segmenter.params().c = 100.0;
                std::vector<Line2> pureLines(lines.size());
                for (int i = 0; i < lines.size(); i++)
                    pureLines[i] = lines[i].component;
                std::tie(segmentedImage, segmentsNum) = segmenter(view.image, pureLines, view.image.cols / 100.0);

                if (1){
                    gui::Visualizer2D::Params vParams;
                    vParams.colorTable = gui::CreateRandomColorTableWithSize(segmentsNum);
                    gui::Visualizer2D(segmentedImage, vParams)
                        << gui::manip2d::Show();
                }

                Save(path, "pre", view, lines, vps, segmentedImage, segmentsNum);
            }
            else{
                Load(path, "pre", view, lines, vps, segmentedImage, segmentsNum);
            }


            if (0){
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


            RLGraph mg;
            RLGraphControls controls;
            RLGraphVars vars;
            std::vector<RegionHandle> rhs;


            if (1){
                int vertVPId = NearestDirectionId(vps, view.camera.upward());
                AppendLines(mg, lines, view.camera, vps);
                rhs = AppendRegions(mg, segmentedImage, view.camera, 0.01, 0.02, 3, 2);

                controls = RLGraphControls(mg, vps);

                AttachGeometricContextConstraints(mg, controls, view.camera, gc,
                    [outdoor, vertVPId](RLGraphComponentControl & c, const Vec<double, 5> & v, double s){
                    int maxlabel = std::max_element(v.val, v.val + 5) - v.val;
                    if (s > 2){
                        if (outdoor){
                            if (maxlabel == 0){ // ground
                                c.orientationClaz = vertVPId;
                            }
                            else if (maxlabel == 1){ // vertical
                                c.orientationClaz = -1;
                                c.orientationNotClaz = vertVPId;
                            }
                            else if (maxlabel == 2){ // clutter
                                c.orientationClaz = -1;
                                c.orientationNotClaz = -1;
                            }
                            else if (maxlabel == 3){ // poros
                                c.used = false;
                            }
                            else if (maxlabel == 4){ // sky
                                c.used = false;
                            }
                        }
                        else{

                        }
                    }
                });

                //AttachPrincipleDirectionConstraints(mg, controls);
                //AttachWallConstriants(mg, controls);

                // set constraint weights
                SetConstraintWeights<LineRelationData>(controls, [&mg](LineRelationHandle h){
                    return std::max(mg.data(h).junctionWeight, 3.0f);
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
                    if (1){
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

                    vars = SolveVariables(mg, controls, true);
                    NormalizeVariables(mg, controls, vars);
                    std::cout << "score = " << Score(mg, controls, vars) << std::endl;

                    LooseOrientationConstraintsOnComponents(mg, controls, vars, 0.2, 0.02, 0.1);
                    if (!AttachAnchorToCenterOfLargestLineIfNoAnchorExists(mg, controls))
                        continue;

                    vars = SolveVariables(mg, controls);
                    NormalizeVariables(mg, controls, vars);
                    /*
                    AttachFloorAndCeilingConstraints(mg, controls, vars, 0.1, 0.6);

                    if (!AttachAnchorToCenterOfLargestRegionIfNoAnchorExists(mg, controls) &&
                    !AttachAnchorToCenterOfLargestLineIfNoAnchorExists(mg, controls))
                    continue;
                    vars = SolveVariables(mg, controls);
                    NormalizeVariables(mg, controls, vars);*/

                    gui::Visualizer vis;
                    Visualize(vis, view, mg, controls, vars);
                    vis.camera(view.camera);
                    vis.renderOptions.cullBackFace = vis.renderOptions.cullFrontFace = false;
                    vis.show(true, true);

                    if (0){
                        LayeredShape3 shape;
                        auto polygons = RegionPolygons(mg, controls, vars);
                        int vertVPId = NearestDirectionId(controls.vanishingPoints);
                        double medianDepth = MedianCenterDepth(mg, controls, vars);
                        Vec3 vertDir = normalize(controls.vanishingPoints[vertVPId]);

                        auto range = experimental::EstimateEffectiveRangeAlongDirection(polygons, vertDir, medianDepth * 0.02, 0.7, -1e5, -1e5);

                        std::vector<Chain3> chains;
                        for (double x = range.first; x <= range.second; x += medianDepth * 0.02){
                            Plane3 cutplane(vertDir * x, vertDir);
                            auto loop = experimental::MakeSectionalPieces(polygons, cutplane);
                            if (loop.empty())
                                continue;
                            chains.push_back(experimental::MakeChain(loop));
                        }

                        for (int i = 0; i < chains.size(); i++){
                            shape.layers.push_back(std::move(chains[i].points));
                        }
                        shape.normal = vertDir;

                        gui::ResourceStore::set("texture", view.image);

                        gui::Visualizer viz;
                        viz.begin(shape).shaderSource(gui::OpenGLShaderSourceDescriptor::XPanorama).resource("texture").end();

                        viz.renderOptions.renderMode = gui::RenderModeFlag::Triangles | gui::RenderModeFlag::Lines;
                        viz.renderOptions.backgroundColor = gui::ColorTag::White;
                        viz.renderOptions.bwColor = 0.0;
                        viz.renderOptions.bwTexColor = 1.0;
                        viz.renderOptions.cullBackFace = false;
                        viz.renderOptions.cullFrontFace = false;
                        viz.camera(PerspectiveCamera(1000, 800, Point2(500, 400),
                            800, Point3(-1, 1, 1), Point3(0, 0, 0), Vec3(0, -1, 0)));

                        viz.show(true, true);
                    }
                }

                //SaveToDisk("./cache/mgp", mg, controls, vars);
            }
            else{
                //LoadFromDisk("./cache/mgp", mg, controls, vars);
            }
        }



        void OutputOrientationMap(const Image & image, const RLGraph & mg,
            const RLGraphControls & controls, const RLGraphVars & vars){
            auto planes = Instances<RegionData>(mg, controls, vars);

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
                Test(name);
            }

        }

    }
}