#include "../../src/core/basic_types.hpp"
#include "../../src/experimental/rl_graph.hpp"
#include "../../src/experimental/tools.hpp"
#include "../../src/misc/matlab_engine.hpp"
#include "../../src/gui/scene.hpp"
#include "../../src/gui/canvas.hpp"
#include "../../src/gui/utility.hpp"

#include "routines.hpp"


namespace panolyz {

    namespace PanoramaIndoor{

        void Run(){

            std::string path;

            using namespace panoramix;
            using namespace core;
            using namespace experimental;

            core::Image3ub image = gui::PickAnImage(PROJECT_TEST_DATA_DIR_STR"/panorama/indoor", &path);
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

                Save(path, "pre", view, cams, lines, vps, segmentedImage, vertVPId);
            }
            else{
                Load(path, "pre", view, cams, lines, vps, segmentedImage, vertVPId);
            }




            std::vector<PerspectiveCamera> hcams;
            std::vector<Weighted<View<PerspectiveCamera, Image5d>>> gcs;
            if (REFRESH){
                // extract gcs
                hcams = CreateHorizontalPerspectiveCameras(view.camera, 16, 500, 400, 300);
                gcs.resize(hcams.size());
                misc::MatlabEngine matlab;
                for (int i = 0; i < hcams.size(); i++){
                    auto pim = view.sampled(hcams[i]);
                    auto pgc = ComputeGeometricContext(pim.image, false, true);
                    gcs[i].component.camera = hcams[i];
                    gcs[i].component.image = pgc;
                    gcs[i].score = abs(normalize(hcams[i].forward()).dot(normalize(view.camera.up())));
                }
                Save(path, "hcamsgcs", hcams, gcs);
            }
            else{
                Load(path, "hcamsgcs", hcams, gcs);                
            }

            Image5d gc;
            gc = Combine(view.camera, gcs).image;
            if (1){
                gui::AsCanvas(gc).show();
            }



            // consider only lines
            if (0){

                RLGraph mg;
                RLGraphControls controls;
                RLGraphVars vars;

                for (int i = 0; i < cams.size(); i++){
                    AppendLines(mg, lines[i], cams[i], vps);
                }

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

                    gui::SceneBuilder vis;
                    Visualize(vis, view, mg, controls, vars);
                    vis.show(true, false, 
                        gui::RenderOptions().camera(PerspectiveCamera(500, 500, Point2(250, 250), 300, Point3(1, 1, 1), Point3(0, 0, 0))));
                }

            }


            RLGraph mg;
            RLGraphControls controls;
            RLGraphVars vars;
            std::vector<RegionHandle> rhs;


            // consider both lines and regions
            if (1){

                for (int i = 0; i < cams.size(); i++){
                    AppendLines(mg, lines[i], cams[i], vps);
                }
                rhs = AppendRegions(mg, segmentedImage, view.camera, 0.03, 0.02, 4, 4, false);

                controls = RLGraphControls(mg, vps);
                AttachPrincipleDirectionConstraints(mg, controls, M_PI / 20.0);
                AttachWallConstriants(mg, controls, M_PI / 10);
                
                auto gcMeanOnRegions = CollectFeatureMeanOnRegions(mg, view.camera, gc);
                auto up = normalize(vps[vertVPId]);
                if (up.dot(- view.camera.up()) < 0){
                    up = -up;
                }
                gui::AsCanvas(Print(mg, segmentedImage, view.camera, rhs, [&up, &mg](RegionHandle rh){
                    return mg.data(rh).normalizedCenter.dot(up) < 0 ? gui::Red : gui::Blue;
                })).show();
                SetComponentControl<RegionData>(controls, 
                    [&mg, &controls, &vps, &gcMeanOnRegions, &view, &up, vertVPId](RegionHandle rh, RLGraphComponentControl & c){
                    auto & gcMean = gcMeanOnRegions[rh];
                    assert(IsFuzzyZero(std::accumulate(std::begin(gcMean.val), std::end(gcMean.val), 0.0) - 1.0, 1e-4));
                    size_t maxid = std::max_element(std::begin(gcMean.val), std::end(gcMean.val)) - gcMean.val;
                    double floorScore = gcMean[ToUnderlying(GeometricContextIndex::FloorOrGround)];
                    if (maxid == ToUnderlying(GeometricContextIndex::FloorOrGround) && mg.data(rh).normalizedCenter.dot(up) < 0){ // lower
                        // assign horizontal constriant
                        c.orientationClaz = vertVPId;
                        c.orientationNotClaz = -1;
                        c.used = true;
                        return;
                    }
                    double wallScore = gcMean[ToUnderlying(GeometricContextIndex::Vertical)];
                    if (wallScore > 0.6){
                        // assign vertical constraint
                        c.orientationClaz = -1;
                        c.orientationNotClaz = vertVPId;
                        c.used = true;
                        return;
                    }
                    double clutterScore = gcMean[ToUnderlying(GeometricContextIndex::ClutterOrPorous)];
                    if (clutterScore > 0.7){
                        c.orientationClaz = c.orientationNotClaz = -1;
                        c.used = false;
                        return;
                    }
                });
                controls.disableAllInvalidConstraints(mg);


                // set constraint weights
                SetNecessaryConstraintWeightedAnchors(mg, controls);

                // cc decompose
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
                    RLGraphVars vars;

                    if (!AttachWeightedAnchorToCenterOfLargestLineIfNoExists(mg, controls))
                        continue;

                    //vars = SolveVariablesWithBoundedAnchors(mg, controls, false, true);
                    ResetToSampledArmorAnchors(mg, controls, 0.05);
                    vars = SolveVariablesWithBoundedAnchors(mg, controls, false, 1000);
                    NormalizeVariables(mg, controls, vars);
                    //SolveVariablesWithBoundedAnchors(mg, controls, vars);
                    std::cout << "score = " << Score(mg, controls, vars) << std::endl;

                    LooseOrientationConstraintsOnComponents(mg, controls, vars, 0.2, 0.01, 0.1);
                    if (!AttachWeightedAnchorToCenterOfLargestLineIfNoExists(mg, controls))
                        continue;

                    vars = SolveVariablesWithBoundedAnchors(mg, controls, false, 1000);
                    NormalizeVariables(mg, controls, vars);

                  /*  AttachFloorAndCeilingConstraints(mg, controls, vars, 0.1, 0.6);

                    if (!AttachWeightedAnchorToCenterOfLargestLineIfNoExists(mg, controls))
                        continue;
                    vars = SolveVariablesWithBoundedAnchors(mg, controls, false, true);
                    NormalizeVariables(mg, controls, vars);*/

                    gui::SceneBuilder vis;
                    Visualize(vis, view, mg, controls, vars);
                    vis.show(true, false,
                        gui::RenderOptions().camera(PerspectiveCamera(500, 500, Point2(250, 250), 300, Point3(1, 1, 1), Point3(0, 0, 0))));

                    {
                        LayeredShape3 shape;
                        auto polygons = RegionPolygons(mg, controls, vars);
                        int vertVPId = NearestDirectionId(controls.vanishingPoints);
                        double medianDepth = MedianCenterDepth(mg, controls, vars);
                        Vec3 vertDir = normalize(controls.vanishingPoints[vertVPId]);

                        auto range = experimental::EstimateEffectiveRangeAlongDirection(polygons, vertDir, medianDepth * 0.02, 0.9, -1e5, -1e5);

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

                        gui::SceneBuilder viz;
                        viz.begin(shape).shaderSource(gui::OpenGLShaderSourceDescriptor::XPanorama).resource("texture").end();
                        viz.show(true, false, gui::RenderOptions()
                            .renderMode(gui::RenderModeFlag::Triangles | gui::RenderModeFlag::Lines)
                            .backgroundColor(gui::White)
                            .bwColor(0.0)
                            .bwTexColor(1.0)
                            .cullBackFace(false)
                            .cullFrontFace(true)
                            .camera(PerspectiveCamera(1000, 800, Point2(500, 400), 800, Point3(-1, 1, 1), Point3(0, 0, 0))));
                    }
                }

                //SaveToDisk("./cache/mgp", mg, controls, vars);
            }
            else{
                //LoadFromDisk("./cache/mgp", mg, controls, vars);
            }

        }
    }
}