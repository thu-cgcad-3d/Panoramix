#include "../../src/core/basic_types.hpp"
#include "../../src/core/containers.hpp"
#include "../../src/experimental/rl_graph_solver.hpp"
#include "../../src/experimental/rl_graph_occlusion.hpp"
#include "../../src/experimental/tools.hpp"
#include "../../src/misc/matlab_api.hpp"
#include "../../src/gui/scene.hpp"
#include "../../src/gui/canvas.hpp"
#include "../../src/gui/utility.hpp"

#include "routines.hpp"


namespace panolyz {

    namespace PanoramaIndoor{

        void Run(){

            std::string path;

            using namespace pano;
            using namespace core;
            using namespace experimental;

            misc::Matlab matlab;

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

            static const bool REFRESH = false;

            if (0 || !Load(path, "pre", view, cams, lines, vps, segmentedImage, vertVPId)) {
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
                std::tie(segmentedImage, segmentsNum) = segmenter(view.image, line3ds, view.camera, DegreesToRadians(1));     
                RemoveThinRegionInSegmentation(segmentedImage, true);
                segmentsNum = DensifySegmentation(segmentedImage, true);
                assert(IsDenseSegmentation(segmentedImage));

                Save(path, "pre", view, cams, lines, vps, segmentedImage, vertVPId);
            }


            if (1) {
                auto ctable = gui::CreateRandomColorTableWithSize(MinMaxValOfImage(segmentedImage).second + 1);
                gui::AsCanvas(ctable(segmentedImage)).show();
                gui::AsCanvas(ctable(segmentedImage)).add(view.image).show();
            }

            
            std::vector<PerspectiveCamera> hcams;
            std::vector<Weighted<View<PerspectiveCamera, Image5d>>> gcs;
            if (REFRESH || !Load(path, "hcamsgcs", hcams, gcs)) {
                // extract gcs
                hcams = CreateHorizontalPerspectiveCameras(view.camera, 16, 500, 600, 300);
                //hcams = CreatePanoContextCameras(view.camera, 500, 400, 300);
                gcs.resize(hcams.size());
                for (int i = 0; i < hcams.size(); i++){
                    auto pim = view.sampled(hcams[i]);
                    auto pgc = ComputeGeometricContext(matlab, pim.image, false, true);
                    gcs[i].component.camera = hcams[i];
                    gcs[i].component.image = pgc;
                    gcs[i].score = abs(1.0 - normalize(hcams[i].forward()).dot(normalize(view.camera.up())));
                }
                Save(path, "hcamsgcs", hcams, gcs);
            }

            Image5d gc;
            if (REFRESH || !Load(path, "gcmerged", gc)) {
                gc = Combine(view.camera, gcs).image;
                Save(path, "gcmerged", gc);
            }

            if (0){
                std::vector<Imaged> gcChannels;
                cv::split(gc, gcChannels);
                gui::AsCanvas(ConvertToImage3d(gc)).show();
                misc::MAT gcMat("./cache/gcMat.mat", misc::MAT::Write);
                gcMat.setVar("gc", gc);                
            }


            // occ
            std::vector<Scored<Chain3>> occbnds;
            if (REFRESH || !Load(path, "occ_panoindoor", occbnds)) {
                for (int i = 0; i < hcams.size(); i++) {
                    auto pim = view.sampled(hcams[i]).image;
                    auto occ = AsDimensionConvertor(hcams[i]).toSpace(DetectOcclusionBoundary(matlab, pim));
                    for (auto & soc : occ) {
                        occbnds.push_back(soc);
                    }
                }
                Save(path, "occ_panoindoor", occbnds);
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
                    auto vars = SolveVariablesWithBoundedAnchors(matlab, mgs[i], controls, true);
                    NormalizeVariables(mg, controls, vars);

                    gui::SceneBuilder vis;
                    Visualize(vis, view, mg, controls, vars);
                    vis.show(true, false, 
                        gui::RenderOptions().camera(PerspectiveCamera(500, 500, Point2(250, 250), 300, Point3(1, 1, 1), Point3(0, 0, 0))));
                }

            }




            // seg topo & t-structs
            SegmentationTopo segtopo;
            std::vector<std::vector<Vec3>> bndSamples;
            std::vector<int> bndClasses;
            std::vector<TStructure> tstructs;
            if (0 || !Load(path, "segtopo", segtopo, bndSamples, bndClasses, tstructs)) {
                segtopo = SegmentationTopo(segmentedImage, true);
                bndSamples = SamplesOnBoundaries(segtopo, view.camera, DegreesToRadians(1));
                bndClasses = ClassifyBoundaries(bndSamples, vps, DegreesToRadians(1.5));
                tstructs = FindTStructuresFuzzy(segtopo, bndSamples, bndClasses, view.camera, vps, DegreesToRadians(10), DegreesToRadians(2));

                // occ from T-struct
                auto occtstructs = GuessOccludedTStructures(vps, segtopo, tstructs, view.camera, gc);
                std::vector<TStructure> occts;
                for (int id : occtstructs) {
                    occts.push_back(tstructs[id]);
                }
                tstructs = std::move(occts);
                Save(path, "segtopo", segtopo, bndSamples, bndClasses, tstructs);
            }


            auto canvas = gui::MakeCanvas(Image3ub(image.size(), Vec3ub()));
            ShowTStructures(canvas.color(gui::White), segtopo, tstructs);
            canvas.show();


            if(0){
                auto canvas = gui::MakeCanvas(Image3ub(image.size(), Vec3ub()));
                gui::ColorTable vpctable = gui::ColorTableDescriptor::RGBGreys;
                for (int i = 0; i < bndSamples.size(); i++) {
                    canvas.color(vpctable[bndClasses[i]]);
                    canvas.add(segtopo.bndpixels[i]);
                    auto & ps = segtopo.bndpixels[i];
                    if (ps.empty()) {
                        continue;
                    }
                    canvas.add(NoteAs(ps[ps.size() / 2], std::to_string(i)));
                }
                canvas.show();
            }





            RLGraph mg;
            std::vector<RegionHandle> rhs;
            std::vector<RegionBoundaryHandle> bhs;

            if (1 || !Load(path, "mg_rhs_bhs", mg, rhs, bhs)) {
                for (int i = 0; i < cams.size(); i++) {
                    AppendLines(mg, lines[i], cams[i], vps);
                }
                std::tie(rhs, bhs) = AppendRegions(mg, segmentedImage, segtopo.bndpixels, segtopo.bnd2segs, 
                    view.camera, 0.03, 0.02, 4, 4, false);                
                //rhs = AppendRegions(mg, segmentedImage, view.camera, 0.03, 0.02, 4, 4, false);
                Save(path, "mg_rhs_bhs", mg, rhs, bhs);
            }


            RLGraphControls controls;
            // initialize controls
            if (1 || !Load(path, "controls", controls)) {

                controls = RLGraphControls(mg, vps);
                AttachPrincipleDirectionConstraints(mg, controls, M_PI / 30.0, false);
                AttachWallConstriants(mg, controls, M_PI / 50);
                ApplyOccludedTStructure(mg, controls, segtopo, bhs, view.camera, tstructs);

                // gc
                auto gcMeanOnRegions = CollectFeatureMeanOnRegions(mg, view.camera, gc);
                auto occOnBoundaries = CollectOcclusionResponseOnBoundaries(mg, occbnds, view.camera);
                auto up = normalize(vps[vertVPId]);
                if (up.dot(-view.camera.up()) < 0) {
                    up = -up;
                }
                SetComponentControl<RegionData>(controls,
                    [&mg, &controls, &vps, &gcMeanOnRegions, &view, &up, vertVPId](RegionHandle rh, RLGraphComponentControl & c) {
                    auto & gcMean = gcMeanOnRegions[rh];
                    double gcMeanSum = std::accumulate(std::begin(gcMean.val), std::end(gcMean.val), 0.0);
                    //if (gcMeanSum == 0.0) {
                    //    c.orientationClaz = -1;
                    //    c.orientationNotClaz = -1;
                    //    c.used = true;
                    //    return;
                    //}
                    size_t maxid = std::max_element(std::begin(gcMean), std::end(gcMean)) - gcMean.val;
                    double floorScore = gcMean[ToUnderlying(GeometricContextIndex::FloorOrGround)];
                    if (maxid == ToUnderlying(GeometricContextIndex::FloorOrGround) && mg.data(rh).normalizedCenter.dot(up) < 0) { // lower
                        // assign horizontal constriant
                        c.orientationClaz = vertVPId;
                        c.orientationNotClaz = -1;
                        c.used = true;
                        return;
                    }
                    double wallScore = gcMean[ToUnderlying(GeometricContextIndex::Vertical)];
                    if (wallScore > 0.6) {
                        // assign vertical constraint
                        c.orientationClaz = -1;
                        c.orientationNotClaz = vertVPId;
                        c.used = true;
                        return;
                    }
                    /*double clutterScore = gcMean[ToUnderlying(GeometricContextIndex::ClutterOrPorous)];
                    if (clutterScore > 0.7) {
                        c.orientationClaz = c.orientationNotClaz = -1;
                        c.used = false;
                        return;
                    }*/
                });


                {
                    auto pim = Print(mg, segmentedImage, view.camera, rhs,
                        [&controls](RegionHandle rh) -> gui::Color {
                        static const gui::ColorTable ctable = gui::ColorTableDescriptor::RGBGreys;
                        if (!controls[rh].used) {
                            return gui::Black;
                        }
                        if (controls[rh].orientationClaz != -1) {
                            return ctable[controls[rh].orientationClaz];
                        }
                        if (controls[rh].orientationNotClaz != -1) {
                            return ctable[controls[rh].orientationNotClaz] * 0.2;
                        }
                        return gui::Transparent;
                    },
                        core::ConstantFunctor<gui::ColorTag>(gui::Transparent),
                        core::ConstantFunctor<gui::ColorTag>(gui::Transparent),
                      2);
                    cv::imshow("occ", Image3f(pim * 0.99f) + Image3f(image / 255.0f * 0.01f));
                    cv::waitKey();
                }


                controls.disableAllInvalidConstraints(mg);

                {
                    auto pim = Print(mg, segmentedImage, view.camera, rhs,
                        core::ConstantFunctor<gui::ColorTag>(gui::Transparent),
                        core::ConstantFunctor<gui::ColorTag>(gui::Transparent),
                        [&occOnBoundaries](RegionBoundaryHandle rrh) {
                        return occOnBoundaries[rrh] ? gui::Yellow : gui::Transparent;
                    }, 2);
                    cv::imshow("occ", Image3f(pim / 2.0f) + Image3f(image / 255.0f / 2.0f));
                    cv::waitKey();
                }


                // set constraint weights
                SetNecessaryConstraintWeightedAnchors(mg, controls);

                Save(path, "controls", controls);
            }



            RLGraphVars vars;

            // consider both lines and regions
            if (1){

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
                    vars = SolveVariablesWithBoundedAnchors(matlab, mg, controls, false, 1000);
                    NormalizeVariables(mg, controls, vars);
                    //SolveVariablesWithBoundedAnchors(mg, controls, vars);
                    std::cout << "score = " << Score(mg, controls, vars) << std::endl;

                    LooseOrientationConstraintsOnComponents(mg, controls, vars, 0.2, 0.01, 0.1);
                    if (!AttachWeightedAnchorToCenterOfLargestLineIfNoExists(mg, controls))
                        continue;

                    vars = SolveVariablesWithBoundedAnchors(matlab, mg, controls, false, 1000);
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

                        auto range = EstimateEffectiveRangeAlongDirection(polygons, vertDir, medianDepth * 0.01, 0.99, -1e5, -1e5);

                        std::vector<Chain3> chains;
                        for (double x = range.first; x <= range.second; x += medianDepth * 0.01){
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