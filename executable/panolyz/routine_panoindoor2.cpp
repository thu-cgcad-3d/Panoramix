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

    namespace PanoramaIndoor2 {

        static const bool refresh = false;

        void Run() {

            using namespace pano;
            using namespace core;
            using namespace experimental;

            misc::Matlab matlab;

            std::vector<std::string> paths;
            std::vector<Image> images = gui::PickImages(PROJECT_TEST_DATA_DIR_STR"/panorama/indoor", &paths);
            assert(paths.size() == images.size());

            for (int k = 0; k < images.size(); k++) {
                const std::string & path = paths[k];
                const Image3ub & original = images[k];
                if (original.empty())
                    continue;


                Image3ub image;
                bool extendedOnTop = false, extendedOnBottom = false;

                static const bool automaticallyRectifyIfNeeded = true;                
                if (refresh || !Load(path, "rectified", image, extendedOnTop, extendedOnBottom)) {
                    image = original.clone();
                    if (!automaticallyRectifyIfNeeded) {
                        gui::MakePanoramaByHand(image, &extendedOnTop, &extendedOnBottom);
                        cv::imshow("rectified", image);
                        cv::waitKey();
                    } else {
                        MakePanorama(image, -1, &extendedOnTop, &extendedOnBottom);
                    }
                    Save(path, "rectified", image, extendedOnTop, extendedOnBottom);
                }


        
                //MakePanorama(image);
                ResizeToHeight(image, 700);

                View<PanoramicCamera, Image3ub> view;

                std::vector<PerspectiveCamera> cams;
                std::vector<std::vector<Classified<Line2>>> lines;
                std::vector<Vec3> vps;
                int vertVPId;

                Imagei segmentedImage;

                if (refresh || !Load(path, "pre", view, cams, lines, vps, segmentedImage, vertVPId)) {
                    view = CreatePanoramicView(image);

                    // collect lines in each view
                    cams = CreateCubicFacedCameras(view.camera, image.rows, image.rows, image.rows * 0.4);
                    lines.resize(cams.size());
                    for (int i = 0; i < cams.size(); i++) {
                        auto pim = view.sampled(cams[i]).image;
                        LineSegmentExtractor lineExtractor;
                        lineExtractor.params().algorithm = LineSegmentExtractor::LSD;
                        auto ls = lineExtractor(pim, 2, 300); // use pyramid
                        lines[i].reserve(ls.size());
                        for (auto & l : ls) {
                            lines[i].push_back(ClassifyAs(l, -1));
                        }
                    }

                    // estimate vp
                    vps = EstimateVanishingPointsAndClassifyLines(cams, lines);
                    vertVPId = NearestDirectionId(vps, Vec3(0, 0, 1));

                    // extract lines from segmentated region boundaries and classify them using estimated vps
                    // make 3d lines
                    std::vector<Line3> line3ds;
                    for (int i = 0; i < cams.size(); i++) {
                        for (auto & l : lines[i]) {
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


                if (refresh) {
                    auto ctable = gui::CreateRandomColorTableWithSize(vps.size());
                    for (int i = 0; i < cams.size(); i++) {
                        auto pim = view.sampled(cams[i]).image;
                        gui::AsCanvas(pim).thickness(3).colorTable(ctable).add(lines[i]).show();
                    }
                }

                if (refresh) {
                    auto ctable = gui::CreateRandomColorTableWithSize(MinMaxValOfImage(segmentedImage).second + 1);
                    //gui::AsCanvas(ctable(segmentedImage)).show();
                    gui::AsCanvas(ctable(segmentedImage)).add(view.image).show();
                }


                std::vector<PerspectiveCamera> hcams;
                std::vector<Weighted<View<PerspectiveCamera, Image5d>>> gcs;
                Image5d gc;
                {
                    static const int hcamNum = 16;
                    static const Sizei hcamScreenSize(500, 500);
                    static const int hcamFocal = 200;

                    std::string hcamsgcsFileName;
                    {
                        std::stringstream ss;
                        ss << "hcamsgcs_" << hcamNum << "_" << hcamScreenSize.width << "_" << hcamScreenSize.height << "_" << hcamFocal;
                        hcamsgcsFileName = ss.str();
                    }
                    if (refresh || !Load(path, hcamsgcsFileName, hcams, gcs)) {
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
                        Save(path, hcamsgcsFileName, hcams, gcs);
                    }
                    if (refresh || !Load(path, "gcmerged", gc)) {
                        gc = Combine(view.camera, gcs).image;
                        Save(path, "gcmerged", gc);
                    }
                }

                static const bool biggerGC = false;
                if (biggerGC) {
                    std::vector<PerspectiveCamera> hcams2;
                    std::vector<Weighted<View<PerspectiveCamera, Image5d>>> gcs2;
                    Image5d gc2;
                    {
                        static const int hcamNum = 16;
                        static const Sizei hcamScreenSize(500, 700);
                        static const int hcamFocal = 200;

                        std::string hcamsgcsFileName;
                        {
                            std::stringstream ss;
                            ss << "hcamsgcs_" << hcamNum << "_" << hcamScreenSize.width << "_" << hcamScreenSize.height << "_" << hcamFocal;
                            hcamsgcsFileName = ss.str();
                        }
                        if (refresh || !Load(path, hcamsgcsFileName, hcams2, gcs2)) {
                            // extract gcs
                            hcams2 = CreateHorizontalPerspectiveCameras(view.camera, hcamNum, hcamScreenSize.width, hcamScreenSize.height, hcamFocal);
                            gcs2.resize(hcams2.size());
                            for (int i = 0; i < hcams2.size(); i++) {
                                auto pim = view.sampled(hcams2[i]);
                                auto pgc = ComputeGeometricContext(matlab, pim.image, false, true);
                                gcs2[i].component.camera = hcams2[i];
                                gcs2[i].component.image = pgc;
                                gcs2[i].score = abs(1.0 - normalize(hcams2[i].forward()).dot(normalize(view.camera.up())));
                            }
                            Save(path, hcamsgcsFileName, hcams2, gcs2);
                        }
                        if (refresh || !Load(path, "gcmerged2", gc2)) {
                            gc2 = Combine(view.camera, gcs2).image;
                            Save(path, "gcmerged2", gc2);
                        }
                    }
                }


                if (1) {
                    std::vector<Imaged> gcChannels;
                    cv::split(gc, gcChannels);
                    gui::AsCanvas(ConvertToImage3d(gc)).show();
                    misc::MAT gcMat("./cache/gcMat.mat", misc::MAT::Write);
                    gcMat.setVar("gc", gc);
                }


                // seg topo
                SegmentationTopo segtopo;
                std::vector<std::vector<Vec3>> bndSamples;
                std::vector<int> bndClasses;
                if (1 || !Load(path, "segtopo2", segtopo, bndSamples, bndClasses)) {
                    segtopo = SegmentationTopo(segmentedImage, true);
                    bndSamples = SamplesOnBoundaries(segtopo, view.camera, DegreesToRadians(1));
                    bndClasses = ClassifyBoundaries(bndSamples, vps, DegreesToRadians(1.5));
                    Save(path, "segtopo2", segtopo, bndSamples, bndClasses);
                }


                if (1) {
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
                    std::tie(rhs, bhs) = AppendRegions2(mg, segmentedImage, segtopo.bndpixels, segtopo.bnd2segs,
                        view.camera, 0.03, 0.02);
                    //std::tie(rhs, bhs) = AppendRegions(mg, segmentedImage, segtopo.bndpixels, segtopo.bnd2segs,
                    //    view.camera, 0.03, 0.02, 1, 1, false); // we have to reserve bnds under lines for reasoning
                    //rhs = AppendRegions(mg, segmentedImage, view.camera, 0.03, 0.02, 4, 4, false);
                    Save(path, "mg_rhs_bhs", mg, rhs, bhs);
                }


                RLGraphControls controls;
                // initialize controls
                if (1 || !Load(path, "controls", controls)) {

                    controls = RLGraphControls(mg, vps);
                    AttachPrincipleDirectionConstraints(mg, controls, M_PI / 40.0, false);
                    AttachWallConstriants(mg, controls, M_PI / 100.0);

                    // top bottom
                    if (extendedOnTop) {
                        int topSegId = segmentedImage(Pixel(0, 0));
                        RegionHandle topRh = rhs[topSegId];
                        if (topRh.valid()) {
                            controls[topRh].used = false;
                        }
                    }
                    if (extendedOnBottom) {
                        int bottomSegId = segmentedImage(Pixel(0, segmentedImage.rows-1));
                        RegionHandle bottomRh = rhs[bottomSegId];
                        if (bottomRh.valid()) {
                            controls[bottomRh].used = false;
                        }
                    }


                    // gc                
                    auto gcMeanOnRegions = CollectFeatureMeanOnRegions(mg, view.camera, gc);
                    //auto occOnBoundaries = CollectOcclusionResponseOnBoundaries(mg, occbnds, view.camera);
                    auto up = normalize(vps[vertVPId]);
                    if (up.dot(-view.camera.up()) < 0) {
                        up = -up;
                    }
                    std::cout << "use gc" << std::endl;

                    SetComponentControl<RegionData>(controls,
                        [&mg, &controls, &vps, &gcMeanOnRegions, &view, &up, vertVPId](RegionHandle rh, RLGraphComponentControl & c) {
                        auto & gcMean = gcMeanOnRegions[rh];
                        double gcMeanSum = std::accumulate(std::begin(gcMean.val), std::end(gcMean.val), 0.0);
                        assert(gcMeanSum <= 1.1);
                        if (gcMeanSum == 0.0) {
                            /*if (mg.data(rh).normalizedCenter.dot(up) < 0) {
                                c.orientationClaz = vertVPId;
                                c.orientationNotClaz = -1;
                                c.used = true;
                                }*/
                            return;
                        }

                        std::vector<size_t> orderedIds(std::distance(std::begin(gcMean), std::end(gcMean)));
                        std::iota(orderedIds.begin(), orderedIds.end(), 0ull);
                        std::sort(orderedIds.begin(), orderedIds.end(), [&gcMean](size_t a, size_t b) {return gcMean[a] > gcMean[b]; });

                        size_t maxid = orderedIds.front();
                        if (maxid == ToUnderlying(GeometricContextIndex::ClutterOrPorous) && gcMean[maxid] < 0.5) {
                            maxid = orderedIds[1];
                        }

                        if (maxid == ToUnderlying(GeometricContextIndex::FloorOrGround)/* && mg.data(rh).normalizedCenter.dot(up) < 0*/) { // lower
                            // assign horizontal constriant
                            c.orientationClaz = vertVPId;
                            c.orientationNotClaz = -1;
                            c.used = true;
                            return;
                        }
                        if (maxid == ToUnderlying(GeometricContextIndex::CeilingOrSky)) {
                            c.orientationClaz = vertVPId;
                            c.orientationNotClaz = -1;
                            c.used = true;
                            return;
                        }
                        double wallScore = gcMean[ToUnderlying(GeometricContextIndex::Vertical)];
                        if (wallScore > 0.5 && mg.data(rh).normalizedCenter.dot(up) < 0) {
                            c.orientationClaz = -1;
                            c.orientationNotClaz = vertVPId;
                            c.used = true;
                            return;
                        }
                    });

                    std::cout << "reasoning occlusions" << std::endl;

                    auto ocontrolsForOcclusionDetection = mg.createComponentTable<RegionData, OrientationControl>();
                    for (auto it = ocontrolsForOcclusionDetection.begin(); it != ocontrolsForOcclusionDetection.end(); ++it) {
                        auto & gcMean = gcMeanOnRegions[it.hd()];
                        double gcMeanSum = std::accumulate(std::begin(gcMean.val), std::end(gcMean.val), 0.0);
                        assert(gcMeanSum <= 1.1);

                        auto rh = it.hd();
                        auto & c = *it;
                        c.used = controls[rh].used;
                        c.orientationClaz = controls[rh].orientationClaz;
                        c.orientationNotClaz = controls[rh].orientationNotClaz;

                        //if (gcMeanSum == 0.0) {
                        //    continue;
                        //}

                        //size_t maxid = std::max_element(std::begin(gcMean), std::end(gcMean)) - gcMean.val;
                        //if (maxid == ToUnderlying(GeometricContextIndex::ClutterOrPorous) && gcMean[maxid] <= 0.5) {
                        //    maxid = ToUnderlying(std::max({
                        //        GeometricContextIndex::FloorOrGround,
                        //        GeometricContextIndex::CeilingOrSky,
                        //        GeometricContextIndex::Vertical
                        //    }, [&gcMean](GeometricContextIndex a, GeometricContextIndex b) {
                        //        return gcMean[ToUnderlying(a)] < gcMean[ToUnderlying(b)];
                        //    }));
                        //}

                        //if (maxid == ToUnderlying(GeometricContextIndex::FloorOrGround) && mg.data(rh).normalizedCenter.dot(up) < 0) { // lower
                        //    // assign horizontal constriant
                        //    c.orientationClaz = vertVPId;
                        //    c.orientationNotClaz = -1;
                        //    c.used = true;
                        //    continue;
                        //}
                        //double wallScore = gcMean[ToUnderlying(GeometricContextIndex::Vertical)];
                        //if (wallScore > 0.7 && mg.data(rh).normalizedCenter.dot(up) < 0) {
                        //    c.orientationClaz = -1;
                        //    c.orientationNotClaz = vertVPId;
                        //    c.used = true;
                        //    continue;
                        //}
                    }
                    //auto occlusions = DetectOcclusions4(mg, ocontrolsForOcclusionDetection, segmentedImage, segtopo,
                    //bndSamples, bndClasses, rhs, bhs, vps, DegreesToRadians(10), DegreesToRadians(1));

                    auto occlusions = DetectOcclusions5(mg, ocontrolsForOcclusionDetection, segmentedImage, segtopo,
                        bndSamples, bndClasses, rhs, bhs, vps, DegreesToRadians(1.5), DegreesToRadians(0.5));

                    std::cout << "applying occlusions" << std::endl;
                    ApplyOcclusions(mg, controls, occlusions, false);

                    DisableTJunctionsInLineRelations(mg, controls, 0.15);


                    {
                        auto pim = Print(mg, segmentedImage, view.camera, rhs,
                            [&ocontrolsForOcclusionDetection](RegionHandle rh) -> gui::Color {
                            static const gui::ColorTable ctable = gui::ColorTableDescriptor::RGBGreys;
                            if (!ocontrolsForOcclusionDetection[rh].used) {
                                return gui::Black;
                            }
                            if (ocontrolsForOcclusionDetection[rh].orientationClaz != -1) {
                                return ctable[ocontrolsForOcclusionDetection[rh].orientationClaz];
                            }
                            if (ocontrolsForOcclusionDetection[rh].orientationNotClaz != -1) {
                                return ctable[ocontrolsForOcclusionDetection[rh].orientationNotClaz] * 0.2;
                            }
                            return gui::White;
                        },
                            core::ConstantFunctor<gui::ColorTag>(gui::Transparent),
                            [&occlusions](RegionBoundaryHandle rrh) {
                            switch (occlusions[rrh]) {
                            case DepthRelation::FirstIsFront:
                                return gui::Green;
                            case DepthRelation::SecondIsFront:
                                return gui::Blue;
                            case DepthRelation::MaybeFolder:
                                return gui::White;
                            default: return gui::Transparent;
                            }
                        }, 2);
                        cv::imshow("occ", Image3f(pim * 0.99f) + Image3f(image / 255.0f * 0.01f));
                        cv::waitKey();
                    }

                    controls.disableAllInvalidConstraints(mg);

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
                            return gui::White;
                        },
                            core::ConstantFunctor<gui::ColorTag>(gui::Transparent),
                            [&controls](RegionBoundaryHandle rrh) {
                            return controls[rrh].used ? gui::Transparent : gui::Yellow;
                        }, 2);
                        cv::imshow("occ2", Image3f(pim * 0.99f) + Image3f(image / 255.0f * 0.01f));
                        cv::waitKey();
                    }

                    // set constraint weights
                    SetNecessaryConstraintWeightedAnchors(mg, controls);
                    //SetFullConstraintWeightedAnchors(mg, controls);
                    //SetMoreConstraintWeightedAnchors(mg, controls, 0.02);

                    Save(path, "controls", controls);
                }



                RLGraphVars vars;

                // consider both lines and regions
                if (1) {

                    // cc decompose
                    auto ccids = MakeHandledTableForAllComponents(mg, -1);
                    int ccnum = ConnectedComponents(mg, controls, ccids, [](const RLGraphConstraintControl & c) {
                        return c.used && c.weightedAnchors.size() > 0;
                    });
                    auto mgs = Decompose(mg, ccids, ccnum);
                    auto cs = Decompose(mg, controls, ccids, ccnum);
                    assert(mgs.size() == cs.size());


                    for (int i = 0; i < ccnum; i++) {
                        auto & mg = mgs[i];
                        auto & controls = cs[i];
                        RLGraphVars vars;

                        if (!AttachWeightedAnchorToCenterOfLargestLineIfNoExists(mg, controls))
                            continue;

                        //vars = SolveVariablesWithBoundedAnchors(mg, controls, false, true);
                        ResetToSampledArmorAnchors(mg, controls, 0.05);

                        ///vars = SolveVariablesWithoutBoundedAnchors(mg, controls, false);
                        vars = SolveVariablesWithBoundedAnchors(matlab, mg, controls, true, 1000);
                        ///vars = SolveVariablesWithBoundedAnchors2(matlab, mg, controls, false, 10);

                        NormalizeVariables(mg, controls, vars);
                        //SolveVariablesWithBoundedAnchors(mg, controls, vars);
                        std::cout << "score = " << Score(mg, controls, vars) << std::endl;

                        LooseOrientationConstraintsOnComponents(mg, controls, vars, 0.2, 0.01, 0.1);
                        if (!AttachWeightedAnchorToCenterOfLargestLineIfNoExists(mg, controls))
                            continue;

                        ///vars = SolveVariablesWithoutBoundedAnchors(mg, controls, false);
                        vars = SolveVariablesWithBoundedAnchors(matlab, mg, controls, true, 1000);
                        ///vars = SolveVariablesWithBoundedAnchors2(matlab, mg, controls, false, 10);

                        NormalizeVariables(mg, controls, vars);



                        gui::SceneBuilder vis;
                        Visualize(vis, view, mg, controls, vars);
                        vis.show(true, false, gui::RenderOptions().cullFrontFace(false).cullBackFace(true).bwColor(0.0).bwTexColor(1.0)
                            .camera(PerspectiveCamera(500, 500, Point2(250, 250), 300, Point3(1, 1, 1), Point3(0, 0, 0))));

                        {
                            LayeredShape3 shape;
                            auto polygons = RegionPolygons(mg, controls, vars);
                            double medianDepth = MedianCenterDepth(mg, controls, vars);
                            Vec3 vertDir = normalize(controls.vanishingPoints[vertVPId]);

                            auto range = EstimateEffectiveRangeAlongDirection(polygons, vertDir, medianDepth * 0.4, 0.3, -1e5, -1e5);

                            std::vector<Chain3> chains;
                            for (double x = range.first; x <= range.second; x += medianDepth * 0.01) {
                                Plane3 cutplane(vertDir * x, vertDir);
                                auto loop = experimental::MakeSectionalPieces(polygons, cutplane);
                                if (loop.empty())
                                    continue;
                                chains.push_back(experimental::MakeChain(loop));
                            }

                            {
                                gui::SceneBuilder vis;
                                vis.installingOptions().discretizeOptions.color = gui::ColorTag::Black;
                                vis.installingOptions().lineWidth = 1.0;
                                vis.begin(chains).shaderSource(gui::OpenGLShaderSourceDescriptor::XLines).lineWidth(1.0).end();
                                Visualize(vis, view, mg, controls, vars);
                                vis.show(true, false,
                                    gui::RenderOptions().camera(PerspectiveCamera(500, 500, Point2(250, 250), 300, Point3(1, 1, 1), Point3(0, 0, 0))));
                            }


                            for (int i = 0; i < chains.size(); i++) {
                                shape.layers.push_back(std::move(chains[i].points));
                            }
                            shape.normal = vertDir;

                            {
                                gui::SceneBuilder viz;
                                gui::ResourceStore::set("texture", image);
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
                    }

                    //SaveToDisk("./cache/mgp", mg, controls, vars);
                } else {
                    //LoadFromDisk("./cache/mgp", mg, controls, vars);
                }

            }

        }
    }
}