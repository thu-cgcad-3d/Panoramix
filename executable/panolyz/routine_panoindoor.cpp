#include "../../src/core/basic_types.hpp"
#include "../../src/experimental/rl_graph.hpp"
#include "../../src/experimental/tools.hpp"
#include "../../src/gui/visualize2d.hpp"
#include "../../src/gui/visualizers.hpp"

#include "routines.hpp"


namespace panolyz {

    ROUTINE_FOR_ALGORITHM(PanoramaIndoor){

        std::string path;
        //path = PROJECT_TEST_DATA_DIR_STR"/panorama/indoor/13.jpg";
        //path = PROJECT_TEST_DATA_DIR_STR"/panorama/indoor/14.jpg";
        //path = PROJECT_TEST_DATA_DIR_STR"/panorama/indoor/x3.jpg";
        //path = PROJECT_TEST_DATA_DIR_STR"/panorama/indoor/45.jpg";
        //path = PROJECT_TEST_DATA_DIR_STR"/panorama/indoor/x2.jpg";
        //path = PROJECT_TEST_DATA_DIR_STR"/panorama/outdoor/univ1.jpg";
        //path = PROJECT_TEST_DATA_DIR_STR"/panorama/indoor/k (9).jpg";// too small
        //path = PROJECT_TEST_DATA_DIR_STR"/panorama/outdoor/yard.jpg";
        //path = PROJECT_TEST_DATA_DIR_STR"/panorama/indoor/k (11).jpg"; // too small
        //path = PROJECT_TEST_DATA_DIR_STR"/panorama/indoor/k (10).jpg"; // too small
        //path = PROJECT_TEST_DATA_DIR_STR"/panorama/indoor/k (7).jpg";
        //path = PROJECT_TEST_DATA_DIR_STR"/panorama/indoor/univlab.jpg";
        //path = PROJECT_TEST_DATA_DIR_STR"/panorama/indoor/univlab2.jpg";
        //path = PROJECT_TEST_DATA_DIR_STR"/panorama/indoor/univlab3.jpg";
        //path = PROJECT_TEST_DATA_DIR_STR"/panorama/indoor/k (2).jpg";
        //path = PROJECT_TEST_DATA_DIR_STR"/panorama/indoor/x5.jpg";
        //path = PROJECT_TEST_DATA_DIR_STR"/panorama/indoor/x6.jpg";
        //path = PROJECT_TEST_DATA_DIR_STR"/panorama/indoor/x7.jpg";
        //path = PROJECT_TEST_DATA_DIR_STR"/panorama/indoor/x8.jpg";
        //path = PROJECT_TEST_DATA_DIR_STR"/panorama/indoor/x9.jpg";
        path = PROJECT_TEST_DATA_DIR_STR"/panorama/indoor/google_chinese.png";


        using namespace panoramix;
        using namespace core;
        using namespace experimental;

        Image image = cv::imread(path);
        ResizeToHeight(image, 700);

        View<PanoramicCamera> view;

        std::vector<PerspectiveCamera> cams;
        std::vector<std::vector<Classified<Line2>>> lines;
        std::vector<Vec3> vps;

        // vp1 vp2 vp3 clutter unknown
        ImageOfType<Vec<double, 5>> gc;
        Imagei gcVotes;

        Imagei segmentedImage;

        if (0){
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
                    gui::Visualizer2D(pim)
                        << gui::manip2d::SetThickness(3)
                        << gui::manip2d::SetColorTable(ctable)
                        << lines[i] << gui::manip2d::Show();
                }
            }

            // extract lines from segmentated region boundaries and classify them using estimated vps
            // make 3d lines
            std::vector<Line3> line3ds;
            for (int i = 0; i < cams.size(); i++){
                for (auto & l : lines[i]){
                    line3ds.emplace_back(normalize(cams[i].spatialDirection(l.component.first)),
                        normalize(cams[i].spatialDirection(l.component.second)));
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
                gui::Visualizer2D(ctable(segmentedImage)) << gui::manip2d::Show();
                gui::Visualizer2D(ctable(segmentedImage)) << view.image << gui::manip2d::Show();
            }

            // extract gcs
            IndoorGeometricContextEstimator gcEstimator;
            std::tie(gc, gcVotes) = gcEstimator(view.image, view.camera, vps);

            Save(path, "all", view, cams, lines, vps, gc, gcVotes, segmentedImage);
        }
        else{
            Load(path, "all", view, cams, lines, vps, gc, gcVotes, segmentedImage);
        }







        // consider only lines
        if (0){

            RLGraph mg;
            RLGraphControls controls;
            RLGraphVars vars;

            for (int i = 0; i < cams.size(); i++){
                AppendLines(mg, lines[i], cams[i], vps);
            }

            controls = MakeControls(mg, vps);
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
                vis.camera(PerspectiveCamera(500, 500, Point2(250, 250), 300, Point3(1, 1, 1), Point3(0, 0, 0)));
                vis.show(true, false);
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
            rhs = AppendRegions(mg, segmentedImage, view.camera, 0.01, 0.02, 3, 2);

            controls = MakeControls(mg, vps);
            AttachPrincipleDirectionConstraints(mg, controls);
            AttachWallConstriants(mg, controls);
            AttachGeometricContextConstraints(mg, controls, view.camera, gc, gcVotes);

            // set constraint weights
            SetConstraintWeights<LineRelationData>(controls, [&mg](LineRelationHandle h){
                return 1.0;// std::max(mg.data(h).junctionWeight, 3.0f);
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
            auto mgs = Decompose(mg, ccids, ccnum);
            auto cs = Decompose(mg, controls, ccids, ccnum);
            assert(mgs.size() == cs.size());


            for (int i = 0; i < ccnum; i++){
                auto & mg = mgs[i];
                auto & controls = cs[i];
                RLGraphVars vars;

                if (!AttachAnchorToCenterOfLargestLineIfNoAnchorExists(mg, controls))
                    continue;

                vars = SolveVariables(mg, controls, false);
                NormalizeVariables(mg, controls, vars);
                std::cout << "score = " << Score(mg, controls, vars) << std::endl;

                LooseOrientationConstraintsOnComponents(mg, controls, vars, 0.2, 0.02, 0.1);
                if (!AttachAnchorToCenterOfLargestLineIfNoAnchorExists(mg, controls))
                    continue;

                vars = SolveVariables(mg, controls);
                NormalizeVariables(mg, controls, vars);

                AttachFloorAndCeilingConstraints(mg, controls, vars, 0.1, 0.6);

                if (!AttachAnchorToCenterOfLargestLineIfNoAnchorExists(mg, controls))
                    continue;
                vars = SolveVariables(mg, controls);
                NormalizeVariables(mg, controls, vars);

                gui::Visualizer vis;
                Visualize(vis, view, mg, controls, vars);
                vis.camera(PerspectiveCamera(500, 500, Point2(250, 250), 300, Point3(1, 1, 1), Point3(0, 0, 0)));
                vis.show(true, false);

                {
                    LayeredShape3 shape;
                    auto polygons = RegionPolygons(mg, controls, vars);
                    int vertVPId = GetVerticalDirectionId(controls.vanishingPoints);
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
                    viz.renderOptions.cullFrontFace = true;
                    viz.camera(PerspectiveCamera(1000, 800, Point2(500, 400),
                        800, Point3(-1, 1, 1), Point3(0, 0, 0)));

                    viz.show(true, false);
                }
            }

            //SaveToDisk("./cache/mgp", mg, controls, vars);
        }
        else{
            //LoadFromDisk("./cache/mgp", mg, controls, vars);
        }


    }
}