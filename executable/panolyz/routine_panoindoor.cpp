#include "../../src/core/basic_types.hpp"
#include "../../src/experimental/rl_graph.hpp"
#include "../../src/gui/visualize2d.hpp"
#include "../../src/gui/visualizers.hpp"

#include "routines.hpp"

namespace panoramix {
    namespace core{
        using namespace experimental;
    }
}


namespace panolyz {

    using namespace panoramix;

    template <class StringT, class ... Ts>
    inline void Save(const std::string & tag, StringT && s, Ts && ... ts){
        core::SaveToDisk("./cache/" + tag + "_" + s, ts...);
    }

    template <class StringT, class ... Ts>
    inline void Load(const std::string & tag, StringT && s, Ts & ... ts){
        core::LoadFromDisk("./cache/" + tag + "_" + s, ts...);
    }

    ROUTINE_FOR_ALGORITHM(PanoramaIndoor_v1){

        core::View<core::PanoramicCamera> view;

        std::vector<core::PerspectiveCamera> cams;
        std::vector<std::vector<core::Classified<core::Line2>>> lines;
        std::vector<core::Vec3> vps;

        // vp1 vp2 vp3 clutter unknown
        core::ImageOfType<core::Vec<double, 5>> gc;
        core::Imagei gcVotes;

        core::Imagei segmentedImage;

        if (0){
            view = core::CreatePanoramicView(image);

            // collect lines in each view
            cams = core::CreateCubicFacedCameras(view.camera, image.rows, image.rows, image.rows * 0.4);
            lines.resize(cams.size());
            for (int i = 0; i < cams.size(); i++){
                auto pim = view.sampled(cams[i]).image;
                core::LineSegmentExtractor lineExtractor;
                lineExtractor.params().algorithm = core::LineSegmentExtractor::LSD;
                auto ls = lineExtractor(pim, 2, 300); // use pyramid
                lines[i].reserve(ls.size());
                for (auto & l : ls){
                    lines[i].push_back(core::ClassifyAs(l, -1));
                }
            }

            // estimate vp
            vps = core::EstimateVanishingPointsAndClassifyLines(cams, lines);
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
            std::vector<core::Line3> line3ds;
            for (int i = 0; i < cams.size(); i++){
                for (auto & l : lines[i]){
                    line3ds.emplace_back(normalize(cams[i].spatialDirection(l.component.first)),
                        normalize(cams[i].spatialDirection(l.component.second)));
                }
            }
            core::SegmentationExtractor segmenter;
            segmenter.params().algorithm = core::SegmentationExtractor::GraphCut;
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
            core::IndoorGeometricContextEstimator gcEstimator;
            std::tie(gc, gcVotes) = gcEstimator(view.image, view.camera, vps);

            Save(tag, "all", view, cams, lines, vps, gc, gcVotes, segmentedImage);
        }
        else{
            Load(tag, "all", view, cams, lines, vps, gc, gcVotes, segmentedImage);
        }


        // consider only lines
        if (0){

            core::RLGraph mg;
            core::RLGraphControls controls;
            core::RLGraphVars vars;

            for(int i = 0; i < cams.size(); i++){
                core::AppendLines(mg, lines[i], cams[i], vps);
            }

            controls = core::MakeControls(mg, vps);
            core::SetConstraintWeights<core::LineRelationData>(controls, [&mg](core::LineRelationHandle h){
                return std::max(mg.data(h).junctionWeight, 1e-1f);
            });

            auto ccids = core::MakeHandledTableForAllComponents(mg, -1);
            int ccnum = core::ConnectedComponents(mg, controls, ccids, [](const core::RLGraphConstraintControl & c){
                return c.used && c.weight > 0;
            });
            auto mgs = core::Decompose(mg, ccids, ccnum);
            auto cs = core::Decompose(mg, controls, ccids, ccnum);
            assert(mgs.size() == cs.size());

            for (int i = 0; i < ccnum; i++){
                auto & mg = mgs[i];
                auto & controls = cs[i];
                core::AttachAnchorToCenterOfLargestLineIfNoAnchorExists(mg, controls);
                auto vars = core::SolveVariables(mgs[i], controls, true);
                core::NormalizeVariables(mg, controls, vars);
                core::Visualize(view, mg, controls, vars);
            }

        }


        core::RLGraph mg;
        core::RLGraphControls controls;
        core::RLGraphVars vars;
        std::vector<core::RegionHandle> rhs;


        // consider both lines and regions
        if (1){

            for (int i = 0; i < cams.size(); i++){
                core::AppendLines(mg, lines[i], cams[i], vps);
            }
            rhs = core::AppendRegions(mg, segmentedImage, view.camera, 0.01, 0.02, 3, 2);

            controls = core::MakeControls(mg, vps);
            core::AttachPrincipleDirectionConstraints(mg, controls);
            core::AttachWallConstriants(mg, controls);
            core::AttachGeometricContextConstraints(mg, controls, view.camera, gc, gcVotes);

            // set constraint weights
            core::SetConstraintWeights<core::LineRelationData>(controls, [&mg](core::LineRelationHandle h){
                return 1.0;// std::max(mg.data(h).junctionWeight, 3.0f);
            });
            core::SetConstraintWeights<core::RegionBoundaryData>(controls, [&mg](core::RegionBoundaryHandle h){
                return std::max(mg.data(h).length / M_PI * 10.0, 1.0);
            });
            core::SetConstraintWeights<core::RegionLineConnectionData>(controls, [&mg](core::RegionLineConnectionHandle h){
                return std::max(mg.data(h).length / M_PI * 10.0, 1.0);
            });

            // cc decompose
            auto ccids = core::MakeHandledTableForAllComponents(mg, -1);
            int ccnum = core::ConnectedComponents(mg, controls, ccids, [](const core::RLGraphConstraintControl & c){
                return c.used && c.weight > 0;
            });
            auto mgs = core::Decompose(mg, ccids, ccnum);
            auto cs = core::Decompose(mg, controls, ccids, ccnum);
            assert(mgs.size() == cs.size());


            for (int i = 0; i < ccnum; i++){
                auto & mg = mgs[i];
                auto & controls = cs[i];
                core::RLGraphVars vars;

                if (!core::AttachAnchorToCenterOfLargestRegionIfNoAnchorExists(mg, controls) &&
                    !core::AttachAnchorToCenterOfLargestLineIfNoAnchorExists(mg, controls))
                    continue;

                vars = core::SolveVariables(mg, controls, false);
                core::NormalizeVariables(mg, controls, vars);
                std::cout << "score = " << core::Score(mg, controls, vars) << std::endl;

                core::LooseOrientationConstraintsOnComponents(mg, controls, vars, 0.2, 0.02, 0.1);
                if (!core::AttachAnchorToCenterOfLargestRegionIfNoAnchorExists(mg, controls) &&
                    !core::AttachAnchorToCenterOfLargestLineIfNoAnchorExists(mg, controls))
                    continue;

                vars = core::SolveVariables(mg, controls);
                core::NormalizeVariables(mg, controls, vars);

                core::AttachFloorAndCeilingConstraints(mg, controls, vars, 0.1, 0.6);

                if (!core::AttachAnchorToCenterOfLargestRegionIfNoAnchorExists(mg, controls) &&
                    !core::AttachAnchorToCenterOfLargestLineIfNoAnchorExists(mg, controls))
                    continue;
                vars = core::SolveVariables(mg, controls);
                core::NormalizeVariables(mg, controls, vars);

                core::Visualize(view, mg, controls, vars);

                {
                    core::LayeredShape3 shape;
                    auto polygons = RegionPolygons(mg, controls, vars);
                    int vertVPId = core::GetVerticalDirectionId(controls.vanishingPoints);
                    double medianDepth = MedianCenterDepth(mg, controls, vars);
                    core::Vec3 vertDir = normalize(controls.vanishingPoints[vertVPId]);

                    auto range = experimental::EstimateEffectiveRangeAlongDirection(polygons, vertDir, medianDepth * 0.02, 0.7, -1e5, -1e5);

                    std::vector<core::Chain3> chains;
                    for (double x = range.first; x <= range.second; x += medianDepth * 0.02){
                        core::Plane3 cutplane(vertDir * x, vertDir);
                        auto loop = experimental::CutRegionLoopAt(polygons, cutplane);
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
                    viz.camera(core::PerspectiveCamera(1000, 800, core::Point2(500, 400),
                        800, core::Point3(-1, 1, 1), core::Point3(0, 0, 0)));

                    viz.show(true, false);
                }
            }

            //core::SaveToDisk("./cache/mgp", mg, controls, vars);
        }
        else{
            //core::LoadFromDisk("./cache/mgp", mg, controls, vars);
        }
}