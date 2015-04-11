#include "../core/algorithms.hpp"
#include "tools.hpp"

namespace panoramix {

    namespace experimental {


        namespace {

            template <class UhClickHandlerFunT, class UhColorizerFunT = core::ConstantFunctor<gui::Color>>
            void VisualizeTemplated(gui::Visualizer & viz, const core::Image & panorama,
                const RLGraph & mg,
                const RLGraphControls & controls, const RLGraphVars & vars,
                UhClickHandlerFunT && uhClicked,
                UhColorizerFunT && uhColorizer) {

                struct ComponentID {
                    int handleID;
                    bool isRegion;
                };

                auto sppCallbackFun = [&controls, &mg, &uhClicked](gui::InteractionID iid,
                    const std::pair<ComponentID, gui::Colored<gui::SpatialProjectedPolygon>> & spp) {
                    std::cout << (spp.first.isRegion ? "Region" : "Line") << spp.first.handleID << std::endl;
                    if (spp.first.isRegion){
                        uhClicked(RegionHandle(spp.first.handleID));
                    }
                    else {
                        uhClicked(LineHandle(spp.first.handleID));
                    }
                };

                gui::ResourceStore::set("texture", panorama);

                viz.renderOptions.bwColor = 1.0;
                viz.renderOptions.bwTexColor = 0.0;
                viz.installingOptions.discretizeOptions.colorTable = gui::ColorTableDescriptor::RGB;
                std::vector<std::pair<ComponentID, gui::Colored<gui::SpatialProjectedPolygon>>> spps;
                std::vector<gui::Colored<core::Line3>> lines;

                for (auto & c : mg.components<RegionData>()){
                    if (!controls[c.topo.hd].used)
                        continue;
                    auto uh = c.topo.hd;
                    auto & region = c.data;
                    gui::SpatialProjectedPolygon spp;
                    // filter corners
                    core::ForeachCompatibleWithLastElement(c.data.normalizedContours.front().begin(), c.data.normalizedContours.front().end(),
                        std::back_inserter(spp.corners),
                        [](const core::Vec3 & a, const core::Vec3 & b) -> bool {
                        return core::AngleBetweenDirections(a, b) > M_PI / 1000.0;
                    });
                    if (spp.corners.size() < 3)
                        continue;

                    spp.projectionCenter = core::Point3(0, 0, 0);
                    spp.plane = Instance(mg, controls, vars, uh);
                    assert(!HasValue(spp.plane, IsInfOrNaN<double>));
                    spps.emplace_back(ComponentID{ uh.id, true }, std::move(gui::ColorAs(spp, uhColorizer(uh))));
                }

                for (auto & c : mg.components<LineData>()){
                    if (!controls[c.topo.hd].used)
                        continue;
                    auto uh = c.topo.hd;
                    auto & line = c.data;
                    lines.push_back(gui::ColorAs(Instance(mg, controls, vars, uh), uhColorizer(uh)));
                }

                viz.begin(spps, sppCallbackFun).shaderSource(gui::OpenGLShaderSourceDescriptor::XPanorama).resource("texture").end();
                viz.installingOptions.discretizeOptions.color = gui::ColorTag::DarkGray;
                viz.installingOptions.lineWidth = 5.0;
                viz.add(lines);

                std::vector<core::Line3> connectionLines;
                for (auto & c : mg.constraints<RegionBoundaryData>()){
                    if (!controls[c.topo.hd].used)
                        continue;
                    auto & samples = c.data.normalizedSampledPoints;
                    auto inst1 = Instance(mg, controls, vars, c.topo.component<0>());
                    auto inst2 = Instance(mg, controls, vars, c.topo.component<1>());
                    for (auto & ss : samples){
                        for (auto & s : ss){
                            double d1 = DepthAt(s, inst1);
                            double d2 = DepthAt(s, inst2);
                            connectionLines.emplace_back(normalize(s) * d1, normalize(s) * d2);
                        }
                    }
                }
                for (auto & c : mg.constraints<RegionLineConnectionData>()){
                    if (!controls[c.topo.hd].used)
                        continue;
                    auto inst1 = Instance(mg, controls, vars, c.topo.component<0>());
                    auto inst2 = Instance(mg, controls, vars, c.topo.component<1>());
                    for (auto & s : c.data.normalizedAnchors){
                        double d1 = DepthAt(s, inst1);
                        double d2 = DepthAt(s, inst2);
                        connectionLines.emplace_back(normalize(s) * d1, normalize(s) * d2);
                    }
                }

                viz.installingOptions.discretizeOptions.color = gui::ColorTag::Black;
                viz.installingOptions.lineWidth = 1.0;
                viz.add(connectionLines);

                viz.renderOptions.renderMode = gui::RenderModeFlag::Triangles | gui::RenderModeFlag::Lines;
                viz.renderOptions.backgroundColor = gui::ColorTag::White;
                viz.renderOptions.bwColor = 0.0;
                viz.renderOptions.bwTexColor = 0.5;

                gui::ResourceStore::clear();

            }


        }


        struct ComponentHandleColorizer {
            const RLGraph & mg;
            const RLGraphControls & controls;
            const RLGraphVars & vars;
            gui::ColorTable colorTableForVPs;
            gui::ColorTable colorTableForRegionAlongOrientations;
            ComponentHandleColorizer(const RLGraph & g, const RLGraphControls & c, const RLGraphVars & v)
                : mg(g), controls(c), vars(v) {
                colorTableForVPs = gui::ColorTable(gui::ColorTableDescriptor::RGB)
                    .appendRandomizedColors(controls.vanishingPoints.size() - 3);
                colorTableForRegionAlongOrientations = gui::CreateGreyColorTableWithSize(controls.vanishingPoints.size());
            }

            inline gui::Color operator()(ComponentHandle<RegionData> rh) const{
                auto & prop = controls[rh];
                if (!prop.used)
                    return gui::ColorTag::Yellow;
                if (prop.orientationClaz >= 0)
                    return colorTableForVPs[prop.orientationClaz];
                if (prop.orientationNotClaz >= 0)
                    return colorTableForRegionAlongOrientations[prop.orientationNotClaz];
                return gui::ColorTag::Yellow;
            }
            inline gui::Color operator()(ComponentHandle<LineData> rh) const {
                return colorTableForVPs[controls[rh].orientationClaz];
            }
        };

        struct ComponentHandleClickHandler {
            const RLGraph & mg;
            const RLGraphControls & controls;
            const RLGraphVars & vars;
            inline bool operator()(ComponentHandle<RegionData> rh) const{
                std::cout << "area: " << mg.data(rh).area
                    << " connected regions: " << mg.topo(rh).constraints<RegionBoundaryData>().size()
                    << " connected lines: " << mg.topo(rh).constraints<RegionLineConnectionData>().size() << std::endl;
                return false;
            }
            inline bool operator()(ComponentHandle<LineData> rh) const {
                return false;
            }
        };


        void Visualize(gui::Visualizer & vis, const View<PanoramicCamera> & texture,
            const RLGraph & mg, const RLGraphControls & controls, const RLGraphVars & vars){
            VisualizeTemplated(vis, texture.image, mg, controls, vars,
                ComponentHandleClickHandler{ mg, controls, vars },
                ComponentHandleColorizer{ mg, controls, vars });
        }

        void Visualize(gui::Visualizer & vis, const View<PerspectiveCamera> & texture,
            const RLGraph & mg, const RLGraphControls & controls, const RLGraphVars & vars) {
            VisualizeTemplated(vis, texture.sampled(PanoramicCamera(500, Point3(0, 0, 0), Point3(1, 0, 0), Vec3(0, 0, 1))).image,
                mg, controls, vars, ComponentHandleClickHandler{ mg, controls, vars },
                ComponentHandleColorizer{ mg, controls, vars });
        }





    }


}