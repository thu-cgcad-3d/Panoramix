#include "../core/algorithms.hpp"
#include "../core/containers.hpp"
#include "tools.hpp"

namespace pano {

    namespace experimental {


        namespace {

            template <class UhClickHandlerFunT, class UhColorizerFunT = core::ConstantFunctor<gui::Color>>
            void VisualizeTemplated(gui::SceneBuilder & viz, const core::Image & panorama,
                const RLGraph & mg,
                const RLGraphControls & controls, const RLGraphVars & vars,
                UhClickHandlerFunT && uhClicked,
                UhColorizerFunT && uhColorizer) {

                struct ComponentID {
                    int handleID;
                    bool isRegion;
                };

                auto sppCallbackFun = [&controls, &mg, &uhClicked](gui::InteractionID iid,
                    const core::Decorated<gui::Colored<gui::SpatialProjectedPolygon>, ComponentID> & spp) {
                    std::cout << (spp.decoration.isRegion ? "Region" : "Line") << spp.decoration.handleID << std::endl;
                    if (spp.decoration.isRegion){
                        uhClicked(RegionHandle(spp.decoration.handleID));
                    }
                    else {
                        uhClicked(LineHandle(spp.decoration.handleID));
                    }
                };

                gui::ResourceStore::set("texture", panorama);

                viz.installingOptions().discretizeOptions.colorTable = gui::ColorTableDescriptor::RGB;
                std::vector<core::Decorated<gui::Colored<gui::SpatialProjectedPolygon>, ComponentID>> spps;
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
                        return core::AngleBetweenDirections(a, b) > 0.0;
                    });
                    if (spp.corners.size() < 3)
                        continue;

                    spp.projectionCenter = core::Point3(0, 0, 0);
                    spp.plane = Instance(mg, controls, vars, uh);
                    if (HasValue(spp.plane, IsInfOrNaN<double>)){
                        WARNNING("inf plane");
                        continue;
                    }
                    spps.push_back(core::DecorateAs(std::move(gui::ColorAs(spp, uhColorizer(uh))), ComponentID{ uh.id, true }));
                }

                for (auto & c : mg.components<LineData>()){
                    if (!controls[c.topo.hd].used)
                        continue;
                    auto uh = c.topo.hd;
                    auto & line = c.data;
                    auto l = Instance(mg, controls, vars, uh);
                    if (HasValue(l, IsInfOrNaN<double>)){
                        WARNNING("inf line");
                        continue;
                    }
                    lines.push_back(gui::ColorAs(l, uhColorizer(uh)));
                }

                viz.begin(spps/*, sppCallbackFun*/).shaderSource(gui::OpenGLShaderSourceDescriptor::XPanorama).resource("texture").end();
                viz.installingOptions().discretizeOptions.color = gui::ColorTag::Magenta;
                viz.installingOptions().lineWidth = 5.0;
                viz.add(lines);

                std::vector<core::Line3> connectionLinesRR;
                std::vector<core::Point3> connectionLineRREnds;
                for (auto & c : mg.constraints<RegionBoundaryData>()){
                    if (!controls[c.topo.hd].used)
                        continue;                  
                    auto & samples = controls[c.topo.hd].weightedAnchors;
                    auto inst1 = Instance(mg, controls, vars, c.topo.component<0>());
                    auto inst2 = Instance(mg, controls, vars, c.topo.component<1>());
                    for (auto & ss : samples){
                        double d1 = DepthAt(ss.component, inst1);
                        double d2 = DepthAt(ss.component, inst2);
                        connectionLinesRR.emplace_back(normalize(ss.component) * d1, normalize(ss.component) * d2);
                        connectionLineRREnds.push_back(connectionLinesRR.back().first);
                        connectionLineRREnds.push_back(connectionLinesRR.back().second);
                    }
                }

                viz.installingOptions().discretizeOptions.color = gui::ColorTag::Black;
                viz.installingOptions().lineWidth = 1.5;
                viz.add(connectionLinesRR);
                viz.installingOptions().pointSize = 3.0;
                viz.add(connectionLineRREnds);


                std::vector<core::Decorated<core::Line3, RegionLineConnectionHandle>> connectionLinesRL;
                std::vector<core::Decorated<core::Point3, RegionLineConnectionHandle>> connectionLineRLEnds;
                for (auto & c : mg.constraints<RegionLineConnectionData>()){
                    if (!controls[c.topo.hd].used)
                        continue;
                    auto inst1 = Instance(mg, controls, vars, c.topo.component<0>());
                    auto inst2 = Instance(mg, controls, vars, c.topo.component<1>());
                    for (auto & s : controls[c.topo.hd].weightedAnchors) {
                        double d1 = DepthAt(s.component, inst1);
                        double d2 = DepthAt(s.component, inst2);
                        connectionLinesRL.push_back(core::DecorateAs(Line3(normalize(s.component) * d1, normalize(s.component) * d2), c.topo.hd));
                        connectionLineRLEnds.push_back(core::DecorateAs(connectionLinesRL.back().component.first, c.topo.hd));
                        connectionLineRLEnds.push_back(core::DecorateAs(connectionLinesRL.back().component.second, c.topo.hd));
                    }
                }


                viz.installingOptions().discretizeOptions.color = gui::ColorTag::Gray;
                viz.installingOptions().lineWidth = 5;
                viz.add(connectionLinesRL, [](gui::InteractionID iid, const core::Decorated<Line3, RegionLineConnectionHandle> & cline) {
                    std::cout << "rlh.id = " << cline.decoration.id << std::endl;
                });
                viz.installingOptions().pointSize = 3.0;
                viz.add(connectionLineRLEnds);

                //gui::ResourceStore::clear();
                viz.show(true, false, gui::RenderOptions().cullFrontFace(false).cullBackFace(true).bwColor(0.0).bwTexColor(1.0)
                    .camera(PerspectiveCamera(500, 500, Point2(250, 250), 300, Point3(1, 1, 1), Point3(0, 0, 0))));
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
                if (rh.invalid()) {
                    return false;
                }
                std::cout << "area: " << mg.data(rh).area
                    << " connected regions: " << mg.topo(rh).constraints<RegionBoundaryData>().size()
                    << " connected lines: " << mg.topo(rh).constraints<RegionLineConnectionData>().size() << std::endl;
                return false;
            }
            inline bool operator()(ComponentHandle<LineData> lh) const {
                if (lh.invalid()) {
                    return false;
                }
                std::cout << "length: " << AngleBetweenDirections(mg.data(lh).line.first, mg.data(lh).line.second) << std::endl;
                return false;
            }
        };


        void Visualize(gui::SceneBuilder & vis, const View<PanoramicCamera> & texture,
            const RLGraph & mg, const RLGraphControls & controls, const RLGraphVars & vars){
            VisualizeTemplated(vis, texture.image, mg, controls, vars,
                ComponentHandleClickHandler{ mg, controls, vars },
                ComponentHandleColorizer{ mg, controls, vars });
        }

        void Visualize(gui::SceneBuilder & vis, const View<PerspectiveCamera> & texture,
            const RLGraph & mg, const RLGraphControls & controls, const RLGraphVars & vars) {
            VisualizeTemplated(vis, texture.sampled(PanoramicCamera(500, Point3(0, 0, 0), Point3(1, 0, 0), Vec3(0, 0, 1))).image,
                mg, controls, vars, ComponentHandleClickHandler{ mg, controls, vars },
                ComponentHandleColorizer{ mg, controls, vars });
        }








    }


}