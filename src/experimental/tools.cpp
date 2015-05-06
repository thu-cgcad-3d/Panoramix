#include "../core/algorithms.hpp"
#include "tools.hpp"

namespace panoramix {

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

                viz.installingOptions().discretizeOptions.colorTable = gui::ColorTableDescriptor::RGB;
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
                    spps.emplace_back(ComponentID{ uh.id, true }, std::move(gui::ColorAs(spp, uhColorizer(uh))));
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

                viz.begin(spps, sppCallbackFun).shaderSource(gui::OpenGLShaderSourceDescriptor::XPanorama).resource("texture").end();
                viz.installingOptions().discretizeOptions.color = gui::ColorTag::DarkGray;
                viz.installingOptions().lineWidth = 5.0;
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

                viz.installingOptions().discretizeOptions.color = gui::ColorTag::Black;
                viz.installingOptions().lineWidth = 1.0;
                viz.add(connectionLines);

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








        std::vector<SectionalPiece> MakeSectionalPieces(const HandledTable<RegionHandle, std::vector<Polygon3>> & polygons,
            const Plane3 & cutplane){

            std::vector<SectionalPiece> segments;
            std::vector<double> startAngles;

            Vec3 x, y;
            std::tie(x, y) = ProposeXYDirectionsFromZDirection(cutplane.normal);
            Point3 original = cutplane.root();

            for (auto it = polygons.begin(); it != polygons.end(); ++it){
                RegionHandle rh = it.hd();

                for (const Polygon3 & polygon : *it){
                    if (polygon.corners.size() <= 2)
                        continue;

                    auto plane = polygon.plane();

                    if (core::IsFuzzyParallel(plane.normal, cutplane.normal, 0.01))
                        continue;

                    Vec3 along = normalize(cutplane.normal.cross(plane.normal));
                    std::vector<Scored<std::pair<Point3, Point2>>> cutpoints;

                    for (int i = 1; i <= polygon.corners.size(); i++){
                        auto & lastP = polygon.corners[i - 1];
                        auto & p = polygon.corners[i % polygon.corners.size()];
                        double lastDist = cutplane.signedDistanceTo(lastP);
                        double dist = cutplane.signedDistanceTo(p);

                        if (lastDist < 0 && dist > 0 || lastDist > 0 && dist < 0){
                            Point3 intersection = (lastP * abs(dist) + p * abs(lastDist)) / (abs(dist) + abs(lastDist));
                            double order = intersection.dot(along);
                            Point2 projOnCutPlane((intersection - original).dot(x), (intersection - original).dot(y));
                            cutpoints.push_back(ScoreAs(std::make_pair(intersection, projOnCutPlane), order));
                        }
                    }

                    if (cutpoints.empty()){
                        continue;
                    }
                    assert(cutpoints.size() % 2 == 0);
                    if (cutpoints.size() % 2 != 0){
                        std::cout << "ODD cutpoints num!!!" << std::endl;
                    }

                    // chain the segments
                    std::sort(cutpoints.begin(), cutpoints.end());
                    for (int i = 0; i < cutpoints.size(); i += 2){
                        SectionalPiece segment;
                        segment.rh = rh;
                        segment.range.first = cutpoints[i].component.first;
                        segment.range.second = cutpoints[i + 1].component.first;
                        segments.push_back(segment);
                        startAngles.push_back(SignedAngleBetweenDirections(Vec2(1, 0), cutpoints[i].component.second));
                    }
                }
            }

            if (segments.empty())
                return segments;

            std::vector<int> ids(segments.size());
            std::iota(ids.begin(), ids.end(), 0);
            std::sort(ids.begin(), ids.end(), [&startAngles](int id1, int id2){return startAngles[id1] < startAngles[id2]; });

            std::vector<SectionalPiece> segs;
            segs.reserve(segments.size());
            for (int id : ids){
                segs.push_back(std::move(segments[id]));
            }
            return segs;
        }



        Chain3 MakeChain(const std::vector<SectionalPiece> & pieces, bool closed){
            Chain3 chain;
            chain.closed = closed;
            chain.points.reserve(pieces.size() * 2);
            for (auto & seg : pieces){
                chain.points.push_back(seg.range.first);
                chain.points.push_back(seg.range.second);
            }

            if (!closed){

            }

            return chain;
        }






        std::pair<double, double> EstimateEffectiveRangeAlongDirection(
            const HandledTable<RegionHandle, std::vector<Polygon3>> & polygons,
            const Vec3 & direction, double stepLen, double minEffectiveAreaRatio,
            double gamma1, double gamma2){

            double minv = std::numeric_limits<double>::max();
            double maxv = std::numeric_limits<double>::lowest();

            auto ndir = normalize(direction);
            for (auto it = polygons.begin(); it != polygons.end(); ++it){
                for (auto & ply : *it){
                    for (auto & p : ply.corners){
                        double pos = p.dot(ndir);
                        if (minv > pos) minv = pos;
                        if (maxv < pos) maxv = pos;
                    }
                }
            }

            assert(maxv - minv > 0);

            std::vector<double> areas;
            std::vector<double> dareas;
            for (double x = minv; x <= maxv; x += stepLen){
                Plane3 cutplane(ndir * x, ndir);
                auto chain = MakeChain(MakeSectionalPieces(polygons, cutplane));
                double area = chain.size() < 3 ? 0.0 : Area(chain);
                assert(area >= 0.0);
                dareas.push_back(area - (areas.empty() ? 0.0 : areas.back()));
                areas.push_back(area);
            }

            std::pair<double, double> range(0, 0);
            if (areas.size() < 2){
                return range;
            }
            dareas[0] = dareas[1];

            size_t maxCutAreaId = std::max_element(areas.begin(), areas.end()) - areas.begin();
            double maxCutArea = areas[maxCutAreaId];

            int start = 0;
            while (start < areas.size() && areas[start] < minEffectiveAreaRatio * maxCutArea) ++start;
            int end = areas.size() - 1;
            while (end >= 0 && areas[end] < minEffectiveAreaRatio * maxCutArea) --end;

            std::cout << "all: " << areas.size() << std::endl;
            std::cout << "start: " << start << std::endl;
            std::cout << "end: " << end << std::endl;

            if (start >= end)
                return range;
            if (start > maxCutAreaId || end < maxCutAreaId)
                return range;

            range.first = start * stepLen + minv;
            range.second = end * stepLen + minv;
            //for (int i = start; i < maxCutAreaId; i++){
            //    if (i != 0 && (dareas[i] - dareas[i - 1]) / maxCutArea * ((maxv - minv) / stepLen) < -gamma1){
            //        range.first = i * stepLen + minv;
            //        break;
            //    }
            //}

            //for (int i = end; i > maxCutAreaId; i--){
            //    if (i != areas.size() - 1 && (dareas[i] - dareas[i + 1]) / maxCutArea * ((maxv - minv) / stepLen) < -gamma2){
            //        range.second = i * stepLen + minv;
            //        break;
            //    }
            //}

            return range;

        }


        void SmoothInstances(const RLGraph & mg, const InstanceTable<RegionData> & planes, const InstanceTable<LineData> & lines){



        }


    }


}