#include "../core/algorithms.hpp"

#include "rl_graph_control.hpp"
#include "pi_graph_vis.hpp"

namespace pano {
    namespace experimental {


        double DepthOfVertexAt(const PIGraph & mg, int vert, const Vec3 & direction, const Point3 & eye = Origin()) {
            auto & v = mg.verts[vert];
            if (v.isSeg) {
                int seg = v.id;
                auto & plane = mg.seg2recPlanes[seg];
                return DepthAt(direction, plane, eye);
            } else {
                int line = v.id;
                auto & l = mg.line2recLines[line];
                return DepthAt(direction, l, eye);
            }
        }


        void VisualizeReconstruction(const std::vector<int> & ccids, const PIGraph & mg,
            const std::function<gui::Color(int vert)> & vertColor,
            const std::function<void(int vert)> & vertClick) {

            auto sppCallbackFun = [&mg, &vertClick](gui::InteractionID iid,
                const core::Decorated<gui::Colored<gui::SpatialProjectedPolygon>, int> & spp) {
                auto & v = mg.verts[spp.decoration];
                std::cout << (v.isSeg ? "Seg" : "Line") << v.id << " is clicked" << std::endl;
                vertClick(spp.decoration);
            };

            gui::ResourceStore::set("texture", mg.view.image);
            for (int ccid : ccids) {

                gui::SceneBuilder viz;
                viz.installingOptions().discretizeOptions.colorTable = gui::ColorTableDescriptor::RGB;
                std::vector<core::Decorated<gui::Colored<gui::SpatialProjectedPolygon>, int>> spps;
                std::vector<gui::Colored<core::Line3>> lines;

                auto & verts = mg.cc2verts.at(ccid);

                for (int vert : verts) {
                    auto & v = mg.verts[vert];
                    if (v.isSeg) {
                        int seg = v.id;
                        if (!mg.seg2control[seg].used)
                            continue;
                        gui::SpatialProjectedPolygon spp;
                        auto & contours = mg.seg2contours[seg];
                        if (contours.empty() || contours.front().empty()) {
                            continue;
                        }
                        // filter corners
                        core::ForeachCompatibleWithLastElement(contours.front().begin(), contours.front().end(),
                            std::back_inserter(spp.corners),
                            [](const core::Vec3 & a, const core::Vec3 & b) -> bool {
                            return core::AngleBetweenDirections(a, b) > 0.0;
                        });
                        if (spp.corners.size() < 3)
                            continue;

                        spp.projectionCenter = core::Point3(0, 0, 0);
                        spp.plane = mg.seg2recPlanes[seg];
                        if (HasValue(spp.plane, IsInfOrNaN<double>)) {
                            WARNNING("inf plane");
                            continue;
                        }
                        spps.push_back(core::DecorateAs(std::move(gui::ColorAs(spp, vertColor(vert))), vert));
                    } else {
                        int line = v.id;
                        auto & l = mg.line2recLines[line];
                        if (HasValue(l, IsInfOrNaN<double>)) {
                            WARNNING("inf line");
                            continue;
                        }
                        lines.push_back(gui::ColorAs(l, vertColor(vert)));
                    }
                }

                viz.begin(spps/*, sppCallbackFun*/).shaderSource(gui::OpenGLShaderSourceDescriptor::XPanorama).resource("texture").end();
                viz.installingOptions().discretizeOptions.color = gui::ColorTag::Magenta;
                viz.installingOptions().lineWidth = 3.0;
                viz.add(lines);

                std::vector<core::Decorated<gui::Colored<core::Line3>, std::pair<int, int>>> connectionLines;
                std::vector<core::Point3> connectionLineEnds;
                for (auto & e : mg.edges) {
                    if (mg.vert2cc[e.vert1] != ccid || mg.vert2cc[e.vert2] != ccid)
                        continue;
                    int snum = mg.verts[e.vert1].isSeg + mg.verts[e.vert2].isSeg;
                    auto & samples = e.anchors;
                    for (auto & ss : samples) {
                        double d1 = DepthOfVertexAt(mg, e.vert1, ss);
                        double d2 = DepthOfVertexAt(mg, e.vert2, ss);
                        Line3 line(normalize(ss) * d1, normalize(ss) * d2);
                        connectionLines.push_back(DecorateAs(gui::ColorAs(line, snum == 0 ? gui::Black : snum == 1 ? gui::Blue : gui::Yellow), std::make_pair(e.vert1, e.vert2)));
                        connectionLineEnds.push_back(connectionLines.back().component.component.first);
                        connectionLineEnds.push_back(connectionLines.back().component.component.second);
                    }
                }

                viz.installingOptions().discretizeOptions.color = gui::ColorTag::Black;
                viz.installingOptions().lineWidth = 5.0;
                viz.begin(connectionLines, [&mg](gui::InteractionID iid,
                    const core::Decorated<gui::Colored<Line3>, std::pair<int, int>> & line) {
                    std::cout << "this is an anchor of edge connecting ";
                    auto & v1 = mg.verts[line.decoration.first];
                    auto & v2 = mg.verts[line.decoration.second];
                    if (v1.isSeg) {
                        std::cout << "seg " << v1.id;
                    } else {
                        std::cout << "line " << v1.id;
                    }
                    std::cout << " and ";
                    if (v2.isSeg) {
                        std::cout << "seg " << v2.id;
                    } else {
                        std::cout << "line " << v2.id;
                    }
                    std::cout << std::endl;
                }).shaderSource(gui::OpenGLShaderSourceDescriptor::XLines).end();
                viz.installingOptions().pointSize = 8.0;
                viz.begin(connectionLineEnds).shaderSource(gui::OpenGLShaderSourceDescriptor::XPoints).end();

                viz.show(true, false, gui::RenderOptions().cullFrontFace(false).cullBackFace(true).bwColor(0.7).bwTexColor(0.3)
                    .camera(PerspectiveCamera(500, 500, Point2(250, 250), 300, Point3(1, 1, 1), Point3(0, 0, 0))));
            }

        }



        void VisualizeLayoutAnnotation(const PILayoutAnnotation & anno) {
            // todo
        }

    }
}