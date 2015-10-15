#include "../core/algorithms.hpp"

#include "rl_graph_control.hpp"
#include "pi_graph_vis.hpp"

namespace pano {
    namespace experimental {


        double DepthOfVertexAt(const PIConstraintGraph & cg, const PIGraph & mg, int ent, 
            const Vec3 & direction, const Point3 & eye = Origin()) {
            auto & v = cg.entities[ent];
            auto & plane = v.supportingPlane.reconstructed;
            return DepthAt(direction, plane, eye);
        }


        void VisualizeReconstruction(const PIConstraintGraph & cg, const PIGraph & mg,
            const std::function<gui::Color(int vert)> & vertColor,
            const std::function<void(int vert)> & vertClick) {

            gui::ResourceStore::set("texture", mg.view.image);

            gui::SceneBuilder viz;
            viz.installingOptions().discretizeOptions.colorTable = gui::ColorTableDescriptor::RGB;
            std::vector<core::Decorated<gui::Colored<gui::SpatialProjectedPolygon>, int>> spps;
            std::vector<core::Decorated<gui::Colored<core::Line3>, int>> lines;

            for (int vert : cg.determinableEnts) {
                auto & v = cg.entities[vert];
                if (v.isSeg()) {
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
                    spp.plane = v.supportingPlane.reconstructed;
                    if (HasValue(spp.plane, IsInfOrNaN<double>)) {
                        WARNNING("inf plane");
                        continue;
                    }
                    spps.push_back(core::DecorateAs(std::move(gui::ColorAs(spp, vertColor(vert))), vert));
                } else {
                    int line = v.id;
                    auto & plane = v.supportingPlane.reconstructed;
                    auto & projectedLine = mg.lines[line].component;
                    Line3 reconstructedLine(IntersectionOfLineAndPlane(Ray3(Origin(), projectedLine.first), plane).position,
                        IntersectionOfLineAndPlane(Ray3(Origin(), projectedLine.second), plane).position);
                    if (HasValue(reconstructedLine, IsInfOrNaN<double>)) {
                        WARNNING("inf line");
                        continue;
                    }
                    lines.push_back(core::DecorateAs(gui::ColorAs(reconstructedLine, vertColor(vert)), vert));
                }
            }

            auto ent2string = [&mg, &cg](int ent) -> std::string {
                auto & e = cg.entities[ent];
                std::stringstream ss;
                if (e.isSeg()) {
                    ss << "seg " << e.id << " dof: " << mg.seg2control[e.id].dof();
                } else {
                    ss << "line " << e.id << " claz: " << mg.lines[e.id].claz;
                }
                return ss.str();
            };
            viz.begin(spps, [&mg, &cg, &vertClick, ent2string](gui::InteractionID iid,
                const core::Decorated<gui::Colored<gui::SpatialProjectedPolygon>, int> & spp) {
                int ent = spp.decoration;
                std::cout << ent2string(ent) << std::endl;
                for (int cons : cg.ent2cons[ent]) {
                    if (Contains(cg.consBetweenDeterminableEnts, cons)) {
                        if (cg.constraints[cons].isConnection()) {
                            std::cout << "connect with " 
                                << ent2string(cg.constraints[cons].ent1 == ent ? cg.constraints[cons].ent2 : cg.constraints[cons].ent1) 
                                << " using " 
                                << cg.constraints[cons].anchors.size() 
                                << " anchors, the weight is " << cg.constraints[cons].weight << std::endl;
                        } else {
                            std::cout << "coplanar with "
                                << ent2string(cg.constraints[cons].ent1 == ent ? cg.constraints[cons].ent2 : cg.constraints[cons].ent1)
                                << ", the weight is " << cg.constraints[cons].weight << std::endl;
                        }
                    }
                }
                vertClick(ent);
            }).shaderSource(gui::OpenGLShaderSourceDescriptor::XPanorama).resource("texture").end();
            viz.installingOptions().discretizeOptions.color = gui::ColorTag::Magenta;
            viz.installingOptions().lineWidth = 4.0;
            viz.add(lines, [&mg, &cg, &vertClick, ent2string](gui::InteractionID iid,
                const core::Decorated<gui::Colored<Line3>, int> & line) {
                int ent = line.decoration;
                std::cout << ent2string(ent) << std::endl;
                for (int cons : cg.ent2cons[ent]) {
                    if (Contains(cg.consBetweenDeterminableEnts, cons)) {
                        if (cg.constraints[cons].isConnection()) {
                            std::cout << "connect with "
                                << ent2string(cg.constraints[cons].ent1 == ent ? cg.constraints[cons].ent2 : cg.constraints[cons].ent1)
                                << " using "
                                << cg.constraints[cons].anchors.size()
                                << " anchors, the weight is " << cg.constraints[cons].weight << std::endl;
                        } else {
                            std::cout << "coplanar with "
                                << ent2string(cg.constraints[cons].ent1 == ent ? cg.constraints[cons].ent2 : cg.constraints[cons].ent1)
                                << ", the weight is " << cg.constraints[cons].weight << std::endl;
                        }
                    }
                }
                vertClick(line.decoration);
            });

            std::vector<core::Decorated<gui::Colored<core::Line3>, std::pair<int, int>>> connectionLines;
            std::vector<core::Point3> connectionLineEnds;
            for (auto & e : cg.constraints) {
                if (!Contains(cg.determinableEnts, e.ent1) || !Contains(cg.determinableEnts, e.ent2))
                    continue;
                int snum = cg.entities[e.ent1].isSeg() + cg.entities[e.ent2].isSeg();
                auto & samples = e.anchors;
                for (auto & ss : samples) {
                    double d1 = DepthOfVertexAt(cg, mg, e.ent1, ss);
                    double d2 = DepthOfVertexAt(cg, mg, e.ent2, ss);
                    Line3 line(normalize(ss) * d1, normalize(ss) * d2);
                    connectionLines.push_back(DecorateAs(gui::ColorAs(line, snum == 0 ? gui::Black : snum == 1 ? gui::Blue : gui::Yellow), std::make_pair(e.ent1, e.ent2)));
                    connectionLineEnds.push_back(connectionLines.back().component.component.first);
                    connectionLineEnds.push_back(connectionLines.back().component.component.second);
                }
            }

            viz.installingOptions().discretizeOptions.color = gui::ColorTag::Black;
            viz.installingOptions().lineWidth = 1.0;
            viz.begin(connectionLines, [&mg, &cg](gui::InteractionID iid,
                const core::Decorated<gui::Colored<Line3>, std::pair<int, int>> & line) {
                std::cout << "this is an anchor of edge connecting ";
                auto & v1 = cg.entities[line.decoration.first];
                auto & v2 = cg.entities[line.decoration.second];
                if (v1.isSeg()) {
                    std::cout << "seg " << v1.id << " dof: " << mg.seg2control[v1.id].dof() << std::endl;
                } else {
                    std::cout << "line " << v1.id << " claz: " << mg.lines[v1.id].claz << std::endl;
                }
                std::cout << " and ";
                if (v2.isSeg()) {
                    std::cout << "seg " << v2.id << " dof: " << mg.seg2control[v2.id].dof() << std::endl;
                } else {
                    std::cout << "line " << v2.id << " claz: " << mg.lines[v2.id].claz << std::endl;
                }
                std::cout << std::endl;
            }).shaderSource(gui::OpenGLShaderSourceDescriptor::XLines).end();
            viz.installingOptions().pointSize = 2.0;
            viz.begin(connectionLineEnds).shaderSource(gui::OpenGLShaderSourceDescriptor::XPoints).end();

            viz.show(true, false, gui::RenderOptions().cullFrontFace(false).cullBackFace(true).bwColor(0.1).bwTexColor(0.9)
                .camera(PerspectiveCamera(500, 500, Point2(250, 250), 300, Point3(1, 1, 1), Point3(0, 0, 0))));

        }




        void VisualizeLayoutAnnotation(const PILayoutAnnotation & anno) {
            
            gui::SceneBuilder viz;
            viz.installingOptions().discretizeOptions.colorTable = gui::ColorTableDescriptor::RGB;
            std::vector<core::Decorated<gui::Colored<gui::SpatialProjectedPolygon>, int>> spps;
            //std::vector<core::Decorated<gui::Colored<Polygon3>, int>> pps;

            gui::ResourceStore::set("texture", anno.rectifiedImage);
            std::vector<Point3> cornerPositions(anno.ncorners());
            std::vector<int> cornerDegrees(anno.ncorners(), 0);
            for (int face = 0; face < anno.nfaces(); face++) {
                for (int c : anno.face2corners[face]) {
                    auto & plane = anno.face2plane[face];
                    auto pos = IntersectionOfLineAndPlane(Ray3(Origin(), anno.corners[c]), plane).position;
                    cornerPositions[c] += pos;
                    cornerDegrees[c] ++;
                }
            }
            for (int c = 0; c < anno.ncorners(); c++) {
                cornerPositions[c] /= cornerDegrees[c];
            }
            for (int face = 0; face < anno.nfaces(); face++) {

                gui::SpatialProjectedPolygon spp;
                spp.projectionCenter = core::Point3(0, 0, 0);
                spp.plane = anno.face2plane[face];                
                for (int c : anno.face2corners[face]) {
                    spp.corners.push_back(normalize(anno.corners[c]));
                }

                static const gui::ColorTable rgbTable = gui::RGB;
                auto & control = anno.face2control[face];

                spps.push_back(core::DecorateAs(std::move(gui::ColorAs(spp, 
                    control.dof() == 1 
                    ? rgbTable[control.orientationClaz] 
                    : (control.dof() == 2 
                        ? rgbTable[control.orientationNotClaz] 
                        : gui::Color(gui::White)))), face));

            /*    Polygon3 pp;
                pp.normal = anno.face2plane[face].normal;
                for (int c : anno.face2corners[face]) {
                    pp.corners.push_back(cornerPositions[c]);
                }

                pps.push_back(core::DecorateAs(std::move(gui::ColorAs(pp,
                    control.dof() == 1
                    ? rgbTable[control.orientationClaz]
                    : (control.dof() == 2
                    ? rgbTable[control.orientationNotClaz]
                    : gui::Color(gui::White)))), face));*/
            }

            viz.begin(spps/*, sppCallbackFun*/).shaderSource(gui::OpenGLShaderSourceDescriptor::XPanorama).resource("texture").end();
            viz.show(true, false, gui::RenderOptions().cullFrontFace(true).cullBackFace(false).bwColor(0.1).bwTexColor(0.9)
                .camera(PerspectiveCamera(500, 500, Point2(250, 250), 300, Point3(1, 1, 1), Point3(0, 0, 0))));

        }


    }
}