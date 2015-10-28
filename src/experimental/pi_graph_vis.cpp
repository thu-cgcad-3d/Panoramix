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


        void VisualizeReconstruction(const PICGDeterminablePart & dp, const PIConstraintGraph & cg, const PIGraph & mg, bool showConnectionLines,
            const std::function<gui::Color(int vert)> & vertColor,
            const std::function<void(int vert)> & vertClick) {

            gui::ResourceStore::set("texture", mg.view.image);

            gui::SceneBuilder viz;
            viz.installingOptions().discretizeOptions.colorTable = gui::ColorTableDescriptor::RGB;
            std::vector<core::Decorated<gui::Colored<gui::SpatialProjectedPolygon>, int>> spps;
            std::vector<core::Decorated<gui::Colored<core::Line3>, int>> lines;

            for (int vert : dp.determinableEnts) {
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
            viz.begin(spps, [&mg, &cg, &dp, &vertClick, ent2string](gui::InteractionID iid,
                const core::Decorated<gui::Colored<gui::SpatialProjectedPolygon>, int> & spp) {
                int ent = spp.decoration;
                std::cout << ent2string(ent) << std::endl;
                for (int cons : cg.ent2cons[ent]) {
                    if (Contains(dp.consBetweenDeterminableEnts, cons)) {
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
            viz.installingOptions().discretizeOptions.color = gui::ColorTag::Black;
            viz.installingOptions().lineWidth = 4.0;
            viz.begin(lines/*, [&mg, &cg, &dp, &vertClick, ent2string](gui::InteractionID iid,
                const core::Decorated<gui::Colored<Line3>, int> & line) {
                int ent = line.decoration;
                std::cout << ent2string(ent) << std::endl;
                for (int cons : cg.ent2cons[ent]) {
                    if (Contains(dp.consBetweenDeterminableEnts, cons)) {
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
            }*/).shaderSource(gui::OpenGLShaderSourceDescriptor::XLines).end();


            if (showConnectionLines) {
                std::vector<core::Decorated<gui::Colored<core::Line3>, std::pair<int, int>>> connectionLines;
                std::vector<core::Point3> connectionLineEnds;
                for (int cons : dp.consBetweenDeterminableEnts) {
                    auto & e = cg.constraints[cons];
                    if (!Contains(dp.determinableEnts, e.ent1) || !Contains(dp.determinableEnts, e.ent2))
                        continue;
                    int snum = cg.entities[e.ent1].isSeg() + cg.entities[e.ent2].isSeg();
                    auto & samples = e.anchors;
                    for (auto & ss : samples) {
                        double d1 = DepthOfVertexAt(cg, mg, e.ent1, ss);
                        double d2 = DepthOfVertexAt(cg, mg, e.ent2, ss);
                        Line3 line(normalize(ss) * d1, normalize(ss) * d2);
                        connectionLines.push_back(DecorateAs(gui::ColorAs(line, snum == 0 ?
                            gui::Black : snum == 1 ? gui::Blue : gui::Yellow), std::make_pair(e.ent1, e.ent2)));
                        connectionLineEnds.push_back(connectionLines.back().component.component.first);
                        connectionLineEnds.push_back(connectionLines.back().component.component.second);
                    }
                }

                viz.installingOptions().discretizeOptions.color = gui::ColorTag::Black;
                viz.installingOptions().lineWidth = 3.0;
                viz.begin(connectionLines, [&mg, &cg, ent2string](gui::InteractionID iid,
                    const core::Decorated<gui::Colored<Line3>, std::pair<int, int>> & line) {
                    std::cout << "this is an anchor of edge connecting ";
                    auto & v1 = cg.entities[line.decoration.first];
                    auto & v2 = cg.entities[line.decoration.second];
                    std::cout << "connecting " << ent2string(line.decoration.first) << " and " << ent2string(line.decoration.second) << std::endl;
                }).shaderSource(gui::OpenGLShaderSourceDescriptor::XLines).end();
                viz.installingOptions().pointSize = 5.0;
                viz.begin(connectionLineEnds).shaderSource(gui::OpenGLShaderSourceDescriptor::XPoints).end();
            }

            viz.show(true, false, gui::RenderOptions().cullFrontFace(false).cullBackFace(true).bwColor(0.1).bwTexColor(0.9)
                .camera(PerspectiveCamera(500, 500, Point2(250, 250), 300, Point3(1, 1, 1), Point3(0, 0, 0))));

        }





        Imaged DepthMap(const PICGDeterminablePart & dp, const PIConstraintGraph & cg, const PIGraph & mg, 
            std::pair<double, double> * validDepthRange) {
            Imaged depths(mg.view.image.size(), 0.0);
            double minv = std::numeric_limits<double>::max();
            double maxv = 0.0;
            for (auto it = depths.begin(); it != depths.end(); ++it) {
                auto pos = it.pos();
                int seg = mg.segs(pos);
                int ent = cg.seg2ent[seg];
                if (ent == -1) {
                    continue;
                }
                if (!Contains(dp.determinableEnts, ent)) {
                    continue;
                }
                auto & plane = cg.entities[ent].supportingPlane.reconstructed;
                Vec3 dir = normalize(mg.view.camera.toSpace(pos));
                double depth = norm(IntersectionOfLineAndPlane(Ray3(Origin(), dir), plane).position);
                if (depth < minv) {
                    minv = depth;
                }
                if (depth > maxv) {
                    maxv = depth;
                }
                *it = depth;
            }
            if (validDepthRange) {
                validDepthRange->first = minv;
                validDepthRange->second = maxv;
            }
            return depths;
        }


        Image3d SurfaceNormalMap(const PICGDeterminablePart & dp, const PIConstraintGraph & cg, const PIGraph & mg) {
            Image3d snm(mg.view.image.size());
            for (auto it = snm.begin(); it != snm.end(); ++it) {
                auto pos = it.pos();
                int seg = mg.segs(pos);
                int ent = cg.seg2ent[seg];
                if (ent == -1) {
                    continue;
                }
                if (!Contains(dp.determinableEnts, ent)) {
                    continue;
                }
                auto & plane = cg.entities[ent].supportingPlane.reconstructed;                
                auto n = normalize(plane.normal);
                *it = n;
            }
            return snm;
        }




        std::vector<Polygon3> CompactModel(const PICGDeterminablePart & dp, const PIConstraintGraph & cg, const PIGraph & mg, double distThres) {

            struct Corner {
                int junction; //
                int bndpiece; // the head of the bndpiece if junction == -1 
                Vec3 dir;
                bool isJunction() const { return junction != -1; }
                bool isHeadOfBndPiece() const { return junction == -1; }
                bool operator < (const Corner & c) const { return std::tie(junction, bndpiece) < std::tie(c.junction, c.bndpiece); }
                bool operator == (const Corner & c) const { return std::tie(junction, bndpiece) == std::tie(c.junction, c.bndpiece); }
            };
            
            std::vector<Corner> corners;
            corners.reserve(mg.njuncs() + mg.nbndPieces());

            // add junction corners
            std::vector<int> junc2corner(mg.njuncs(), -1);
            for (int junc = 0; junc < mg.njuncs(); junc++) {
                Corner c;
                c.junction = junc;
                c.bndpiece = -1;
                c.dir = normalize(mg.junc2positions[junc]);
                corners.push_back(c);
                junc2corner[junc] = corners.size() - 1;
            }

            // add bndpiece head corners
            std::vector<int> bndPiece2corner(mg.nbndPieces(), -1);
            for (int bnd = 0; bnd < mg.nbnds(); bnd++) {
                auto & bps = mg.bnd2bndPieces[bnd];
                assert(!bps.empty());
                for (int k = 1; k < bps.size(); k++) {
                    int bp = bps[k];
                    Corner c;
                    c.junction = -1;
                    c.bndpiece = bp;
                    c.dir = normalize(mg.bndPiece2dirs[bp].front());
                    corners.push_back(c);
                    bndPiece2corner[bp] = corners.size() - 1;
                }
            }

            const int ncorners = corners.size();

            // collect corner2segs
            std::vector<std::vector<int>> corner2segs(ncorners);
            for (int corner = 0; corner < ncorners; corner++) {
                auto & c = corners[corner];
                if (c.isJunction()) {
                    std::set<int> relatedSegs;
                    for (int bnd : mg.junc2bnds[c.junction]) {
                        relatedSegs.insert(mg.bnd2segs[bnd].first);
                        relatedSegs.insert(mg.bnd2segs[bnd].second);
                    }
                    corner2segs[corner] = std::vector<int>(relatedSegs.begin(), relatedSegs.end());
                } else {
                    int bnd = mg.bndPiece2bnd[c.bndpiece];
                    corner2segs[corner] = {
                        mg.bnd2segs[bnd].first,
                        mg.bnd2segs[bnd].second
                    };
                }
            }

            // collect seg2corners
            std::vector<std::vector<int>> seg2corners(mg.nsegs);
            std::vector<std::vector<bool>> seg2bpCornersReversed(mg.nsegs);
            for (int seg = 0; seg < mg.nsegs; seg++) {
                std::set<int> bndsVisited;

                for (int firstBnd : mg.seg2bnds[seg]) {
                    if (Contains(bndsVisited, firstBnd)) {
                        continue;
                    }
                    
                    bool firstBndReversed = mg.bnd2segs[firstBnd].second != seg;

                    std::vector<int> orderedBnds = { firstBnd };
                    std::vector<bool> bndReversed = { firstBndReversed };

                    int fromJunc = mg.bnd2juncs[firstBnd].first;
                    int toJunc = mg.bnd2juncs[firstBnd].second;
                    if (firstBndReversed) {
                        std::swap(fromJunc, toJunc);
                    }

                    std::vector<int> orderedJuncs = { fromJunc, toJunc };
                    // the order should be [junc0, bnd0, junc1, bnd1, ..., ] -> [junc0 ,...

                    // lastBndRightHand

                    while (orderedJuncs.back() != orderedJuncs.front()) {
                        auto & bndCands = mg.junc2bnds[orderedJuncs.back()];
                        //int bndSelected = -1;
                        //int nextJuncSelected = -1;
                        bool foundBnd = false;
                        for (int bnd : bndCands) {
                            if (mg.bnd2segs[bnd].first != seg && mg.bnd2segs[bnd].second != seg) {
                                continue;
                            }
                            if (Contains(orderedBnds, bnd)) {
                                continue;
                            }
                            orderedBnds.push_back(bnd);
                            bool reversed = mg.bnd2segs[bnd].second != seg;
                            bndReversed.push_back(reversed);
                            int nextJunc = reversed ? mg.bnd2juncs[bnd].first : mg.bnd2juncs[bnd].second;
                            orderedJuncs.push_back(nextJunc);
                            foundBnd = true;
                            break;
                        }
                        assert(foundBnd);
                        if (!foundBnd) {
                            break;
                        }
                    }

                    std::vector<int> cs;
                    std::vector<bool> creversed;
                    for (int i = 0; i < orderedBnds.size(); i++) {
                        int junc = orderedJuncs[i];
                        int bnd = orderedBnds[i];
                        bool rev = bndReversed[i];

                        assert(junc2corner[junc] != -1);
                        cs.push_back(junc2corner[junc]);
                        creversed.push_back(false);
                        auto & bps = mg.bnd2bndPieces[bnd];
                        if (!rev) {
                            for (int k = 1; k < bps.size(); k++) {
                                int bp = bps[k];
                                assert(bndPiece2corner[bp] != -1);
                                cs.push_back(bndPiece2corner[bp]);
                                creversed.push_back(false);
                            }
                        } else {
                            for (int k = bps.size() - 1; k > 0; k--) {
                                int bp = bps[k];
                                assert(bndPiece2corner[bp] != -1);
                                cs.push_back(bndPiece2corner[bp]);
                                creversed.push_back(true);
                            }
                        }
                    }

                    if (cs.size() > seg2corners[seg].size()) {
                        for (int i = 0; i < orderedBnds.size(); i++) {
                            int bnd = orderedBnds[i];
                            bndsVisited.insert(bnd);
                        }
                        seg2corners[seg] = std::move(cs);
                        seg2bpCornersReversed[seg] = std::move(creversed);
                    }
                }
            }


            std::vector<std::map<int, double>> corner2segDepths(ncorners);
            for (int corner = 0; corner < ncorners; corner++) {
                auto & segs = corner2segs[corner];
                auto & segDepths = corner2segDepths[corner];
                for (int seg : segs) {
                    int ent = cg.seg2ent[seg];
                    if (ent == -1) {
                        continue;
                    }
                    if (!Contains(dp.determinableEnts, ent)) {
                        continue;
                    }
                    auto & plane = cg.entities[ent].supportingPlane.reconstructed;
                    double depth = norm(IntersectionOfLineAndPlane(Ray3(Origin(), corners[corner].dir), plane).position);
                    segDepths[seg] = depth;
                }

                if (segDepths.empty()) {
                    continue;
                }

                // cluster these depths
                std::vector<double> groupedDepths;
                std::map<int, int> seg2group;
                for (auto & segAndDepth : segDepths) {
                    int seg = segAndDepth.first;
                    double depth = segAndDepth.second;
                    double minDepthDiff = distThres;
                    for (int g = 0; g < groupedDepths.size(); g++) { // find the group with min depth diff
                        double depthDiff = abs(depth - groupedDepths[g]);
                        if (depthDiff < minDepthDiff) {
                            seg2group[seg] = g;
                            minDepthDiff = depthDiff;
                        }
                    }
                    if (!Contains(seg2group, seg)) { // if not group is found, add a new group
                        groupedDepths.push_back(depth);
                        seg2group[seg] = groupedDepths.size() - 1;
                    }
                }

                // update each seg depth using group depth
                for (auto & segAndDepth : segDepths) {
                    segAndDepth.second = groupedDepths[seg2group.at(segAndDepth.first)];
                }
            }



            // fix the surfaces
            std::vector<bool> seg2determined(mg.nsegs, false);
            for (int seg = 0; seg < mg.nsegs; seg++) {
                int ent = cg.seg2ent[seg];
                if (ent == -1) {
                    continue;
                }
                if (!Contains(dp.determinableEnts, ent)) {
                    continue;
                }
                seg2determined[seg] = true;
            }

            std::vector<Polygon3> seg2compact(mg.nsegs);
            for (int seg = 0; seg < mg.nsegs; seg++) {
                auto & compact = seg2compact[seg];
                if (seg2determined[seg]) {
                    auto & cs = seg2corners[seg];
                    compact.corners.resize(cs.size());
                    for (int k = 0; k < cs.size(); k++) {
                        int corner = cs[k];
                        auto dir = normalize(corners[corner].dir);
                        double depth = corner2segDepths[corner].at(seg);
                        compact.corners[k] = dir * depth;
                    }
                    compact.normal = cg.entities[cg.seg2ent[seg]].supportingPlane.reconstructed.normal;
                } else {
                    auto & cs = seg2corners[seg];                    
                    std::vector<std::vector<Point3>> candContourTable;
                    for (int k = 0; k < cs.size(); k++) {
                        int corner = cs[k];
                        auto dir = normalize(corners[corner].dir);
                        std::vector<Point3> candContour;
                        for (auto & dcand : corner2segDepths[corner]) {
                            double depth = dcand.second;
                            candContour.push_back(dir * depth);
                        }
                        if (!candContour.empty()) {
                            candContourTable.push_back(std::move(candContour));
                        }
                    }
                    if (candContourTable.size() >= 3) {
                        std::vector<int> candContourIds(candContourTable.size());
                        std::iota(candContourIds.begin(), candContourIds.end(), 0);
                        // make cand contour with fewer cands forward
                        std::sort(candContourIds.begin(), candContourIds.end(), [&candContourTable](int a, int b) {
                            return candContourTable[a].size() < candContourTable[b].size();
                        });
                        // iterate all candidate planes to find the best that fit the most
                        Plane3 bestPlane;
                        double bestScore = -1.0;
                        for (auto & p1 : candContourTable[candContourIds[0]]) {
                            for (auto & p2 : candContourTable[candContourIds[1]]) {
                                for (auto & p3 : candContourTable[candContourIds[2]]) {
                                    auto candPlane = Plane3From3Points(p1, p2, p3);
                                    if (candPlane.distanceTo(Origin()) < distThres) {
                                        continue;
                                    }
                                    double candPlaneScore = 0.0;
                                    // is this candPlane fit other contour cands?
                                    for (int i = 3; i < candContourIds.size(); i++) {
                                        for (auto & pRest : candContourTable[candContourIds[i]]) {
                                            if (DistanceFromPointToPlane(pRest, candPlane).first < distThres) {
                                                candPlaneScore += 1.0;
                                            }
                                        }
                                    }
                                    if (std::any_of(mg.vps.begin(), mg.vps.end(), [&candPlane](const Vec3 & vp) {
                                        return AngleBetweenUndirectedVectors(vp, candPlane.normal) < DegreesToRadians(5);
                                    })) {
                                        candPlaneScore + 0.5;
                                    }
                                    if (candPlaneScore > bestScore) {
                                        bestScore = candPlaneScore;
                                        bestPlane = candPlane;
                                    }
                                }
                            }
                        }
                        // find a best plane
                        if (bestPlane.normal != Origin()) {
                            compact.corners.resize(cs.size());
                            compact.normal = bestPlane.normal;
                            for (int k = 0; k < cs.size(); k++) {
                                int corner = cs[k];
                                auto dir = normalize(corners[corner].dir);
                                // find most fit depth in corner2segDepths
                                double depth = -1.0;
                                double minDist = std::numeric_limits<double>::infinity();
                                for (auto & depthCand : corner2segDepths[corner]) {
                                    Point3 p = dir * depthCand.second;
                                    double dist = bestPlane.distanceTo(p);
                                    if (dist < minDist) {
                                        depth = depthCand.second;
                                        minDist = dist;
                                    }
                                }
                                // corner2segDepths[corner] is empty, then compute depth directly
                                if (depth > 0) {
                                    compact.corners[k] = dir * depth;
                                } else {
                                    compact.corners[k] = IntersectionOfLineAndPlane(Ray3(Origin(), dir), bestPlane).position;
                                }
                            }
                            if (compact.normal.dot(compact.corners.front()) < 0) {
                                compact.normal = - compact.normal;
                            }
                        }
                    }
                }
            }
         
            return seg2compact;
        }






        void VisualizeReconstructionCompact(const PICGDeterminablePart & dp, const PIConstraintGraph & cg, const PIGraph & mg, bool doModel) {
            gui::ResourceStore::set("texture", mg.view.image);

            auto compactPolygons = CompactModel(dp, cg, mg, 0.1);
            auto e = std::remove_if(compactPolygons.begin(), compactPolygons.end(),
                [](const Polygon3 & poly) {return poly.normal == Origin() || poly.corners.size() <= 2; });

            compactPolygons.erase(e, compactPolygons.end());

            gui::SceneBuilder viz;
            viz.installingOptions().discretizeOptions.colorTable = gui::ColorTableDescriptor::RGB;
            viz.begin(compactPolygons/*, sppCallbackFun*/).shaderSource(gui::OpenGLShaderSourceDescriptor::XPanorama).resource("texture").end();


            std::vector<core::Decorated<gui::Colored<core::Line3>, int>> lines;
            for (int vert : dp.determinableEnts) {
                auto & v = cg.entities[vert];
                if (v.isLine()) {
                    int line = v.id;
                    auto & plane = v.supportingPlane.reconstructed;
                    auto & projectedLine = mg.lines[line].component;
                    Line3 reconstructedLine(IntersectionOfLineAndPlane(Ray3(Origin(), projectedLine.first), plane).position,
                        IntersectionOfLineAndPlane(Ray3(Origin(), projectedLine.second), plane).position);
                    if (HasValue(reconstructedLine, IsInfOrNaN<double>)) {
                        WARNNING("inf line");
                        continue;
                    }
                    lines.push_back(core::DecorateAs(gui::ColorAs(reconstructedLine, gui::Black), vert));
                }
            }
            viz.begin(lines).shaderSource(gui::OpenGLShaderSourceDescriptor::XLines).end();


            viz.show(doModel, false, gui::RenderOptions().cullFrontFace(true).cullBackFace(false).bwColor(0.1).bwTexColor(0.9)
                .camera(PerspectiveCamera(500, 500, Point2(250, 250), 300, Point3(1, 1, 1), Point3(0, 0, 0))));
        }




        void VisualizeLayoutAnnotation(const PILayoutAnnotation & anno) {

            if (anno.nfaces() == 0) {
                return;
            }
            
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