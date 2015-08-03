#include "../core/algorithms.hpp"
#include "../core/containers.hpp"
#include "../core/cameras.hpp"
#include "../core/utility.hpp"
#include "../core/clock.hpp"
#include "../core/homo_graph.hpp"
#include "../core/single_view.hpp"

#include "rlgraph.hpp"

namespace pano {
    namespace experimental {



        Polygon3 RData::asPolygon(const Plane3 & plane) const {
            auto & outerContour = normalizedContours.front();
            std::vector<Point3> projectedContur(outerContour.size());
            for (int i = 0; i < outerContour.size(); i++) {
                projectedContur[i] = IntersectionOfLineAndPlane(Ray3(Origin(), outerContour[i]), plane).position;
            }
            return Polygon3(std::move(projectedContur), plane.normal);
        }


        RLGraphBuilder::Params::Params() : 
            intersectionAngleThresholdForLL(0.04),
            incidenceParaAngleThresholdForLL(0.1),
            incidenceVertAngleThresholdForLL(0.02),
            incidenceParaAngleThresholdForLLAcrossViews(0.15),
            incidenceVertAngleThresholdForLLAcrossViews(0.03),
            samplingStepAngleOnRR(2 / 500.0), 
            samplingStepAngleOnRL(2 / 500.0), 
            samplingStepAngleOnOcc(1 / 500.0),
            samplingStepSizeOnRR(2), 
            samplingStepSizeOnRL(2),
            samplingStepSizeOnOcc(1){
        }


        namespace {

            // lines graph in 2d
            struct LineData2D {
                Line2 line;
                int claz;
                std::vector<double> vpScores;
            };
            struct LineRelationData2D {
                Point2 relationCenter;
                bool isIntersection;
                LLData::ManhattanJunctionType mjType;
            };
            using LinesGraph2D = HomogeneousGraph02<LineData2D, LineRelationData2D>;
            using LineHandle2D = HandleOfTypeAtLevel<LinesGraph2D, 0>;
            using LineRelationHandle2D = HandleOfTypeAtLevel<LinesGraph2D, 1>;


            LLData::ManhattanJunctionType GetManhattanJunctionType(const Mat<float, 3, 2> & v) {
                double junctionWeight = 0.0;
                // Y
                double Y = 0.0;
                for (int s = 0; s < 2; s++) {
                    Y += v(0, s) * v(1, s) * v(2, s) * DiracDelta(v(0, 1 - s) + v(1, 1 - s) + v(2, 1 - s), 1e-4);
                }

                // W
                double W = 0.0;
                for (int i = 0; i < 3; i++) {
                    for (int j = 0; j < 3; j++) {
                        if (i == j)
                            continue;
                        int k = 3 - i - j;
                        for (int s = 0; s < 2; s++) {
                            W += v(i, s) * v(j, 1 - s) * v(k, 1 - s) * DiracDelta(v(i, 1 - s) + v(j, s) + v(k, s), 1e-4);
                        }
                    }
                }

                // K
                double K = 0.0;
                for (int i = 0; i < 3; i++) {
                    for (int j = 0; j < 3; j++) {
                        if (i == j)
                            continue;
                        int k = 3 - i - j;
                        K += v(i, 0) * v(i, 1) * v(j, 0) * v(k, 1) * DiracDelta(v(j, 1) + v(k, 0), 1e-4);
                        K += v(i, 0) * v(i, 1) * v(j, 1) * v(k, 0) * DiracDelta(v(j, 0) + v(k, 1), 1e-4);
                    }
                }

                // compute X junction
                double X = 0.0;
                for (int i = 0; i < 3; i++) {
                    for (int j = 0; j < 3; j++) {
                        if (i == j)
                            continue;
                        int k = 3 - i - j;
                        X += v(i, 0) * v(i, 1) * v(j, 0) * v(j, 1) * DiracDelta(v(k, 0) + v(k, 1), 1e-4);
                    }
                }

                // compute T junction
                double T = 0.0;
                for (int i = 0; i < 3; i++) {
                    for (int j = 0; j < 3; j++) {
                        if (i == j)
                            continue;
                        int k = 3 - i - j;
                        T += v(i, 0) * v(i, 1) * v(j, 0) * DiracDelta(v(j, 1) + v(k, 0) + v(k, 1), 1e-4);
                        T += v(i, 0) * v(i, 1) * v(j, 1) * DiracDelta(v(j, 0) + v(k, 0) + v(k, 1), 1e-4);
                    }
                }

                // compute L junction
                double L = 0.0;
                for (int i = 0; i < 3; i++) {
                    for (int j = 0; j < 3; j++) {
                        if (i == j)
                            continue;
                        int k = 3 - i - j;
                        for (int a = 0; a < 2; a++) {
                            int nota = 1 - a;
                            for (int b = 0; b < 2; b++) {
                                int notb = 1 - b;
                                L += v(i, a) * v(j, b) * DiracDelta(v(i, nota) + v(j, notb) + v(k, 0) + v(k, 1), 1e-4);
                            }
                        }
                    }
                }

                //std::cout << " Y-" << Y << " W-" << W << " K-" << K << 
                //    " X-" << X << " T-" << T << " L-" << L << std::endl; 
                static const double threshold = 0;
                if (Y > threshold) {
                    return LLData::Y;
                } else if (W > threshold) {
                    return LLData::W;
                } else if (L > threshold) {
                    return LLData::L;
                } else if (K > threshold) {
                    return LLData::K;
                } else if (X > threshold) {
                    return LLData::X;
                } else if (T > threshold) {
                    return LLData::T;
                } else {
                    return LLData::Hex;
                }
            }



            LinesGraph2D CreateLinesGraph2D(const std::vector<Classified<Line2>> & lines,
                const std::vector<HPoint2> & vps,
                const DenseMatd & lineVPScores,
                double intersectionDistanceThreshold,
                double incidenceDistanceAlongDirectionThreshold,
                double incidenceDistanceVerticalDirectionThreshold,
                bool includeIncidencesBetweenUnclassifiedLines = false) {

                /// TODO include all line relations is NOT STABLE!

                LinesGraph2D graph;

                // insert lines
                std::vector<LineHandle2D> handles;
                handles.reserve(lines.size());
                graph.internalElements<0>().reserve(lines.size());

                for (int i = 0; i < lines.size(); i++) {
                    auto & line = lines[i];
                    LineData2D ld;
                    ld.line = line.component;
                    ld.claz = line.claz;
                    lineVPScores.row(i).copyTo(ld.vpScores);
                    assert(ld.claz == -1 || 
                        std::max_element(ld.vpScores.begin(), ld.vpScores.end()) - ld.vpScores.begin() == ld.claz);
                    handles.push_back(graph.add(std::move(ld)));
                }

                // construct incidence/intersection relations
                auto & linesData = graph.internalElements<0>();
                for (int i = 0; i < linesData.size(); i++) {
                    auto & linei = linesData[i].data.line;
                    int clazi = linesData[i].data.claz;

                    for (int j = i + 1; j < linesData.size(); j++) {
                        auto & linej = linesData[j].data.line;
                        int clazj = linesData[j].data.claz;

                        auto nearest = DistanceBetweenTwoLines(linei, linej);
                        double d = nearest.first;

                        if (clazi == clazj && clazi >= 0) { // incidences for classified lines
                            auto conCenter = (nearest.second.first.position + nearest.second.second.position) / 2.0;
                            auto conDir = (nearest.second.first.position - nearest.second.second.position);

                            auto & vp = vps[clazi];

                            if (Distance(vp.value(), conCenter) < intersectionDistanceThreshold)
                                continue;

                            auto dir = normalize((vp - HPoint2(conCenter)).numerator);
                            double dAlong = abs(conDir.dot(dir));
                            double dVert = sqrt(Square(norm(conDir)) - dAlong*dAlong);

                            if (dAlong < incidenceDistanceAlongDirectionThreshold &&
                                dVert < incidenceDistanceVerticalDirectionThreshold) { // incidence
                                LineRelationData2D lrd;
                                lrd.isIntersection = false;
                                lrd.relationCenter = conCenter;

                                if (HasValue(lrd.relationCenter, IsInfOrNaN<double>))
                                    continue;

                                graph.add<1>({ handles[i], handles[j] }, std::move(lrd));
                            }
                        } else if (clazi != clazj && clazi >= 0 && clazj >= 0) { // intersections for classified lines
                            if (d < intersectionDistanceThreshold) {
                                auto conCenter = HPointFromVector(GetCoeffs(linei.ray())
                                    .cross(GetCoeffs(linej.ray()))).value();

                                assert(conCenter != Point2(0.0, 0.0));

                                if (DistanceFromPointToLine(conCenter, linei).first > intersectionDistanceThreshold * 4 ||
                                    DistanceFromPointToLine(conCenter, linej).first > intersectionDistanceThreshold * 4)
                                    continue;

                                LineRelationData2D lrd;
                                lrd.isIntersection = true;
                                lrd.relationCenter = conCenter;

                                if (HasValue(lrd.relationCenter, IsInfOrNaN<double>))
                                    continue;

                                graph.add<1>({ handles[i], handles[j] }, std::move(lrd));
                            }
                        } else if (includeIncidencesBetweenUnclassifiedLines && clazi == clazj && clazi == -1) { // incidence for unclassified lines
                            if (d >= incidenceDistanceAlongDirectionThreshold / 3.0)
                                continue;

                            double angle = AngleBetweenUndirectedVectors(linei.direction(), linej.direction());
                            if (angle >= DegreesToRadians(2))
                                continue;

                            double angle2 = std::max(AngleBetweenUndirectedVectors(linei.direction(),
                                nearest.second.first.position - nearest.second.second.position),
                                AngleBetweenUndirectedVectors(linej.direction(),
                                nearest.second.first.position - nearest.second.second.position));
                            if (angle2 >= DegreesToRadians(3))
                                continue;

                            LineRelationData2D lrd;
                            lrd.isIntersection = false;
                            lrd.relationCenter = (nearest.second.first.position + nearest.second.second.position) / 2.0;

                            if (HasValue(lrd.relationCenter, IsInfOrNaN<double>))
                                continue;

                            graph.add<1>({ handles[i], handles[j] }, std::move(lrd));
                        }

                    }
                }

                {
                    static const double angleThreshold = M_PI / 32;
                    static const double sigma = 0.1;

                    enum LineVotingDirection : int {
                        TowardsVanishingPoint = 0,
                        TowardsOppositeOfVanishingPoint = 1
                    };

                    Box2 linesBBox = BoundingBoxOfContainer(lines);
                    linesBBox.expand(linesBBox.outerSphere().radius / 2.0);

                    // compute junction weights
                    for (auto & lr : graph.elements<1>()) {
                        auto & lrd = lr.data;
                        if (!lrd.isIntersection) {
                            lrd.mjType = LLData::I;
                        } else {
                            Mat<float, 3, 2> votingData;
                            std::fill(std::begin(votingData.val), std::end(votingData.val), 0);

                            for (auto & ld : graph.elements<0>()) {
                                auto & line = ld.data.line;
                                int claz = ld.data.claz;
                                if (claz == -1 || claz >= 3)
                                    continue;

                                auto & vp = vps[claz];
                                Point2 center = line.center();

                                Vec2 center2vp = vp.value() - center;
                                Vec2 center2pos = lrd.relationCenter - center;

                                if (norm(center2pos) <= 1)
                                    continue;

                                double angle = AngleBetweenDirections(center2pos, center2vp);
                                double angleSmall = angle > M_PI_2 ? (M_PI - angle) : angle;
                                if (IsInfOrNaN(angleSmall))
                                    continue;

                                assert(angleSmall >= 0 && angleSmall <= M_PI_2);

                                double angleScore =
                                    exp(-(angleSmall / angleThreshold) * (angleSmall / angleThreshold) / sigma / sigma / 2);

                                auto proj = ProjectionOfPointOnLine(lrd.relationCenter, line);
                                double projRatio = BoundBetween(proj.ratio, 0.0, 1.0);

                                if (AngleBetweenDirections(center2vp, line.direction()) < M_PI_2) { // first-second-vp
                                    votingData(claz, TowardsVanishingPoint) += angleScore * line.length() * (1 - projRatio);
                                    votingData(claz, TowardsOppositeOfVanishingPoint) += angleScore * line.length() * projRatio;
                                } else { // vp-first-second
                                    votingData(claz, TowardsOppositeOfVanishingPoint) += angleScore * line.length() * (1 - projRatio);
                                    votingData(claz, TowardsVanishingPoint) += angleScore * line.length() * projRatio;
                                }
                            }


                            auto p = ToPixel(lrd.relationCenter);
                            if (core::IsBetween(lrd.relationCenter[0], 0, linesBBox.size()[0] - 1) &&
                                core::IsBetween(lrd.relationCenter[1], 0, linesBBox.size()[1] - 1)) {
                                lrd.mjType = GetManhattanJunctionType(votingData);
                            } else {
                                lrd.mjType = LLData::Outsided;
                            }
                        }
                    }
                }
                return graph;

            }




            void AppendLines(RLGraph & mg, const std::vector<Classified<Line2>> & lineSegments, const DenseMatd & lineVPScore,
                const PerspectiveCamera & cam,
                const std::vector<Vec3> & vps,
                double intersectionAngleThreshold /*= 0.04*/,
                double incidenceAngleAlongDirectionThreshold /*= 0.1*/,
                double incidenceAngleVerticalDirectionThreshold /*= 0.02*/,
                double interViewIncidenceAngleAlongDirectionThreshold /*= 0.15*/, // for new line-line incidence recognition
                double interViewIncidenceAngleVerticalDirectionThreshold /*= 0.03*/) {

                if (lineSegments.empty())
                    return;

                assert(mg.internalComponents<RData>().empty() && "Regions must be added AFTER all lines!!!");

                std::vector<HPoint2> vps2d(vps.size());
                for (int i = 0; i < vps.size(); i++) {
                    vps2d[i] = cam.toScreenInHPoint(vps[i]);
                }
                LinesGraph2D graph2d = CreateLinesGraph2D(lineSegments, vps2d, lineVPScore,
                    cam.focal() * intersectionAngleThreshold,
                    cam.focal() * incidenceAngleAlongDirectionThreshold,
                    cam.focal() * incidenceAngleVerticalDirectionThreshold, false);

                mg.internalComponents<LData>().reserve(mg.internalComponents<LData>().size() +
                    graph2d.internalElements<0>().size());
                mg.internalConstraints<LLData>().reserve(mg.internalConstraints<LLData>().size() +
                    graph2d.internalElements<1>().size());

                assert(graph2d.isDense());

                std::unordered_set<LHandle> newLineHandles;
                std::unordered_map<LineHandle2D, LHandle> lh2dToLh;

                for (auto & l2d : graph2d.elements<0>()) {
                    auto & ld2d = l2d.data;
                    LData ld;
                    ld.initialClaz = ld2d.claz;
                    ld.normalizedLine.first = normalize(cam.toSpace(ld2d.line.first));
                    ld.normalizedLine.second = normalize(cam.toSpace(ld2d.line.second));
                    ld.vpScores = ld2d.vpScores;
                    lh2dToLh[l2d.topo.hd] = mg.addComponent(std::move(ld));
                    newLineHandles.insert(lh2dToLh[l2d.topo.hd]);
                }

                for (auto & r2d : graph2d.elements<1>()) {
                    auto & rd2d = r2d.data;
                    LLData rd;
                    rd.mjType = rd2d.mjType;
                    rd.normalizedRelationCenter = normalize(cam.toSpace(rd2d.relationCenter));
                    mg.addConstraint(std::move(rd), lh2dToLh.at(r2d.topo.lowers.front()), lh2dToLh.at(r2d.topo.lowers.back()));
                }

                // line incidence relations between old lines and new lines
                std::map<std::pair<LHandle, LHandle>, Vec3> interViewLineIncidences;

                // build rtree for lines
                auto lookupLineNormal = [&mg](const LHandle & li) -> Box3 {
                    auto normal = mg.data(li).normalizedLine.first.cross(mg.data(li).normalizedLine.second);
                    Box3 b = BoundingBox(normalize(normal));
                    static const double s = 0.2;
                    b.minCorner = b.minCorner - Vec3(s, s, s);
                    b.maxCorner = b.maxCorner + Vec3(s, s, s);
                    return b;
                };

                RTreeWrapper<LHandle, decltype(lookupLineNormal)> linesRTree(lookupLineNormal);
                for (auto & l : mg.components<LData>()) {
                    linesRTree.insert(l.topo.hd);
                }

                // recognize incidence constraints between lines of different views
                for (auto & l : mg.components<LData>()) {
                    auto lh = l.topo.hd;
                    if (Contains(newLineHandles, lh))
                        continue; // stick old lh, find new lh
                    auto & lineData = l.data;
                    linesRTree.search(lookupLineNormal(lh),
                        [interViewIncidenceAngleAlongDirectionThreshold, interViewIncidenceAngleVerticalDirectionThreshold,
                        &l, &mg, &interViewLineIncidences, &newLineHandles](const LHandle & relatedLh) -> bool {
                        if (!Contains(newLineHandles, relatedLh))
                            return true;

                        assert(l.topo.hd < relatedLh);

                        auto & line1 = l.data.normalizedLine;
                        auto & line2 = mg.data(relatedLh).normalizedLine;

                        if (l.data.initialClaz != mg.data(relatedLh).initialClaz) // only incidence relations are recognized here
                            return true;

                        auto normal1 = normalize(line1.first.cross(line1.second));
                        auto normal2 = normalize(line2.first.cross(line2.second));

                        if (std::min(std::abs(AngleBetweenDirections(normal1, normal2)),
                            std::abs(AngleBetweenDirections(normal1, -normal2))) <
                            interViewIncidenceAngleVerticalDirectionThreshold) {

                            auto nearest = DistanceBetweenTwoLines(normalize(line1), normalize(line2));
                            if (AngleBetweenDirections(nearest.second.first.position, nearest.second.second.position) >
                                interViewIncidenceAngleAlongDirectionThreshold) // ignore too far-away relations
                                return true;

                            auto relationCenter = (nearest.second.first.position + nearest.second.second.position) / 2.0;
                            relationCenter /= norm(relationCenter);

                            interViewLineIncidences[std::make_pair(l.topo.hd, relatedLh)] = relationCenter;
                        }
                        return true;
                    });
                }

                // add to mixed graph
                for (auto & ivli : interViewLineIncidences) {
                    LLData lrd;
                    lrd.mjType = LLData::IAcrossView;
                    lrd.normalizedRelationCenter = ivli.second;
                    mg.addConstraint(std::move(lrd), ivli.first.first, ivli.first.second);
                }

            }


        }


        namespace {

            template <class CameraT>
            std::vector<RHandle> AppendRegionsTemplate(
                RLGraph & mg,
                const Imagei & segmentedRegions, const Image5d & gc, const std::vector<Chain3> & occ,
                const CameraT & cam,
                double samplingStepAngleOnBoundary,
                double samplingStepAngleOnLine,
                double samplingStepAngleOnOcc,
                int samplerSizeOnBoundary,
                int samplerSizeOnLine, // the neighborhood radius for each sampler line point 
                int samplerSizeOnOcc
                ) 
            {

                std::vector<std::vector<std::vector<Vec3>>> ncontours;
                std::vector<Vec3> ncenters;
                std::vector<double> areas;
                // regions too small are removed and for them regIdOld2New[regId] = -1
                auto regIdOld2New = ComputeSpatialRegionProperties(segmentedRegions, cam, &ncontours, &ncenters, &areas);
                assert(AllSameInContainer({ ncontours.size(), ncenters.size(), areas.size() }));

                auto & rdatatable = mg.internalComponents<RData>();
                rdatatable.reserve(rdatatable.size() + regIdOld2New.size());
                auto gcView = MakeView(gc, cam);

                // add RData
                std::vector<RHandle> regId2Rh(regIdOld2New.size());
                for (int i = 0; i < regIdOld2New.size(); i++) {
                    int newRegId = regIdOld2New[i];
                    if (newRegId == -1)
                        continue;
                    RData rd;
                    rd.area = areas[newRegId];
                    rd.normalizedCenter = ncenters[newRegId];
                    rd.normalizedContours = std::move(ncontours[newRegId]);
                    rd.gcMean = MeanInMask(gcView, core::PerfectRegionMaskView(rd.normalizedContours, rd.normalizedCenter, 80.0));
                    auto rh = mg.addComponent(std::move(rd));
                    regId2Rh[i] = rh;
                }

                


                int regionNum = areas.size();

                // add RLs and RRLs
                std::map<std::tuple<RHandle, LHandle>, RLData> rlCons;
                std::map<std::tuple<RHandle, RHandle, LHandle>, RRLData> rrlCons;

                auto dc = AsDimensionConvertor(cam);

                for (auto & ld : mg.components<LData>()) {
                    auto & line = ld.data.normalizedLine;
                    Line2 line2 = dc.toScreen(line);
                    Vec2 vertToLineDir = normalize(PerpendicularDirection(line2.direction())); // vertical to this line
                    double stepOnImage = std::max(cam.focal() * samplingStepAngleOnLine, 0.5);

                    double spanAngle = AngleBetweenDirections(line.first, line.second);
                    int stepNum = static_cast<int>(std::ceil(spanAngle / samplingStepAngleOnLine));
                    if (stepNum == 0)
                        continue;

                    for (int step = 0; step <= stepNum; step++) {
                        double angle = step * samplingStepAngleOnLine;
                        Vec3 sample = RotateDirection(line.first, line.second, angle);
                        if (!cam.isVisibleOnScreen(sample))
                            continue;
                        Point2 sampleP = cam.toScreen(sample);
                        Pixel originalP = ToPixel(sampleP);
                        // collect neighbors
                        std::set<RHandle> connectedRhs;
                        for (int x = originalP.x - samplerSizeOnLine; x <= originalP.x + samplerSizeOnLine; x++) {
                            for (int y = originalP.y - samplerSizeOnLine; y <= originalP.y + samplerSizeOnLine; y++) {
                                Pixel p(x, y);
                                if (!Contains(segmentedRegions, p))
                                    continue;
                                auto rh = regId2Rh[segmentedRegions(p)];
                                if (rh.invalid())
                                    continue;
                                connectedRhs.insert(rh);
                            }
                        }
                        if (connectedRhs.size() == 1) {
                            // a RL here
                            RHandle rh = *connectedRhs.begin();
                            auto & rd = rlCons[std::make_tuple(rh, ld.topo.hd)];
                            rd.normalizedAnchors.push_back(normalize(sample));
                        } else if (connectedRhs.size() == 2) {
                            // a RRL here
                            auto it = connectedRhs.begin();
                            RHandle rh1 = *it;
                            ++it;
                            RHandle rh2 = *it;
                            if (rh2 < rh1)
                                std::swap(rh1, rh2);
                            auto & rrd = rrlCons[std::make_tuple(rh1, rh2, ld.topo.hd)];
                            rrd.normalizedAnchors.push_back(normalize(sample));
                        } else {
                            // the line lines on a boundary junction, we just ignore the sample point here :)
                        }
                    }
                }

                // add rl
                for (auto & rlc : rlCons) {
                    if (rlc.second.normalizedAnchors.size() <= 1)
                        continue;
                    rlc.second.length =
                        AngleBetweenDirections(rlc.second.normalizedAnchors.front(), rlc.second.normalizedAnchors.back());
                    mg.addConstraint(std::move(rlc.second), Get<0>(rlc.first), Get<1>(rlc.first));
                }
                // add rrl
                for (auto & rrlc : rrlCons) {
                    if (rrlc.second.normalizedAnchors.size() <= 1)
                        continue;
                    rrlc.second.length =
                        AngleBetweenDirections(rrlc.second.normalizedAnchors.front(), rrlc.second.normalizedAnchors.back());
                    mg.addConstraint(std::move(rrlc.second), Get<0>(rrlc.first), Get<1>(rrlc.first), Get<2>(rrlc.first));
                }



                // add RRs
                // todo 
                //THERE_ARE_BUGS_HERE("is panoramic image ok?");
                std::map<std::pair<int, int>, std::vector<std::vector<Pixel>>> boundaryEdges =
                    FindRegionBoundaries(segmentedRegions, samplerSizeOnBoundary, false);

                for (auto & bep : boundaryEdges) {
                    auto edges = bep.second;

                    RHandle rh1 = regId2Rh[bep.first.first];
                    RHandle rh2 = regId2Rh[bep.first.second];
                    if (rh1.invalid() || rh2.invalid())
                        continue;

                    RRData rrd;
                    rrd.normalizedEdges.resize(edges.size());
                    for (int k = 0; k < edges.size(); k++) {
                        rrd.normalizedEdges[k].reserve(edges[k].size());
                        for (auto & p : edges[k]) {
                            rrd.normalizedEdges[k].push_back(normalize(cam.toSpace(p)));
                        }
                    }

                    rrd.length = 0;
                    for (auto & e : rrd.normalizedEdges) {
                        assert(!e.empty() && "edges should never be empty!");
                        // get edge point projections
                        for (int i = 0; i < e.size() - 1; i++) {
                            rrd.length += AngleBetweenDirections(e[i], e[i + 1]);
                        }
                    }

                    rrd.normalizedSampledPoints.resize(rrd.normalizedEdges.size());
                    for (int k = 0; k < rrd.normalizedEdges.size(); k++) {
                        auto & edge = rrd.normalizedEdges[k];
                        if (edge.empty()) {
                            continue;
                        }
                        std::vector<Vec3> points = { edge.front() };
                        for (auto & edgeP : edge) {
                            double remainedAngle = AngleBetweenDirections(points.back(), edgeP);
                            while (remainedAngle >= samplingStepAngleOnBoundary) {
                                points.push_back(normalize(RotateDirection(points.back(), edgeP, samplingStepAngleOnBoundary)));
                                remainedAngle -= samplingStepAngleOnBoundary;
                            }
                        }
                        rrd.normalizedSampledPoints[k] = std::move(points);
                    }

                    rrd.occDetectionResult = RRData::Connected;
                    mg.addConstraint(std::move(rrd), rh1, rh2);
                }

                
                HandledTable<RRHandle, std::vector<std::vector<int>>> rrSamplePointsTouched(mg.internalConstraints<RRData>().size());
                for (auto & rr : mg.constraints<RRData>()) {
                    auto & touched = rrSamplePointsTouched[rr.topo.hd];
                    touched.resize(rr.data.normalizedSampledPoints.size());
                    for (int i = 0; i < touched.size(); i++) {
                        touched[i].resize(rr.data.normalizedSampledPoints[i].size(), 0);
                    }
                }

                RTree<Box3, std::tuple<RRHandle, int, int>> rrSamplingPointsTree;
                for (auto & rr : mg.constraints<RRData>()) {
                    for (int i = 0; i < rr.data.normalizedSampledPoints.size(); i++) {
                        auto & ps = rr.data.normalizedSampledPoints[i];
                        for (int j = 0; j < ps.size(); j++) {
                            auto & p = ps[j];
                            rrSamplingPointsTree.insert(BoundingBox(p).expand(2 * samplingStepAngleOnOcc / cam.focal()), std::make_tuple(rr.topo.hd, i, j));
                        }
                    }
                }

                for (int i = 0; i < occ.size(); i++) {
                    auto & points = occ[i].points;
                    std::vector<Vec3> pts; // sampling pts
                    if (points.empty()) {
                        continue;
                    }
                    auto last = points.front();
                    for (int j = 1; j < points.size(); j++) {
                        double dist = AngleBetweenDirections(last, points[j]);
                        while (dist >= samplingStepAngleOnOcc) {
                            auto p = RotateDirection(last, points[j], samplingStepAngleOnOcc);
                            pts.push_back(normalize(p));
                            dist -= samplingStepAngleOnOcc;
                            last = p;
                        }
                    }
                    // for debugging
                    if(false){
                        double originalLen = 0.0;
                        for (int j = 1; j < points.size(); j++) {
                            double dist = AngleBetweenDirections(points[j-1], points[j]);
                            originalLen += dist;
                        }
                        double len = 0.0;
                        for (int j = 1; j < pts.size(); j++) {
                            double dist = AngleBetweenDirections(pts[j - 1], pts[j]);
                            len += dist;
                        }
                        std::cout << "original len: " << originalLen << " len:" << len << std::endl;
                    }
                    for (auto & p : pts) {
                        rrSamplingPointsTree.search(BoundingBox(p).expand(2 * samplerSizeOnOcc / cam.focal()),
                            [&p, &rrSamplePointsTouched, &mg, samplerSizeOnOcc, &cam](const std::tuple<RRHandle, int, int> & spId) {
                            const Vec3 & pos = mg.data(Get<0>(spId)).normalizedSampledPoints[Get<1>(spId)][Get<2>(spId)];
                            if (Distance(p, pos) <= samplerSizeOnOcc / cam.focal()) {
                                rrSamplePointsTouched[Get<0>(spId)][Get<1>(spId)][Get<2>(spId)] = 1;
                            }
                            return true;
                        });                        
                    }
                }
                for (auto & rd : mg.constraints<RRData>()) {
                    int touchedCount = 0;
                    int allCount = 0;
                    for (auto & ts : rrSamplePointsTouched[rd.topo.hd]) {
                        for (int t : ts) {
                            allCount++;
                            if (t) {
                                touchedCount++;
                            }
                        }
                    }
                    if (touchedCount > 0.9 * allCount) {
                        rd.data.occDetectionResult = RRData::NotConnected;
                    }
                }



                // add RRRs and RRRRs
                for (auto && bhp : ExtractBoundaryJunctions(segmentedRegions)) {
                    assert(bhp.first.size() == 3 || bhp.first.size() == 4);
                    Vec3 dir = normalize(cam.toSpace(bhp.second));
                    if (bhp.first.size() == 3) {
                        RHandle rhs[] = {
                            regId2Rh[bhp.first[0]],
                            regId2Rh[bhp.first[1]],
                            regId2Rh[bhp.first[2]]
                        };
                        if (std::any_of(std::begin(rhs), std::end(rhs), [](RHandle rh) {return rh.invalid(); }))
                            continue;
                        RRRData rrrd;
                        rrrd.normalizedCenter = dir;
                        mg.addConstraint(std::move(rrrd), rhs[0], rhs[1], rhs[2]);
                    } else {
                        RHandle rhs[] = {
                            regId2Rh[bhp.first[0]],
                            regId2Rh[bhp.first[1]],
                            regId2Rh[bhp.first[2]],
                            regId2Rh[bhp.first[3]]
                        };
                        if (std::any_of(std::begin(rhs), std::end(rhs), [](RHandle rh) {return rh.invalid(); }))
                            continue;
                        RRRRData rrrrd;
                        rrrrd.normalizedCenter = dir;
                        mg.addConstraint(std::move(rrrrd), rhs[0], rhs[1], rhs[2], rhs[3]);
                    }
                }

                return regId2Rh;

            }


            template <class CameraT>
            RLGraph BuildRLGraph(const RLGraphBuilder::Params & params, const std::vector<Vec3> & vps,
                const std::vector<std::vector<Classified<Line2>>> & lines, const std::vector<DenseMatd> & lineVPScores,
                const std::vector<PerspectiveCamera> & cams, const Imagei & segmentedRegions, const Image5d & gc,
                const std::vector<Chain3> & occ, const CameraT & cam, std::vector<RHandle> * segId2Rhs /*= nullptr*/) {
                
                RLGraph g;
                assert(AllSameInContainer({ lines.size(), lineVPScores.size(), cams.size() }));
                for (int i = 0; i < lines.size(); i++) {
                    AppendLines(g, lines[i], lineVPScores[i], cams[i], vps,
                        params.intersectionAngleThresholdForLL,
                        params.incidenceParaAngleThresholdForLL,
                        params.incidenceVertAngleThresholdForLL,
                        params.incidenceParaAngleThresholdForLLAcrossViews,
                        params.incidenceVertAngleThresholdForLLAcrossViews);
                }
                auto rhs = AppendRegionsTemplate(g, segmentedRegions, gc, occ, cam,
                    params.samplingStepAngleOnRR,
                    params.samplingStepAngleOnRL,
                    params.samplingStepAngleOnOcc,
                    params.samplingStepSizeOnRR,
                    params.samplingStepSizeOnRL,
                    params.samplingStepSizeOnOcc);
                if (segId2Rhs) {
                    *segId2Rhs = std::move(rhs);
                }
                return g;

            }
        }




        RLGraph RLGraphBuilder::operator()(const std::vector<Vec3> & vps,
            const std::vector<std::vector<Classified<Line2>>> & lines, const std::vector<DenseMatd> & lineVPScores, 
            const std::vector<PerspectiveCamera> & cams, const Imagei & segmentedRegions, const Image5d & gc,
            const std::vector<Chain3> & occ, const PerspectiveCamera & cam, std::vector<RHandle> * segId2Rhs /*= nullptr*/) const {
            return BuildRLGraph(_params, vps, lines, lineVPScores, cams, segmentedRegions, gc, occ, cam, segId2Rhs);
        }

        RLGraph RLGraphBuilder::operator()(const std::vector<Vec3> & vps,
            const std::vector<std::vector<Classified<Line2>>> & lines, const std::vector<DenseMatd> & lineVPScores, 
            const std::vector<PerspectiveCamera> & cams, const Imagei & segmentedRegions, const Image5d & gc,
            const std::vector<Chain3> & occ, const PanoramicCamera & cam, std::vector<RHandle> * segId2Rhs /*= nullptr*/) const {
            return BuildRLGraph(_params, vps, lines, lineVPScores, cams, segmentedRegions, gc, occ, cam, segId2Rhs);
        }




    }
}