
extern "C" {
    #include <mosek.h>
}
#include <Eigen/Dense>
#include <Eigen/Sparse>

#include "algorithms.hpp"
#include "containers.hpp"
#include "matlab.hpp"
#include "mixed_graph.hpp"

#include "../vis/visualizers.hpp"
#include "matlab.hpp"

namespace panoramix {
    namespace core {

        double MGUnaryVariable::rawDepth() const{
            double sqSum = 0.0;
            for (auto v : variables){
                sqSum += v * v;
            }
            return 1.0 / sqrt(sqSum);
        }

        Plane3 MGUnaryVariable::interpretAsPlane() const {
            return Plane3FromEquation(variables[0], variables[1], variables[2]);
        }

        Line3 MGUnaryVariable::interpretAsLine(const MGUnary & line, const std::vector<Vec3> & vps) const {
            if (line.lineClaz >= 0){
                InfiniteLine3 infLine(line.normalizedCenter / variables[0], vps[line.lineClaz]);
                return Line3(DistanceBetweenTwoLines(InfiniteLine3(Point3(0, 0, 0), line.normalizedCorners.front()), infLine).second.second,
                    DistanceBetweenTwoLines(InfiniteLine3(Point3(0, 0, 0), line.normalizedCorners.back()), infLine).second.second);
            }
            else{
            /*           | sin(theta) | | p | | q |
                len:---------------------------------   => 1/len: [(1/q)sin(phi) - (1/p)sin(phi-theta)] / sin(theta)
                    | p sin(phi) - q sin(phi - theta) |
            */
                // variables[0] -> 1/p
                // variables[1] -> 1/q
                return Line3(line.normalizedCorners.front() / variables[0], line.normalizedCorners.back() / variables[1]);
            }
        }

        std::vector<double> MGUnaryVariable::variableCoeffsForInverseDepthAtDirection(const Vec3 & direction,
            const MGUnary & unary, const std::vector<Vec3> & vps) const{
            if (unary.type == MGUnary::Region){
                assert(variables.size() == 3);
                // depth = 1.0 / (ax + by + cz) where (x, y, z) = direction, (a, b, c) = variables
                // -> 1.0/depth = ax + by + cz
                return std::vector<double>{direction[0], direction[1], direction[2]};
            }
            else /*if (unary.type == MGUnary::Line)*/
            if (unary.lineClaz >= 0){
                assert(variables.size() == 1);
                const auto & line = unary;
                InfiniteLine3 infLine(line.normalizedCenter, vps[line.lineClaz]);
                // variable is 1.0/centerDepth
                // corresponding coeff is 1.0/depthRatio
                // so that 1.0/depth = 1.0/centerDepth * 1.0/depthRatio -> depth = centerDepth * depthRatio
                double depthRatio = norm(DistanceBetweenTwoLines(InfiniteLine3(Point3(0, 0, 0), direction), infLine).second.first);
                return std::vector<double>{1.0 / depthRatio};
            }
            else /*if(unary.lineClaz == -1)*/{
                assert(variables.size() == 2);
                const auto & line = unary;
                double theta = AngleBetweenDirections(line.normalizedCorners.front(), line.normalizedCorners.back());
                double phi = AngleBetweenDirections(line.normalizedCorners.front(), direction);
                /*           | sin(theta) | | p | | q |
                    len:---------------------------------   => 1/len: [(1/q)sin(phi) - (1/p)sin(phi-theta)] / sin(theta)
                       | p sin(phi) - q sin(phi - theta) |
                */
                // variables[0] -> 1/p
                // variables[1] -> 1/q
                double coeffFor1_p = -sin(phi - theta) / sin(theta);
                double coeffFor1_q = sin(phi) / sin(theta);
                assert(!IsInfOrNaN(coeffFor1_p) && !IsInfOrNaN(coeffFor1_q));
                return std::vector<double>{coeffFor1_p, coeffFor1_q};
            }
        }

        double MGUnaryVariable::inverseDepthAtDirection(const Vec3 & direction,
            const MGUnary & unary, const std::vector<Vec3> & vps) const {
            if (unary.type == MGUnary::Region){
                assert(variables.size() == 3);
                // depth = 1.0 / (ax + by + cz) where (x, y, z) = direction, (a, b, c) = variables
                // -> 1.0/depth = ax + by + cz
                return variables[0] * direction[0] + variables[1] * direction[1] + variables[2] * direction[2];
            }
            else /*if (unary.type == MGUnary::Line)*/
            if(unary.lineClaz >= 0){
                assert(variables.size() == 1);
                const auto & line = unary;
                InfiniteLine3 infLine(line.normalizedCenter, vps[line.lineClaz]);
                // variable is 1.0/centerDepth
                // corresponding coeff is 1.0/depthRatio
                // so that 1.0/depth = 1.0/centerDepth * 1.0/depthRatio -> depth = centerDepth * depthRatio
                double depthRatio = norm(DistanceBetweenTwoLines(InfiniteLine3(Point3(0, 0, 0), direction), infLine).second.first);
                return variables[0] / depthRatio; 
            }
            else /*if(unary.lineClaz == -1)*/ {
                assert(variables.size() == 2);
                const auto & line = unary;
                double theta = AngleBetweenDirections(line.normalizedCorners.front(), line.normalizedCorners.back());
                double phi = AngleBetweenDirections(line.normalizedCorners.front(), direction);
                /*           | sin(theta) | | p | | q |
                len:---------------------------------   => 1/len: [(1/q)sin(phi) - (1/p)sin(phi-theta)] / sin(theta)
                | p sin(phi) - q sin(phi - theta) |
                */
                // variables[0] -> 1/p
                // variables[1] -> 1/q
                double coeffFor1_p = -sin(phi - theta) / sin(theta);
                double coeffFor1_q = sin(phi) / sin(theta);
                assert(!IsInfOrNaN(coeffFor1_p) && !IsInfOrNaN(coeffFor1_q));
                return variables[0] * coeffFor1_p + variables[1] * coeffFor1_q;
            }
        }

        double MGUnaryVariable::depthAtCenter(const MGUnary & unary, const std::vector<Vec3> & vps) const {
            return 1.0 / inverseDepthAtDirection(unary.normalizedCenter, unary, vps);
        }


        MixedGraph BuildMixedGraph(const std::vector<View<PerspectiveCamera>> & views,
            const std::vector<RegionsGraph> & regionsGraphs, const std::vector<LinesGraph> & linesGraphs,
            const ComponentIndexHashMap<std::pair<RegionIndex, RegionIndex>, double> & regionOverlappingsAcrossViews,
            const ComponentIndexHashMap<std::pair<LineIndex, LineIndex>, Vec3> & lineIncidencesAcrossViews,
            const std::vector<std::map<std::pair<RegionHandle, LineHandle>, std::vector<Point2>>> & regionLineConnections,
            const std::vector<Vec3> & vps,
            MGUnaryVarTable & unaryVars, MGBinaryVarTable& binaryVars,
            double initialDepth){

            MixedGraph mg;
            ComponentIndexHashMap<RegionIndex, MGUnaryHandle> ri2mgh;
            ComponentIndexHashMap<LineIndex, MGUnaryHandle> li2mgh;
            std::unordered_map<MGUnaryHandle, RegionIndex> mgh2ri;
            std::unordered_map<MGUnaryHandle, LineIndex> mgh2li;

            // add components in each view
            for (int i = 0; i < views.size(); i++){
                auto & cam = views[i].camera;
                // regions
                for (auto & rd : regionsGraphs[i].elements<0>()){
                    auto ri = RegionIndex{ i, rd.topo.hd };
                    std::vector<Vec3> normalizedContour;
                    for (auto & ps : rd.data.contours){
                        for (auto & p : ps){
                            normalizedContour.push_back(normalize(cam.spatialDirection(p)));
                        }
                    }
                    if (normalizedContour.size() <= 2){
                        continue;
                    }
                    auto center = normalize(cam.spatialDirection(rd.data.center));
                    int initialClaz = std::distance(vps.begin(), std::min_element(vps.begin(), vps.end(),
                        [&center](const Vec3 & vp1, const Vec3 & vp2) -> bool {
                        return AngleBetweenUndirectedVectors(center, vp1) <
                            AngleBetweenUndirectedVectors(center, vp2);
                    }));
                    ri2mgh[ri] = mg.add(MGUnary{
                        MGUnary::Region,
                        std::move(normalizedContour),
                        center,
                        -1
                    });
                    mgh2ri[ri2mgh[ri]] = ri;
                    int sign = center.dot(vps[initialClaz]) < 0 ? -1 : 1;
                    unaryVars[ri2mgh[ri]].variables = { 
                        sign * vps[initialClaz][0],
                        sign * vps[initialClaz][1], 
                        sign * vps[initialClaz][2] 
                    };
                    unaryVars[ri2mgh[ri]].fixed = false;
                }
                // lines
                for (auto & ld : linesGraphs[i].elements<0>()){
                    auto li = LineIndex{ i, ld.topo.hd };
                    li2mgh[li] = mg.add(MGUnary{
                        MGUnary::Line,
                        {
                            normalize(cam.spatialDirection(ld.data.line.component.first)),
                            normalize(cam.spatialDirection(ld.data.line.component.second))
                        },
                        normalize(cam.spatialDirection(ld.data.line.component.center())),
                        ld.data.line.claz
                    });
                    mgh2li[li2mgh[li]] = li;
                    unaryVars[li2mgh[li]].variables = ld.data.line.claz == -1 ? 
                        std::vector<double>{ 1.0, 1.0 } : std::vector<double>{ 1.0 };
                    unaryVars[li2mgh[li]].fixed = false;
                }
                // region-region in each view
                for (auto & bd : regionsGraphs[i].elements<1>()){
                    MGBinary rr;
                    rr.type = MGBinary::RegionRegionConnection;
                    rr.weight = 1.0;
                    rr.normalizedAnchors.reserve(bd.data.sampledPoints.size());
                    for (auto & ps : bd.data.sampledPoints){
                        for (auto & p : ps){
                            rr.normalizedAnchors.push_back(normalize(cam.spatialDirection(p)));
                        }
                    }
                    auto r1 = RegionIndex{ i, bd.topo.lowers.front() };
                    auto r2 = RegionIndex{ i, bd.topo.lowers.back() };
                    if (!Contains(ri2mgh, r1) || !Contains(ri2mgh, r2))
                        continue;
                    mg.add<1>({ ri2mgh[r1], ri2mgh[r2] }, std::move(rr));
                }
                // region-line
                for (auto & regionLine : regionLineConnections[i]){
                    MGBinary rl;
                    rl.type = MGBinary::RegionLineConnection;
                    rl.weight = 1.0;
                    rl.normalizedAnchors.reserve(regionLine.second.size());
                    for (auto & p : regionLine.second){
                        rl.normalizedAnchors.push_back(normalize(cam.spatialDirection(p)));
                    }
                    auto ri = RegionIndex{ i, regionLine.first.first };
                    auto li = RegionIndex{ i, regionLine.first.second };
                    if (!Contains(ri2mgh, ri))
                        continue;
                    mg.add<1>({ ri2mgh[ri], li2mgh[li] }, std::move(rl));
                }
                // line-line
                for (auto & rd : linesGraphs[i].elements<1>()){
                    auto l1 = LineIndex{ i, rd.topo.lowers.front() };
                    auto l2 = LineIndex{ i, rd.topo.lowers.back() };
                    if (rd.data.type == LineRelationData::Intersection){
                        MGBinary llinter;
                        llinter.type = MGBinary::LineLineIntersection;
                        llinter.weight = rd.data.junctionWeight * 10;
                        llinter.normalizedAnchors = { normalize(cam.spatialDirection(rd.data.relationCenter)) };
                        mg.add<1>({ li2mgh[l1], li2mgh[l2] }, std::move(llinter));
                    }
                    else if (rd.data.type == LineRelationData::Incidence){
                        MGBinary llincid;
                        llincid.type = MGBinary::LineLineIncidence;
                        llincid.weight = rd.data.junctionWeight * 10;
                        llincid.normalizedAnchors = { normalize(cam.spatialDirection(rd.data.relationCenter)) };
                        mg.add<1>({ li2mgh[l1], li2mgh[l2] }, std::move(llincid));
                    }
                }
            }
            // add cross view constraints
            // region-region overlappings
            for (auto & regionOverlapping : regionOverlappingsAcrossViews){
                if (regionOverlapping.second < 0.2)
                    continue;
                MGBinary rro;
                rro.type = MGBinary::RegionRegionOverlapping;
                rro.weight = 100;
                auto & r1 = regionOverlapping.first.first;
                auto & r2 = regionOverlapping.first.second;

                if (!Contains(ri2mgh, r1) || !Contains(ri2mgh, r2))
                    continue;

                // get samples
                Vec3 z = mg.data(ri2mgh[r1]).normalizedCenter + mg.data(ri2mgh[r2]).normalizedCenter;
                Vec3 x, y;
                std::tie(x, y) = ProposeXYDirectionsFromZDirection(z);
                double minx = std::numeric_limits<double>::max(),
                    miny = std::numeric_limits<double>::max();
                double maxx = std::numeric_limits<double>::lowest(),
                    maxy = std::numeric_limits<double>::lowest();
                rro.normalizedAnchors.resize(4, z);
                for (auto & a : mg.data(ri2mgh[r1]).normalizedCorners){
                    double dx = a.dot(x), dy = a.dot(y);
                    if (dx < minx){
                        rro.normalizedAnchors[0] = a;
                        minx = dx;
                    }
                    else if (dx > maxx){
                        rro.normalizedAnchors[1] = a;
                        maxx = dx;
                    }
                    if (dy < miny){
                        rro.normalizedAnchors[2] = a;
                        miny = dy;
                    }
                    else if (dy > maxy){
                        rro.normalizedAnchors[3] = a;
                        maxy = dy;
                    }
                }
                for (auto & a : mg.data(ri2mgh[r2]).normalizedCorners){
                    double dx = a.dot(x), dy = a.dot(y);
                    if (dx < minx){
                        rro.normalizedAnchors[0] = a;
                        minx = dx;
                    }
                    else if (dx > maxx){
                        rro.normalizedAnchors[1] = a;
                        maxx = dx;
                    }
                    if (dy < miny){
                        rro.normalizedAnchors[2] = a;
                        miny = dy;
                    }
                    else if (dy > maxy){
                        rro.normalizedAnchors[3] = a;
                        maxy = dy;
                    }
                }
                mg.add<1>({ ri2mgh[r1], ri2mgh[r2] }, std::move(rro));
            }
            // line-line incidencs
            for (auto & lineIncidence : lineIncidencesAcrossViews){
                MGBinary llincid;
                llincid.type = MGBinary::LineLineIncidence;
                llincid.weight = IncidenceJunctionWeight(true) * 10;
                llincid.normalizedAnchors = { normalize(lineIncidence.second) };
                auto & l1 = lineIncidence.first.first;
                auto & l2 = lineIncidence.first.second;
                mg.add<1>({ li2mgh[l1], li2mgh[l2] }, std::move(llincid));
            }

            // compute importance ratios
            std::vector<double> unaryWeightSums(mg.internalElements<0>().size(), 0.0);
            for (auto & b : mg.internalElements<1>()){
                unaryWeightSums[b.topo.lowers.front().id] += b.data.weight * b.data.normalizedAnchors.size();
                unaryWeightSums[b.topo.lowers.back().id] += b.data.weight * b.data.normalizedAnchors.size();
            }
            for (auto & b : mg.internalElements<1>()){
                b.data.importanceRatioInRelatedUnaries.front() = b.data.weight * b.data.normalizedAnchors.size() /
                    unaryWeightSums[b.topo.lowers.front().id];
                b.data.importanceRatioInRelatedUnaries.back() = b.data.weight * b.data.normalizedAnchors.size() /
                    unaryWeightSums[b.topo.lowers.back().id];
            }

#ifdef _DEBUG
            for (auto & u : mg.elements<0>()){
                double importanceRatioSum = 0.0;
                for (auto & bh : u.topo.uppers){
                    importanceRatioSum +=
                        mg.data(bh).importanceRatioInRelatedUnaries[u.topo.hd == mg.topo(bh).lowers[0] ? 0 : 1];
                }
                assert(FuzzyEquals(importanceRatioSum, 1.0, 0.01));
            }
#endif
            for (auto & b : mg.internalElements<1>()){
                binaryVars[b.topo.hd].enabled = true;
            }

            return mg;
        }


        MixedGraph BuildMixedGraph(const std::vector<View<PerspectiveCamera>> & views,
            std::vector<Vec3> & vps,
            MGUnaryVarTable & unaryVars, MGBinaryVarTable& binaryVars,

            double initialDepth,
            const core::LineSegmentExtractor & lineseger,
            double intersectionDistanceThreshold,
            double incidenceDistanceAlongDirectionThreshold,
            double incidenceDistanceVerticalDirectionThreshold,

            const core::SegmentationExtractor & segmenter,
            double samplingStepLengthOnBoundary ,
            double samplingStepLengthOnLines,
            int dilationSize,

            double interViewIncidenceAngleAlongDirectionThreshold,
            double interViewIncidenceAngleVerticalDirectionThreshold){

            std::vector<core::LinesGraph> linesGraphs;
            core::EstimateVanishingPointsAndBuildLinesGraphs(views, vps, linesGraphs, lineseger,
                intersectionDistanceThreshold, incidenceDistanceAlongDirectionThreshold, incidenceDistanceVerticalDirectionThreshold, 
                false, true);

            std::vector<core::Imagei> segmentedRegionsArray;
            std::vector<core::RegionsGraph> regionsGraphs;

            for (int i = 0; i < views.size(); i++){
                auto & v = views[i];
                std::vector<core::Line2> lines;
                for (auto & ld : linesGraphs[i].elements<0>()){
                    auto line = ld.data.line.component;
                    line.first -= core::normalize(line.direction()) * 15.0;
                    line.second += core::normalize(line.direction()) * 15.0;
                    lines.push_back(line);
                }
                auto segmentedRegions = segmenter(v.image, lines).first;
                int samplePointsOnBoundariesSum = 0;
                auto regions = core::CreateRegionsGraph(segmentedRegions, samplingStepLengthOnBoundary, dilationSize);
                for (auto & r : regions.elements<1>()){
                    for (auto & ps : r.data.sampledPoints){
                        samplePointsOnBoundariesSum += ps.size();
                    }
                }
                regionsGraphs.push_back(std::move(regions));
                segmentedRegionsArray.push_back(segmentedRegions);
            }

            std::vector<std::map<std::pair<core::RegionHandle, core::LineHandle>, std::vector<core::Point2>>>
                regionLineConnectionsArray(views.size());
            for (int i = 0; i < views.size(); i++){
                regionLineConnectionsArray[i] = 
                    core::RecognizeRegionLineConnections(segmentedRegionsArray[i], linesGraphs[i], samplingStepLengthOnLines);
            }

            auto regionOverlappingsAcrossViews =
                core::RecognizeRegionOverlappingsAcrossViews(views, regionsGraphs);
            auto lineIncidencesAcrossViews =
                core::RecognizeLineIncidencesAcrossViews(views, linesGraphs,
                interViewIncidenceAngleAlongDirectionThreshold, interViewIncidenceAngleVerticalDirectionThreshold);

            return BuildMixedGraph(views, regionsGraphs, linesGraphs,
                regionOverlappingsAcrossViews, lineIncidencesAcrossViews, regionLineConnectionsArray, vps,
                unaryVars, binaryVars, initialDepth);

        }




        bool BinaryHandlesAreValidInPatch(const MixedGraph & mg, const MGPatch & patch) {
            for (auto & bhv : patch.bhs){
                auto bh = bhv.first;
                if (!Contains(patch.uhs, mg.topo(bh).lowers.front()) ||
                    !Contains(patch.uhs, mg.topo(bh).lowers.back()))
                    return false;
            }
            return true;
        }

        bool UnariesAreConnectedInPatch(const MixedGraph & mg, const MGPatch & patch){
            std::unordered_map<MGUnaryHandle, bool> visited;
            visited.reserve(patch.uhs.size());
            for (auto & uhv : patch.uhs){
                visited[uhv.first] = false;
            }
            std::vector<MGUnaryHandle> uhs;
            uhs.reserve(patch.uhs.size());
            for (auto & uhv : patch.uhs){
                uhs.push_back(uhv.first);
            }
            BreadthFirstSearch(uhs.begin(), uhs.end(), [&mg, &patch](const MGUnaryHandle & uh){
                std::vector<MGUnaryHandle> neighbors;
                for (auto & bh : mg.topo(uh).uppers){
                    if (Contains(patch.bhs, bh)){
                        auto anotherUh = mg.topo(bh).lowers.front();
                        if (anotherUh == uh)
                            anotherUh = mg.topo(bh).lowers.back();
                        neighbors.push_back(anotherUh);
                    }
                }
                return neighbors;
            }, [&visited](MGUnaryHandle uh){
                visited[uh] = true;
                return true;
            });
            for (auto & v : visited){
                if (!v.second)
                    return false;
            }
            return true;
        }


        MGPatch MakePatchOnBinary(const MixedGraph & mg, const MGBinaryHandle & bh,
            const MGUnaryVarTable & unaryVars, const MGBinaryVarTable & binaryVars){
            MGPatch patch;
            patch.bhs[bh] = binaryVars.at(bh);
            auto uhs = mg.topo(bh).lowers;
            patch.uhs[uhs[0]] = unaryVars.at(uhs[0]);
            patch.uhs[uhs[1]] = unaryVars.at(uhs[1]);
            assert(BinaryHandlesAreValidInPatch(mg, patch));
            assert(UnariesAreConnectedInPatch(mg, patch));
            return patch;
        }


        MGPatch MakeStarPatchAroundUnary(const MixedGraph & mg, const MGUnaryHandle & uh,
            const MGUnaryVarTable & unaryVars, const MGBinaryVarTable & binaryVars){
            MGPatch patch;
            patch.uhs[uh] = unaryVars.at(uh);
            for (auto & bh : mg.topo(uh).uppers){
                auto anotherUh = mg.topo(bh).lowers.front();
                if (anotherUh == uh)
                    anotherUh = mg.topo(bh).lowers.back();
                patch.bhs[bh] = binaryVars.at(bh);
                patch.uhs[anotherUh] = unaryVars.at(uh);
            }
            assert(BinaryHandlesAreValidInPatch(mg, patch));
            assert(UnariesAreConnectedInPatch(mg, patch));
            return patch;
        }


        double AnchorDistanceSumOnBinaryOfPatch(const MixedGraph & mg, const MGBinaryHandle & bh, 
            const MGPatch & patch, const std::vector<Vec3> & vps){
            assert(Contains(patch, bh));
            auto uh1 = mg.topo(bh).lowers.front();
            auto uh2 = mg.topo(bh).lowers.back();
            double distanceSum = 0.0;
            for (auto & a : mg.data(bh).normalizedAnchors){
                distanceSum += abs(1.0 / patch.uhs.at(uh1).inverseDepthAtDirection(a, mg.data(uh1), vps) - 
                    1.0 / patch.uhs.at(uh2).inverseDepthAtDirection(a, mg.data(uh2), vps));
            }
            return distanceSum;
        }

        double BinaryDistanceOfPatch(const MixedGraph & mg, const MGBinaryHandle & bh,
            const MGPatch & patch, const std::vector<Vec3> & vps){
            assert(Contains(patch, bh));
            auto uh1 = mg.topo(bh).lowers.front();
            auto uh2 = mg.topo(bh).lowers.back();
            double distanceSum = 0.0;
            for (auto & a : mg.data(bh).normalizedAnchors){
                distanceSum += abs(1.0 / patch.uhs.at(uh1).inverseDepthAtDirection(a, mg.data(uh1), vps) -
                    1.0 / patch.uhs.at(uh2).inverseDepthAtDirection(a, mg.data(uh2), vps));
                assert(!core::IsInfOrNaN(distanceSum));
            }
            return distanceSum / mg.data(bh).normalizedAnchors.size();
        }

        double AverageBinaryDistanceOfPatch(const MixedGraph & mg,
            const MGPatch & patch, const std::vector<Vec3> & vps){
            double distanceSum = 0.0;
            for (auto & bhv : patch.bhs){
                distanceSum += BinaryDistanceOfPatch(mg, bhv.first, patch, vps);
            }
            return distanceSum / patch.bhs.size();
        }

        double AverageUnaryCenterDepthOfPatch(const MixedGraph & mg, const MGPatch & patch, const std::vector<Vec3> & vps){
            double depthSum = 0.0;
            for (auto & uhv : patch.uhs){
                depthSum += uhv.second.depthAtCenter(mg.data(uhv.first), vps);
            }
            return depthSum / patch.uhs.size();
        }

        double AverageRawDepthOfPatch(const MGPatch & patch){
            double depthSum = 0.0;
            for (auto & uhv : patch.uhs){
                depthSum += uhv.second.rawDepth();
            }
            return depthSum / patch.uhs.size();
        }


        void ScalePatch(MGPatch & patch, double scale){
            for (auto & uhv : patch.uhs)
            for (auto & v : uhv.second.variables){
                v /= scale;
            }
        }



        std::vector<MGPatch> SplitMixedGraphIntoPatches(const MixedGraph & mg,
            const MGUnaryVarTable & unaryVars, const MGBinaryVarTable & binaryVars){

            std::unordered_map<MGUnaryHandle, int> ccids;
            std::vector<MGUnaryHandle> uhs;
            uhs.reserve(mg.internalElements<0>().size());
            for (auto & u : mg.elements<0>()){
                uhs.push_back(u.topo.hd);
            }
            int ccNum = core::ConnectedComponents(uhs.begin(), uhs.end(), [&mg](MGUnaryHandle uh){
                std::vector<MGUnaryHandle> neighbors;
                for (auto & bh : mg.topo(uh).uppers){
                    auto anotherUh = mg.topo(bh).lowers.front();
                    if (anotherUh == uh)
                        anotherUh = mg.topo(bh).lowers.back();
                    neighbors.push_back(anotherUh);
                }
                return neighbors;
            }, [&ccids](MGUnaryHandle uh, int ccid){
                ccids[uh] = ccid;
            });

            std::vector<MGPatch> patches(ccNum);
            for (auto & uhccid : ccids){
                auto uh = uhccid.first;
                int ccid = uhccid.second;
                patches[ccid].uhs[uh] = unaryVars.at(uh);
            }
            for (auto & b : mg.elements<1>()){
                auto & uhs = mg.topo(b.topo.hd).lowers;
                assert(ccids[uhs.front()] == ccids[uhs.back()]);
                patches[ccids[uhs.front()]].bhs[b.topo.hd] = binaryVars.at(b.topo.hd);
            }

            for (auto & p : patches){
                assert(UnariesAreConnectedInPatch(mg, p));
                assert(BinaryHandlesAreValidInPatch(mg, p));
            }

            return patches;
        }


        std::vector<MGPatch> SplitPatch(const MixedGraph & mg, const MGPatch & patch,
            std::function<bool(MGBinaryHandle bh)> useBh){

            std::vector<MGUnaryHandle> uhs;
            uhs.reserve(patch.uhs.size());
            for (auto & uhv : patch.uhs){
                uhs.push_back(uhv.first);
            }

            std::unordered_map<MGUnaryHandle, int> ccids;
            int ccNum = core::ConnectedComponents(uhs.begin(), uhs.end(), [&mg, &patch, &useBh](MGUnaryHandle uh){
                std::vector<MGUnaryHandle> neighbors;
                for (auto & bh : mg.topo(uh).uppers){
                    if (!Contains(patch.bhs, bh))
                        continue;
                    if (!useBh(bh))
                        continue;
                    auto anotherUh = mg.topo(bh).lowers.front();
                    if (anotherUh == uh)
                        anotherUh = mg.topo(bh).lowers.back();
                    neighbors.push_back(anotherUh);
                }
                return neighbors;
            }, [&ccids](MGUnaryHandle uh, int ccid){
                ccids[uh] = ccid;
            });

            std::vector<MGPatch> patches(ccNum);
            for (auto & uhccid : ccids){
                auto uh = uhccid.first;
                int ccid = uhccid.second;
                patches[ccid].uhs[uh] = patch.uhs.at(uh);
            }
            for (auto & b : mg.elements<1>()){
                auto & uhs = mg.topo(b.topo.hd).lowers;
                if (ccids[uhs.front()] == ccids[uhs.back()]){
                    patches[ccids[uhs.front()]].bhs[b.topo.hd] = patch.bhs.at(b.topo.hd);
                }
            }

            for (auto & p : patches){
                assert(UnariesAreConnectedInPatch(mg, p));
                assert(BinaryHandlesAreValidInPatch(mg, p));
            }

            return patches;
        }




        MGPatch MinimumSpanningTreePatch(const MixedGraph & mg, const MGPatch & patch,
            std::function<bool(MGBinaryHandle, MGBinaryHandle)> compareBh){

            assert(UnariesAreConnectedInPatch(mg, patch));
            assert(BinaryHandlesAreValidInPatch(mg, patch));

            MGPatch mst;
            mst.uhs = patch.uhs;

            std::vector<MGUnaryHandle> uhs;
            uhs.reserve(patch.uhs.size());
            for (auto & uhv : patch.uhs){
                uhs.push_back(uhv.first);
            }
            std::vector<MGBinaryHandle> bhs;
            bhs.reserve(patch.bhs.size());
            for (auto & bhv : patch.bhs){
                bhs.push_back(bhv.first);
            }

            std::vector<MGBinaryHandle> bhsReserved;

            core::MinimumSpanningTree(uhs.begin(), uhs.end(), bhs.begin(), bhs.end(),
                std::back_inserter(bhsReserved),
                [&mg](const MGBinaryHandle & bh){
                return mg.topo(bh).lowers;
            }, compareBh);

            for (auto & bh : bhsReserved){
                mst.bhs[bh] = patch.bhs.at(bh);
            }

            assert(UnariesAreConnectedInPatch(mg, mst));
            assert(BinaryHandlesAreValidInPatch(mg, mst));

            return mst;
        }




        namespace {



            template <class UhColorizerFunT = core::ConstantFunctor<vis::Color>>
            void ManuallyOptimizeMixedGraph(const core::Image & panorama,
                const core::MixedGraph & mg,
                core::MGPatch & patch,
                const std::vector<core::Vec3> & vps,
                UhColorizerFunT uhColorizer = UhColorizerFunT(vis::ColorTag::White),
                bool optimizeInEachIteration = false) {

                bool modified = true;
                auto sppCallbackFun = [&patch, &vps, &mg, &modified](vis::InteractionID iid,
                    const std::pair<core::MGUnaryHandle, vis::Colored<vis::SpatialProjectedPolygon>> & spp) {
                    std::cout << "uh: " << spp.first.id << std::endl;
                   /* if (iid == vis::InteractionID::PressSpace){
                        std::cout << "space pressed!" << std::endl;
                        int & claz = patch.uhs[spp.first].claz;
                        claz = (claz + 1) % vps.size();
                        std::cout << "current orientation is : " << vps[claz] << std::endl;
                        modified = true;
                    }*/
                };

                while (modified){

                    vis::ResourceStore::set("texture", panorama);

                    modified = false;
                    if (optimizeInEachIteration){
                        core::MGPatchDepthsOptimizer(mg, patch, vps, false, core::MGPatchDepthsOptimizer::EigenSparseQR)
                            .optimize();
                    }
                    patch /= core::AverageRawDepthOfPatch(patch);

                    vis::Visualizer viz("mixed graph optimizable");
                    viz.renderOptions.bwColor = 1.0;
                    viz.renderOptions.bwTexColor = 0.0;
                    viz.installingOptions.discretizeOptions.colorTable = vis::ColorTableDescriptor::RGB;
                    std::vector<std::pair<core::MGUnaryHandle, vis::Colored<vis::SpatialProjectedPolygon>>> spps;
                    std::vector<vis::Colored<core::Line3>> lines;

                    for (auto & uhv : patch.uhs){
                        auto uh = uhv.first;
                        auto & v = mg.data(uh);
                        if (v.type == core::MGUnary::Region){
                            auto & region = v;
                            vis::SpatialProjectedPolygon spp;
                            // filter corners
                            core::ForeachCompatibleWithLastElement(region.normalizedCorners.begin(), region.normalizedCorners.end(),
                                std::back_inserter(spp.corners),
                                [](const core::Vec3 & a, const core::Vec3 & b) -> bool {
                                return core::AngleBetweenDirections(a, b) > M_PI / 100.0;
                            });
                            if (spp.corners.size() < 3)
                                continue;

                            spp.projectionCenter = core::Point3(0, 0, 0);
                            spp.plane = uhv.second.interpretAsPlane();
                            spps.emplace_back(uh, std::move(vis::ColorAs(spp, uhColorizer(uh))));
                        }
                        else if (v.type == core::MGUnary::Line){
                            auto & line = v;
                            lines.push_back(vis::ColorAs(uhv.second.interpretAsLine(mg.data(uh), vps), uhColorizer(uh)));
                        }
                    }

                    viz.begin(spps, sppCallbackFun).shaderSource(vis::OpenGLShaderSourceDescriptor::XPanorama).resource("texture").end();
                    viz.installingOptions.lineWidth = 4.0;
                    viz.add(lines);

                   /* std::vector<core::Line3> connectionLines;
                    for (auto & bhv : patch.bhs){
                        auto bh = bhv.first;
                        auto & v = bhv.second;
                        auto & samples = mg.data(bh).normalizedAnchors;
                        for (int i = 0; i < samples.size(); i++){
                            connectionLines.emplace_back(normalize(samples[i]) * v.sampleDepthsOnRelatedUnaries.front()[i],
                                normalize(samples[i]) * v.sampleDepthsOnRelatedUnaries.back()[i]);
                        }
                    }*/

                    viz.installingOptions.discretizeOptions.color = vis::ColorTag::DarkGray;
                    viz.installingOptions.lineWidth = 2.0;
                    viz.renderOptions.renderMode = vis::RenderModeFlag::Triangles | vis::RenderModeFlag::Lines;
                    //viz.add(connectionLines);
                    viz.camera(core::PerspectiveCamera(800, 800, 500, { 1.0, 1.0, -1.0 }, { 0.0, 0.0, 0.0 }, { 0.0, 0.0, -1.0 }));
                    viz.renderOptions.backgroundColor = vis::ColorTag::Black;
                    viz.show(true, false);

                    vis::ResourceStore::clear();
                }

            }


            template <class T>
            inline std::vector<T> SelectSubset(const std::vector<T> & data, int numLimit){
                assert(numLimit > 0);
                if (data.size() <= numLimit)
                    return data;
                if (numLimit == 1)
                    return std::vector<T>{data[data.size() / 2]};
                if (numLimit == 2)
                    return std::vector<T>{data.front(), data.back()};
                int step = (data.size() + numLimit - 1) / numLimit;
                std::vector<T> filtered;
                filtered.reserve(numLimit);
                for (int i = 0; i < data.size(); i += step){
                    filtered.push_back(data[i]);
                }
                if (filtered.size() == numLimit - 1)
                    filtered.push_back(data.back());
                return filtered;
            }

            std::vector<Vec3> NecessaryAnchorsForBinary(const MixedGraph & mg, MGBinaryHandle bh){
                auto & b = mg.data(bh);
                if (b.type == MGBinary::LineLineIncidence || b.type == MGBinary::LineLineIntersection){
                    assert(b.normalizedAnchors.size() == 1);
                    return b.normalizedAnchors;
                }
                if (b.type == MGBinary::RegionRegionOverlapping){
                    assert(b.normalizedAnchors.size() >= 3);
                    return { b.normalizedAnchors[0], b.normalizedAnchors[1], b.normalizedAnchors[2] };
                }
                if (b.type == MGBinary::RegionLineConnection){
                    return { b.normalizedAnchors.front(), b.normalizedAnchors.back() };
                }
                if (b.type == MGBinary::RegionRegionConnection){
                    Vec3 alignDir = normalize(b.normalizedAnchors.front().cross(b.normalizedAnchors.back()));
                    double maxDotProd = 0.0;
                    int maxOffsetedAnchorId = -1;
                    for (int k = 1; k < b.normalizedAnchors.size() - 1; k++){
                        double dotProd = abs(alignDir.dot(b.normalizedAnchors[k]));
                        if (dotProd > maxDotProd){
                            maxDotProd = dotProd;
                            maxOffsetedAnchorId = k;
                        }
                    }
                    if (maxDotProd > 5e-1){
                        return{ b.normalizedAnchors.front(), b.normalizedAnchors[maxOffsetedAnchorId], b.normalizedAnchors.back() };
                    }
                    else {
                        return{ b.normalizedAnchors.front(), b.normalizedAnchors.back() };
                    }
                }
                else{
                    return{};
                }
            }


            




            template <class SparseMatElementT>
            void FormulateConstraintsAsMatrices(const MixedGraph & mg, MGPatch & patch,
                const std::vector<Vec3> & vanishingPoints,
                std::unordered_map<MGUnaryHandle, int> & uh2varStartPosition,
                std::unordered_map<MGBinaryHandle, int> & bh2consStartPosition,
                std::unordered_map<MGBinaryHandle, std::vector<Vec3>> & appliedBinaryAnchors,
                int & varNum, int & consNum,
                std::vector<SparseMatElementT> & Atriplets,
                std::vector<SparseMatElementT> & Wtriplets,
                Eigen::VectorXd & B,
                bool addAnchor = true){

                assert(BinaryHandlesAreValidInPatch(mg, patch));
                assert(UnariesAreConnectedInPatch(mg, patch));

                varNum = 0;
                bool hasFixedUnary = false;
                for (auto & uhv : patch.uhs){
                    if (uhv.second.fixed){
                        hasFixedUnary = true;
                        continue;
                    }
                    uh2varStartPosition[uhv.first] = varNum;
                    varNum += uhv.second.variables.size();
                }

                consNum = 0;

                if (addAnchor){
                    if (!hasFixedUnary){
                        consNum++;
                    }
                }
                for (auto & bhv : patch.bhs){
                    if (!bhv.second.enabled)
                        continue;
                    auto bh = bhv.first;
                    auto uh1 = mg.topo(bh).lowers.front();
                    auto uh2 = mg.topo(bh).lowers.back();
                    auto & u1 = mg.data(uh1);
                    auto & u2 = mg.data(uh2);

                    bool u1IsFixed = !Contains(uh2varStartPosition, uh1);
                    bool u2IsFixed = !Contains(uh2varStartPosition, uh2);
                    if (u1IsFixed && u2IsFixed)
                        continue;

                    bh2consStartPosition[bhv.first] = consNum;
                    appliedBinaryAnchors[bhv.first] = NecessaryAnchorsForBinary(mg, bh);
                    consNum += appliedBinaryAnchors[bhv.first].size();
                }

                //A.resize(consNum, varNum);
                //W.resize(consNum, consNum);
                B.resize(consNum);

                Atriplets.reserve(consNum * 6);
                Wtriplets.reserve(consNum);

                // write equations
                int eid = 0;
                if (addAnchor){
                    if (!hasFixedUnary){ // the anchor constraint
                        MGUnaryHandle uh = patch.uhs.begin()->first;
                        auto & uhVar = patch.uhs.begin()->second;
                        int uhVarNum = uhVar.variables.size();
                        Vec3 uhCenter = mg.data(uh).normalizedCenter;
                        auto uhVarCoeffsAtCenter = uhVar.variableCoeffsForInverseDepthAtDirection(uhCenter, mg.data(uh), vanishingPoints);
                        assert(uhVarCoeffsAtCenter.size() == uhVar.variables.size());
                        int uhVarStartPosition = uh2varStartPosition.at(uh);
                        for (int i = 0; i < uhVarCoeffsAtCenter.size(); i++){
                            //A.insert(eid, uhVarStartPosition + i) = uhVarCoeffsAtCenter[i];
                            Atriplets.emplace_back(eid, uhVarStartPosition + i, uhVarCoeffsAtCenter[i]);
                        }
                        B(eid) = 1.0;
                        //W.insert(eid, eid) = 1.0;
                        Wtriplets.emplace_back(eid, eid, 1.0);
                        eid++;
                    }
                }
                for (auto & bhv : patch.bhs){
                    if (!bhv.second.enabled)
                        continue;

                    auto & bh = bhv.first;
                    auto uh1 = mg.topo(bh).lowers.front();
                    auto uh2 = mg.topo(bh).lowers.back();
                    auto & u1 = mg.data(uh1);
                    auto & u2 = mg.data(uh2);

                    bool u1IsFixed = !Contains(uh2varStartPosition, uh1);
                    bool u2IsFixed = !Contains(uh2varStartPosition, uh2);

                    if (u1IsFixed && u2IsFixed){
                        continue;
                    }

                    int u1VarStartPosition = u1IsFixed ? -1 : uh2varStartPosition.at(uh1);
                    auto & u1Var = patch.uhs.at(uh1);
                    int u1VarNum = u1Var.variables.size();

                    int u2VarStartPosition = u2IsFixed ? -1 : uh2varStartPosition.at(uh2);
                    auto & u2Var = patch.uhs.at(uh2);
                    int u2VarNum = u2Var.variables.size();

                    for (auto & a : appliedBinaryAnchors.at(bh)){

                        B(eid) = 0.0;
                        assert(mg.data(bh).weight >= 0.0);
                        //W.insert(eid, eid) = mg.data(bh).weight;
                        Wtriplets.emplace_back(eid, eid, mg.data(bh).weight);

                        if (u1IsFixed){
                            double inverseDepthAtA = u1Var.inverseDepthAtDirection(a, u1, vanishingPoints);
                            B(eid) = -inverseDepthAtA;
                        }
                        else{
                            auto u1VarCoeffs = u1Var.variableCoeffsForInverseDepthAtDirection(a, u1, vanishingPoints);
                            assert(u1VarCoeffs.size() == u1VarNum);
                            for (int i = 0; i < u1VarCoeffs.size(); i++){
                                //A.insert(eid, u1VarStartPosition + i) = u1VarCoeffs[i]; // pos
                                Atriplets.emplace_back(eid, u1VarStartPosition + i, u1VarCoeffs[i]);
                            }
                        }

                        if (u2IsFixed){
                            double inverseDepthAtA = u2Var.inverseDepthAtDirection(a, u2, vanishingPoints);
                            B(eid) = inverseDepthAtA;
                        }
                        else{
                            auto u2VarCoeffs = u2Var.variableCoeffsForInverseDepthAtDirection(a, u2, vanishingPoints);
                            assert(u2VarCoeffs.size() == u2VarNum);
                            for (int i = 0; i < u2VarCoeffs.size(); i++){
                                //A.insert(eid, u2VarStartPosition + i) = -u2VarCoeffs[i]; // neg
                                Atriplets.emplace_back(eid, u2VarStartPosition + i, - u2VarCoeffs[i]);
                            }
                        }

                        eid++;
                    }
                }
                assert(eid == consNum);

            }



            



            struct MGPatchDepthsOptimizerInternalBase {
                virtual void initialize(const MixedGraph & mg, MGPatch & patch,
                    const std::vector<Vec3> & vanishingPoints, bool useWeights) = 0;
                virtual bool optimize(const MixedGraph & mg, MGPatch & patch,
                    const std::vector<Vec3> & vanishingPoints) = 0;
                virtual void finalize() = 0;
            };

            static struct Mosek {
                MSKenv_t env;
                inline Mosek() { MSK_makeenv(&env, nullptr); }
                inline ~Mosek() { MSK_deleteenv(&env); }
            } mosek;
            static void MSKAPI printstr(void *handle,
                MSKCONST char str[]){
                printf("%s", str);
            }

            struct MGPatchDepthsOptimizerInternalMosekSimplified : MGPatchDepthsOptimizerInternalBase {

                MSKtask_t task;
                int varNum, consNum;
                int inverseDepthVarNum, slackNum;
                std::unordered_map<MGUnaryHandle, int> uh2varStartPosition;
                std::unordered_map<MGBinaryHandle, int> bh2anchorStartPosition;
                std::unordered_map<MGBinaryHandle, std::vector<Vec3>> appliedBinaryAnchors;

                std::unordered_map<MGUnaryHandle, bool> uhFixed;
                std::map<int, std::vector<MGUnaryHandle>> uhCCs;
                std::unordered_map<MGUnaryHandle, int> uhCCIds;

                virtual void initialize(const MixedGraph & mg, MGPatch & patch,
                    const std::vector<Vec3> & vanishingPoints, bool useWeights){

                    assert(BinaryHandlesAreValidInPatch(mg, patch));
                    assert(UnariesAreConnectedInPatch(mg, patch));
                    
                    for (auto & uhv : patch.uhs){
                        uhFixed[uhv.first] = uhv.second.fixed;
                    }

                    // prepare registering constraints
                    // find cc with more than 3 necessary anchors
                    for (auto & bhv : patch.bhs){
                        if (!bhv.second.enabled)
                            continue;
                        auto bh = bhv.first;
                        auto uh1 = mg.topo(bh).lowers.front();
                        auto uh2 = mg.topo(bh).lowers.back();
                        auto & u1 = mg.data(uh1);
                        auto & u2 = mg.data(uh2);

                        bool u1IsFixed = uhFixed.at(uh1);
                        bool u2IsFixed = uhFixed.at(uh2);
                        if (u1IsFixed && u2IsFixed)
                            continue;

                        appliedBinaryAnchors[bhv.first] = NecessaryAnchorsForBinary(mg, bh);
                    }


                    {
                        std::vector<MGUnaryHandle> uhs;
                        uhs.reserve(patch.uhs.size());
                        for (auto & uhv : patch.uhs){
                            uhs.push_back(uhv.first);
                        }
                        core::ConnectedComponents(uhs.begin(), uhs.end(),
                            [&mg, &patch, this](MGUnaryHandle uh) -> std::vector<MGUnaryHandle> {
                            std::vector<MGUnaryHandle> neighbors;
                            for (auto & bh : mg.topo(uh).uppers){
                                if (!Contains(appliedBinaryAnchors, bh))
                                    continue;
                                if (appliedBinaryAnchors.at(bh).size() < 3) // use only STRONG connections!!
                                    continue;
                                MGUnaryHandle anotherUh = mg.topo(bh).lowers[0];
                                if (anotherUh == uh)
                                    anotherUh = mg.topo(bh).lowers[1];
                                if (!Contains(patch, uh))
                                    continue;
                                neighbors.push_back(anotherUh);
                            }
                            return neighbors;
                        }, [this](MGUnaryHandle uh, int ccId){
                            uhCCs[ccId].push_back(uh);
                            uhCCIds[uh] = ccId;
                        });
                    }


                    // spread the fix status and assign the fixed variable
                    for (auto & cc : uhCCs) {
                        auto & uhs = cc.second;
                        // collect fixed uhs in this cc
                        Vec3 abc(0, 0, 0);
                        int fixedNum = 0;
                        for (auto uh : uhs){
                            auto & uhVar = patch.uhs.at(uh);
                            if (uhVar.fixed){
                                Vec3 thisAbc(uhVar.variables[0], uhVar.variables[1], uhVar.variables[2]);
                                if (Distance(abc / fixedNum, thisAbc) > 1e-2){
                                    std::cerr << "variables of fixed unaries do not match!!!!" << std::endl;
                                }
                                abc += thisAbc;
                                fixedNum++;
                            }
                        }
                        if (fixedNum == 0)
                            continue;
                        abc /= fixedNum; // the averaged variables
                        for (auto uh : uhs){
                            assert(mg.data(uh).type == MGUnary::Region);
                            auto & uhVar = patch.uhs.at(uh);
                            if (!uhVar.fixed){
                                uhVar.variables = { abc[0], abc[1], abc[2] };
                            }
                            uhFixed[uh] = true;
                        }
                    }

                    // register non-fixed variables
                    // register same variable for all uhs in one CC
                    varNum = 0;
                    for (auto & uhv : patch.uhs){
                        if (uhFixed.at(uhv.first)){
                            continue;
                        }
                        int ccId = uhCCIds.at(uhv.first);
                        auto firstUhInThisCC = uhCCs.at(ccId).front();
                        if (firstUhInThisCC == uhv.first){ // this is the first uh
                            uh2varStartPosition[uhv.first] = varNum;
                            varNum += uhv.second.variables.size();
                        }
                        else{
                            uh2varStartPosition[uhv.first] = 
                                uh2varStartPosition.at(firstUhInThisCC); // share the variables!!
                        }
                    }

                    { // erase STRONG connections since they are used in CCs
                        for (auto it = appliedBinaryAnchors.cbegin(); it != appliedBinaryAnchors.cend();){
                            auto nextIt = std::next(it);
                            if (uhCCIds.at(mg.topo(it->first).lowers[0]) == uhCCIds.at(mg.topo(it->first).lowers[1])){
                                appliedBinaryAnchors.erase(it);
                            }
                            it = nextIt;
                        }
                    }

                    // register constraints
                    consNum = 0;
                    for (auto & bha : appliedBinaryAnchors){
                        bh2anchorStartPosition[bha.first] = consNum;
                        consNum += bha.second.size();
                    }


                    slackNum = consNum;
                    inverseDepthVarNum = varNum;

                    varNum += slackNum; // real vars and slacks
                    consNum *= 2;   // [equation = 0] < slackVar;  
                                    // -[equation = 0] < slackVar
                    
                    task = nullptr;
                    auto & env = mosek.env;

                    MSK_maketask(env, consNum, varNum, &task);
                    MSK_linkfunctotaskstream(task, MSK_STREAM_LOG, NULL, printstr);

                    MSK_appendcons(task, consNum);
                    MSK_appendvars(task, varNum);

                    assert(slackNum + inverseDepthVarNum == varNum);
                    assert(slackNum * 2 == consNum);

                    // set weights
                    {
                        //for (int i = 0; i < inverseDepthVarNum; i++)
                        //    MSK_putcj(task, i, 0.0);
                        int slackVarId = inverseDepthVarNum;
                        for (auto & bha : appliedBinaryAnchors){
                            auto & bd = mg.data(bha.first);
                            assert(bd.weight >= 0.0);
                            for (int k = 0; k < bha.second.size(); k++){
                                MSK_putcj(task, slackVarId, useWeights ? bd.weight : 1.0);
                                slackVarId++;
                            }
                        }
                        assert(slackVarId == slackNum + inverseDepthVarNum);
                        assert(slackVarId == varNum);
                    }

                    // bounds for vars
                    {
                        int varId = 0;
                        for (; varId < inverseDepthVarNum; varId++){ // for depths
                            MSK_putvarbound(task, varId, MSK_BK_LO, 1.0, +MSK_INFINITY);
                        }
                        //for (; varId < varNum; varId++){ // for slacks
                        //    MSK_putvarbound(task, varId, MSK_BK_LO, 1e-9, +MSK_INFINITY);
                        //}
                    }


                    { // bounds for constraints without any fixed uhs
                        int consId = 0;
                        for (; consId < consNum; consId++){
                            // [invdepth1 * invratio1 - invdepth2 * invratio2 - slackVar < 0]; 
                            // [- invdepth1 * invratio1 + invdepth2 * invratio2 - slackVar < 0];
                            MSK_putconbound(task, consId, MSK_BK_UP, -MSK_INFINITY, 0.0);
                        }
                    }


                    { // fill constraints related to each vars
                        for (auto & uhv : patch.uhs){ // depths
                            auto uh = uhv.first;
                            auto & relatedBhs = mg.topo(uh).uppers;
                            int relatedAnchorsNum = 0;
                            for (auto & bh : relatedBhs){
                                if (!Contains(appliedBinaryAnchors, bh))
                                    continue;
                                relatedAnchorsNum += appliedBinaryAnchors.at(bh).size();
                            }
                            int relatedConsNum = relatedAnchorsNum * 2;

                            if (uhFixed.at(uh)){
                                for (auto & bh : relatedBhs){
                                    if (!Contains(appliedBinaryAnchors, bh))
                                        continue;
                                    assert(!uhFixed.at(mg.topo(bh).lowers[0]) || !uhFixed.at(mg.topo(bh).lowers[1]));
                                    int firstAnchorPosition = bh2anchorStartPosition.at(bh);
                                    auto & samples = appliedBinaryAnchors.at(bh);
                                    for (int k = 0; k < samples.size(); k++){
                                        double curInverseDepth = uhv.second.inverseDepthAtDirection(samples[k], mg.data(uh), vanishingPoints);
                                        assert(!core::IsInfOrNaN(curInverseDepth));
                                        int consId1 = (firstAnchorPosition + k) * 2; // one for [invdepth1 * invratio1 - invdepth2 * invratio2 - slackVar < 0]; 
                                        int consId2 = (firstAnchorPosition + k) * 2 + 1; // another for [- invdepth1 * invratio1 + invdepth2 * invratio2 - slackVar < 0];
                                        bool isOnLeftSide = uh == mg.topo(bh).lowers.front();
                                        assert(consId1 < consNum && consId2 < consNum);
                                        if (isOnLeftSide){
                                            MSK_putconbound(task, consId1, MSK_BK_UP, -MSK_INFINITY, - curInverseDepth);
                                            MSK_putconbound(task, consId2, MSK_BK_UP, -MSK_INFINITY, + curInverseDepth);
                                        }
                                        else{
                                            MSK_putconbound(task, consId1, MSK_BK_UP, -MSK_INFINITY, + curInverseDepth);
                                            MSK_putconbound(task, consId2, MSK_BK_UP, -MSK_INFINITY, - curInverseDepth);
                                        }
                                    }
                                }
                            }
                            else{
                                int uhVarStartPosition = uh2varStartPosition.at(uh);
                                int uhVarNum = uhv.second.variables.size();

                                std::vector<MSKint32t> consIds;
                                std::vector<std::vector<MSKrealt>> consValues(uhVarNum);
                                consIds.reserve(relatedConsNum);
                                for (auto & cvalues : consValues)
                                    cvalues.reserve(relatedConsNum);

                                for (auto & bh : relatedBhs){
                                    if (!Contains(appliedBinaryAnchors, bh))
                                        continue;
                                    assert(!uhFixed.at(mg.topo(bh).lowers[0]) || !uhFixed.at(mg.topo(bh).lowers[1]));

                                    int firstAnchorPosition = bh2anchorStartPosition.at(bh);
                                    auto & samples = appliedBinaryAnchors.at(bh);
                                    assert(samples.size() == 1 || samples.size() == 2);

                                    for (int k = 0; k < samples.size(); k++){
                                        assert(firstAnchorPosition + k < slackNum);
                                        consIds.push_back((firstAnchorPosition + k) * 2); // one for [invdepth1 * invratio1 - invdepth2 * invratio2 - slackVar < 0]; 
                                        consIds.push_back((firstAnchorPosition + k) * 2 + 1); // another for [- invdepth1 * invratio1 + invdepth2 * invratio2 - slackVar < 0];

                                        auto & a = samples[k];
                                        auto uhVarCoeffs = uhv.second.variableCoeffsForInverseDepthAtDirection(a, mg.data(uh), vanishingPoints);
                                        assert(uhVarCoeffs.size() == uhVarNum);

                                        bool isOnLeftSide = uh == mg.topo(bh).lowers.front();
                                        if (isOnLeftSide){ // as depth1
                                            for (int m = 0; m < uhVarNum; m++){
                                                assert(!core::IsInfOrNaN(uhVarCoeffs[m]));
                                                consValues[m].push_back(uhVarCoeffs[m]);
                                                consValues[m].push_back(-uhVarCoeffs[m]);
                                            }
                                        }
                                        else{
                                            for (int m = 0; m < uhVarNum; m++){
                                                assert(!core::IsInfOrNaN(uhVarCoeffs[m]));
                                                consValues[m].push_back(-uhVarCoeffs[m]);
                                                consValues[m].push_back(uhVarCoeffs[m]);
                                            }
                                        }
                                    }
                                }

                                for (int consId : consIds)
                                    assert(consId < consNum);

                                for (int m = 0; m < uhVarNum; m++){
                                    assert(uhVarStartPosition + m < inverseDepthVarNum);
                                    MSK_putacol(task, uhVarStartPosition + m, relatedConsNum, consIds.data(), consValues[m].data());
                                }
                            }
                        }


                        for (int varId = inverseDepthVarNum; varId < inverseDepthVarNum + slackNum; varId++){ // slack vars
                            MSKint32t consIds[] = {
                                (varId - inverseDepthVarNum) * 2, // one for [depth1 * ratio1 - depth2 * ratio2 - slackVar < 0]; 
                                (varId - inverseDepthVarNum) * 2 + 1  // another for [- depth1 * ratio1 + depth2 * ratio2 - slackVar < 0];
                            };
                            MSKrealt consValues[] = { -1.0, -1.0 };
                            MSK_putacol(task, varId, 2, consIds, consValues);
                        }
                    }
                }

                virtual bool optimize(const MixedGraph & mg, MGPatch & patch,
                    const std::vector<Vec3> & vanishingPoints){                    

                    MSK_putobjsense(task, MSK_OBJECTIVE_SENSE_MINIMIZE);
                    MSKrescodee trmcode;
                    MSK_optimizetrm(task, &trmcode);
                    MSK_solutionsummary(task, MSK_STREAM_LOG);

                    MSKsolstae solsta;
                    MSK_getsolsta(task, MSK_SOL_BAS, &solsta);

                    //Tock(tick);

                    switch (solsta){
                    case MSK_SOL_STA_OPTIMAL:
                    case MSK_SOL_STA_NEAR_OPTIMAL:
                    {
                                                     double *xx = new double[varNum];
                                                     MSK_getxx(task, MSK_SOL_BAS, xx);
                                                     /*std::cout << "results: ";
                                                     for (int i = 0; i < varNum; i++)
                                                     std::cout << xx[i] << ' ';
                                                     std::cout << std::endl;*/
                                                     //printf("Optimal primal solution\n");
                                                     for (auto & uhv : patch.uhs){ // get resulted depths
                                                         if (uhFixed.at(uhv.first)){
                                                             continue;
                                                         }
                                                         int uhVarStartPosition = uh2varStartPosition.at(uhv.first);
                                                         for (int k = 0; k < uhv.second.variables.size(); k++){
                                                             uhv.second.variables[k] = xx[uhVarStartPosition + k];
                                                         }
                                                     }
                                                     delete[] xx;
                                                     return true;
                    }
                    case MSK_SOL_STA_DUAL_INFEAS_CER:
                    case MSK_SOL_STA_PRIM_INFEAS_CER:
                    case MSK_SOL_STA_NEAR_DUAL_INFEAS_CER:
                    case MSK_SOL_STA_NEAR_PRIM_INFEAS_CER:
                        printf("Primal or dual infeasibility certificate found.\n");
                        return false;
                    case MSK_SOL_STA_UNKNOWN:
                    {
                                                char symname[MSK_MAX_STR_LEN];
                                                char desc[MSK_MAX_STR_LEN];
                                                MSK_getcodedesc(trmcode,
                                                    symname,
                                                    desc);

                                                printf("The solution status is unknown.\n");
                                                printf("The optimizer terminitated with code: %s\n", symname);
                                                return false;
                    }
                    default:
                        printf("Other solution status.\n");
                        return false;
                    }

                }
                virtual void finalize() {
                    MSK_deletetask(&task);
                }
            };


            struct MGPatchDepthsOptimizerInternalMosek : MGPatchDepthsOptimizerInternalBase {
                
                MSKtask_t task;
                int varNum, consNum;
                int inverseDepthVarNum, slackNum;
                std::unordered_map<MGUnaryHandle, int> uh2varStartPosition;
                std::unordered_map<MGBinaryHandle, int> bh2anchorStartPosition;
                std::unordered_map<MGBinaryHandle, std::vector<Vec3>> appliedBinaryAnchors;
                
                virtual void initialize(const MixedGraph & mg, MGPatch & patch,
                    const std::vector<Vec3> & vanishingPoints, bool useWeights) {

                    assert(BinaryHandlesAreValidInPatch(mg, patch));
                    assert(UnariesAreConnectedInPatch(mg, patch));

                    // prepare registering constraints
                    // find cc with more than 3 necessary anchors
                    for (auto & bhv : patch.bhs){
                        if (!bhv.second.enabled)
                            continue;
                        auto bh = bhv.first;
                        auto uh1 = mg.topo(bh).lowers.front();
                        auto uh2 = mg.topo(bh).lowers.back();
                        auto & u1 = mg.data(uh1);
                        auto & u2 = mg.data(uh2);

                        //bool u1IsFixed = patch.uhs.at(uh1).fixed;
                        //bool u2IsFixed = patch.uhs.at(uh2).fixed;
                        //if (u1IsFixed && u2IsFixed)
                        //    continue;

                        auto anchors = NecessaryAnchorsForBinary(mg, bh);
                        appliedBinaryAnchors[bhv.first] = SelectSubset(anchors, 2);
                    }

                    // register non-fixed variables
                    // register same variable for all uhs in one CC
                    varNum = 0;
                    for (auto & uhv : patch.uhs){
                        //if (uhv.second.fixed){
                        //    continue;
                        //}
                        uh2varStartPosition[uhv.first] = varNum;
                        varNum += uhv.second.variables.size();
                    }

                    // register constraints
                    consNum = 0;
                    for (auto & bha : appliedBinaryAnchors){
                        bh2anchorStartPosition[bha.first] = consNum;
                        consNum += bha.second.size();
                    }


                    slackNum = consNum;
                    inverseDepthVarNum = varNum;

                    varNum += slackNum; // real vars and slacks
                    consNum *= 2;   // [equation = 0] < slackVar;  
                    // -[equation = 0] < slackVar

                    task = nullptr;
                    auto & env = mosek.env;

                    MSK_maketask(env, consNum, varNum, &task);
                    MSK_linkfunctotaskstream(task, MSK_STREAM_LOG, NULL, printstr);

                    MSK_appendcons(task, consNum);
                    MSK_appendvars(task, varNum);

                    assert(slackNum + inverseDepthVarNum == varNum);
                    assert(slackNum * 2 == consNum);

                    // set weights
                    {
                        int slackVarId = inverseDepthVarNum;
                        for (auto & bha : appliedBinaryAnchors){
                            auto & bd = mg.data(bha.first);
                            assert(bd.weight >= 0.0);
                            for (int k = 0; k < bha.second.size(); k++){
                                MSK_putcj(task, slackVarId, useWeights ? bd.weight : 1.0);
                                slackVarId++;
                            }
                        }
                        assert(slackVarId == slackNum + inverseDepthVarNum);
                    }

                    // bounds for vars
                    {
                        int varId = 0;
                        for (; varId < inverseDepthVarNum; varId++){ // for depths
                            MSK_putvarbound(task, varId, MSK_BK_LO, 1.0, +MSK_INFINITY);
                        }
                    }


                    { // bounds for constraints without any fixed uhs
                        int consId = 0;
                        for (; consId < consNum; consId++){
                            // [invdepth1 * invratio1 - invdepth2 * invratio2 - slackVar < 0]; 
                            // [- invdepth1 * invratio1 + invdepth2 * invratio2 - slackVar < 0];
                            MSK_putconbound(task, consId, MSK_BK_UP, -MSK_INFINITY, 0.0);
                        }
                    }


                    { // fill constraints related to each vars
                        for (auto & uhv : patch.uhs){ // depths
                            auto uh = uhv.first;
                            auto & relatedBhs = mg.topo(uh).uppers;
                            int relatedAnchorsNum = 0;
                            for (auto & bh : relatedBhs){
                                if (!Contains(appliedBinaryAnchors, bh))
                                    continue;
                                relatedAnchorsNum += appliedBinaryAnchors.at(bh).size();
                            }
                            int relatedConsNum = relatedAnchorsNum * 2;

                            //if (uhv.second.fixed){
                            //    for (auto & bh : relatedBhs){
                            //        if (!Contains(appliedBinaryAnchors, bh))
                            //            continue;

                            //        int firstAnchorPosition = bh2anchorStartPosition.at(bh);
                            //        auto & samples = appliedBinaryAnchors.at(bh);
                            //        for (int k = 0; k < samples.size(); k++){
                            //            double curInverseDepth = uhv.second.inverseDepthAtDirection(samples[k], mg.data(uh), vanishingPoints);
                            //            assert(!core::IsInfOrNaN(curInverseDepth));
                            //            int consId1 = (firstAnchorPosition + k) * 2; // one for [invdepth1 * invratio1 - invdepth2 * invratio2 - slackVar < 0]; 
                            //            int consId2 = (firstAnchorPosition + k) * 2 + 1; // another for [- invdepth1 * invratio1 + invdepth2 * invratio2 - slackVar < 0];
                            //            bool isOnLeftSide = uh == mg.topo(bh).lowers.front();
                            //            assert(consId1 < consNum && consId2 < consNum);
                            //            if (isOnLeftSide){
                            //                MSK_putconbound(task, consId1, MSK_BK_UP, -MSK_INFINITY, -curInverseDepth);
                            //                MSK_putconbound(task, consId2, MSK_BK_UP, -MSK_INFINITY, +curInverseDepth);
                            //            }
                            //            else{
                            //                MSK_putconbound(task, consId1, MSK_BK_UP, -MSK_INFINITY, +curInverseDepth);
                            //                MSK_putconbound(task, consId2, MSK_BK_UP, -MSK_INFINITY, -curInverseDepth);
                            //            }
                            //        }
                            //    }
                            //}
                            //else
                            {
                                int uhVarStartPosition = uh2varStartPosition.at(uh);
                                int uhVarNum = uhv.second.variables.size();

                                std::vector<MSKint32t> consIds;
                                std::vector<std::vector<MSKrealt>> consValues(uhVarNum);
                                consIds.reserve(relatedConsNum);
                                for (auto & cvalues : consValues)
                                    cvalues.reserve(relatedConsNum);

                                for (auto & bh : relatedBhs){
                                    if (!Contains(appliedBinaryAnchors, bh))
                                        continue;

                                    int firstAnchorPosition = bh2anchorStartPosition.at(bh);
                                    auto & samples = appliedBinaryAnchors.at(bh);
                                    assert(samples.size() == 1 || samples.size() == 2);

                                    for (int k = 0; k < samples.size(); k++){
                                        assert(firstAnchorPosition + k < slackNum);
                                        consIds.push_back((firstAnchorPosition + k) * 2); // one for [invdepth1 * invratio1 - invdepth2 * invratio2 - slackVar < 0]; 
                                        consIds.push_back((firstAnchorPosition + k) * 2 + 1); // another for [- invdepth1 * invratio1 + invdepth2 * invratio2 - slackVar < 0];

                                        auto & a = samples[k];
                                        auto uhVarCoeffs = uhv.second.variableCoeffsForInverseDepthAtDirection(a, mg.data(uh), vanishingPoints);
                                        assert(uhVarCoeffs.size() == uhVarNum);
                                        assert(!HasValue(uhVarCoeffs, core::IsInfOrNaN<double>));

                                        bool isOnLeftSide = uh == mg.topo(bh).lowers.front();
                                        if (isOnLeftSide){ // as depth1
                                            for (int m = 0; m < uhVarNum; m++){
                                                assert(!core::IsInfOrNaN(uhVarCoeffs[m]));
                                                consValues[m].push_back(uhVarCoeffs[m]);
                                                consValues[m].push_back(-uhVarCoeffs[m]);
                                            }
                                        }
                                        else{
                                            for (int m = 0; m < uhVarNum; m++){
                                                assert(!core::IsInfOrNaN(uhVarCoeffs[m]));
                                                consValues[m].push_back(-uhVarCoeffs[m]);
                                                consValues[m].push_back(uhVarCoeffs[m]);
                                            }
                                        }
                                    }
                                }

                                for (int consId : consIds)
                                    assert(consId < consNum);

                                for (int m = 0; m < uhVarNum; m++){
                                    assert(uhVarStartPosition + m < inverseDepthVarNum);
                                    MSK_putacol(task, uhVarStartPosition + m, relatedConsNum, consIds.data(), consValues[m].data());
                                }
                            }
                        }


                        for (int varId = inverseDepthVarNum; varId < inverseDepthVarNum + slackNum; varId++){ // slack vars
                            MSKint32t consIds[] = {
                                (varId - inverseDepthVarNum) * 2, // one for [depth1 * ratio1 - depth2 * ratio2 - slackVar < 0]; 
                                (varId - inverseDepthVarNum) * 2 + 1  // another for [- depth1 * ratio1 + depth2 * ratio2 - slackVar < 0];
                            };
                            MSKrealt consValues[] = { -1.0, -1.0 };
                            MSK_putacol(task, varId, 2, consIds, consValues);
                        }
                    }

                }

                virtual bool optimize(const MixedGraph & mg, MGPatch & patch,
                    const std::vector<Vec3> & vanishingPoints) {

                    MSK_putobjsense(task, MSK_OBJECTIVE_SENSE_MINIMIZE);
                    MSKrescodee trmcode;
                    MSK_optimizetrm(task, &trmcode);
                    MSK_solutionsummary(task, MSK_STREAM_LOG);

                    MSKsolstae solsta;
                    MSK_getsolsta(task, MSK_SOL_BAS, &solsta);

                    //Tock(tick);

                    switch (solsta){
                    case MSK_SOL_STA_OPTIMAL:
                    case MSK_SOL_STA_NEAR_OPTIMAL:
                    {
                                                     double *xx = new double[varNum];
                                                     MSK_getxx(task, MSK_SOL_BAS, xx);
                                                     /*std::cout << "results: ";
                                                     for (int i = 0; i < varNum; i++)
                                                     std::cout << xx[i] << ' ';
                                                     std::cout << std::endl;*/
                                                     //printf("Optimal primal solution\n");
                                                     for (auto & uhv : patch.uhs){ // get resulted depths
                                                        /* if (uhv.second.fixed){
                                                             continue;
                                                         }*/
                                                         int uhVarStartPosition = uh2varStartPosition.at(uhv.first);
                                                         for (int k = 0; k < uhv.second.variables.size(); k++){
                                                             uhv.second.variables[k] = xx[uhVarStartPosition + k];
                                                         }
                                                     }
                                                     delete[] xx;
                                                     return true;
                    }
                    case MSK_SOL_STA_DUAL_INFEAS_CER:
                    case MSK_SOL_STA_PRIM_INFEAS_CER:
                    case MSK_SOL_STA_NEAR_DUAL_INFEAS_CER:
                    case MSK_SOL_STA_NEAR_PRIM_INFEAS_CER:
                        printf("Primal or dual infeasibility certificate found.\n");
                        return false;
                    case MSK_SOL_STA_UNKNOWN:
                    {
                                                char symname[MSK_MAX_STR_LEN];
                                                char desc[MSK_MAX_STR_LEN];
                                                MSK_getcodedesc(trmcode,
                                                    symname,
                                                    desc);

                                                printf("The solution status is unknown.\n");
                                                printf("The optimizer terminitated with code: %s\n", symname);
                                                return false;
                    }
                    default:
                        printf("Other solution status.\n");
                        return false;
                    }

                }

                virtual void finalize() {
                    MSK_deletetask(&task);
                }
            };



            struct MGPatchDepthsOptimizerInternalEigenSparseQRSimplified {

                int varNum, consNum;
                int inverseDepthVarNum, slackNum;
                std::unordered_map<MGUnaryHandle, int> uh2varStartPosition;
                std::unordered_map<MGBinaryHandle, int> bh2anchorStartPosition;
                std::unordered_map<MGBinaryHandle, std::vector<Vec3>> appliedBinaryAnchors;

                std::unordered_map<MGUnaryHandle, bool> uhFixed;
                std::map<int, std::vector<MGUnaryHandle>> uhCCs;
                std::unordered_map<MGUnaryHandle, int> uhCCIds;

                Eigen::SparseMatrix<double> A, W;
                Eigen::VectorXd B;
                bool useWeights;


                virtual void initialize(const MixedGraph & mg, MGPatch & patch,
                    const std::vector<Vec3> & vanishingPoints, bool useWeights) {
                    assert(BinaryHandlesAreValidInPatch(mg, patch));
                    assert(UnariesAreConnectedInPatch(mg, patch));

                    for (auto & uhv : patch.uhs){
                        uhFixed[uhv.first] = uhv.second.fixed;
                    }

                    // prepare registering constraints
                    // find cc with more than 3 necessary anchors
                    for (auto & bhv : patch.bhs){
                        if (!bhv.second.enabled)
                            continue;
                        auto bh = bhv.first;
                        auto uh1 = mg.topo(bh).lowers.front();
                        auto uh2 = mg.topo(bh).lowers.back();
                        auto & u1 = mg.data(uh1);
                        auto & u2 = mg.data(uh2);

                        bool u1IsFixed = uhFixed.at(uh1);
                        bool u2IsFixed = uhFixed.at(uh2);
                        if (u1IsFixed && u2IsFixed)
                            continue;

                        appliedBinaryAnchors[bhv.first] = NecessaryAnchorsForBinary(mg, bh);
                    }

                    int ccNum = 0;
                    {
                        std::vector<MGUnaryHandle> uhs;
                        uhs.reserve(patch.uhs.size());
                        for (auto & uhv : patch.uhs){
                            uhs.push_back(uhv.first);
                        }
                        ccNum = core::ConnectedComponents(uhs.begin(), uhs.end(),
                            [&mg, &patch, this](MGUnaryHandle uh) -> std::vector<MGUnaryHandle> {
                            std::vector<MGUnaryHandle> neighbors;
                            for (auto & bh : mg.topo(uh).uppers){
                                if (!Contains(appliedBinaryAnchors, bh))
                                    continue;
                                if (appliedBinaryAnchors.at(bh).size() < 3) // use only STRONG connections!!
                                    continue;
                                //assert(mg.data(bh).type == MGBinary::RegionRegionOverlapping);
                                MGUnaryHandle anotherUh = mg.topo(bh).lowers[0];
                                if (anotherUh == uh)
                                    anotherUh = mg.topo(bh).lowers[1];
                                if (!Contains(patch, uh))
                                    continue;
                                neighbors.push_back(anotherUh);
                            }
                            return neighbors;
                        }, [this](MGUnaryHandle uh, int ccId){
                            uhCCs[ccId].push_back(uh);
                            uhCCIds[uh] = ccId;
                        });
                    }

                    auto ctable = vis::CreateRandomColorTableWithSize(ccNum);

                    {
                        core::Image panorama;
                        core::LoadFromDisk("./cache/test_view.View.SampleViews.panorama", panorama);
                        ManuallyOptimizeMixedGraph(panorama, mg, patch, vanishingPoints, [this, &ctable](MGUnaryHandle uh){
                            return ctable[uhCCIds.at(uh)];
                        });
                    }

                    // spread the fix status and assign the fixed variable
                    for (auto & cc : uhCCs) {
                        auto & uhs = cc.second;
                        // collect fixed uhs in this cc
                        Vec3 abc(0, 0, 0);
                        int fixedNum = 0;
                        for (auto uh : uhs){
                            auto & uhVar = patch.uhs.at(uh);
                            if (uhVar.fixed){
                                Vec3 thisAbc(uhVar.variables[0], uhVar.variables[1], uhVar.variables[2]);
                                if (Distance(abc / fixedNum, thisAbc) > 1e-2){
                                    std::cerr << "variables of fixed unaries do not match!!!!" << std::endl;
                                }
                                abc += thisAbc;
                                fixedNum++;
                            }
                        }
                        if (fixedNum == 0)
                            continue;
                        abc /= fixedNum; // the averaged variables
                        for (auto uh : uhs){
                            assert(mg.data(uh).type == MGUnary::Region);
                            auto & uhVar = patch.uhs.at(uh);
                            if (!uhVar.fixed){
                                uhVar.variables = { abc[0], abc[1], abc[2] };
                            }
                            uhFixed[uh] = true;
                        }
                    }


                    // register non-fixed variables
                    // register same variable for all uhs in one CC
                    varNum = 0;
                    bool hasFixedUnary = false;
                    for (auto & uhv : patch.uhs){
                        if (uhFixed.at(uhv.first)){
                            hasFixedUnary = true;
                            continue;
                        }
                        int ccId = uhCCIds.at(uhv.first);
                        auto firstUhInThisCC = uhCCs.at(ccId).front();
                        if (firstUhInThisCC == uhv.first){ // this is the first uh
                            uh2varStartPosition[uhv.first] = varNum;
                            varNum += uhv.second.variables.size();
                        }
                        else{
                            uh2varStartPosition[uhv.first] =
                                uh2varStartPosition.at(firstUhInThisCC); // share the variables!!
                        }
                    }

                    { // erase STRONG connections since they are used in CCs
                        for (auto it = appliedBinaryAnchors.cbegin(); it != appliedBinaryAnchors.cend();){
                            auto nextIt = std::next(it);
                            if (uhCCIds.at(mg.topo(it->first).lowers[0]) == uhCCIds.at(mg.topo(it->first).lowers[1])){
                                appliedBinaryAnchors.erase(it);
                            }
                            it = nextIt;
                        }
                    }

                    // register constraints
                    consNum = 0;
                    if (!hasFixedUnary){
                        consNum++;
                    }
                    for (auto & bha : appliedBinaryAnchors){
                        bh2anchorStartPosition[bha.first] = consNum;
                        consNum += bha.second.size();
                    }

                    A.resize(consNum, varNum);
                    W.resize(consNum, consNum);
                    B.resize(consNum);

                    A.reserve(consNum * 6);
                    W.reserve(consNum);

                    // write equations
                    int eid = 0;
                    if (!hasFixedUnary){ // the anchor constraint
                        MGUnaryHandle uh = patch.uhs.begin()->first;
                        auto & uhVar = patch.uhs.begin()->second;
                        int uhVarNum = uhVar.variables.size();
                        Vec3 uhCenter = mg.data(uh).normalizedCenter;
                        auto uhVarCoeffsAtCenter = uhVar.variableCoeffsForInverseDepthAtDirection(uhCenter, mg.data(uh), vanishingPoints);
                        assert(uhVarCoeffsAtCenter.size() == uhVar.variables.size());
                        int uhVarStartPosition = uh2varStartPosition.at(uh);
                        for (int i = 0; i < uhVarCoeffsAtCenter.size(); i++){
                            A.insert(eid, uhVarStartPosition + i) = uhVarCoeffsAtCenter[i];
                        }
                        B(eid) = 1.0;
                        W.insert(eid, eid) = 1.0;
                        eid++;
                    }
                    for (auto & bha : appliedBinaryAnchors){

                        auto & bh = bha.first;
                        auto uh1 = mg.topo(bh).lowers.front();
                        auto uh2 = mg.topo(bh).lowers.back();
                        auto & u1 = mg.data(uh1);
                        auto & u2 = mg.data(uh2);

                        bool u1IsFixed = uhFixed.at(uh1);
                        bool u2IsFixed = uhFixed.at(uh2);

                        assert(!u1IsFixed || !u2IsFixed);

                        int u1VarStartPosition = u1IsFixed ? -1 : uh2varStartPosition.at(uh1);
                        auto & u1Var = patch.uhs.at(uh1);
                        int u1VarNum = u1Var.variables.size();

                        int u2VarStartPosition = u2IsFixed ? -1 : uh2varStartPosition.at(uh2);
                        auto & u2Var = patch.uhs.at(uh2);
                        int u2VarNum = u2Var.variables.size();

                        bool u1IsFirstInItsCC = uhCCs.at(uhCCIds.at(uh1)).front() == uh1;
                        bool u2IsFirstInItsCC = uhCCs.at(uhCCIds.at(uh2)).front() == uh2;

                        for (auto & a : bha.second){

                            B(eid) = 0.0;
                            assert(mg.data(bh).weight >= 0.0);
                            W.insert(eid, eid) = mg.data(bh).weight;

                            if (u1IsFixed){
                                double inverseDepthAtA = u1Var.inverseDepthAtDirection(a, u1, vanishingPoints);
                                B(eid) = -inverseDepthAtA;
                            }
                            else{
                                //if (u1IsFirstInItsCC){
                                    auto u1VarCoeffs = u1Var.variableCoeffsForInverseDepthAtDirection(a, u1, vanishingPoints);
                                    assert(u1VarCoeffs.size() == u1VarNum);
                                    for (int i = 0; i < u1VarCoeffs.size(); i++){
                                        A.insert(eid, u1VarStartPosition + i) = u1VarCoeffs[i]; // pos
                                    }
                                //}
                            }

                            if (u2IsFixed){
                                double inverseDepthAtA = u2Var.inverseDepthAtDirection(a, u2, vanishingPoints);
                                B(eid) = inverseDepthAtA;
                            }
                            else{
                                //if (u2IsFirstInItsCC){
                                    auto u2VarCoeffs = u2Var.variableCoeffsForInverseDepthAtDirection(a, u2, vanishingPoints);
                                    assert(u2VarCoeffs.size() == u2VarNum);
                                    for (int i = 0; i < u2VarCoeffs.size(); i++){
                                        A.insert(eid, u2VarStartPosition + i) = -u2VarCoeffs[i]; // neg
                                    }
                               // }
                            }

                            eid++;
                        }
                    }
                    assert(eid == consNum);

                }

                virtual bool optimize(const MixedGraph & mg, MGPatch & patch,
                    const std::vector<Vec3> & vanishingPoints) {
                    using namespace Eigen;

                    SparseQR<Eigen::SparseMatrix<double>, COLAMDOrdering<int>> solver;
                    static_assert(!(Eigen::SparseMatrix<double>::IsRowMajor), "COLAMDOrdering only supports column major");
                    Eigen::SparseMatrix<double> WA = W * A;
                    A.makeCompressed();
                    WA.makeCompressed();

                    solver.compute(useWeights ? WA : A);

                    if (solver.info() != Success) {
                        assert(0);
                        std::cout << "computation error" << std::endl;
                        return false;
                    }

                    VectorXd WB = W * B;
                    VectorXd X = solver.solve(useWeights ? WB : B);
                    if (solver.info() != Success) {
                        assert(0);
                        std::cout << "solving error" << std::endl;
                        return false;
                    }

                    for (auto & uhv : patch.uhs){
                        if (!Contains(uh2varStartPosition, uhv.first))
                            continue;
                        int uhStartPosition = uh2varStartPosition.at(uhv.first);
                        for (int i = 0; i < uhv.second.variables.size(); i++){
                            uhv.second.variables[i] = X(uhStartPosition + i);
                        }
                    }

                    return true;

                }

                virtual void finalize() { }
            };


            struct MGPatchDepthsOptimizerInternalEigen : MGPatchDepthsOptimizerInternalBase {
                Eigen::SparseMatrix<double> A, W;
                Eigen::VectorXd B;
                bool useWeights;

                std::unordered_map<MGUnaryHandle, int> uh2varStartPosition;
                std::unordered_map<MGBinaryHandle, int> bh2consStartPosition;

                std::unordered_map<MGBinaryHandle, std::vector<Vec3>> appliedBinaryAnchors;

                virtual void initialize(const MixedGraph & mg, MGPatch & patch,
                    const std::vector<Vec3> & vanishingPoints, bool useWeights) override {

                    int varNum, consNum;
                    std::vector<Eigen::Triplet<double>> Atriplets, Wtriplets;
                    FormulateConstraintsAsMatrices(mg, patch, vanishingPoints, 
                        uh2varStartPosition, bh2consStartPosition, appliedBinaryAnchors, varNum, consNum, Atriplets, Wtriplets, B);
                    A.resize(consNum, varNum);
                    W.resize(consNum, consNum);
                    A.setFromTriplets(Atriplets.begin(), Atriplets.end());
                    W.setFromTriplets(Wtriplets.begin(), Wtriplets.end());

                }

                virtual bool optimize(const MixedGraph & mg, MGPatch & patch,
                    const std::vector<Vec3> & vanishingPoints)  override {

                    using namespace Eigen;

                    SparseQR<Eigen::SparseMatrix<double>, COLAMDOrdering<int>> solver;
                    static_assert(!(Eigen::SparseMatrix<double>::IsRowMajor), "COLAMDOrdering only supports column major");
                    Eigen::SparseMatrix<double> WA = W * A;
                    A.makeCompressed();
                    WA.makeCompressed();

                    solver.compute(useWeights ? WA : A);

                    if (solver.info() != Success) {
                        assert(0);
                        std::cout << "computation error" << std::endl;
                        return false;
                    }

                    VectorXd WB = W * B;
                    VectorXd X = solver.solve(useWeights ? WB : B);
                    if (solver.info() != Success) {
                        assert(0);
                        std::cout << "solving error" << std::endl;
                        return false;
                    }

                    for (auto & uhv : patch.uhs){
                        if (!Contains(uh2varStartPosition, uhv.first))
                            continue;
                        int uhStartPosition = uh2varStartPosition.at(uhv.first);
                        for (int i = 0; i < uhv.second.variables.size(); i++){
                            uhv.second.variables[i] = X(uhStartPosition + i);
                        }
                    }

                    return true;
                }
                virtual void finalize() {}
            };



            struct MGPatchDepthsOptimizerInternalMATLABCVX : MGPatchDepthsOptimizerInternalBase {
                SparseMat<double> A, W;
                Eigen::VectorXd B;
                bool useWeights;

                std::unordered_map<MGUnaryHandle, int> uh2varStartPosition;
                std::unordered_map<MGBinaryHandle, int> bh2consStartPosition;

                std::unordered_map<MGBinaryHandle, std::vector<Vec3>> appliedBinaryAnchors;

                virtual void initialize(const MixedGraph & mg, MGPatch & patch,
                    const std::vector<Vec3> & vanishingPoints, bool useWeights) override {

                    int varNum, consNum;
                    std::vector<SparseMatElement<double>> Atriplets, Wtriplets;
                    FormulateConstraintsAsMatrices(mg, patch, vanishingPoints,
                        uh2varStartPosition, bh2consStartPosition, appliedBinaryAnchors, varNum, consNum, Atriplets, Wtriplets, B);
                    A = MakeSparseMatFromElements(consNum, varNum, Atriplets.begin(), Atriplets.end());
                    W = MakeSparseMatFromElements(consNum, consNum, Wtriplets.begin(), Wtriplets.end());
                }

                virtual bool optimize(const MixedGraph & mg, MGPatch & patch,
                    const std::vector<Vec3> & vanishingPoints)  override {

                    using namespace Eigen;

                    core::Matlab::RunScript("clear");
                    
                    core::Matlab::PutVariable("A", A);
                    core::Matlab::PutVariable("W", W);
                    std::vector<double> Bv(B.size());
                    std::copy_n(B.data(), B.size(), Bv.begin());
                    core::Matlab::PutVariable("B", Bv);
                    core::Matlab::RunScript("save tempfile");

                    core::Matlab matlab;
                    matlab
                        << "n = size(A, 2);"
                        << "m = size(A, 1);"
                        << "cvx_setup;"
                        << "cvx_begin"
                        << "   variable X(n)"
                        //<< "   variable slack(m)"
                        << "   minimize(norm((A * X - B')))"
                        //<< "   subject to"
                        //<< "      ones(n, 1) <= X"
                        //<< "       <= slack"
                        //<< "      B' - A * X <= slack"
                        //<< "      zeros(m, 1) <= slack"
                        << "cvx_end"
                        << "X = X';";
                    
                    std::vector<double> X;
                    core::Matlab::GetVariable("X", X, false);

                    for (auto & uhv : patch.uhs){
                        if (!Contains(uh2varStartPosition, uhv.first))
                            continue;
                        int uhStartPosition = uh2varStartPosition.at(uhv.first);
                        for (int i = 0; i < uhv.second.variables.size(); i++){
                            uhv.second.variables[i] = X[uhStartPosition + i];
                        }
                    }

                    return true;
                }
                virtual void finalize() {}

            };



            struct MGPatchDepthsOptimizerInternalMATLABCVX_v2 : MGPatchDepthsOptimizerInternalBase {
                Eigen::SparseMatrix<double> unaryFixedToAnchorDepths;
                Eigen::SparseMatrix<double> unaryLeftToAnchorDepths, unaryRightToAnchorDepths;
                Eigen::SparseMatrix<double> anchorDepthsWeights;
                Eigen::VectorXd Rhs;
                bool useWeights;

                std::unordered_map<MGUnaryHandle, int> uh2varStartPosition;
                std::unordered_map<MGBinaryHandle, int> bh2consStartPosition;

                std::unordered_map<MGBinaryHandle, std::vector<Vec3>> appliedBinaryAnchors;
                
                virtual void initialize(const MixedGraph & mg, MGPatch & patch,
                    const std::vector<Vec3> & vanishingPoints, bool useWeights) override {

                    bool addAnchor = false;

                    assert(BinaryHandlesAreValidInPatch(mg, patch));
                    assert(UnariesAreConnectedInPatch(mg, patch));

                    int varNum = 0;
                    bool hasFixedUnary = false;
                    for (auto & uhv : patch.uhs){
                        if (uhv.second.fixed){
                            hasFixedUnary = true;
                            continue;
                        }
                        uh2varStartPosition[uhv.first] = varNum;
                        varNum += uhv.second.variables.size();
                    }

                    int consNum = 0;

                    if (addAnchor){
                        if (!hasFixedUnary){
                            consNum++;
                        }
                    }
                    for (auto & bhv : patch.bhs){
                        if (!bhv.second.enabled)
                            continue;
                        auto bh = bhv.first;
                        auto uh1 = mg.topo(bh).lowers.front();
                        auto uh2 = mg.topo(bh).lowers.back();
                        auto & u1 = mg.data(uh1);
                        auto & u2 = mg.data(uh2);

                        bool u1IsFixed = !Contains(uh2varStartPosition, uh1);
                        bool u2IsFixed = !Contains(uh2varStartPosition, uh2);
                        if (u1IsFixed && u2IsFixed)
                            continue;

                        bh2consStartPosition[bhv.first] = consNum;
                        appliedBinaryAnchors[bhv.first] = NecessaryAnchorsForBinary(mg, bh);
                        consNum += appliedBinaryAnchors[bhv.first].size();
                    }

                    unaryLeftToAnchorDepths.resize(consNum, varNum);
                    unaryRightToAnchorDepths.resize(consNum, varNum);
                    anchorDepthsWeights.resize(consNum, consNum);
                    Rhs.resize(consNum);

                    unaryLeftToAnchorDepths.reserve(consNum * 3);
                    unaryRightToAnchorDepths.reserve(consNum * 3);
                    anchorDepthsWeights.reserve(consNum);

                    // write equations
                    int eid = 0;
                    if (addAnchor){
                        if (!hasFixedUnary){ // the anchor constraint
                            MGUnaryHandle uh = patch.uhs.begin()->first;
                            auto & uhVar = patch.uhs.begin()->second;
                            int uhVarNum = uhVar.variables.size();
                            Vec3 uhCenter = mg.data(uh).normalizedCenter;
                            auto uhVarCoeffsAtCenter = uhVar.variableCoeffsForInverseDepthAtDirection(uhCenter, mg.data(uh), vanishingPoints);
                            assert(uhVarCoeffsAtCenter.size() == uhVar.variables.size());
                            int uhVarStartPosition = uh2varStartPosition.at(uh);
                            for (int i = 0; i < uhVarCoeffsAtCenter.size(); i++){
                                unaryLeftToAnchorDepths.insert(eid, uhVarStartPosition + i) = uhVarCoeffsAtCenter[i];
                            }
                            Rhs(eid) = 1.0;
                            anchorDepthsWeights.insert(eid, eid) = 1.0;
                            eid++;
                        }
                    }
                    for (auto & bhv : patch.bhs){
                        if (!bhv.second.enabled)
                            continue;

                        auto & bh = bhv.first;
                        auto uh1 = mg.topo(bh).lowers.front();
                        auto uh2 = mg.topo(bh).lowers.back();
                        auto & u1 = mg.data(uh1);
                        auto & u2 = mg.data(uh2);

                        bool u1IsFixed = !Contains(uh2varStartPosition, uh1);
                        bool u2IsFixed = !Contains(uh2varStartPosition, uh2);

                        if (u1IsFixed && u2IsFixed){
                            continue;
                        }

                        int u1VarStartPosition = u1IsFixed ? -1 : uh2varStartPosition.at(uh1);
                        auto & u1Var = patch.uhs.at(uh1);
                        int u1VarNum = u1Var.variables.size();

                        int u2VarStartPosition = u2IsFixed ? -1 : uh2varStartPosition.at(uh2);
                        auto & u2Var = patch.uhs.at(uh2);
                        int u2VarNum = u2Var.variables.size();

                        for (auto & a : appliedBinaryAnchors.at(bh)){

                            Rhs(eid) = 0.0;
                            assert(mg.data(bh).weight >= 0.0);
                            anchorDepthsWeights.insert(eid, eid) = mg.data(bh).weight;

                            if (u1IsFixed){
                                double inverseDepthAtA = u1Var.inverseDepthAtDirection(a, u1, vanishingPoints);
                                Rhs(eid) = -inverseDepthAtA;
                            }
                            else{
                                auto u1VarCoeffs = u1Var.variableCoeffsForInverseDepthAtDirection(a, u1, vanishingPoints);
                                assert(u1VarCoeffs.size() == u1VarNum);
                                for (int i = 0; i < u1VarCoeffs.size(); i++){
                                    unaryLeftToAnchorDepths.insert(eid, u1VarStartPosition + i) = u1VarCoeffs[i]; // pos
                                }
                            }

                            if (u2IsFixed){
                                double inverseDepthAtA = u2Var.inverseDepthAtDirection(a, u2, vanishingPoints);
                                Rhs(eid) = inverseDepthAtA;
                            }
                            else{
                                auto u2VarCoeffs = u2Var.variableCoeffsForInverseDepthAtDirection(a, u2, vanishingPoints);
                                assert(u2VarCoeffs.size() == u2VarNum);
                                for (int i = 0; i < u2VarCoeffs.size(); i++){
                                    unaryRightToAnchorDepths.insert(eid, u2VarStartPosition + i) = u2VarCoeffs[i]; // neg
                                }
                            }

                            eid++;
                        }
                    }
                    assert(eid == consNum);

                }

                static inline void putSparseMatrixIntoMatlab(const std::string & name, const Eigen::SparseMatrix<double> & m) {
                    assert(!m.isCompressed());
                    std::vector<double> i, j;
                    std::vector<double> s;
                    i.reserve(m.nonZeros());
                    j.reserve(m.nonZeros());
                    s.reserve(m.nonZeros());
                    for (int r = 0; r < m.rows(); r++){
                        for (int c = 0; c < m.cols(); c++){
                            if (m.coeff(r, c) != 0){
                                i.push_back(r + 1);
                                j.push_back(c + 1);
                                s.push_back(m.coeff(r, c));
                            }
                        }
                    }
                    core::Matlab::PutVariable(name + "_i", i);
                    core::Matlab::PutVariable(name + "_j", j);
                    core::Matlab::PutVariable(name + "_s", s);
                    core::Matlab::PutVariable(name + "_m", m.rows());
                    core::Matlab::PutVariable(name + "_n", m.cols());
                    core::Matlab::RunScript(name + "=" +
                        "sparse(" + name + "_i," + name + "_j," + name + "_s," + name + "_m," + name + "_n" + ");");
                    core::Matlab::RunScript("clear " + name + "_i " + name + "_j " + name + "_s " + name + "_m " + name + "_n;");
                }

                virtual bool optimize(const MixedGraph & mg, MGPatch & patch,
                    const std::vector<Vec3> & vanishingPoints)  override {

                    //using namespace Eigen;

                    //core::Matlab::RunScript("clear");

                    //putSparseMatrixIntoMatlab("A", A);
                    //putSparseMatrixIntoMatlab("W", W);
                    //std::vector<double> Bv(B.size());
                    //std::copy_n(B.data(), B.size(), Bv.begin());
                    //core::Matlab::PutVariable("B", Bv);
                    //core::Matlab::RunScript("save tempfile");

                    //core::Matlab matlab;
                    //matlab
                    //    << "n = size(A, 2);"
                    //    << "m = size(A, 1);"
                    //    << "cvx_setup;"
                    //    << "cvx_begin"
                    //    << "   variable X(n)"
                    //    << "   variable slack(m)"
                    //    << "   minimize(norm((A * X - B')))"
                    //    //<< "   subject to"
                    //    //<< "      ones(n, 1) <= X"
                    //    //<< "       <= slack"
                    //    //<< "      B' - A * X <= slack"
                    //    //<< "      zeros(m, 1) <= slack"
                    //    << "cvx_end"
                    //    << "X = X';";

                    //std::vector<double> X;
                    //core::Matlab::GetVariable("X", X, false);

                    //for (auto & uhv : patch.uhs){
                    //    if (!Contains(uh2varStartPosition, uhv.first))
                    //        continue;
                    //    int uhStartPosition = uh2varStartPosition.at(uhv.first);
                    //    for (int i = 0; i < uhv.second.variables.size(); i++){
                    //        uhv.second.variables[i] = X[uhStartPosition + i];
                    //    }
                    //}

                    return true;
                }
                virtual void finalize() {}
            };
        }


        MGPatchDepthsOptimizer::MGPatchDepthsOptimizer(const MixedGraph & mg, MGPatch & patch,
            const std::vector<Vec3> & vanishingPoints, bool useWeights, Algorithm at)
            : _mg(mg), _patch(patch), _vanishingPoints(vanishingPoints), _at(at){

            if (_at == Algorithm::MosekLinearProgramming){
                _internal = new MGPatchDepthsOptimizerInternalMosek;
            }
            else if (_at == Algorithm::MosekLinearProgrammingSimplified) {
                _internal = new MGPatchDepthsOptimizerInternalMosekSimplified;
            }
            else if (_at == Algorithm::EigenSparseQR){
                _internal = new MGPatchDepthsOptimizerInternalEigen;
            }
            else if (_at == Algorithm::EigenSparseQRSimplified) {
                _internal = new MGPatchDepthsOptimizerInternalEigenSparseQRSimplified;
            }
            else if (_at == Algorithm::MATLAB_CVX) {
                _internal = new MGPatchDepthsOptimizerInternalMATLABCVX;
            }
            else if (_at == Algorithm::MATLAB_CVX_v2){
                _internal = new MGPatchDepthsOptimizerInternalMATLABCVX_v2;
            }

            auto internalData = static_cast<MGPatchDepthsOptimizerInternalBase*>(_internal);
            internalData->initialize(mg, patch, vanishingPoints, useWeights);
        }

        MGPatchDepthsOptimizer::~MGPatchDepthsOptimizer(){
            if (_internal){
                auto internalData = static_cast<MGPatchDepthsOptimizerInternalBase*>(_internal);
                internalData->finalize();
                delete internalData;
            }
        }


        bool MGPatchDepthsOptimizer::optimize() {
            return static_cast<MGPatchDepthsOptimizerInternalBase*>(_internal)->optimize(_mg, _patch, _vanishingPoints);
        }



   	}
}