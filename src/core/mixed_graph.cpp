
extern "C" {
    #include <mosek.h>
}
#include <Eigen/Dense>
#include <Eigen/Sparse>

#include "algorithms.hpp"
#include "containers.hpp"
#include "mixed_graph.hpp"

namespace panoramix {
    namespace core {
     
        Plane3 PlaneOfMGUnary(const MGUnary & region, const std::vector<Vec3> & vps, const MGUnaryVariable & var) {
            assert(region.type == MGUnary::Region);
            return Plane3(region.normalizedCenter * var.depthOfCenter, vps[var.claz]);
        }

        Line3 LineOfMGUnary(const MGUnary & line, const std::vector<Vec3> & vps, const MGUnaryVariable & var) {
            assert(line.type == MGUnary::Line);
            InfiniteLine3 infLine(line.normalizedCenter * var.depthOfCenter, vps[var.claz]);
            return Line3(DistanceBetweenTwoLines(InfiniteLine3(Point3(0, 0, 0), line.normalizedCorners.front()), infLine).second.second,
                DistanceBetweenTwoLines(InfiniteLine3(Point3(0, 0, 0), line.normalizedCorners.back()), infLine).second.second);
        }

        double DepthRatioOnMGUnary(const Vec3 & direction, const MGUnary & unary, const std::vector<Vec3> & vps, int claz) {
            if (unary.type == MGUnary::Region){
                const auto & region = unary;
                assert(!region.normalizedCorners.empty());
                Plane3 plane(region.normalizedCenter, vps[claz]);
                return norm(IntersectionOfLineAndPlane(InfiniteLine3(Point3(0, 0, 0), direction), plane).position);
            }
            else if (unary.type == MGUnary::Line){
                const auto & line = unary;
                InfiniteLine3 infLine(line.normalizedCenter, vps[claz]);
                return norm(DistanceBetweenTwoLines(InfiniteLine3(Point3(0, 0, 0), direction), infLine).second.first);
            }
        }

        Point3 LocationOnMGUnary(const Vec3 & direction, const MGUnary & unary, const std::vector<Vec3> & vps,
            const MGUnaryVariable & var){
            if (unary.type == MGUnary::Region){
                const auto & region = unary;
                assert(!region.normalizedCorners.empty());
                Plane3 plane(region.normalizedCenter * var.depthOfCenter, vps[var.claz]);
                return IntersectionOfLineAndPlane(InfiniteLine3(Point3(0, 0, 0), direction), plane).position;
            }
            else if (unary.type == MGUnary::Line){
                const auto & line = unary;
                InfiniteLine3 infLine(line.normalizedCenter * var.depthOfCenter, vps[var.claz]);
                return DistanceBetweenTwoLines(InfiniteLine3(Point3(0, 0, 0), direction), infLine).second.first;
            }
        }




        void InitializeUnaryVarDepths(MGUnaryVarTable & unaryVars, double depth){
            for (auto & uhv : unaryVars){
                uhv.second.depthOfCenter = depth;
            }
        }

        void UpdateBinaryVars(const MixedGraph & mg, const std::vector<Vec3> & vps,
            const MGUnaryVarTable & unaryVars, MGBinaryVarTable & binaryVars){
            for (auto & bv : binaryVars){
                auto & uhs = mg.topo(bv.first).lowers;
                auto & normalizedAnchors = mg.data(bv.first).normalizedAnchors;
                for (int i = 0; i < normalizedAnchors.size(); i++){
                    bv.second.sampleDepthsOnRelatedUnaries.front()[i] = unaryVars.at(uhs.front()).depthOfCenter *
                        DepthRatioOnMGUnary(normalizedAnchors[i], mg.data(uhs.front()), vps, unaryVars.at(uhs.front()).claz);
                    bv.second.sampleDepthsOnRelatedUnaries.back()[i] = unaryVars.at(uhs.back()).depthOfCenter *
                        DepthRatioOnMGUnary(normalizedAnchors[i], mg.data(uhs.back()), vps, unaryVars.at(uhs.back()).claz);
                }
            }
        }


        double FeasibilityOfBinary(const MGBinary & b, int uv1Claz, int uv2Claz,
            const std::vector<Vec3> & vps){

            if (b.normalizedAnchors.size() == 1)
                return 1.0;

            if (b.type == MGBinary::RegionLineConnection){
                return 1.0 - abs(normalize(vps[uv1Claz])
                    .dot(normalize(vps[uv2Claz])));
            }
            if (b.type == MGBinary::RegionRegionOverlapping){
                return abs(normalize(vps[uv1Claz])
                    .dot(normalize(vps[uv2Claz])));
            }
            assert(b.type == MGBinary::RegionRegionConnection);
            if (uv1Claz == uv2Claz)
                return 1.0;

            Vec3 alignDir = normalize(vps[uv1Claz].cross(vps[uv2Claz]));
            double maxDotProd = 0.0;
            for (int k = 1; k < b.normalizedAnchors.size(); k++){
                Vec3 rotateNorm = b.normalizedAnchors[k - 1].cross(b.normalizedAnchors[k]);
                double dotProd = abs(alignDir.dot(rotateNorm));
                if (dotProd > maxDotProd)
                    maxDotProd = dotProd;
            }
            return 1.0 - maxDotProd;

        }


        //double FeasibilityOfBinary(const MixedGraph & mg, const MGBinaryHandle & bh,
        //    const MGUnaryVarTable & unaryVars, const std::vector<Vec3> & vps){
        //    if (mg.data(bh).normalizedAnchors.size() == 1)
        //        return 1.0;
        //    auto & uhs = mg.topo(bh).lowers;
        //    if (mg.data(bh).type == MGBinary::RegionLineConnection){
        //        return 1.0 - abs(normalize(vps[unaryVars.at(uhs.front()).claz])
        //            .dot(normalize(vps[unaryVars.at(uhs.back()).claz])));
        //    }
        //    if (mg.data(bh).type == MGBinary::RegionRegionOverlapping){
        //        return abs(normalize(vps[unaryVars.at(uhs.front()).claz])
        //            .dot(normalize(vps[unaryVars.at(uhs.back()).claz])));
        //    }
        //    assert(mg.data(bh).type == MGBinary::RegionRegionConnection);
        //    if (unaryVars.at(uhs.front()).claz == unaryVars.at(uhs.back()).claz)
        //        return 1.0;
        //    
        //    Vec3 alignDir = normalize(vps[unaryVars.at(uhs.front()).claz].cross(vps[unaryVars.at(uhs.back()).claz]));
        //    double maxDotProd = 0.0;
        //    for (int k = 1; k < mg.data(bh).normalizedAnchors.size(); k++){
        //        Vec3 rotateNorm = mg.data(bh).normalizedAnchors[k-1].cross(mg.data(bh).normalizedAnchors[k]);
        //        double dotProd = abs(alignDir.dot(rotateNorm));
        //        if (dotProd > maxDotProd)
        //            maxDotProd = dotProd;
        //    }
        //    return 1.0 - maxDotProd;
        //}


        MixedGraph BuildMixedGraph(const std::vector<View<PerspectiveCamera>> & views,
            const std::vector<RegionsGraph> & regionsGraphs, const std::vector<LinesGraph> & linesGraphs,
            const ComponentIndexHashMap<std::pair<RegionIndex, RegionIndex>, double> & regionOverlappingsAcrossViews,
            const ComponentIndexHashMap<std::pair<LineIndex, LineIndex>, Vec3> & lineIncidencesAcrossViews,
            const std::vector<std::map<std::pair<RegionHandle, LineHandle>, std::vector<Point2>>> & regionLineConnections,
            const std::vector<Vec3> & vps,
            MGUnaryVarTable & unaryVars, MGBinaryVarTable & binaryVars,
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
                        center
                    });
                    mgh2ri[ri2mgh[ri]] = ri;
                    /* if (ri2mgh[ri].id == 1307){
                    GetVersion();
                    }*/
                    unaryVars[ri2mgh[ri]].claz = initialClaz;
                    unaryVars[ri2mgh[ri]].depthOfCenter = initialDepth;
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
                        normalize(cam.spatialDirection(ld.data.line.component.center()))
                    });
                    mgh2li[li2mgh[li]] = li;
                    unaryVars[li2mgh[li]].claz = ld.data.line.claz;
                    unaryVars[li2mgh[li]].depthOfCenter = initialDepth;
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

            // make overlapped region claz consistent
            static const bool makeOverlappedRegionClassConsistent = false;
            if (makeOverlappedRegionClassConsistent){
                std::vector<MGUnaryHandle> uhs;
                for (auto & v : mg.elements<0>()){
                    if (v.data.type == MGUnary::Region){
                        uhs.push_back(v.topo.hd);
                    }
                }
                std::map<int, std::vector<MGUnaryHandle>> ccUhs;
                core::ConnectedComponents(uhs.begin(), uhs.end(), [&mg](const MGUnaryHandle & uh) {
                    std::vector<MGUnaryHandle> neighbors;
                    for (auto & bh : mg.topo(uh).uppers){
                        if (mg.data(bh).type == MGBinary::RegionRegionOverlapping){
                            auto anotherUh = mg.topo(bh).lowers.front();
                            if (anotherUh == uh)
                                anotherUh = mg.topo(bh).lowers.back();
                            neighbors.push_back(anotherUh);
                        }
                    }
                    return neighbors;
                }, [&ccUhs](const MGUnaryHandle & uh, int ccid){
                    ccUhs[ccid].push_back(uh);
                });
                for (auto & ccUh : ccUhs){
                    auto & uhs = ccUh.second;
                    Vec3 center(0, 0, 0);
                    for (auto uh : uhs){
                        if (uh.id == 1307){
                            GetVersion();
                        }
                        center += (mg.data(uh).normalizedCenter * GetData(mgh2ri[uh], regionsGraphs).area);
                    }
                    center /= norm(center);
                    int initialClaz = std::distance(vps.begin(), std::min_element(vps.begin(), vps.end(),
                        [&center](const Vec3 & vp1, const Vec3 & vp2) -> bool {
                        return AngleBetweenUndirectedVectors(center, vp1) <
                            AngleBetweenUndirectedVectors(center, vp2);
                    }));
                    for (auto uh : uhs){
                        unaryVars[uh].claz = initialClaz;
                    }
                }
            }

            binaryVars.clear();
            for (auto & b : mg.elements<1>()){
                binaryVars[b.topo.hd].sampleDepthsOnRelatedUnaries.front().resize(b.data.normalizedAnchors.size());
                binaryVars[b.topo.hd].sampleDepthsOnRelatedUnaries.back().resize(b.data.normalizedAnchors.size());
            }
            UpdateBinaryVars(mg, vps, unaryVars, binaryVars);

            return mg;
        }



        int MarkConnectedComponentIds(const MixedGraph & mg, std::unordered_map<MGUnaryHandle, int> & ccids) {
            ccids.reserve(mg.internalElements<0>().size());
            std::vector<MGUnaryHandle> uhs;
            uhs.reserve(mg.internalElements<0>().size());
            for (auto & v : mg.elements<0>()){
                uhs.push_back(v.topo.hd);
            }
            return core::ConnectedComponents(uhs.begin(), uhs.end(),
                [&mg](const MGUnaryHandle & uh){
                std::vector<MGUnaryHandle> neighborUhs;
                for (auto & bh : mg.topo(uh).uppers){
                    auto anotherUh = mg.topo(bh).lowers.front();
                    if (anotherUh == uh)
                        anotherUh = mg.topo(bh).lowers.back();
                    neighborUhs.push_back(anotherUh);
                }
                return neighborUhs;
            }, [&ccids](const MGUnaryHandle & uh, int ccid){
                ccids[uh] = ccid;
            });
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
            DepthFirstSearch(uhs.begin(), uhs.end(), [&mg, &patch](const MGUnaryHandle & uh){
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

        double FeasibilityOfPatch(const MixedGraph & mg, const MGPatch & patch, const std::vector<Vec3> & vps){
            if (!BinaryHandlesAreValidInPatch(mg, patch) || !UnariesAreConnectedInPatch(mg, patch))
                return 0.0;
            double minConsistency = 1.0;
            for (auto & bhv : patch.bhs){
                double c = FeasibilityOfBinary(mg, bhv.first, patch.uhs, vps);
                if (c < minConsistency)
                    c = minConsistency;
            }
            return minConsistency;
        }

        double AnchorDistanceSumOnBinaryOfPatch(const MGBinaryHandle & bh, const MGPatch & patch){
            assert(Contains(patch, bh));
            auto & sampleDepths = patch.bhs.at(bh).sampleDepthsOnRelatedUnaries;
            double distanceSum = 0.0;
            for (int i = 0; i < sampleDepths[0].size(); i++){
                distanceSum += abs(sampleDepths[0][i] - sampleDepths[1][i]);
            }
            return distanceSum;
        }

        double AnchorDistanceSumOnBinary(const MixedGraph & mg, const MGBinaryHandle & bh, const MGPatch & patch, const std::vector<Vec3> & vps){
            auto uhs = mg.topo(bh).lowers;
            assert(Contains(patch, uhs.front()) && Contains(patch, uhs.back()));
            auto & nanchors = mg.data(bh).normalizedAnchors;
            double distanceSum = 0.0;
            for (auto & na : nanchors){
                double d1 = patch.uhs.at(uhs.front()).depthOfCenter *
                    DepthRatioOnMGUnary(na, mg.data(uhs.front()), vps, patch.uhs.at(uhs.front()).claz);
                double d2 = patch.uhs.at(uhs.back()).depthOfCenter *
                    DepthRatioOnMGUnary(na, mg.data(uhs.back()), vps, patch.uhs.at(uhs.back()).claz);
                distanceSum += abs(d1 - d2);
            }
            return distanceSum;
        }

        double BinaryDistanceOfPatch(const MGBinaryHandle & bh, const MGPatch & patch){
            return AnchorDistanceSumOnBinaryOfPatch(bh, patch) / patch.bhs.at(bh).sampleDepthsOnRelatedUnaries[0].size();
        }

        double AverageBinaryDistanceOfPatch(const MGPatch & patch, int pow){
            double distanceSum = 0.0;
            int distanceCount = 0;
            for (auto & bhv : patch.bhs){
                for (int i = 0; i < bhv.second.sampleDepthsOnRelatedUnaries.front().size(); i++){
                    distanceSum += std::pow(abs(bhv.second.sampleDepthsOnRelatedUnaries.front()[i] -
                        bhv.second.sampleDepthsOnRelatedUnaries.back()[i]), pow);
                }
                distanceCount += bhv.second.sampleDepthsOnRelatedUnaries.front().size();
            }
            return distanceSum / distanceCount;
        }

        double AverageDepthOfPatch(const MGPatch & patch){
            double depthSum = 0.0;
            for (auto & uhv : patch.uhs){
                depthSum += uhv.second.depthOfCenter;
            }
            return depthSum / patch.uhs.size();
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


        bool IsTreePatch(const MixedGraph & mg, const MGPatch & patch){
            if (!(BinaryHandlesAreValidInPatch(mg, patch) && UnariesAreConnectedInPatch(mg, patch)))
                return false;
            return patch.uhs.size() == patch.bhs.size() + 1;
        }

        void CompletePatch(const MixedGraph & mg, MGPatch & patch, const MGBinaryVarTable & binaryVars){
            for (auto & uhv : patch.uhs){
                auto & uh = uhv.first;
                for (auto & bh : mg.topo(uh).uppers){
                    if (Contains(patch.uhs, mg.topo(bh).lowers.front()) &&
                        Contains(patch.uhs, mg.topo(bh).lowers.back()) &&
                        !Contains(patch.bhs, bh)){
                        patch.bhs[bh] = binaryVars.at(bh);
                    }
                }
            }
        }


        void ScalePatch(MGPatch & patch, double scale){
            for (auto & uhv : patch.uhs)
                uhv.second.depthOfCenter *= scale;
            for (auto & bhv : patch.bhs){
                for (auto & arr : bhv.second.sampleDepthsOnRelatedUnaries){
                    for (double & d : arr){
                        d *= scale;
                    }
                }
            }
        }






        void CommitPatchToVariableTable(const MGPatch & patch,
            MGUnaryVarTable & unaryVars, MGBinaryVarTable & binaryVars){
            for (auto & uhv : patch.uhs){
                unaryVars[uhv.first] = uhv.second;
            }
            for (auto & bhv : patch.bhs){
                binaryVars[bhv.first] = bhv.second;
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



        static const int g_NumOfAnchorsUsedForEachBinary = 1;

        struct MGPatchDepthsOptimizerInternalBase {
            virtual void initialize(const MixedGraph & mg, MGPatch & patch,
            const std::vector<Vec3> & vanishingPoints, bool useWeights) = 0;
            virtual void setDepthBounds(MGPatch & patch, double depthLb, double depthUb) = 0;
            virtual void setDepthsAllGreaterThan(MGPatch & patch, double lob) = 0;
            virtual void setUnaryClass(const MixedGraph & mg, MGPatch & patch,
                const std::vector<Vec3> & vanishingPoints,
                const MGUnaryHandle & uh, int claz) = 0;
            virtual bool optimize(const MixedGraph & mg, MGPatch & patch,
                const std::vector<Vec3> & vanishingPoints) = 0;
            virtual void finalize() = 0;
        };



        static struct MosekGlobal {
            MSKenv_t env;
            inline MosekGlobal() { MSK_makeenv(&env, nullptr); }
            inline ~MosekGlobal() { MSK_deleteenv(&env); }
        } mosekGlobal;

        static const bool g_DepthBoundAsConstraint = false;
        static const bool g_BoundSlackPositive = false;

        static void MSKAPI printstr(void *handle,
            MSKCONST char str[]){
            printf("%s", str);
        }


        struct MGPatchDepthsOptimizerInternalMosek : MGPatchDepthsOptimizerInternalBase {
            MSKtask_t task;

            std::unordered_map<MGUnaryHandle, int> uhPositions;
            std::unordered_map<MGBinaryHandle, std::pair<std::vector<Vec3>, int>> bhAnchors;

            int varNum;
            int consNum;

            virtual void initialize(const MixedGraph & mg, MGPatch & patch,
                const std::vector<Vec3> & vanishingPoints, bool useWeights) override;
            virtual void setDepthBounds(MGPatch & patch, double depthLb, double depthUb) override;
            virtual void setDepthsAllGreaterThan(MGPatch & patch, double lob) override;
            virtual void setUnaryClass(const MixedGraph & mg, MGPatch & patch,
                const std::vector<Vec3> & vanishingPoints,
                const MGUnaryHandle & uh, int claz) override;
            virtual bool optimize(const MixedGraph & mg, MGPatch & patch,
                const std::vector<Vec3> & vanishingPoints) override;
            virtual void finalize() override;
        };

        void MGPatchDepthsOptimizerInternalMosek::initialize(const MixedGraph & mg, MGPatch & patch,
            const std::vector<Vec3> & vanishingPoints, bool useWeights){
            auto internalData = this;

            auto & task = internalData->task;

            //auto tick = Tick("preparing optimization");

            assert(BinaryHandlesAreValidInPatch(mg, patch));
            int depthNum = patch.uhs.size();

            assert(UnariesAreConnectedInPatch(mg, patch));

            // temp data  
            auto & uhPositions = internalData->uhPositions;
            auto & bhAnchors = internalData->bhAnchors;

            uhPositions.reserve(patch.uhs.size());
            int unaryNum = 0;
            for (auto & uhv : patch.uhs){
                uhPositions[uhv.first] = unaryNum;
                unaryNum++;
            }

            bhAnchors.reserve(patch.bhs.size());

            // count sample connections
            int samplesNum = 0;
            for (auto & bhv : patch.bhs){
                bhAnchors[bhv.first].first = SelectSubset(mg.data(bhv.first).normalizedAnchors, g_NumOfAnchorsUsedForEachBinary);
                bhAnchors[bhv.first].second = samplesNum;
                samplesNum += bhAnchors[bhv.first].first.size();
            }

            int slackVarNum = samplesNum;
            auto & varNum = internalData->varNum;
            varNum = depthNum + slackVarNum; // depths, slackVars
            auto & consNum = internalData->consNum;
            consNum = samplesNum * 2;   // depth1 * ratio1 - depth2 * ratio2 < slackVar;  
            // depth2 * ratio2 - depth1 * ratio1 < slackVar

            //make as bounds:   depth = defaultValue
            //                  depth \in [depthLb, depthUb]
            if (g_DepthBoundAsConstraint){
                consNum += depthNum;
            }

            //env = nullptr;
            task = nullptr;
            auto & env = mosekGlobal.env;

            MSK_maketask(env, consNum, varNum, &task);
            //MSK_linkfunctotaskstream(task, MSK_STREAM_LOG, NULL, printstr);

            MSK_appendcons(task, consNum);
            MSK_appendvars(task, varNum);

            // set weights
            int slackVarId = 0;
            for (auto & bhv : patch.bhs){
                auto & bd = mg.data(bhv.first);
                assert(bd.weight >= 0.0);
                for (int k = 0; k < bhAnchors[bhv.first].first.size(); k++){
                    MSK_putcj(task, slackVarId, useWeights ? bd.weight : 1.0);
                    slackVarId++;
                }
            }

            // bounds for vars
            {
                int varId = 0;
                if (!g_DepthBoundAsConstraint){
                    for (; varId < depthNum; varId++){ // for depths
                        MSK_putvarbound(task, varId, MSK_BK_LO, 1.0, +MSK_INFINITY);
                    }
                }
                if (g_BoundSlackPositive){
                    for (; varId < depthNum + slackVarNum; varId++){ // for slack vars 
                        MSK_putvarbound(task, varId, MSK_BK_LO, 0.0, +MSK_INFINITY);
                    }
                }
            }

            // fill constraints related to each vars
            {
                int varId = 0;
                for (auto & uhv : patch.uhs){ // depths
                    auto & uh = uhv.first;

                    auto & relatedBhs = mg.topo(uh).uppers;
                    int relatedAnchorsNum = 0;
                    for (auto & bh : relatedBhs){
                        if (!Contains(patch.bhs, bh))
                            continue;
                        relatedAnchorsNum += bhAnchors[bh].first.size();
                    }

                    int relatedConsNum = relatedAnchorsNum * 2;
                    if (g_DepthBoundAsConstraint){
                        relatedConsNum += 1;
                    }

                    std::vector<MSKint32t> consIds;
                    std::vector<MSKrealt> consValues;

                    consIds.reserve(relatedConsNum);
                    consValues.reserve(relatedConsNum);

                    for (auto & bh : relatedBhs){
                        if (!Contains(patch.bhs, bh))
                            continue;

                        int firstAnchorPosition = bhAnchors[bh].second;
                        auto & samples = bhAnchors[bh].first;
                        for (int k = 0; k < samples.size(); k++){
                            consIds.push_back((firstAnchorPosition + k) * 2); // one for [depth1 * ratio1 - depth2 * ratio2 - slackVar < 0]; 
                            consIds.push_back((firstAnchorPosition + k) * 2 + 1); // another for [- depth1 * ratio1 + depth2 * ratio2 - slackVar < 0];

                            auto & a = samples[k];
                            double ratio = DepthRatioOnMGUnary(a, mg.data(uh), vanishingPoints, patch.uhs.at(uh).claz);
                            bool isOnLeftSide = uh == mg.topo(bh).lowers.front();
                            if (isOnLeftSide){ // as depth1
                                consValues.push_back(ratio);
                                consValues.push_back(-ratio);
                            }
                            else{
                                consValues.push_back(-ratio);
                                consValues.push_back(ratio);
                            }
                        }
                    }

                    if (g_DepthBoundAsConstraint){ // add depth bound
                        consIds.push_back(samplesNum * 2 + varId);
                        consValues.push_back(1.0); // (1.0) * depth > 1.0;
                    }

                    assert(varId == uhPositions[uh]);
                    MSK_putacol(task, uhPositions[uh], relatedConsNum, consIds.data(), consValues.data());
                    varId++;
                }

                for (; varId < depthNum + slackVarNum; varId++){ // slack vars
                    MSKint32t consIds[] = {
                        (varId - depthNum) * 2, // one for [depth1 * ratio1 - depth2 * ratio2 - slackVar < 0]; 
                        (varId - depthNum) * 2 + 1  // another for [- depth1 * ratio1 + depth2 * ratio2 - slackVar < 0];
                    };
                    MSKrealt consValues[] = { -1.0, -1.0 };
                    MSK_putacol(task, varId, 2, consIds, consValues);
                }
            }

            // bounds for constraints
            int consId = 0;
            for (; consId < samplesNum * 2; consId++){
                // all [depth1 * ratio1 - depth2 * ratio2 - slackVar < 0];
                MSK_putconbound(task, consId, MSK_BK_UP, -MSK_INFINITY, 0.0);
            }
            if (g_DepthBoundAsConstraint){
                for (; consId < consNum; consId++){
                    // depth > 1.0
                    MSK_putconbound(task, consId, MSK_BK_LO, 1.0, +MSK_INFINITY);
                }
            }
            else{
                assert(consId == consNum);
            }


        }

        void MGPatchDepthsOptimizerInternalMosek::setDepthBounds(MGPatch & patch, double depthLb, double depthUb) {
            for (int varId = 0; varId < patch.uhs.size(); varId++){ // for depths
                MSK_putvarbound(task,
                    varId, MSK_BK_RA, depthLb, depthUb);
            }
        }

        void MGPatchDepthsOptimizerInternalMosek::setDepthsAllGreaterThan(MGPatch & patch, double lob) {
            for (int varId = 0; varId < patch.uhs.size(); varId++){ // for depths
                MSK_putvarbound(task,
                    varId, MSK_BK_LO, lob, +MSK_INFINITY);
            }
        }

        void MGPatchDepthsOptimizerInternalMosek::setUnaryClass(const MixedGraph & mg, MGPatch & patch,
            const std::vector<Vec3> & vanishingPoints,
            const MGUnaryHandle & uh, int claz) {

            auto & relatedBhs = mg.topo(uh).uppers;
            int relatedAnchorsNum = 0;
            for (auto & bh : relatedBhs){
                if (!Contains(patch.bhs, bh))
                    continue;
                relatedAnchorsNum += bhAnchors[bh].first.size();
            }
            int relatedConsNum = relatedAnchorsNum * 2;

            std::vector<MSKint32t> consIds;
            std::vector<MSKrealt> consValues;

            consIds.reserve(relatedConsNum);
            consValues.reserve(relatedConsNum);

            for (auto & bh : relatedBhs){
                if (!Contains(patch.bhs, bh))
                    continue;
                int firstAnchorPosition = bhAnchors.at(bh).second;
                auto & samples = bhAnchors.at(bh).first;
                for (int k = 0; k < samples.size(); k++){
                    consIds.push_back((firstAnchorPosition + k) * 2); // one for [depth1 * ratio1 - depth2 * ratio2 - slackVar < 0]; 
                    consIds.push_back((firstAnchorPosition + k) * 2 + 1); // another for [- depth1 * ratio1 + depth2 * ratio2 - slackVar < 0];

                    auto & a = samples[k];
                    double ratio = DepthRatioOnMGUnary(a, mg.data(uh), vanishingPoints, patch.uhs.at(uh).claz);
                    bool isOnLeftSide = uh == mg.topo(bh).lowers.front();
                    if (isOnLeftSide){ // as depth1
                        consValues.push_back(ratio);
                        consValues.push_back(-ratio);
                    }
                    else{
                        consValues.push_back(-ratio);
                        consValues.push_back(ratio);
                    }
                }
            }

            MSK_putacol(task, uhPositions.at(uh), relatedConsNum, consIds.data(), consValues.data());
        }

        bool MGPatchDepthsOptimizerInternalMosek::optimize(const MixedGraph & mg, MGPatch & patch,
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
                                                 uhv.second.depthOfCenter = xx[uhPositions.at(uhv.first)];
                                             }
                                             delete[] xx;
                                             UpdateBinaryVars(mg, vanishingPoints, patch.uhs, patch.bhs);
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

        void MGPatchDepthsOptimizerInternalMosek::finalize(){
            MSK_deletetask(&task);
        }




        struct MGPatchDepthsOptimizerInternalEigen : MGPatchDepthsOptimizerInternalBase {
            Eigen::SparseMatrix<double> A, W;
            Eigen::VectorXd B;
            bool useWeights;
            std::unordered_map<MGUnaryHandle, int> uhPositions;
            std::unordered_map<MGBinaryHandle, std::vector<Vec3>> bhAnchors;

            virtual void initialize(const MixedGraph & mg, MGPatch & patch,
                const std::vector<Vec3> & vanishingPoints, bool useWeights) override;
            virtual void setDepthBounds(MGPatch & patch, double depthLb, double depthUb) override;
            virtual void setDepthsAllGreaterThan(MGPatch & patch, double lob) override;
            virtual void setUnaryClass(const MixedGraph & mg, MGPatch & patch,
                const std::vector<Vec3> & vanishingPoints,
                const MGUnaryHandle & uh, int claz) override;
            virtual bool optimize(const MixedGraph & mg, MGPatch & patch,
                const std::vector<Vec3> & vanishingPoints) override;
            virtual void finalize() override;
        };

        void MGPatchDepthsOptimizerInternalEigen::initialize(const MixedGraph & mg, MGPatch & patch,
            const std::vector<Vec3> & vanishingPoints, bool useWeights) {

            this->useWeights = useWeights;

            assert(BinaryHandlesAreValidInPatch(mg, patch));
            assert(UnariesAreConnectedInPatch(mg, patch));

            int uhid = 0;
            for (auto & uhv : patch.uhs){
                uhPositions[uhv.first] = uhid++;
            }

            int varNum = patch.uhs.size();
            int consNum = 0;

            consNum++;
            for (auto & bhv : patch.bhs){
                bhAnchors[bhv.first] = SelectSubset(mg.data(bhv.first).normalizedAnchors, g_NumOfAnchorsUsedForEachBinary);
                consNum += bhAnchors.at(bhv.first).size();
            }

            A.resize(consNum, varNum);
            W.resize(consNum, consNum);
            B.resize(consNum);

            // write equations
            int eid = 0;
            A.insert(eid, 0) = 1.0;
            B(eid) = 1.0;
            W.insert(eid, eid) = 1.0;
            eid++;

            for (auto & bhv : patch.bhs){
                auto & bh = bhv.first;
                for (auto & a : bhAnchors.at(bh)){
                    auto uh1 = mg.topo(bh).lowers.front();
                    auto uh2 = mg.topo(bh).lowers.back();
                    auto & u1 = mg.data(uh1);
                    auto & u2 = mg.data(uh2);
                    double ratio1 = DepthRatioOnMGUnary(a, u1, vanishingPoints, patch.uhs.at(uh1).claz);
                    double ratio2 = DepthRatioOnMGUnary(a, u2, vanishingPoints, patch.uhs.at(uh2).claz);
                    A.insert(eid, uhPositions.at(uh1)) = ratio1;
                    A.insert(eid, uhPositions.at(uh2)) = -ratio2;
                    B(eid) = 0.0;
                    assert(mg.data(bh).weight >= 0.0);
                    W.insert(eid, eid) = mg.data(bh).weight;
                    eid++;
                }
            }

            assert(eid == consNum);
        }

        void MGPatchDepthsOptimizerInternalEigen::setDepthBounds(MGPatch & patch, double depthLb, double depthUb) {
            NOT_IMPLEMENTED_YET();
        }

        void MGPatchDepthsOptimizerInternalEigen::setDepthsAllGreaterThan(MGPatch & patch, double lob) {
            NOT_IMPLEMENTED_YET();
        }

        void MGPatchDepthsOptimizerInternalEigen::setUnaryClass(const MixedGraph & mg, MGPatch & patch,
            const std::vector<Vec3> & vanishingPoints,
            const MGUnaryHandle & uh, int claz) {
            int eid = 1;
            for (auto & bhv : patch.bhs){
                auto & bh = bhv.first;
                auto uh1 = mg.topo(bh).lowers.front();
                auto uh2 = mg.topo(bh).lowers.back();
                if (uh1 == uh || uh2 == uh){
                    for (auto & a : bhAnchors.at(bh)){
                        auto & u1 = mg.data(uh1);
                        auto & u2 = mg.data(uh2);
                        double ratio1 = DepthRatioOnMGUnary(a, u1, vanishingPoints, claz);
                        double ratio2 = DepthRatioOnMGUnary(a, u2, vanishingPoints, claz);
                        A.insert(eid, uhPositions.at(uh1)) = ratio1;
                        A.insert(eid, uhPositions.at(uh2)) = -ratio2;
                        B(eid) = 0.0;
                        W.insert(eid, eid) = mg.data(bh).weight;
                    }
                }
                eid += bhAnchors.at(bh).size();
            }
        }

        bool MGPatchDepthsOptimizerInternalEigen::optimize(const MixedGraph & mg, MGPatch & patch,
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
                uhv.second.depthOfCenter = X(uhPositions.at(uhv.first));
            }

            UpdateBinaryVars(mg, vanishingPoints, patch.uhs, patch.bhs);

            return true;

        }

        void MGPatchDepthsOptimizerInternalEigen::finalize() {}






        MGPatchDepthsOptimizer::MGPatchDepthsOptimizer(const MixedGraph & mg, MGPatch & patch,
            const std::vector<Vec3> & vanishingPoints, bool useWeights, AlgorithmType at)
            : _mg(mg), _patch(patch), _vanishingPoints(vanishingPoints), _at(at){

            if (_at == AlgorithmType::MosekLinearProgramming){
                _internal = new MGPatchDepthsOptimizerInternalMosek;
            }
            else if (_at == AlgorithmType::EigenSparseQR){
                _internal = new MGPatchDepthsOptimizerInternalEigen;
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

        void MGPatchDepthsOptimizer::setDepthBounds(double depthLb, double depthUb){
            static_cast<MGPatchDepthsOptimizerInternalBase*>(_internal)->setDepthBounds(_patch, depthLb, depthUb);
        }

        void MGPatchDepthsOptimizer::setDepthsAllGreaterThan(double lob){
            static_cast<MGPatchDepthsOptimizerInternalBase*>(_internal)->setDepthsAllGreaterThan(_patch, lob);
        }

        void MGPatchDepthsOptimizer::setUnaryClass(const MGUnaryHandle & uh, int claz){
            static_cast<MGPatchDepthsOptimizerInternalBase*>(_internal)->setUnaryClass(_mg, _patch, _vanishingPoints, uh, claz);
        }

        bool MGPatchDepthsOptimizer::optimize() {
            return static_cast<MGPatchDepthsOptimizerInternalBase*>(_internal)->optimize(_mg, _patch, _vanishingPoints);
        }








        void AddUnaryAndRelatedBinariesToPatch(const MixedGraph & mg, MGPatch & patch, const MGUnaryHandle & uh,
            MGUnaryVarTable & unaryVars, MGBinaryVarTable & binaryVars){
            if (Contains(patch, uh))
                return;
            patch.uhs[uh] = unaryVars.at(uh);
            int relatedBhCount = 0;
            for (auto & bh : mg.topo(uh).uppers){
                auto anotherUh = mg.topo(bh).lowers.front();
                if (anotherUh == uh)
                    anotherUh = mg.topo(bh).lowers.back();
                if (Contains(patch, anotherUh)){
                    relatedBhCount++;
                    patch.bhs[bh] = binaryVars.at(bh);
                }
            }
            assert(relatedBhCount > 0);
        }

        // returns a confidence for each newly added unary
        // the confidence should consider: 
        //  //1. the type of the unary : region/line [this can be expressed in 4.]
        //  2. the center direction of this unary in space (consistent with vps)
        //  3. the consistency(as feasibility of binaries) with current patch unaries
        //  4. the importance ratio of feasible binaries

        namespace  {

            struct UOption {
                MGUnaryHandle uh;
                int claz;
            };

            inline bool operator == (const UOption & a, const UOption & b){
                return a.uh == b.uh && a.claz == b.claz;
            }
            struct UOptionHasher {
                inline uint64_t operator()(const UOption & o) const {
                    return (o.uh.id << 2) + o.claz;
                }
            };
            using UOptionQueue = MaxHeap<UOption, double, std::unordered_map<UOption, int, UOptionHasher>>;

            struct UDependency{
                std::unordered_map<MGBinaryHandle, double> bhs;
                double precomputedDepth;
            };

            std::vector<Scored<UDependency>> DependenciesOfOption(const MixedGraph & mg,
                const MGPatch & patch, const UOption & option, const std::vector<Vec3> & vps,
                double depthThreshold = 0.02){
                assert(!Contains(patch, option.uh));

                std::vector<Scored<UDependency>> dependencies;
                for (auto & bh : mg.topo(option.uh).uppers){
                    auto insider = mg.topo(bh).lowers[0];
                    if (insider == option.uh)
                        insider = mg.topo(bh).lowers[1];
                    if (!Contains(patch, insider))
                        continue;

                    auto & anchor = mg.data(bh).normalizedAnchors.front();
                    double depth = DepthRatioOnMGUnary(anchor, mg.data(insider), vps, patch.uhs.at(insider).claz)
                        * patch.uhs.at(insider).depthOfCenter
                        / DepthRatioOnMGUnary(anchor, mg.data(option.uh), vps, option.claz);

                    bool hasDuplicates = false;
                    double feas = FeasibilityOfBinary(mg.data(bh), patch.uhs.at(insider).claz, option.claz, vps);
                    feas = 0.5 - Pitfall(1.0 - feas, 0.2); // -0.5 ~ +0.5

                    double scoreOnThisBh = mg.data(bh).importanceRatioInRelatedUnaries[option.uh == mg.topo(bh).lowers[0] ? 0 : 1]
                        * feas;
                    assert(IsBetween(scoreOnThisBh, -0.6, +0.6));
                    for (auto & ds : dependencies){
                        if (abs(ds.component.precomputedDepth - depth) < depthThreshold){
                            ds.score += scoreOnThisBh;
                            ds.component.bhs[bh] = scoreOnThisBh;
                            hasDuplicates = true;
                        }
                    }
                    if (!hasDuplicates){
                        dependencies.push_back(ScoreAs(UDependency{ { std::make_pair(bh, scoreOnThisBh) }, depth }, scoreOnThisBh));
                    }
                }

                return dependencies;
            }

            void UpdateOptionsOfAnOutsiderUnary(const MixedGraph & mg, const MGUnaryHandle & outsider,
                UOptionQueue & Q,
                std::unordered_map<UOption, UDependency, UOptionHasher> & optionDependencies,
                const MGPatch & patch, const MGUnaryVarTable & unaryVars, const std::vector<Vec3> & vps){

                if (Contains(patch, outsider))
                    return;

                // candidate options for this uh
                std::vector<Scored<UOption>> outOptionsWithVPScores;
                if (mg.data(outsider).type == MGUnary::Region){ // if is a Region
                    outOptionsWithVPScores.resize(vps.size());
                    const Vec3 & ncenter = mg.data(outsider).normalizedCenter;
                    for (int i = 0; i < vps.size(); i++){
                        double maxAngle = 0.0;
                        for (auto & ncorner : mg.data(outsider).normalizedCorners){
                            double angle = AngleBetweenUndirectedVectors(ncorner, vps[i]);
                            if (angle > maxAngle){
                                maxAngle = angle;
                            }
                        }
                        double score = 1.0;
                        if (maxAngle >= M_PI_2 * 0.8){
                            score *= Square((maxAngle - M_PI_2) / (M_PI_2 * 0.2));
                        }
                        outOptionsWithVPScores[i] = ScoreAs(UOption{ outsider, i }, score);
                    }
                }
                else if (mg.data(outsider).type == MGUnary::Line){ // if is a Line
                    outOptionsWithVPScores = { ScoreAs(UOption{ outsider, unaryVars.at(outsider).claz }, 0.1) };
                }

                // find best dependencies with current patch
                for (Scored<UOption> & outOptionWithVPScore : outOptionsWithVPScores){
                    auto & outOption = outOptionWithVPScore.component;
                    // get the depth with the highest score
                    auto dependencies = DependenciesOfOption(mg, patch, outOption, vps, 0.02);

                    ///// [[commented out the Entropy calculation]]
                    // get entropy of all choice scores
                    //std::vector<double> scores;
                    //scores.reserve(dependencies.size());
                    //double scoresSum = 0.0;
                    //for (auto & dep : dependencies){
                    //    if (dep.score >= 0.0){
                    //        scores.push_back(dep.score);
                    //        scoresSum += scores.back();
                    //    }
                    //}
                    //for (auto & score : scores){
                    //    score /= scoresSum;
                    //}
                    //double scoreEntropy = 1.0 - exp(-EntropyOfContainer(scores, 1.0)); // 0 ~ 1
                    //assert(IsBetween(scoreEntropy, -1e-4, 1.0 + 1e-4));

                    auto maxDSIter = std::max_element(dependencies.begin(), dependencies.end());
                    double maxScore = maxDSIter->score;
                    double bhsSize = maxDSIter->component.bhs.size();
                    assert(IsBetween(maxScore,
                        -bhsSize,
                        bhsSize));

                    // register the depth
                    optionDependencies[outOption] = std::move(maxDSIter->component);
                    // update the score
                    double outOptionScore = (maxDSIter->score/* - scoreEntropy*/) * outOptionWithVPScore.score;
                    if (std::isnan(outOptionScore)){
                        outOptionScore = std::numeric_limits<double>::lowest();
                        std::cerr << "nan option score!" << std::endl;
                    }

                    if (!Contains(Q, outOption)){
                        Q.push(outOption, outOptionScore);
                    }
                    else{
                        Q.setScore(outOption, outOptionScore);
                    }
                }

            }

        }

        void ExtandPatch(const MixedGraph & mg, MGPatch & patch,
            const MGUnaryVarTable & unaryVars, const MGBinaryVarTable & binaryVars,
            const std::vector<Vec3> & vps,
            HandledTable<MGUnaryHandle, double> & uhScores,
            HandledTable<MGBinaryHandle, double> & bhScores,
            std::vector<MGUnaryHandle> & uhsOrder,
            double scoreThreshold,
            const std::function<bool(MGUnaryHandle)> & uhIsValid){

            std::unordered_map<UOption, UDependency, UOptionHasher> optionDependencies;

            UOptionQueue Q; // all adjacent uh options
            for (auto & uhv : patch.uhs){
                auto uh = uhv.first;
                for (auto & bh : mg.topo(uh).uppers){
                    auto outsider = mg.topo(bh).lowers[0];
                    if (outsider == uh)
                        outsider = mg.topo(bh).lowers[1];
                    if (Contains(patch, outsider))
                        continue;
                    if (!uhIsValid(outsider))
                        continue;
                    UpdateOptionsOfAnOutsiderUnary(mg, outsider, Q, optionDependencies, patch, unaryVars, vps);
                }
            }

            std::unordered_set<UOption, UOptionHasher> blockedOptions;

            while (true){

                //std::cout << "Q score " << Q.topScore() << "  Q size " << Q.size() << std::endl;
                if (!(!Q.empty() && Q.topScore() >= scoreThreshold))
                    break;

                // pop best option
                UOption o = Q.top();
                double score = Q.topScore();
                Q.pop();

                if (Contains(patch, o.uh))
                    continue;

                // install option to patch
                patch.uhs[o.uh] = MGUnaryVariable{ o.claz, optionDependencies.at(o).precomputedDepth };
                for (auto & bhs : optionDependencies.at(o).bhs){
                    assert(!Contains(patch, bhs.first));
                    patch.bhs[bhs.first] = binaryVars.at(bhs.first);
                    bhScores[bhs.first] = bhs.second;
                }
                uhScores[o.uh] = score;
                uhsOrder.push_back(o.uh);

                //std::cout << patch.uhs.size() << "-th uh installed" << std::endl;

                patch.updateBinaryVars(mg, vps);
                //MGPatchDepthsOptimizer(mg, patch, vps, false, MGPatchDepthsOptimizer::EigenSparseQR).optimize();                

                // update neighbor uhs
                for (auto & bh : mg.topo(o.uh).uppers){
                    auto outsider = mg.topo(bh).lowers[0];
                    if (outsider == o.uh)
                        outsider = mg.topo(bh).lowers[1];
                    if (Contains(patch, outsider))
                        continue;
                    if (!uhIsValid(outsider))
                        continue;
                    UpdateOptionsOfAnOutsiderUnary(mg, outsider, Q, optionDependencies, patch, unaryVars, vps);
                }

            }

            patch.updateBinaryVars(mg, vps);

        }
                
    }
}
 