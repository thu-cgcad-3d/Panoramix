
#include "../core/containers.hpp"
#include "../ml/factor_graph.hpp"

#include "rl_graph_occlusion.hpp"

namespace pano {
    namespace experimental {

        using namespace core;



        namespace {

            void ShrinkLine(Line3 & line, double ratio){
                double angle = AngleBetweenDirections(line.first, line.second);
                double dangle = angle * ratio;
                line = {
                    RotateDirection(line.first, line.second, dangle),
                    RotateDirection(line.second, line.first, dangle)
                };
            }

        }


        HandledTable<LineRelationHandle, DepthRelationGuess> GuessLineDepthRelation(const RLGraph & mg, double maxDistance){

            auto dr = mg.createConstraintTable<LineRelationData>(DepthRelationGuess::Unknown);

            // find T-junctions in line relations
            for (auto & lr : mg.constraints<LineRelationData>()){
                auto & lh1 = lr.topo.component<0>();
                auto & lh2 = lr.topo.component<1>();
                auto & ld1 = mg.data(lh1);
                auto & ld2 = mg.data(lh2);
                if (ld1.initialClaz != -1 && ld1.initialClaz == ld2.initialClaz){
                    // incidence // connected?
                    // todo
                    dr[lr.topo.hd] = DepthRelationGuess::MaybeConnected;
                    continue;
                }
                if (ld1.initialClaz == -1 && ld2.initialClaz == -1){
                    continue;
                }
                if (DistanceBetweenTwoLines(normalize(ld1.line), normalize(ld2.line)).first > maxDistance)
                    continue;

                Line3 sline1 = ld1.line; ShrinkLine(sline1, 0.3);
                Line3 sline2 = ld2.line; ShrinkLine(sline2, 0.3);

                Vec3 eq1 = ld1.line.first.cross(ld1.line.second);
                Vec3 eq2 = ld2.line.first.cross(ld2.line.second);
                Vec3 inter = normalize(eq1.cross(eq2));
                bool interInLine1 = (inter.cross(ld1.line.first)).dot(ld1.line.second.cross(inter)) >= 0;
                bool interInLine2 = (inter.cross(ld2.line.first)).dot(ld2.line.second.cross(inter)) >= 0;
                if (interInLine1 && interInLine2){
                    // todo
                    dr[lr.topo.hd] = DepthRelationGuess::MaybeConnected;
                    continue;
                }
                if (interInLine1 && !interInLine2){
                    // |-
                    // todo
                    dr[lr.topo.hd] = DepthRelationGuess::FirstMaybeCloser;
                    continue;
                }
                if (interInLine2 && !interInLine1){
                    // -|
                    // todo
                    dr[lr.topo.hd] = DepthRelationGuess::SecondMaybeCloser;
                    continue;
                }
            }

            return dr;

        }



        HandledTable<LineHandle, LineDetachStatus> GuessLineDetachStatus(
            const RLGraph & mg, const HandledTable<LineRelationHandle, DepthRelationGuess> & ldr){

            auto result = mg.createComponentTable<LineData>(LineDetachStatus{ false, false });
            for (auto & l : mg.components<LineData>()){
                if (l.data.initialClaz == -1){
                    continue;
                }
                Vec3 rightdir = normalize(l.data.line.first.cross(l.data.line.second));

                int rightOccCandNum = 0;
                int rightConNum = 0;
                int rightAll = 0;
                int leftOccCandNum = 0;
                int leftConNum = 0;
                int leftAll = 0;

                for (auto & lr : l.topo.constraints<LineRelationData>()){
                    const DepthRelationGuess & dr = ldr[lr];
                    auto another = mg.topo(lr).component<0>();
                    if (another == l.topo.hd){
                        another = mg.topo(lr).component<1>();
                    }
                    auto & anotherld = mg.data(another);                    

                    bool isOnRight = (anotherld.line.center() - l.data.line.center()).dot(rightdir) > 0;
                    if (mg.topo(lr).component<0>() == l.topo.hd && dr == DepthRelationGuess::FirstMaybeCloser ||
                        mg.topo(lr).component<1>() == l.topo.hd && dr == DepthRelationGuess::SecondMaybeCloser){
                        if (isOnRight){
                            rightOccCandNum++;
                        }
                        else{
                            leftOccCandNum++;
                        }
                    }

                    if (isOnRight){
                        rightAll++;
                    }
                    else{
                        leftAll++;
                    }

                    if (dr == DepthRelationGuess::MaybeConnected && mg.data(lr).type == LineRelationData::Intersection){
                        if (isOnRight){
                            rightConNum++;
                        }
                        else{
                            leftConNum++;
                        }

                    }
                }

                for (auto & rl : l.topo.constraints<RegionLineConnectionData>()){
                    bool isOnRight = (mg.data(mg.topo(rl).component<0>()).normalizedCenter - l.data.line.center()).dot(rightdir) > 0;
                    if (isOnRight){
                        rightAll++;
                    }
                    else{
                        leftAll++;
                    }
                }

                result[l.topo.hd].lineRightMayDetach = rightOccCandNum > 0 && rightConNum == 0 && 
                    leftAll > 0 &&
                    leftOccCandNum < rightOccCandNum;
                result[l.topo.hd].lineLeftMayDetach = leftOccCandNum > 0 && leftConNum == 0 &&
                    rightAll > 0 &&
                    rightOccCandNum < leftOccCandNum;
            }

            return result;

        }


        bool MayOccludes(const RLGraph & mg, const HandledTable<LineHandle, LineDetachStatus> & lrds,
            LineHandle lh,
            LineRelationHandle lrh){
            auto another = mg.topo(lrh).component<0>();
            if (another == lh){
                another = mg.topo(lrh).component<1>();
            }
            else{
                assert(mg.topo(lrh).component<1>() == lh);
            }
            auto & anotherld = mg.data(another);
            Vec3 rightdir = normalize(mg.data(lh).line.first.cross(mg.data(lh).line.second));
            bool isOnRight = (anotherld.line.center() - mg.data(lh).line.center()).dot(rightdir) > 0;
            if (isOnRight && lrds[lh].lineRightMayDetach)
                return true;
            if (!isOnRight && lrds[lh].lineLeftMayDetach)
                return true;
            return false;
        }

        bool MayOccludes(const RLGraph & mg, const HandledTable<LineHandle, LineDetachStatus> & lrds,
            LineHandle lh,
            RegionLineConnectionHandle rlh){
            auto another = mg.topo(rlh).component<0>();
            auto & anotherld = mg.data(another);
            Vec3 rightdir = normalize(mg.data(lh).line.first.cross(mg.data(lh).line.second));
            bool isOnRight = (anotherld.normalizedCenter - mg.data(lh).line.center()).dot(rightdir) > 0;
            if (isOnRight && lrds[lh].lineRightMayDetach)
                return true;
            if (!isOnRight && lrds[lh].lineLeftMayDetach)
                return true;
            return false;
        }











        namespace {

            bool AllAlong(const std::vector<std::vector<Vec3>> & pts, const Vec3 & from, const Vec3 & to, double angleThres) {
                auto n = normalize(from.cross(to));
                return std::all_of(pts.begin(), pts.end(), [&n, &angleThres](const std::vector<Vec3> & ps) {
                    return std::all_of(ps.begin(), ps.end(), [&n, &angleThres](const Vec3 & p) {
                        return abs(M_PI_2 - AngleBetweenDirections(n, p)) < angleThres;
                    });
                });
            }

            bool AllAlong(const std::vector<Vec3> & pts, const Vec3 & from, const Vec3 & to, double angleThres) {
                auto n = normalize(from.cross(to));
                return std::all_of(pts.begin(), pts.end(), [&n, angleThres](const Vec3 & p) {
                    return abs(M_PI_2 - AngleBetweenDirections(n, p)) < angleThres;
                });
            }

        }

        bool operator < (const Pixel & a, const Pixel & b) {
            return std::tie(a.x, a.y) < std::tie(b.x, b.y);
        }



        std::vector<std::vector<Vec3>> SamplesOnBoundaries(const SegmentationTopo & segtopo, const PanoramicCamera & cam, 
            double sampleStep /*= DegreesToRadians(1)*/) {

            std::vector<std::vector<Vec3>> ss(segtopo.nboundaries());

            for (int i = 0; i < segtopo.nboundaries(); i++) {
                auto & samples = ss[i];
                std::vector<Vec3> edge;
                edge.reserve(segtopo.bndpixels[i].size());
                for (auto & p : segtopo.bndpixels[i]) {
                    edge.push_back(normalize(cam.toSpace(p)));
                }
                samples = { edge.front() };
                for (auto & p : edge) {
                    double remainedAngle = AngleBetweenDirections(samples.back(), p);
                    while (remainedAngle >= sampleStep) {
                        samples.push_back(normalize(RotateDirection(samples.back(), p, sampleStep)));
                        remainedAngle -= sampleStep;
                    }
                }
            }

            return ss;
        }


        std::vector<int> ClassifyBoundaries(const std::vector<std::vector<Vec3>> & bndSamples,
            const std::vector<Vec3> & vps, double angleThres) {

            std::vector<int> cs(bndSamples.size(), -1);
            for (int i = 0; i < bndSamples.size(); i++) {
                if (bndSamples[i].size() <= 2) {
                    continue;
                }
                std::set<int> alongVPIds;
                for (int j = 0; j < vps.size(); j++) {
                    if (AllAlong(bndSamples[i], bndSamples[i].front(), vps[j], angleThres)) {
                        alongVPIds.insert(j);
                    }
                }
                if (alongVPIds.size() == 1) {
                    cs[i] = *alongVPIds.begin();
                }
            }
            return cs;

        }


        HandledTable<RegionBoundaryHandle, DepthRelation> DetectOcclusions(
            const RLGraph & mg, const RLGraphControls & controls,
            const SegmentationTopo & segtopo, const std::vector<int> & bndclasses, 
            const std::vector<RegionHandle> & rhs, const std::vector<RegionBoundaryHandle> & bhs,
            const std::vector<Vec3> & vps) {

            auto occlusions = mg.createConstraintTable<RegionBoundaryData>(DepthRelation::Connected);

            assert(segtopo.nsegs() == rhs.size());
            assert(segtopo.nboundaries() == bhs.size());
            for (int bndid = 0; bndid < segtopo.nboundaries(); bndid++) {
                auto & segids = segtopo.bnd2segs[bndid];
                auto rh1 = rhs[segids.first];
                auto rh2 = rhs[segids.second];
                if (rh1.invalid() || rh2.invalid()) {
                    continue;
                }
                auto bh = bhs[bndid];
                if (bh.invalid()) {
                    continue;
                }

                if (!controls[rh1].used || !controls[rh2].used) {
                    continue;
                }
                int bndclass = bndclasses[bndid];

                // 
                if (controls[rh1].orientationClaz != -1 && controls[rh2].orientationClaz != -1
                    && controls[rh1].orientationClaz != controls[rh2].orientationClaz
                    && (bndclass == -1 
                    || !IsFuzzyPerpendicular(vps[bndclass], vps[controls[rh1].orientationClaz]) 
                    || !IsFuzzyPerpendicular(vps[bndclass], vps[controls[rh2].orientationClaz]))) {
                    // break this
                    occlusions[bh] = DepthRelation::Disconnected;
                    auto & rhvp1 = vps[controls[rh1].orientationClaz];
                    auto & rhvp2 = vps[controls[rh2].orientationClaz];
                    if (bndclass != -1) {
                        auto & bndvp = vps[bndclass];
                        if (IsFuzzyPerpendicular(rhvp1, bndvp) && !IsFuzzyPerpendicular(rhvp2, bndvp)) {
                            occlusions[bh] = DepthRelation::FirstIsFront;
                        } else if (!IsFuzzyPerpendicular(rhvp1, bndvp) && IsFuzzyPerpendicular(rhvp2, bndvp)) {
                            occlusions[bh] = DepthRelation::SecondIsFront;
                        }
                    }
                } else if (controls[rh1].orientationClaz != -1 && 
                    controls[rh2].orientationClaz == -1 && controls[rh2].orientationNotClaz == controls[rh1].orientationClaz
                    && (bndclass == -1
                    || !IsFuzzyPerpendicular(vps[bndclass], vps[controls[rh1].orientationClaz]))) {
                    // break this
                    occlusions[bh] = DepthRelation::Disconnected;
                    if (bndclass == controls[rh1].orientationClaz) {
                        occlusions[bh] = DepthRelation::SecondIsFront;
                    }
                } else if (controls[rh2].orientationClaz != -1 &&
                    controls[rh1].orientationClaz == -1 && controls[rh1].orientationNotClaz == controls[rh2].orientationClaz
                    && (bndclass == -1
                    || !IsFuzzyPerpendicular(vps[bndclass], vps[controls[rh2].orientationClaz]))) {
                    // break this
                    occlusions[bh] = DepthRelation::Disconnected;
                    if (bndclass == controls[rh2].orientationClaz) {
                        occlusions[bh] = DepthRelation::FirstIsFront;
                    }
                }
            }

            return occlusions;

        }



        namespace {

            bool NoOrientationControl(const RLGraphComponentControl & c) {
                return c.orientationClaz == -1 && c.orientationNotClaz == -1;
            }

            bool NoOrientationControl(const OrientationControl & c) {
                return c.orientationClaz == -1 && c.orientationNotClaz == -1;
            }

        }



        HandledTable<RegionBoundaryHandle, DepthRelation> DetectOcclusions2(
            const RLGraph & mg, const RLGraphControls & controls,
            const Imagei & segs,
            const SegmentationTopo & segtopo, const std::vector<std::vector<Vec3>> & bndsamples, const std::vector<int> & bndclasses,
            const std::vector<RegionHandle> & rhs, const std::vector<RegionBoundaryHandle> & bhs,
            const std::vector<Vec3> & vps, double angleDistThres) {

            auto occlusions = mg.createConstraintTable<RegionBoundaryData>(DepthRelation::Connected);
            assert(segtopo.nsegs() == rhs.size());
            assert(segtopo.nboundaries() == bhs.size());


            // find closest bnd sample pairs to locate edge bnds along narrow segs
            struct BndSample {
                int bndid;
                int sampleid;
                bool operator == (const BndSample & bs) const { return std::tie(bndid, sampleid) == std::tie(bs.bndid, bs.sampleid); }
                bool operator < (const BndSample & bs) const { return std::tie(bndid, sampleid) < std::tie(bs.bndid, bs.sampleid); }
            };
            
            RTreeMap<Vec3, BndSample> samples;
            for (int bndid = 0; bndid < bndsamples.size(); bndid++) {
                auto & segids = segtopo.bnd2segs[bndid];
                auto rh1 = rhs[segids.first];
                auto rh2 = rhs[segids.second];
                if (rh1.invalid() || rh2.invalid()) {
                    continue;
                }
                if (std::tie(controls[rh1].orientationClaz, controls[rh1].orientationNotClaz) ==
                    std::tie(controls[rh2].orientationClaz, controls[rh2].orientationNotClaz)) {
                    continue;
                }
                if (!NoOrientationControl(controls[rh1]) && !NoOrientationControl(controls[rh2])) {
                    continue; // only consider detached edges
                }
                for (int sampleid = 0; sampleid < bndsamples[bndid].size(); sampleid++) {
                    samples.insert(std::make_pair(bndsamples[bndid][sampleid], BndSample{ bndid, sampleid }));
                }
            }

            std::map<BndSample, BndSample> closestBs;
            for (int bndid = 0; bndid < bndsamples.size(); bndid++) {
                for (int sampleid = 0; sampleid < bndsamples[bndid].size(); sampleid++) {

                    auto & segids = segtopo.bnd2segs[bndid];
                    auto rh1 = rhs[segids.first];
                    auto rh2 = rhs[segids.second];
                    if (rh1.invalid() || rh2.invalid()) {
                        continue;
                    }
                    if (std::tie(controls[rh1].orientationClaz, controls[rh1].orientationNotClaz) ==
                        std::tie(controls[rh2].orientationClaz, controls[rh2].orientationNotClaz)) {
                        continue;
                    }
                    if (!NoOrientationControl(controls[rh1]) && !NoOrientationControl(controls[rh2])) {
                        continue; // only consider detached edges
                    }

                    auto controledSide = controls[rh1];
                    if (NoOrientationControl(controledSide)) {
                        controledSide = controls[rh2];
                    }
                    assert(!NoOrientationControl(controledSide));

                    auto & dir = bndsamples[bndid][sampleid];

                    BndSample closest = { -1, -1 };
                    double minAngle = angleDistThres;
                    
                    samples.search(BoundingBox(normalize(dir)).expand(angleDistThres * 2),
                        [&closest, &minAngle, &segtopo, &rhs, &controls, &controledSide, &dir](const std::pair<Vec3, BndSample> & bs) {

                        auto & segids = segtopo.bnd2segs[bs.second.bndid];
                        auto rh1 = rhs[segids.first];
                        auto rh2 = rhs[segids.second];

                        auto thisControledSide = controls[rh1];
                        if (NoOrientationControl(thisControledSide)) {
                            thisControledSide = controls[rh2];
                        }
                        assert(!NoOrientationControl(thisControledSide));

                        if (std::tie(thisControledSide.orientationClaz, thisControledSide.orientationNotClaz) ==
                            std::tie(controledSide.orientationClaz, controledSide.orientationNotClaz)) {
                            return true;
                        }

                        double angle = AngleBetweenDirections(dir, bs.first);
                        if (angle < minAngle) {
                            closest = bs.second;
                            minAngle = angle;
                        }

                        return true;
                    });

                    if (closest.bndid != -1) {
                        closestBs[BndSample{ bndid, sampleid }] = closest;
                    }
                }
            }

            std::vector<std::pair<BndSample, BndSample>> closestPairs;
            for (auto & p : closestBs) {
                if (closestBs.at(p.second) == p.first) {
                    closestPairs.emplace_back(p.first, p.second);
                }
            }

            std::map<std::pair<int, int>, double> bndPairCorrespRatios;
            for (auto & samplePair : closestPairs) {
                bndPairCorrespRatios[std::make_pair(samplePair.first.bndid, samplePair.second.bndid)] += 
                    1.0 / bndsamples[samplePair.first.bndid].size();
                bndPairCorrespRatios[std::make_pair(samplePair.second.bndid, samplePair.first.bndid)] +=
                    1.0 / bndsamples[samplePair.second.bndid].size();
            }

            // segidpair -> [bndid -> ratio]
            std::map<std::pair<int, int>, std::map<int, double>> segpair2bnds;
            for (int bndid = 0; bndid < bndsamples.size(); bndid++) {
                auto & segids = segtopo.bnd2segs[bndid];
                auto rh1 = rhs[segids.first];
                auto rh2 = rhs[segids.second];
                if (rh1.invalid() || rh2.invalid()) {
                    continue;
                }
                if (std::tie(controls[rh1].orientationClaz, controls[rh1].orientationNotClaz) ==
                    std::tie(controls[rh2].orientationClaz, controls[rh2].orientationNotClaz)) {
                    continue;
                }
                if (!NoOrientationControl(controls[rh1]) && !NoOrientationControl(controls[rh2])) {
                    auto orderedSegIds = segids.first < segids.second ? segids : std::make_pair(segids.second, segids.first);
                    segpair2bnds[orderedSegIds][bndid] = 1.0;
                }
            }
            for (auto & p : closestPairs) {
                int bndid1 = p.first.bndid;
                int bndid2 = p.second.bndid;
                assert(bndid1 != p.second.bndid);

                int constrainedSegId1 = segtopo.bnd2segs[bndid1].first;
                if (NoOrientationControl(controls[rhs[segtopo.bnd2segs[bndid1].first]])) {
                    constrainedSegId1 = segtopo.bnd2segs[bndid1].second;
                }

                int constrainedSegId2 = segtopo.bnd2segs[bndid2].first;
                if (NoOrientationControl(controls[rhs[segtopo.bnd2segs[bndid2].first]])) {
                    constrainedSegId2 = segtopo.bnd2segs[bndid2].second;
                }

                assert(!NoOrientationControl(controls[rhs[constrainedSegId1]]) &&
                    !NoOrientationControl(controls[rhs[constrainedSegId2]]));

                auto & bnds =
                    segpair2bnds[constrainedSegId1 < constrainedSegId2 ? 
                    std::make_pair(constrainedSegId1, constrainedSegId2) : 
                    std::make_pair(constrainedSegId2, constrainedSegId2)];

                bnds[bndid1] += 1.0 / bndsamples[bndid1].size();
                bnds[bndid2] += 1.0 / bndsamples[bndid2].size();
            }

            for (auto & segpairbnd : segpair2bnds) {
                const auto & segids = segpairbnd.first;
                auto & bnds = segpairbnd.second;

                auto rh1 = rhs[segids.first];
                auto rh2 = rhs[segids.second];
                if (rh1.invalid() || rh2.invalid()) {
                    continue;
                }

                for (auto & bndwithscore : bnds) {
                    int bndid = bndwithscore.first;

                    double occupationRatio = bndwithscore.second;
                    if (occupationRatio < 0.5) {
                        continue;
                    }

                    auto bh = bhs[bndid];
                    if (bh.invalid()) {
                        continue;
                    }

                    if (!controls[rh1].used || !controls[rh2].used) {
                        continue;
                    }
                    int bndclass = bndclasses[bndid];


                    if (controls[rh1].orientationClaz != -1 && controls[rh2].orientationClaz != -1
                        && controls[rh1].orientationClaz != controls[rh2].orientationClaz
                        && bndclass != -1
                        && IsFuzzyPerpendicular(vps[bndclass], vps[controls[rh1].orientationClaz])
                        && IsFuzzyPerpendicular(vps[bndclass], vps[controls[rh2].orientationClaz])) {
                        // break this
                        occlusions[bh] = DepthRelation::MaybeFolder;   
                    } else if (controls[rh1].orientationClaz != -1 && controls[rh2].orientationClaz != -1
                        && controls[rh1].orientationClaz != controls[rh2].orientationClaz
                        && (bndclass == -1
                        || !IsFuzzyPerpendicular(vps[bndclass], vps[controls[rh1].orientationClaz])
                        || !IsFuzzyPerpendicular(vps[bndclass], vps[controls[rh2].orientationClaz]))) {
                        // break this
                        occlusions[bh] = DepthRelation::Disconnected;
                        auto & rhvp1 = vps[controls[rh1].orientationClaz];
                        auto & rhvp2 = vps[controls[rh2].orientationClaz];
                        if (bndclass != -1) {
                            auto & bndvp = vps[bndclass];
                            if (IsFuzzyPerpendicular(rhvp1, bndvp) && !IsFuzzyPerpendicular(rhvp2, bndvp)) {
                                occlusions[bh] = DepthRelation::FirstIsFront;
                            } else if (!IsFuzzyPerpendicular(rhvp1, bndvp) && IsFuzzyPerpendicular(rhvp2, bndvp)) {
                                occlusions[bh] = DepthRelation::SecondIsFront;
                            }
                        }
                    } else if (controls[rh1].orientationClaz != -1 &&
                        controls[rh2].orientationClaz == -1 && controls[rh2].orientationNotClaz == controls[rh1].orientationClaz
                        && (bndclass == -1
                        || !IsFuzzyPerpendicular(vps[bndclass], vps[controls[rh1].orientationClaz]))) {
                        // break this
                        occlusions[bh] = DepthRelation::Disconnected;
                        if (bndclass == controls[rh1].orientationClaz) {
                            occlusions[bh] = DepthRelation::SecondIsFront;
                        }
                    } else if (controls[rh2].orientationClaz != -1 &&
                        controls[rh1].orientationClaz == -1 && controls[rh1].orientationNotClaz == controls[rh2].orientationClaz
                        && (bndclass == -1
                        || !IsFuzzyPerpendicular(vps[bndclass], vps[controls[rh2].orientationClaz]))) {
                        // break this
                        occlusions[bh] = DepthRelation::Disconnected;
                        if (bndclass == controls[rh2].orientationClaz) {
                            occlusions[bh] = DepthRelation::FirstIsFront;
                        }
                    } else if (controls[rh1].orientationClaz != -1 &&
                        controls[rh2].orientationClaz == -1 && controls[rh2].orientationNotClaz == controls[rh1].orientationClaz
                        && bndclass != -1
                        && IsFuzzyPerpendicular(vps[bndclass], vps[controls[rh1].orientationClaz])) {
                        occlusions[bh] = DepthRelation::MaybeFolder;
                    } else if (controls[rh2].orientationClaz != -1 &&
                        controls[rh1].orientationClaz == -1 && controls[rh1].orientationNotClaz == controls[rh2].orientationClaz
                        && bndclass != -1
                        && IsFuzzyPerpendicular(vps[bndclass], vps[controls[rh2].orientationClaz])) {
                        occlusions[bh] = DepthRelation::MaybeFolder;
                    }
                }
            }

            return occlusions;

        }


        namespace {
            template <class T>
            inline std::pair<T, T> MakeOrderedPair(const T & a, const T & b) {
                return a < b ? std::make_pair(a, b) : std::make_pair(b, a);
            }
        }


        HandledTable<RegionBoundaryHandle, DepthRelation> DetectOcclusions3(
            const RLGraph & mg, const RLGraphControls & controls,
            const Imagei & segs,
            const SegmentationTopo & segtopo, const std::vector<std::vector<Vec3>> & bndsamples, const std::vector<int> & bndclasses,
            const std::vector<RegionHandle> & rhs, const std::vector<RegionBoundaryHandle> & bhs,
            const std::vector<Vec3> & vps, double angleDistThres, double angleSampleStepOnLine) {

            // register all sample points           
            struct BndSample {
                int bndid;
                int sampleid;
                bool operator == (const BndSample & bs) const { return std::tie(bndid, sampleid) == std::tie(bs.bndid, bs.sampleid); }
                bool operator < (const BndSample & bs) const { return std::tie(bndid, sampleid) < std::tie(bs.bndid, bs.sampleid); }
            };
            RTreeMap<Vec3, BndSample> samples;
            for (int bndid = 0; bndid < bndsamples.size(); bndid++) {
                auto & segids = segtopo.bnd2segs[bndid]; 
                auto rh1 = rhs[segids.first];
                auto rh2 = rhs[segids.second];
                if (rh1.invalid() || rh2.invalid()) {
                    continue;
                }
                for (int sampleid = 0; sampleid < bndsamples[bndid].size(); sampleid++) {
                    samples.insert(std::make_pair(bndsamples[bndid][sampleid], BndSample{ bndid, sampleid }));
                }
            }

            // find closest bnd sample pairs to locate edge bnds along narrow segs
            // segidpair -> [bndid -> ratio]
            std::map<std::pair<int, int>, std::map<int, double>> segpair2bnds;
            std::map<std::pair<int, int>, int> bndPairCorrespSampleNums;
            {
                std::map<BndSample, BndSample> closestBs;
                for (int bndid = 0; bndid < bndsamples.size(); bndid++) {
                    for (int sampleid = 0; sampleid < bndsamples[bndid].size(); sampleid++) {
                        auto & segids = segtopo.bnd2segs[bndid];
                        auto rh1 = rhs[segids.first];
                        auto rh2 = rhs[segids.second];
                        if (rh1.invalid() || rh2.invalid()) {
                            continue;
                        }
                        if (std::tie(controls[rh1].orientationClaz, controls[rh1].orientationNotClaz) ==
                            std::tie(controls[rh2].orientationClaz, controls[rh2].orientationNotClaz)) {
                            continue;
                        }
                        if (!NoOrientationControl(controls[rh1]) && !NoOrientationControl(controls[rh2])) {
                            continue; // only consider detached edges
                        }

                        auto controledSide = controls[rh1];
                        if (NoOrientationControl(controledSide)) {
                            controledSide = controls[rh2];
                        }
                        assert(!NoOrientationControl(controledSide));

                        auto & dir = bndsamples[bndid][sampleid];

                        BndSample closest = { -1, -1 };
                        double minAngle = angleDistThres;

                        samples.search(BoundingBox(normalize(dir)).expand(angleDistThres * 2),
                            [&closest, &minAngle, &segtopo, &rhs, &controls, &controledSide, &dir](const std::pair<Vec3, BndSample> & bs) {

                            auto & segids = segtopo.bnd2segs[bs.second.bndid];
                            auto rh1 = rhs[segids.first];
                            auto rh2 = rhs[segids.second];

                            if (std::tie(controls[rh1].orientationClaz, controls[rh1].orientationNotClaz) ==
                                std::tie(controls[rh2].orientationClaz, controls[rh2].orientationNotClaz)) {
                                return true;
                            }
                            if (!NoOrientationControl(controls[rh1]) && !NoOrientationControl(controls[rh2])) {
                                return true; // only consider detached edges
                            }

                            auto thisControledSide = controls[rh1];
                            if (NoOrientationControl(thisControledSide)) {
                                thisControledSide = controls[rh2];
                            }
                            assert(!NoOrientationControl(thisControledSide));

                            if (std::tie(thisControledSide.orientationClaz, thisControledSide.orientationNotClaz) ==
                                std::tie(controledSide.orientationClaz, controledSide.orientationNotClaz)) {
                                return true;
                            }

                            double angle = AngleBetweenDirections(dir, bs.first);
                            if (angle < minAngle) {
                                closest = bs.second;
                                minAngle = angle;
                            }

                            return true;
                        });

                        if (closest.bndid != -1) {
                            closestBs[BndSample{ bndid, sampleid }] = closest;
                        }
                    }
                }

                std::vector<std::pair<BndSample, BndSample>> closestPairs;
                for (auto & p : closestBs) {
                    if (closestBs.at(p.second) == p.first) {
                        closestPairs.emplace_back(p.first, p.second);
                    }
                }
                for (auto & samplePair : closestPairs) {
                    bndPairCorrespSampleNums[MakeOrderedPair(samplePair.first.bndid, samplePair.second.bndid)] ++;
                }


                for (int bndid = 0; bndid < bndsamples.size(); bndid++) {
                    auto & segids = segtopo.bnd2segs[bndid];
                    auto rh1 = rhs[segids.first];
                    auto rh2 = rhs[segids.second];
                    if (rh1.invalid() || rh2.invalid()) {
                        continue;
                    }
                    if (std::tie(controls[rh1].orientationClaz, controls[rh1].orientationNotClaz) ==
                        std::tie(controls[rh2].orientationClaz, controls[rh2].orientationNotClaz)) {
                        continue;
                    }
                    if (!NoOrientationControl(controls[rh1]) && !NoOrientationControl(controls[rh2])) {
                        auto orderedSegIds = segids.first < segids.second ? segids : std::make_pair(segids.second, segids.first);
                        segpair2bnds[orderedSegIds][bndid] = 1.0;
                    }
                }
                for (auto & p : closestPairs) {
                    int bndid1 = p.first.bndid;
                    int bndid2 = p.second.bndid;
                    assert(bndid1 != p.second.bndid);

                    int constrainedSegId1 = segtopo.bnd2segs[bndid1].first;
                    if (NoOrientationControl(controls[rhs[segtopo.bnd2segs[bndid1].first]])) {
                        constrainedSegId1 = segtopo.bnd2segs[bndid1].second;
                    }

                    int constrainedSegId2 = segtopo.bnd2segs[bndid2].first;
                    if (NoOrientationControl(controls[rhs[segtopo.bnd2segs[bndid2].first]])) {
                        constrainedSegId2 = segtopo.bnd2segs[bndid2].second;
                    }

                    assert(!NoOrientationControl(controls[rhs[constrainedSegId1]]) &&
                        !NoOrientationControl(controls[rhs[constrainedSegId2]]));

                    auto & bnds =
                        segpair2bnds[constrainedSegId1 < constrainedSegId2 ?
                        std::make_pair(constrainedSegId1, constrainedSegId2) :
                        std::make_pair(constrainedSegId2, constrainedSegId2)];

                    bnds[bndid1] += 1.0 / bndsamples[bndid1].size();
                    bnds[bndid2] += 1.0 / bndsamples[bndid2].size();
                }
            }

            // get first guess
            auto firstGuess = mg.createConstraintTable<RegionBoundaryData>(DepthRelation::Unknown);
            for (auto & segpairbnd : segpair2bnds) {
                const auto & segids = segpairbnd.first;
                auto & bnds = segpairbnd.second;
                
                auto rh1 = rhs[segids.first];
                auto rh2 = rhs[segids.second];
                if (rh1.invalid() || rh2.invalid()) {
                    continue;
                }
                for (auto & bndwithscore : bnds) {
                    int bndid = bndwithscore.first;

                    double occupationRatio = bndwithscore.second;
                    if (occupationRatio < 0.5) {
                        continue;
                    }

                    auto bh = bhs[bndid];
                    if (bh.invalid()) {
                        continue;
                    }

                    if (!controls[rh1].used || !controls[rh2].used) {
                        continue;
                    }
                    int bndclass = bndclasses[bndid];

                    if (controls[rh1].orientationClaz != -1 && controls[rh2].orientationClaz != -1
                        && controls[rh1].orientationClaz != controls[rh2].orientationClaz
                        && (bndclass == -1
                        || !IsFuzzyPerpendicular(vps[bndclass], vps[controls[rh1].orientationClaz])
                        || !IsFuzzyPerpendicular(vps[bndclass], vps[controls[rh2].orientationClaz]))) {
                        // break this
                        firstGuess[bh] = DepthRelation::Disconnected;
                        auto & rhvp1 = vps[controls[rh1].orientationClaz];
                        auto & rhvp2 = vps[controls[rh2].orientationClaz];
                        if (bndclass != -1) {
                            auto & bndvp = vps[bndclass];
                            if (IsFuzzyPerpendicular(rhvp1, bndvp) && !IsFuzzyPerpendicular(rhvp2, bndvp)) {
                                firstGuess[bh] = DepthRelation::FirstIsFront;
                            } else if (!IsFuzzyPerpendicular(rhvp1, bndvp) && IsFuzzyPerpendicular(rhvp2, bndvp)) {
                                firstGuess[bh] = DepthRelation::SecondIsFront;
                            }
                        }
                    } else if (controls[rh1].orientationClaz != -1 &&
                        controls[rh2].orientationClaz == -1 && controls[rh2].orientationNotClaz == controls[rh1].orientationClaz
                        && (bndclass == -1
                        || !IsFuzzyPerpendicular(vps[bndclass], vps[controls[rh1].orientationClaz]))) {
                        // break this
                        firstGuess[bh] = DepthRelation::Disconnected;
                        if (bndclass == controls[rh1].orientationClaz) {
                            firstGuess[bh] = DepthRelation::SecondIsFront;
                        }
                    } else if (controls[rh2].orientationClaz != -1 &&
                        controls[rh1].orientationClaz == -1 && controls[rh1].orientationNotClaz == controls[rh2].orientationClaz
                        && (bndclass == -1
                        || !IsFuzzyPerpendicular(vps[bndclass], vps[controls[rh2].orientationClaz]))) {
                        // break this
                        firstGuess[bh] = DepthRelation::Disconnected;
                        if (bndclass == controls[rh2].orientationClaz) {
                            firstGuess[bh] = DepthRelation::FirstIsFront;
                        }
                    } else if (controls[rh1].orientationClaz != -1 && controls[rh2].orientationClaz != -1 && 
                        controls[rh1].orientationClaz == controls[rh2].orientationClaz) {
                        firstGuess[bh] = DepthRelation::Connected;
                    } else if (controls[rh1].orientationClaz == -1 && controls[rh1].orientationNotClaz != -1
                        && controls[rh2].orientationClaz == -1 && controls[rh2].orientationNotClaz != -1
                        && controls[rh1].orientationNotClaz == controls[rh2].orientationNotClaz
                        && bndclass != -1 && bndclass != controls[rh1].orientationNotClaz) {
                        firstGuess[bh] = DepthRelation::Connected;
                    }
                }
            }


            // get continuous bnds
            //std::set<std::pair<int, int>> continousBnds;
            //{
            //    for (int bndid = 0; bndid < bndsamples.size(); bndid++) {
            //        if (!bhs[bndid].valid()) {
            //            continue;
            //        }
            //        if (bndclasses[bndid] == -1) {
            //            continue;
            //        }
            //        for (int sampleid = 0; sampleid < bndsamples[bndid].size(); sampleid++) {
            //            auto & segids = segtopo.bnd2segs[bndid];
            //            auto rh1 = rhs[segids.first];
            //            auto rh2 = rhs[segids.second];
            //            if (rh1.invalid() || rh2.invalid()) {
            //                continue;
            //            }
            //            auto & dir = bndsamples[bndid][sampleid];

            //            samples.search(BoundingBox(normalize(dir)).expand(angleDistThres * 2),
            //                [bndid, &bndsamples, &dir, angleDistThres, &continousBnds, &bndclasses, &vps](const std::pair<Vec3, BndSample> & bs) {
            //                if (bndclasses[bs.second.bndid] != bndclasses[bndid]) {
            //                    return true;
            //                }
            //                if (bs.second.bndid == bndid) {
            //                    return true;
            //                }
            //                auto & thisDir = bndsamples[bs.second.bndid][bs.second.sampleid];

            //                if (AngleBetweenDirections(dir, thisDir) < angleDistThres * 2 && 
            //                    AllAlong(std::vector<std::vector<Vec3>>{bndsamples[bndid], bndsamples[bs.second.bndid]}, dir, vps[bndclasses[bndid]], angleDistThres)) {
            //                    continousBnds.insert(MakeOrderedPair(bndid, bs.second.bndid));
            //                }
            //                return true;
            //            });
            //        }
            //    }
            //}

            auto bndsAttachedToLines = mg.createComponentTable<LineData, std::vector<int>>();
            {               
                for (auto & l : mg.components<LineData>()) {
                    std::vector<std::set<int>> bndSamplesNearLines(segtopo.nboundaries());
                    auto & line = l.data.line;
                    double angle = AngleBetweenDirections(line.first, line.second);
                    if (angle < DegreesToRadians(1)) {
                        continue;
                    }
                    double angleStep = angleSampleStepOnLine;
                    for (double a = 0.0; a <= angle; a += angleStep) {
                        Vec3 dir = RotateDirection(line.first, line.second, a);
                        BndSample nearest = { -1, -1 };
                        double minAngleDist = M_PI;

                        // find nearest bndSample
                        samples.search(BoundingBox(normalize(dir)).expand(angleDistThres * 2),
                            [&bndsamples, &dir, &nearest, &minAngleDist](const std::pair<Vec3, BndSample> & bs) {
                            auto & thisDir = bndsamples[bs.second.bndid][bs.second.sampleid];
                            double angleDist = AngleBetweenDirections(thisDir, dir);
                            if (minAngleDist > angleDist) {
                                minAngleDist = angleDist;
                                nearest = bs.second;
                            }
                            return true;
                        });

                        if (nearest.bndid == -1) {
                            continue;
                        }
                        bndSamplesNearLines[nearest.bndid].insert(nearest.sampleid);
                    }

                    auto & bndsAttached = bndsAttachedToLines[l.topo.hd];
                    for (int bndid = 0; bndid < segtopo.nboundaries(); bndid++) {
                        if (bndSamplesNearLines[bndid].size() > bndsamples[bndid].size() * 0.6) {
                            bndsAttached.push_back(bndid);
                        }
                    }
                }
            }


            std::cout << "building factor graph" << std::endl;

            // factor graph for bnds
            ml::FactorGraph fg;
            static const double spnum2weightRatio = 1.0 / 100.0;

            // add bnd as vars
            std::vector<ml::FactorGraph::VarHandle> vhs(segtopo.nboundaries());
            for (int bndid = 0; bndid < segtopo.nboundaries(); bndid++) {
                if (bhs[bndid].invalid()) {
                    continue;
                }
                vhs[bndid] = fg.addVar(fg.addVarCategory(3, bndsamples[bndid].size() / 10.0)); // first front, second front, connected
            }

            // first guess factor
            for (int bndid = 0; bndid < segtopo.nboundaries(); bndid++) {
                if (vhs[bndid].invalid()) {
                    continue;
                }
                auto guess = firstGuess[bhs[bndid]];
                size_t spnum = bndsamples[bndid].size();
                int unaryFC = fg.addFactorCategory(ml::FactorGraph::FactorCategory{ [guess, spnum](const int * varlabels, size_t nvar,
                    ml::FactorGraph::FactorCategoryId fcid, void * givenData) -> double {
                    double weight = spnum * spnum2weightRatio;
                    assert(nvar == 1);
                    DepthRelation label = static_cast<DepthRelation>(varlabels[0]);
                    switch (guess) {
                    case DepthRelation::Unknown: return label == DepthRelation::Connected ? 0.0 : 1.0 * weight;  // always prefer connected
                    case DepthRelation::Connected: return label == DepthRelation::Connected ? 0.0 : 3.0 * weight;
                    case DepthRelation::Disconnected: 
                        return Contains({ 
                            DepthRelation::Disconnected, DepthRelation::FirstIsFront, DepthRelation::SecondIsFront 
                        }, label) ? 0.0 : 10.0 * weight;
                    case DepthRelation::FirstIsFront: return label == DepthRelation::FirstIsFront ? 0.0 : 20.0 * weight;
                    case DepthRelation::SecondIsFront: return label == DepthRelation::SecondIsFront ? 0.0 : 20.0 * weight;
                    default: SHOULD_NEVER_BE_CALLED("this label should never occur in first guess!");
                    }
                }, 10.0 });
                fg.addFactor({ vhs[bndid] }, unaryFC);
            }

            // overlap factor
            for (auto & bndPair : bndPairCorrespSampleNums) {
                if (bhs[bndPair.first.first].invalid() || bhs[bndPair.first.second].invalid()) {
                    continue;
                }
                int overlapSPNum = bndPair.second;
                double weight = spnum2weightRatio * overlapSPNum;
                int overlapFC = fg.addFactorCategory(ml::FactorGraph::FactorCategory{ [weight](const int * varlabels, size_t nvar,
                    ml::FactorGraph::FactorCategoryId fcid, void * givenData) -> double {
                    assert(nvar == 2);
                    DepthRelation label1 = static_cast<DepthRelation>(varlabels[0]);
                    DepthRelation label2 = static_cast<DepthRelation>(varlabels[1]);
                    if (label1 != DepthRelation::Connected && label2 != DepthRelation::Connected) {
                        return 10.0 * weight;
                    }
                    return 0.0;
                }, 1.0 });
                fg.addFactor({ vhs[bndPair.first.first], vhs[bndPair.first.second] }, overlapFC);
            }

            // junction validity factor
            for (int jid = 0; jid < segtopo.njunctions(); jid++) {
                auto & bndids = segtopo.junc2bnds[jid];
                if (std::any_of(bndids.begin(), bndids.end(), [&bhs](int bndid) {return bhs[bndid].invalid(); })) {
                    continue;
                }
                double weight = 0.0;
                for (int bndid : bndids) {
                    weight += spnum2weightRatio * bndsamples[bndid].size();
                }
                int juncValidFC = fg.addFactorCategory(ml::FactorGraph::FactorCategory{ [&bndids, weight](const int * varlabels, size_t nvar,
                    ml::FactorGraph::FactorCategoryId fcid, void * givenData) -> double {
                    assert(nvar == bndids.size());
                    int disconNum = 0;
                    for (int i = 0; i < nvar; i++) {
                        if (DepthRelation(varlabels[i]) != DepthRelation::Connected) {
                            disconNum++;
                        }
                    }
                    if (disconNum >= 3) {
                        return 10.0 * weight;
                    }
                    return 0.0;
                }, 1.0 });
                std::vector<ml::FactorGraph::VarHandle> relatedVHs(bndids.size());
                for (int i = 0; i < bndids.size(); i++) {
                    relatedVHs[i] = vhs[bndids[i]];
                }
                fg.addFactor(relatedVHs.begin(), relatedVHs.end(), juncValidFC);
            }

            // continuity factor
            //for (auto & contBndPair : continousBnds) {
            //    // todo judge bnd label side: right left
            //    auto & segs1 = segtopo.bnd2segs[contBndPair.first];
            //    Vec3 bnd1SegPairDir = normalize(mg.data(rhs[segs1.first]).normalizedCenter.cross(mg.data(rhs[segs1.second]).normalizedCenter));
            //    auto & segs2 = segtopo.bnd2segs[contBndPair.second];
            //    Vec3 bnd2SegPairDir = normalize(mg.data(rhs[segs2.first]).normalizedCenter.cross(mg.data(rhs[segs2.second]).normalizedCenter));
            //    bool sameSegPairDir = bnd1SegPairDir.dot(bnd2SegPairDir) > 0;

            //    double weight = (bndsamples[contBndPair.first].size() + bndsamples[contBndPair.second].size()) * spnum2weightRatio;
            //    int bndContFC = fg.addFactorCategory(ml::FactorGraph::FactorCategory{ [weight, sameSegPairDir](const int * varlabels, size_t nvar,
            //        ml::FactorGraph::FactorCategoryId fcid, void * givenData) -> double {
            //        assert(nvar == 2);                    
            //        DepthRelation rel1 = (DepthRelation)varlabels[0];
            //        DepthRelation rel2 = (DepthRelation)varlabels[1];
            //        if (rel1 == DepthRelation::Connected && rel2 == DepthRelation::Connected ||
            //            rel1 == DepthRelation::Disconnected && rel2 == DepthRelation::Disconnected ||
            //            rel1 == DepthRelation::Disconnected && Contains({DepthRelation::FirstIsFront, DepthRelation::SecondIsFront}, rel2) ||
            //            rel2 == DepthRelation::Disconnected && Contains({DepthRelation::FirstIsFront, DepthRelation::SecondIsFront}, rel1)) {
            //            return 0.0;
            //        }                    
            //        if (sameSegPairDir) {
            //            if (rel1 == DepthRelation::FirstIsFront && rel2 == DepthRelation::FirstIsFront ||
            //                rel1 == DepthRelation::SecondIsFront && rel2 == DepthRelation::SecondIsFront) {
            //                return 0.0;
            //            }
            //        } else {
            //            if (rel1 == DepthRelation::FirstIsFront && rel2 == DepthRelation::SecondIsFront ||
            //                rel1 == DepthRelation::SecondIsFront && rel2 == DepthRelation::FirstIsFront) {
            //                return 0.0;
            //            }
            //        }
            //        return 2.0 * weight;
            //    }, 1.0 });
            //}

            //std::map<int, std::vector<int>> bndSegPairDirSameWithFirstTable;
            //for (auto & bnds : bndsAttachedToLines) {
            //    if (bnds.size() <= 1) {
            //        continue;
            //    }

            //    int firstBnd = bnds.front();
            //    auto & firstSegs = segtopo.bnd2segs[firstBnd];
            //    Vec3 firstBndSegPairDir = normalize(mg.data(rhs[firstSegs.first]).normalizedCenter.cross(mg.data(rhs[firstSegs.second]).normalizedCenter));

            //    std::vector<int> bndSegPairDirSameWithFirst(bnds.size());
            //    for (int i = 0; i < bnds.size(); i++) {
            //        int thisBnd = bnds[i];
            //        auto & thisSegs = segtopo.bnd2segs[thisBnd];
            //        Vec3 thisBndSegPairDir = normalize(mg.data(rhs[thisSegs.first]).normalizedCenter.cross(mg.data(rhs[thisSegs.second]).normalizedCenter));
            //        bndSegPairDirSameWithFirst[i] = (thisBndSegPairDir.dot(firstBndSegPairDir) > 0);
            //    }

            //    double weight = 0.0;
            //    for (int bndid : bnds) {
            //        weight += bndsamples[bndid].size() * spnum2weightRatio;
            //    }
            //    int bndContFC = fg.addFactorCategory(ml::FactorGraph::FactorCategory{ [weight, &bndSegPairDirSameWithFirstTable](const int * varlabels, size_t nvar,
            //        ml::FactorGraph::FactorCategoryId fcid, void * givenData) -> double {
            //        auto & bndSegPairDirSameWithFirst = bndSegPairDirSameWithFirstTable.at(fcid);
            //        assert(nvar == bndSegPairDirSameWithFirst.size());
            //        enum BndContStatus { 
            //            SegPairDirSameWithFirstTableMeansFirstIsFront, 
            //            SegPairDirSameWithFirstTableMeansSecondIsFront, 
            //            JustDisconnected,
            //            AllConnected,
            //            Undefined
            //        };
            //        BndContStatus status = Undefined;
            //        bool hasViolations = false;
            //        for (int i = 0; i < nvar; i++) {
            //            auto rel = (DepthRelation)varlabels[i];
            //            if (rel == DepthRelation::Connected) {
            //                if (status == Undefined) {
            //                    status = AllConnected;
            //                    continue;
            //                }
            //                if (status == AllConnected) {
            //                    continue;
            //                }
            //                hasViolations = true;
            //                break;
            //            } else if (rel == DepthRelation::Disconnected) {
            //                if (status == Undefined) {
            //                    status = JustDisconnected;
            //                    continue;
            //                }
            //                if (status == AllConnected) {
            //                    hasViolations = true;
            //                    break;
            //                }                            
            //            } else if (rel == DepthRelation::FirstIsFront) {
            //                int segPairDirSame = bndSegPairDirSameWithFirst[i];
            //                if (status == Undefined) {
            //                    status = segPairDirSame ? SegPairDirSameWithFirstTableMeansFirstIsFront : SegPairDirSameWithFirstTableMeansSecondIsFront;
            //                    continue;
            //                }
            //                if (status == AllConnected) {
            //                    hasViolations = true;
            //                    break;
            //                }
            //                if (status == JustDisconnected) {
            //                    status = segPairDirSame ? SegPairDirSameWithFirstTableMeansFirstIsFront : SegPairDirSameWithFirstTableMeansSecondIsFront;
            //                    continue;
            //                }
            //                if (status == SegPairDirSameWithFirstTableMeansSecondIsFront) {
            //                    if (segPairDirSame) {
            //                        hasViolations = true;
            //                        break;
            //                    }
            //                }
            //                if (status == SegPairDirSameWithFirstTableMeansFirstIsFront) {
            //                    if (!segPairDirSame) {
            //                        hasViolations = true;
            //                        break;
            //                    }
            //                }
            //            } else if (rel == DepthRelation::SecondIsFront) {
            //                int segPairDirSame = bndSegPairDirSameWithFirst[i];
            //                if (status == Undefined) {
            //                    status = segPairDirSame ? SegPairDirSameWithFirstTableMeansSecondIsFront : SegPairDirSameWithFirstTableMeansFirstIsFront;
            //                    continue;
            //                }
            //                if (status == AllConnected) {
            //                    hasViolations = true;
            //                    break;
            //                }
            //                if (status == JustDisconnected) {
            //                    status = segPairDirSame ? SegPairDirSameWithFirstTableMeansSecondIsFront : SegPairDirSameWithFirstTableMeansFirstIsFront;
            //                    continue;
            //                }
            //                if (status == SegPairDirSameWithFirstTableMeansSecondIsFront) {
            //                    if (!segPairDirSame) {
            //                        hasViolations = true;
            //                        break;
            //                    }
            //                }
            //                if (status == SegPairDirSameWithFirstTableMeansFirstIsFront) {
            //                    if (segPairDirSame) {
            //                        hasViolations = true;
            //                        break;
            //                    }
            //                }
            //            }
            //        }
            //        return hasViolations ? weight * 5.0 : 0.0;
            //    }, 1.0 });

            //    bndSegPairDirSameWithFirstTable[bndContFC] = std::move(bndSegPairDirSameWithFirst);

            //    std::vector<ml::FactorGraph::VarHandle> relatedVHs(bnds.size());
            //    for (int i = 0; i < bnds.size(); i++) {
            //        relatedVHs[i] = vhs[bnds[i]];
            //    }
            //    fg.addFactor(relatedVHs.begin(), relatedVHs.end(), bndContFC);
            //}


            for (auto & bnds : bndsAttachedToLines) {
                if (bnds.size() <= 1) {
                    continue;
                }
                for (int i = 0; i < bnds.size(); i++) {
                    for (int j = i + 1; (j < i + 9 && j < bnds.size()); j++) {
                        int bnd1 = bnds[i];
                        int bnd2 = bnds[j];
                        if (bhs[bnd1].invalid() || bhs[bnd2].invalid()) {
                            continue;
                        }

                        auto & segs1 = segtopo.bnd2segs[bnd1];
                        auto & bh1 = bhs[bnd1];
                        auto & bh2 = bhs[bnd2];
                        Vec3 bnd1SegPairDir = normalize(mg.data(mg.topo(bh1).component<0>()).normalizedCenter.cross(mg.data(mg.topo(bh1).component<1>()).normalizedCenter));
                        auto & segs2 = segtopo.bnd2segs[bnd2];
                        Vec3 bnd2SegPairDir = normalize(mg.data(mg.topo(bh2).component<0>()).normalizedCenter.cross(mg.data(mg.topo(bh2).component<1>()).normalizedCenter));
                        bool sameSegPairDir = bnd1SegPairDir.dot(bnd2SegPairDir) > 0;

                        double weight = (bndsamples[bnd1].size() + bndsamples[bnd2].size()) * spnum2weightRatio;
                        int bndContFC = fg.addFactorCategory(ml::FactorGraph::FactorCategory{ [weight, sameSegPairDir](const int * varlabels, size_t nvar,
                            ml::FactorGraph::FactorCategoryId fcid, void * givenData) -> double {
                            assert(nvar == 2);
                            DepthRelation rel1 = (DepthRelation)varlabels[0];
                            DepthRelation rel2 = (DepthRelation)varlabels[1];
                            if (rel1 == DepthRelation::Connected && rel2 == DepthRelation::Connected ||
                                rel1 == DepthRelation::Disconnected && rel2 == DepthRelation::Disconnected ||
                                rel1 == DepthRelation::Disconnected && Contains({ DepthRelation::FirstIsFront, DepthRelation::SecondIsFront }, rel2) ||
                                rel2 == DepthRelation::Disconnected && Contains({ DepthRelation::FirstIsFront, DepthRelation::SecondIsFront }, rel1)) {
                                return 0.0;
                            }
                            if (sameSegPairDir) {
                                if (rel1 == DepthRelation::FirstIsFront && rel2 == DepthRelation::FirstIsFront ||
                                    rel1 == DepthRelation::SecondIsFront && rel2 == DepthRelation::SecondIsFront) {
                                    return 0.0;
                                }
                            } else {
                                if (rel1 == DepthRelation::FirstIsFront && rel2 == DepthRelation::SecondIsFront ||
                                    rel1 == DepthRelation::SecondIsFront && rel2 == DepthRelation::FirstIsFront) {
                                    return 0.0;
                                }
                            }
                            return 5.0 * weight;
                        }, 1.0 });
                        fg.addFactor({ vhs[bnd1], vhs[bnd2] }, bndContFC);
                    }
                }
            }



            std::cout << "solving factor graph" << std::endl;

            ml::FactorGraph::ResultTable resultLabels;
            fg.solve(100, 10, [&resultLabels](int epoch, double energy, double denergy, const ml::FactorGraph::ResultTable & results) -> bool {
                std::cout << "epoch: " << epoch << "\t energy: " << energy << std::endl;
                if (denergy < 0.1) {
                    resultLabels = results;
                    return true;
                }
                return false;
            });
            
            auto occlusions = mg.createConstraintTable<RegionBoundaryData>(DepthRelation::Connected);
            for (int bndid = 0; bndid < segtopo.nboundaries(); bndid++) {
                if (vhs[bndid].invalid()) {
                    continue;
                }
                auto & bh = bhs[bndid];
                occlusions[bh] = (DepthRelation)resultLabels[vhs[bndid]];
            }

            return occlusions;
        }




        HandledTable<RegionBoundaryHandle, DepthRelation> DetectOcclusions4(
            const RLGraph & mg, const HandledTable<RegionHandle, OrientationControl> & controls,
            const Imagei & segs,
            const SegmentationTopo & segtopo, const std::vector<std::vector<Vec3>> & bndsamples, const std::vector<int> & bndclasses,
            const std::vector<RegionHandle> & rhs, const std::vector<RegionBoundaryHandle> & bhs,
            const std::vector<Vec3> & vps, double angleDistThres, double angleSampleStepOnLine) {

            // register all sample points           
            struct BndSample {
                int bndid;
                int sampleid;
                bool operator == (const BndSample & bs) const { return std::tie(bndid, sampleid) == std::tie(bs.bndid, bs.sampleid); }
                bool operator < (const BndSample & bs) const { return std::tie(bndid, sampleid) < std::tie(bs.bndid, bs.sampleid); }
            };
            RTreeMap<Vec3, BndSample> samples;
            for (int bndid = 0; bndid < bndsamples.size(); bndid++) {
                auto & segids = segtopo.bnd2segs[bndid];
                auto rh1 = rhs[segids.first];
                auto rh2 = rhs[segids.second];
                if (rh1.invalid() || rh2.invalid()) {
                    continue;
                }
                for (int sampleid = 0; sampleid < bndsamples[bndid].size(); sampleid++) {
                    samples.insert(std::make_pair(bndsamples[bndid][sampleid], BndSample{ bndid, sampleid }));
                }
            }

            // find closest bnd sample pairs to locate edge bnds along narrow segs
            // segidpair -> [bndid -> ratio]
            std::map<std::pair<int, int>, std::map<int, double>> segpair2bnds;
            std::map<std::pair<int, int>, int> bndPairCorrespSampleNums;
            {
                std::map<BndSample, BndSample> closestBs;
                for (int bndid = 0; bndid < bndsamples.size(); bndid++) {
                    for (int sampleid = 0; sampleid < bndsamples[bndid].size(); sampleid++) {
                        auto & segids = segtopo.bnd2segs[bndid];
                        auto rh1 = rhs[segids.first];
                        auto rh2 = rhs[segids.second];
                        if (rh1.invalid() || rh2.invalid()) {
                            continue;
                        }
                        if (std::tie(controls[rh1].orientationClaz, controls[rh1].orientationNotClaz) ==
                            std::tie(controls[rh2].orientationClaz, controls[rh2].orientationNotClaz)) {
                            continue;
                        }
                        if (!NoOrientationControl(controls[rh1]) && !NoOrientationControl(controls[rh2])) {
                            continue; // only consider detached edges
                        }

                        auto controledSide = controls[rh1];
                        if (NoOrientationControl(controledSide)) {
                            controledSide = controls[rh2];
                        }
                        assert(!NoOrientationControl(controledSide));

                        auto & dir = bndsamples[bndid][sampleid];

                        BndSample closest = { -1, -1 };
                        double minAngle = angleDistThres;

                        samples.search(BoundingBox(normalize(dir)).expand(angleDistThres * 2),
                            [&closest, &minAngle, &segtopo, &rhs, &controls, &controledSide, &dir](const std::pair<Vec3, BndSample> & bs) {

                            auto & segids = segtopo.bnd2segs[bs.second.bndid];
                            auto rh1 = rhs[segids.first];
                            auto rh2 = rhs[segids.second];

                            if (std::tie(controls[rh1].orientationClaz, controls[rh1].orientationNotClaz) ==
                                std::tie(controls[rh2].orientationClaz, controls[rh2].orientationNotClaz)) {
                                return true;
                            }
                            if (!NoOrientationControl(controls[rh1]) && !NoOrientationControl(controls[rh2])) {
                                return true; // only consider detached edges
                            }

                            auto thisControledSide = controls[rh1];
                            if (NoOrientationControl(thisControledSide)) {
                                thisControledSide = controls[rh2];
                            }
                            assert(!NoOrientationControl(thisControledSide));

                            if (std::tie(thisControledSide.orientationClaz, thisControledSide.orientationNotClaz) ==
                                std::tie(controledSide.orientationClaz, controledSide.orientationNotClaz)) {
                                return true;
                            }

                            double angle = AngleBetweenDirections(dir, bs.first);
                            if (angle < minAngle) {
                                closest = bs.second;
                                minAngle = angle;
                            }

                            return true;
                        });

                        if (closest.bndid != -1) {
                            closestBs[BndSample{ bndid, sampleid }] = closest;
                        }
                    }
                }

                std::vector<std::pair<BndSample, BndSample>> closestPairs;
                for (auto & p : closestBs) {
                    if (closestBs.at(p.second) == p.first) {
                        closestPairs.emplace_back(p.first, p.second);
                    }
                }
                for (auto & samplePair : closestPairs) {
                    bndPairCorrespSampleNums[MakeOrderedPair(samplePair.first.bndid, samplePair.second.bndid)] ++;
                }


                for (int bndid = 0; bndid < bndsamples.size(); bndid++) {
                    auto & segids = segtopo.bnd2segs[bndid];
                    auto rh1 = rhs[segids.first];
                    auto rh2 = rhs[segids.second];
                    if (rh1.invalid() || rh2.invalid()) {
                        continue;
                    }
                    if (std::tie(controls[rh1].orientationClaz, controls[rh1].orientationNotClaz) ==
                        std::tie(controls[rh2].orientationClaz, controls[rh2].orientationNotClaz)) {
                        continue;
                    }
                    if (!NoOrientationControl(controls[rh1]) && !NoOrientationControl(controls[rh2])) {
                        auto orderedSegIds = segids.first < segids.second ? segids : std::make_pair(segids.second, segids.first);
                        segpair2bnds[orderedSegIds][bndid] = 1.0;
                    }
                }
                for (auto & p : closestPairs) {
                    int bndid1 = p.first.bndid;
                    int bndid2 = p.second.bndid;
                    assert(bndid1 != p.second.bndid);

                    int constrainedSegId1 = segtopo.bnd2segs[bndid1].first;
                    if (NoOrientationControl(controls[rhs[segtopo.bnd2segs[bndid1].first]])) {
                        constrainedSegId1 = segtopo.bnd2segs[bndid1].second;
                    }

                    int constrainedSegId2 = segtopo.bnd2segs[bndid2].first;
                    if (NoOrientationControl(controls[rhs[segtopo.bnd2segs[bndid2].first]])) {
                        constrainedSegId2 = segtopo.bnd2segs[bndid2].second;
                    }

                    assert(!NoOrientationControl(controls[rhs[constrainedSegId1]]) &&
                        !NoOrientationControl(controls[rhs[constrainedSegId2]]));

                    auto & bnds =
                        segpair2bnds[constrainedSegId1 < constrainedSegId2 ?
                        std::make_pair(constrainedSegId1, constrainedSegId2) :
                        std::make_pair(constrainedSegId2, constrainedSegId2)];

                    bnds[bndid1] += 1.0 / bndsamples[bndid1].size();
                    bnds[bndid2] += 1.0 / bndsamples[bndid2].size();
                }
            }

            // get first guess
            auto firstGuess = mg.createConstraintTable<RegionBoundaryData>(DepthRelation::Unknown);
            for (auto & segpairbnd : segpair2bnds) {
                const auto & segids = segpairbnd.first;
                auto & bnds = segpairbnd.second;

                auto rh1 = rhs[segids.first];
                auto rh2 = rhs[segids.second];
                if (rh1.invalid() || rh2.invalid()) {
                    continue;
                }
                for (auto & bndwithscore : bnds) {
                    int bndid = bndwithscore.first;

                    double occupationRatio = bndwithscore.second;
                    if (occupationRatio < 0.5) {
                        continue;
                    }

                    auto bh = bhs[bndid];
                    if (bh.invalid()) {
                        continue;
                    }

                    if (!controls[rh1].used || !controls[rh2].used) {
                        continue;
                    }

                    int bndclass = bndclasses[bndid];

                    if (controls[rh1].orientationClaz != -1 && controls[rh2].orientationClaz != -1
                        && controls[rh1].orientationClaz != controls[rh2].orientationClaz
                        && (bndclass == -1
                        || !IsFuzzyPerpendicular(vps[bndclass], vps[controls[rh1].orientationClaz])
                        || !IsFuzzyPerpendicular(vps[bndclass], vps[controls[rh2].orientationClaz]))) {
                        // break this
                        firstGuess[bh] = DepthRelation::Disconnected;
                        auto & rhvp1 = vps[controls[rh1].orientationClaz];
                        auto & rhvp2 = vps[controls[rh2].orientationClaz];
                        if (bndclass != -1) {
                            auto & bndvp = vps[bndclass];
                            if (IsFuzzyPerpendicular(rhvp1, bndvp) && !IsFuzzyPerpendicular(rhvp2, bndvp)) {
                                firstGuess[bh] = DepthRelation::FirstIsFront;
                            } else if (!IsFuzzyPerpendicular(rhvp1, bndvp) && IsFuzzyPerpendicular(rhvp2, bndvp)) {
                                firstGuess[bh] = DepthRelation::SecondIsFront;
                            }
                        }
                    } else if (controls[rh1].orientationClaz != -1 &&
                        controls[rh2].orientationClaz == -1 && controls[rh2].orientationNotClaz == controls[rh1].orientationClaz
                        && (bndclass == -1
                        || !IsFuzzyPerpendicular(vps[bndclass], vps[controls[rh1].orientationClaz]))) {
                        // break this
                        firstGuess[bh] = DepthRelation::Disconnected;
                        if (bndclass == controls[rh1].orientationClaz) {
                            firstGuess[bh] = DepthRelation::SecondIsFront;
                        }
                    } else if (controls[rh2].orientationClaz != -1 &&
                        controls[rh1].orientationClaz == -1 && controls[rh1].orientationNotClaz == controls[rh2].orientationClaz
                        && (bndclass == -1
                        || !IsFuzzyPerpendicular(vps[bndclass], vps[controls[rh2].orientationClaz]))) {
                        // break this
                        firstGuess[bh] = DepthRelation::Disconnected;
                        if (bndclass == controls[rh2].orientationClaz) {
                            firstGuess[bh] = DepthRelation::FirstIsFront;
                        }
                    } else if (controls[rh1].orientationClaz != -1 && controls[rh2].orientationClaz != -1 &&
                        controls[rh1].orientationClaz == controls[rh2].orientationClaz) {
                        firstGuess[bh] = DepthRelation::Connected;
                    } else if (controls[rh1].orientationClaz == -1 && controls[rh1].orientationNotClaz != -1
                        && controls[rh2].orientationClaz == -1 && controls[rh2].orientationNotClaz != -1
                        && controls[rh1].orientationNotClaz == controls[rh2].orientationNotClaz
                        && bndclass != -1 && bndclass != controls[rh1].orientationNotClaz) {
                        firstGuess[bh] = DepthRelation::Connected;
                    }
                }
            }

            auto bndsAttachedToLines = mg.createComponentTable<LineData, std::vector<int>>();
            {
                for (auto & l : mg.components<LineData>()) {
                    std::vector<std::set<int>> bndSamplesNearLines(segtopo.nboundaries());
                    auto & line = l.data.line;
                    double angle = AngleBetweenDirections(line.first, line.second);
                    if (angle < DegreesToRadians(1)) {
                        continue;
                    }
                    double angleStep = angleSampleStepOnLine;
                    for (double a = 0.0; a <= angle; a += angleStep) {
                        Vec3 dir = RotateDirection(line.first, line.second, a);
                        BndSample nearest = { -1, -1 };
                        double minAngleDist = M_PI;

                        // find nearest bndSample
                        samples.search(BoundingBox(normalize(dir)).expand(angleDistThres * 2),
                            [&bndsamples, &dir, &nearest, &minAngleDist](const std::pair<Vec3, BndSample> & bs) {
                            auto & thisDir = bndsamples[bs.second.bndid][bs.second.sampleid];
                            double angleDist = AngleBetweenDirections(thisDir, dir);
                            if (minAngleDist > angleDist) {
                                minAngleDist = angleDist;
                                nearest = bs.second;
                            }
                            return true;
                        });

                        if (nearest.bndid == -1) {
                            continue;
                        }
                        bndSamplesNearLines[nearest.bndid].insert(nearest.sampleid);
                    }

                    auto & bndsAttached = bndsAttachedToLines[l.topo.hd];
                    for (int bndid = 0; bndid < segtopo.nboundaries(); bndid++) {
                        if (bndSamplesNearLines[bndid].size() > bndsamples[bndid].size() * 0.6) {
                            bndsAttached.push_back(bndid);
                        }
                    }
                }
            }


            std::cout << "building factor graph" << std::endl;

            // factor graph for bnds
            ml::FactorGraph fg;
            static const double spnum2weightRatio = 1.0 / 100.0;

            // add bnd as vars
            std::vector<ml::FactorGraph::VarHandle> vhs(segtopo.nboundaries());
            for (int bndid = 0; bndid < segtopo.nboundaries(); bndid++) {
                if (bhs[bndid].invalid()) {
                    continue;
                }
                vhs[bndid] = fg.addVar(fg.addVarCategory(3, bndsamples[bndid].size() / 10.0)); // first front, second front, connected
            }

            // first guess factor
            for (int bndid = 0; bndid < segtopo.nboundaries(); bndid++) {
                if (vhs[bndid].invalid()) {
                    continue;
                }
                auto guess = firstGuess[bhs[bndid]];
                size_t spnum = bndsamples[bndid].size();
                int unaryFC = fg.addFactorCategory(ml::FactorGraph::FactorCategory{ [guess, spnum](const int * varlabels, size_t nvar,
                    ml::FactorGraph::FactorCategoryId fcid, void * givenData) -> double {
                    double weight = spnum * spnum2weightRatio;
                    assert(nvar == 1);
                    DepthRelation label = static_cast<DepthRelation>(varlabels[0]);
                    switch (guess) {
                    case DepthRelation::Unknown: return label == DepthRelation::Connected ? 0.0 : 1.0 * weight;  // always prefer connected
                    case DepthRelation::Connected: return label == DepthRelation::Connected ? 0.0 : 3.0 * weight;
                    case DepthRelation::Disconnected:
                        return Contains({
                            DepthRelation::Disconnected, DepthRelation::FirstIsFront, DepthRelation::SecondIsFront
                        }, label) ? 0.0 : 10.0 * weight;
                    case DepthRelation::FirstIsFront: return label == DepthRelation::FirstIsFront ? 0.0 : 20.0 * weight;
                    case DepthRelation::SecondIsFront: return label == DepthRelation::SecondIsFront ? 0.0 : 20.0 * weight;
                    default: SHOULD_NEVER_BE_CALLED("this label should never occur in first guess!");
                    }
                }, 10.0 });
                fg.addFactor({ vhs[bndid] }, unaryFC);
            }

            // overlap factor
            for (auto & bndPair : bndPairCorrespSampleNums) {
                if (bhs[bndPair.first.first].invalid() || bhs[bndPair.first.second].invalid()) {
                    continue;
                }
                int overlapSPNum = bndPair.second;
                double weight = spnum2weightRatio * overlapSPNum;
                int overlapFC = fg.addFactorCategory(ml::FactorGraph::FactorCategory{ [weight](const int * varlabels, size_t nvar,
                    ml::FactorGraph::FactorCategoryId fcid, void * givenData) -> double {
                    assert(nvar == 2);
                    DepthRelation label1 = static_cast<DepthRelation>(varlabels[0]);
                    DepthRelation label2 = static_cast<DepthRelation>(varlabels[1]);
                    if (label1 != DepthRelation::Connected && label2 != DepthRelation::Connected) {
                        return 10.0 * weight;
                    }
                    return 0.0;
                }, 1.0 });
                fg.addFactor({ vhs[bndPair.first.first], vhs[bndPair.first.second] }, overlapFC);
            }

            // junction validity factor
            for (int jid = 0; jid < segtopo.njunctions(); jid++) {
                std::set<int> bndids(segtopo.junc2bnds[jid].begin(), segtopo.junc2bnds[jid].end());
                if (bndids.size() < 3) {
                    continue;
                }
                if (std::any_of(bndids.begin(), bndids.end(), [&bhs](int bndid) {return bhs[bndid].invalid(); })) {
                    continue;
                }
                double weight = 0.0;
                for (int bndid : bndids) {
                    weight += spnum2weightRatio * bndsamples[bndid].size();
                }
                int juncValidFC = fg.addFactorCategory(ml::FactorGraph::FactorCategory{ [weight](const int * varlabels, size_t nvar,
                    ml::FactorGraph::FactorCategoryId fcid, void * givenData) -> double {
                    assert(nvar >= 3);
                    int disconNum = 0;
                    for (int i = 0; i < nvar; i++) {
                        if (DepthRelation(varlabels[i]) != DepthRelation::Connected) {
                            disconNum++;
                        }
                    }
                    if (disconNum >= 3) {
                        return 10.0 * weight;
                    }
                    return 0.0;
                }, 1.0 });
                std::vector<ml::FactorGraph::VarHandle> relatedVHs(bndids.size());
                auto it = bndids.begin();
                for (int i = 0; i < bndids.size(); i++) {
                    relatedVHs[i] = vhs[*it];
                    ++it;
                }
                fg.addFactor(relatedVHs.begin(), relatedVHs.end(), juncValidFC);
            }

            for (auto & bnds : bndsAttachedToLines) {
                if (bnds.size() <= 1) {
                    continue;
                }
                for (int i = 0; i < bnds.size(); i++) {
                    for (int j = i + 1; (j < i + 9 && j < bnds.size()); j++) {
                        int bnd1 = bnds[i];
                        int bnd2 = bnds[j];
                        if (bhs[bnd1].invalid() || bhs[bnd2].invalid()) {
                            continue;
                        }

                        auto & segs1 = segtopo.bnd2segs[bnd1];
                        auto & bh1 = bhs[bnd1];
                        auto & bh2 = bhs[bnd2];
                        Vec3 bnd1SegPairDir = normalize(mg.data(mg.topo(bh1).component<0>()).normalizedCenter.cross(mg.data(mg.topo(bh1).component<1>()).normalizedCenter));
                        auto & segs2 = segtopo.bnd2segs[bnd2];
                        Vec3 bnd2SegPairDir = normalize(mg.data(mg.topo(bh2).component<0>()).normalizedCenter.cross(mg.data(mg.topo(bh2).component<1>()).normalizedCenter));
                        bool sameSegPairDir = bnd1SegPairDir.dot(bnd2SegPairDir) > 0;

                        double weight = (bndsamples[bnd1].size() + bndsamples[bnd2].size()) * spnum2weightRatio;
                        int bndContFC = fg.addFactorCategory(ml::FactorGraph::FactorCategory{ [weight, sameSegPairDir](const int * varlabels, size_t nvar,
                            ml::FactorGraph::FactorCategoryId fcid, void * givenData) -> double {
                            assert(nvar == 2);
                            DepthRelation rel1 = (DepthRelation)varlabels[0];
                            DepthRelation rel2 = (DepthRelation)varlabels[1];
                            if (rel1 == DepthRelation::Connected && rel2 == DepthRelation::Connected ||
                                rel1 == DepthRelation::Disconnected && rel2 == DepthRelation::Disconnected ||
                                rel1 == DepthRelation::Disconnected && Contains({ DepthRelation::FirstIsFront, DepthRelation::SecondIsFront }, rel2) ||
                                rel2 == DepthRelation::Disconnected && Contains({ DepthRelation::FirstIsFront, DepthRelation::SecondIsFront }, rel1)) {
                                return 0.0;
                            }
                            if (sameSegPairDir) {
                                if (rel1 == DepthRelation::FirstIsFront && rel2 == DepthRelation::FirstIsFront ||
                                    rel1 == DepthRelation::SecondIsFront && rel2 == DepthRelation::SecondIsFront) {
                                    return 0.0;
                                }
                            } else {
                                if (rel1 == DepthRelation::FirstIsFront && rel2 == DepthRelation::SecondIsFront ||
                                    rel1 == DepthRelation::SecondIsFront && rel2 == DepthRelation::FirstIsFront) {
                                    return 0.0;
                                }
                            }
                            return 5.0 * weight;
                        }, 1.0 });
                        fg.addFactor({ vhs[bnd1], vhs[bnd2] }, bndContFC);
                    }
                }
            }



            std::cout << "solving factor graph" << std::endl;

            ml::FactorGraph::ResultTable resultLabels;
            fg.solve(100, 10, [&resultLabels](int epoch, double energy, double denergy, const ml::FactorGraph::ResultTable & results) -> bool {
                std::cout << "epoch: " << epoch << "\t energy: " << energy << std::endl;
                if (energy == 0.0) {
                    resultLabels = results;
                    return false;
                }
                if (denergy < 0.1) {
                    resultLabels = results;
                    return true;
                }
                return false;
            });

            auto occlusions = mg.createConstraintTable<RegionBoundaryData>(DepthRelation::Connected);
            for (int bndid = 0; bndid < segtopo.nboundaries(); bndid++) {
                if (vhs[bndid].invalid()) {
                    continue;
                }
                auto & bh = bhs[bndid];
                occlusions[bh] = (DepthRelation)resultLabels[vhs[bndid]];
            }

            return occlusions;

        }




        void ApplyOcclusions(const RLGraph & mg, RLGraphControls & controls,
            const HandledTable<RegionBoundaryHandle, DepthRelation> & occlusions) {

            for (auto it = occlusions.begin(); it != occlusions.end(); ++it) {

                auto bh = it.hd();
                DepthRelation relation = *it;

                if (true) {
                    size_t spnum = 0;
                    for (auto & sps : mg.data(bh).normalizedSampledPoints) {
                        spnum += sps.size();
                    }
                    if (spnum <= 1) {
                        controls[bh].used = false;
                    }
                }

                if (relation == DepthRelation::Connected || relation == DepthRelation::MaybeFolder) {
                    continue;
                }
                if (bh.invalid()) {
                    continue;
                }
                controls[bh].used = false;

                // the front side region
                auto rh1 = mg.topo(bh).component<0>();
                auto rh2 = mg.topo(bh).component<1>();

                // related lh ?
                for (auto & l : mg.components<LineData>()) {
                    bool hasRh1 = false, hasRh2 = false;
                    RegionLineConnectionHandle rlh1, rlh2;
                    for (RegionLineConnectionHandle rlh : l.topo.constraints<RegionLineConnectionData>()) {
                        if (mg.topo(rlh).component<0>() == rh1) {
                            hasRh1 = true;
                            rlh1 = rlh;
                        }
                        if (mg.topo(rlh).component<0>() == rh2) {
                            hasRh2 = true;
                            rlh2 = rlh;
                        }
                        if (hasRh1 && hasRh2) {
                            break;
                        }
                    }
                    if (hasRh1 && hasRh2) {
                        if (relation == DepthRelation::FirstIsFront) {
                            controls[rlh2].used = false;
                            controls[rlh1].used = true;
                        } else if(relation == DepthRelation::SecondIsFront) {
                            controls[rlh1].used = false;
                            controls[rlh2].used = true;
                        } /*else if (relation == DepthRelation::Disconnected) {
                            controls[rlh1].used = false;
                            controls[rlh2].used = false;
                        } */else {
                            SHOULD_NEVER_BE_CALLED();
                        }

                        // spread to all rhs related to this lh
                        if (true) {
                            Vec3 lineNormal = normalize(l.data.line.first.cross(l.data.line.second));
                            double dot1 = mg.data(rh1).normalizedCenter.dot(lineNormal);
                            double dot2 = mg.data(rh2).normalizedCenter.dot(lineNormal);
                            if (dot1 * dot2 >= 0) {
                                continue; // confused case...
                            }
                            for (RegionLineConnectionHandle rlh : l.topo.constraints<RegionLineConnectionData>()) {
                                double dotHere = mg.data(mg.topo(rlh).component<0>()).normalizedCenter.dot(lineNormal);
                                bool tendToRh1 = (dot1 * dotHere >= 0);
                                if (tendToRh1 && relation != DepthRelation::FirstIsFront || !tendToRh1 && relation == DepthRelation::FirstIsFront) {
                                    controls[rlh].used = false;
                                }
                            }
                        }
                    }
                }
            }

        }




        std::vector<TStructure> FindTStructures(const SegmentationTopo & segtopo,
            const PanoramicCamera & cam, const std::vector<Vec3> & vps, 
            double minSpanAngle, double angleThres){

            std::vector<TStructure> tstructs;
            
            struct BndData {
                std::vector<Vec3> edge;
                std::vector<Vec3> sampledEdge;
            };
            std::vector<BndData> bndData(segtopo.nboundaries());
            // compute bnd data
            static const double samplingStepAngleOnBoundary = 0.01;
            for (int i = 0; i < segtopo.nboundaries(); i++) {
                auto & bd = bndData[i];
                bd.edge.reserve(segtopo.bndpixels[i].size());
                for (auto & p : segtopo.bndpixels[i]) {
                    bd.edge.push_back(normalize(cam.toSpace(p)));
                }
                bd.sampledEdge = { bd.edge.front() };
                for (auto & p : bd.edge) {
                    double remainedAngle = AngleBetweenDirections(bd.sampledEdge.back(), p);
                    while (remainedAngle >= samplingStepAngleOnBoundary) {
                        bd.sampledEdge.push_back(normalize(RotateDirection(bd.sampledEdge.back(), p, samplingStepAngleOnBoundary)));
                        remainedAngle -= samplingStepAngleOnBoundary;
                    }
                }
            }

            
            // for each junction
            for (int jid = 0; jid < segtopo.njunctions(); jid++) {
                auto bndids = segtopo.junc2bnds[jid];
                Vec3 center = normalize(cam.toSpace(segtopo.juncpositions[jid]));

                if (bndids.size() != 3) {
                    continue;
                }
  
                std::vector<std::vector<int>> candvpids(bndids.size());
                for (int i = 0; i < bndids.size(); i++) {
                    auto & bndid = bndids[i];
                    
                    // do they follow any vp direction ?
                    for (int j = 0; j < vps.size(); j++) {
                        auto & vp = vps[j];
                        if (AllAlong(bndData[bndid].sampledEdge, center, vp, angleThres)) {
                            candvpids[i].push_back(j);
                        }
                    }
                }

                if (std::any_of(candvpids.begin(), candvpids.end(), 
                    [](const std::vector<int> & vpids) {return vpids.empty(); })) {
                    continue;
                }

                // find a T structure
                std::vector<int> vpids;
                for (int i1 : candvpids[0]) {
                    if (!vpids.empty())
                        break;
                    for (int i2 : candvpids[1]) {
                        if (!vpids.empty())
                            break;
                        for (int i3 : candvpids[2]) {
                            if (std::set<int>({i1, i2, i3}).size() == 2) {
                                vpids = { i1, i2, i3 };
                                break;
                            }
                        }
                    }
                }

                if (vpids.empty()) {
                    continue;
                }

                // extend the bnds
                std::vector<std::vector<int>> extendedJuncIds(bndids.size());
                std::vector<std::vector<int>> extendedBndIds(bndids.size());
                for (int i = 0; i < bndids.size(); i++) {
                    int curjunc = jid;
                    int curbnd = bndids[i];
                    const auto & vp = vps[vpids[i]];

                    while (true) {
                        extendedJuncIds[i].push_back(curjunc);
                        extendedBndIds[i].push_back(curbnd);

                        // next junc
                        int nextjunc = segtopo.bnd2juncs[curbnd].first;
                        if (nextjunc == curjunc) {
                            nextjunc = segtopo.bnd2juncs[curbnd].second;
                        }

                        // next bnd
                        bool hasnext = false;
                        for (int nextcandbnd : segtopo.junc2bnds[nextjunc]) {
                            if (Contains(extendedBndIds[i], nextcandbnd))
                                continue;
                            auto & bd = bndData[nextcandbnd];
                            if (AllAlong(bd.sampledEdge, center, vp, angleThres)) {
                                curjunc = nextjunc;
                                curbnd = nextcandbnd;
                                hasnext = true;
                                break;
                            }
                        }
                        if (!hasnext) {
                            break;
                        }
                    }
                }


                // reorder bndids as [short, long, lont]
                if (vpids[0] == vpids[1]) {
                    extendedJuncIds = { std::move(extendedJuncIds[2]), std::move(extendedJuncIds[0]), std::move(extendedJuncIds[1]) };
                    extendedBndIds = { std::move(extendedBndIds[2]), std::move(extendedBndIds[0]), std::move(extendedBndIds[1]) };
                    bndids = { bndids[2], bndids[0], bndids[1] };
                    vpids = { vpids[2], vpids[0], vpids[1] };
                } else if (vpids[1] == vpids[2]) {
                    // remain the same
                } else if (vpids[2] == vpids[0]) {
                    extendedJuncIds = { std::move(extendedJuncIds[1]), std::move(extendedJuncIds[0]), std::move(extendedJuncIds[2]) };
                    extendedBndIds = { std::move(extendedBndIds[1]), std::move(extendedBndIds[0]), std::move(extendedBndIds[2]) };
                    bndids = { bndids[1], bndids[0], bndids[2] };
                    vpids = { vpids[1], vpids[0], vpids[2] };
                }

                assert(vpids[1] == vpids[2]);
                assert((vpids[0] ^ vpids[1] ^ vpids[2]) == vpids[0]);

                TStructure ts;
                ts.centerJunctionId = jid;
                ts.shortBndIds = std::move(extendedBndIds[0]);
                ts.longBndIds[0] = std::move(extendedBndIds[1]);
                ts.longBndIds[1] = std::move(extendedBndIds[2]);
                ts.shortVPId = vpids[0];
                ts.longVPId = vpids[1];

                ts.shortEndJunctionId = extendedJuncIds[0].back();
                ts.longEndJunctionIds[0] = extendedJuncIds[1].back();
                ts.longEndJunctionIds[1] = extendedJuncIds[2].back();

                if (std::set<int>({ ts.shortEndJunctionId, ts.longEndJunctionIds[0], ts.longEndJunctionIds[1], ts.centerJunctionId }).size() < 4) {
                    continue;
                }

                {
                    Vec3 center = normalize(cam.toSpace(segtopo.juncpositions[ts.centerJunctionId]));
                    Vec3 longEnds[] = {
                        normalize(cam.toSpace(segtopo.juncpositions[ts.longEndJunctionIds[0]])),
                        normalize(cam.toSpace(segtopo.juncpositions[ts.longEndJunctionIds[1]]))
                    };
                    Vec3 shortEnd = normalize(cam.toSpace(segtopo.juncpositions[ts.shortEndJunctionId]));

                    double angleSpan = std::max(AngleBetweenDirections(shortEnd, longEnds[0]), AngleBetweenDirections(shortEnd, longEnds[1]));
                    if (angleSpan < minSpanAngle) // too small!
                        continue;
                }

                tstructs.push_back(std::move(ts));
            }

            return tstructs;
        }




        std::vector<TStructure> FindTStructuresFuzzy(const SegmentationTopo & segtopo, 
            const std::vector<std::vector<Vec3>> & bndsamples,
            const std::vector<int> & bndclasses,
            const PanoramicCamera & cam, const std::vector<Vec3> & vps, 
            double minSpanAngle, double angleThres /*= DegreesToRadians(1.5)*/) {
            double testSpanAngle = minSpanAngle / 2.0;

            std::vector<TStructure> tstructs;

            // for each junction
            for (int jid = 0; jid < segtopo.njunctions(); jid++) {
                auto & jcenter = segtopo.juncpositions[jid];
                Vec3 jcenterdir = normalize(cam.toSpace(jcenter));
                auto bndids = segtopo.junc2bnds[jid];

                struct Branch {
                    std::vector<int> bnds;
                    std::vector<Vec3> samples; // not all
                    int endJuncId;
                };
                std::map<int, std::vector<Branch>> goodBranches;

                struct BndJunc {
                    int bndid;
                    int juncid; // the geometric order is: [juncid -> bndid] -> [next juncid -> next bndid] -> ...
                    double juncAngleDist;
                };
                std::deque<BndJunc> stack;
                std::map<int, int> bndPrevPosInStack;

                for (int bndid : bndids) {
                    stack.push_front(BndJunc{ bndid, jid, 0.0 });
                    bndPrevPosInStack[bndid] = -1;
                }

                // dfs to get good bnds
                while (!stack.empty()) {

                    BndJunc curbj = stack.front();
                    stack.pop_front();

                    // get next BndJuncs
                    // next juncid
                    int nextjuncid = segtopo.bnd2juncs[curbj.bndid].first;
                    if (nextjuncid == curbj.juncid) {
                        nextjuncid = segtopo.bnd2juncs[curbj.bndid].second;
                    }
                    // next juncAngleDist
                    double nextJuncAngleDist = 
                        AngleBetweenDirections(jcenterdir, normalize(cam.toSpace(segtopo.juncpositions[nextjuncid])));
                    if (nextJuncAngleDist <= curbj.juncAngleDist) {
                        // cur bnd turns back! it's invalid!
                        // todo
                        continue;
                    }

                    if (nextJuncAngleDist > testSpanAngle) {
                        // it's sufficiently large now
                        // get the path using .prevPosInStack
                        int allSamplePointCount = 0;
                        int classifiedSamplePointCount = 0;
                        std::vector<int> bndsOnPath;
                        std::vector<int> bndsClassifiedOnPath;
                        int bndiditer = curbj.bndid;
                        assert(Contains(bndPrevPosInStack, bndiditer));
                        while (bndiditer != -1) {
                            allSamplePointCount += bndsamples[bndiditer].size();
                            if (bndclasses[bndiditer] != -1) {
                                classifiedSamplePointCount += bndsamples[bndiditer].size();
                                bndsClassifiedOnPath.push_back(bndiditer);
                            }
                            bndsOnPath.push_back(bndiditer);
                            bndiditer = bndPrevPosInStack.at(bndiditer);
                        }
                        std::reverse(bndsOnPath.begin(), bndsOnPath.end());
                        std::reverse(bndsClassifiedOnPath.begin(), bndsClassifiedOnPath.end());

                        if (bndsClassifiedOnPath.empty()) {
                            continue;
                        }
                        if (classifiedSamplePointCount < allSamplePointCount * 0.8) {
                            continue;
                        }

                        assert(segtopo.bnd2juncs[bndsOnPath.front()].first == jid ||
                            segtopo.bnd2juncs[bndsOnPath.front()].second == jid);

                        std::vector<std::vector<Vec3>> samples(vps.size());
                        for (int b : bndsClassifiedOnPath) {
                            auto & ss = samples[bndclasses[b]];
                            ss.insert(ss.end(), bndsamples[b].begin(), bndsamples[b].end());
                        } 
                        for (int i = 0; i < vps.size(); i++) {
                            if (samples[i].size() < allSamplePointCount * 0.5) {
                                continue;
                            }
                            if (AllAlong(samples[i], jcenterdir, vps[i], angleThres)) {
                                goodBranches[i].push_back(Branch{ std::move(bndsOnPath), std::move(samples[i]), nextjuncid });
                                break;
                            }
                        }

                        continue;
                    }

                    // next bndid
                    auto & nextBndCands = segtopo.junc2bnds[nextjuncid];
                    for (int nextBndCand : nextBndCands) {
                        if (nextBndCand == curbj.bndid) {
                            continue;
                        }
                        int nextjuncidcand = segtopo.bnd2juncs[nextBndCand].first;
                        if (nextjuncidcand == nextjuncid) {
                            nextjuncidcand = segtopo.bnd2juncs[nextBndCand].second;
                        }
                        double ad = AngleBetweenDirections(jcenterdir, normalize(cam.toSpace(segtopo.juncpositions[nextjuncidcand])));
                        if (ad <= curbj.juncAngleDist) {
                            continue;
                        }
                        stack.push_front(BndJunc{ nextBndCand, nextjuncid, nextJuncAngleDist });
                        bndPrevPosInStack[nextBndCand] = curbj.bndid;
                    }

                }

                if (goodBranches.size() != 2) { // the T structure
                    continue;
                }


                // extend the bnds
                for (auto & branches : goodBranches) {
                    int vpid = branches.first;
                    for (auto & branch : branches.second) {
                        assert(segtopo.bnd2juncs[branch.bnds.front()].first == jid || segtopo.bnd2juncs[branch.bnds.front()].second == jid);
                        int endBnd = branch.bnds.back();
                        auto & endJuncCands = segtopo.bnd2juncs[endBnd];
                        double endJuncCandAngleDists[2] = {
                            AngleBetweenDirections(jcenterdir, cam.toSpace(segtopo.juncpositions[endJuncCands.first])),
                            AngleBetweenDirections(jcenterdir, cam.toSpace(segtopo.juncpositions[endJuncCands.second]))
                        };
                        int endJunc = endJuncCandAngleDists[0] < endJuncCandAngleDists[1] ? endJuncCands.second : endJuncCands.first;

                        // try extending this branch
                        int lastBnd = endBnd;
                        int lastJunc = endJunc;
                        double lastJuncAngleDist = std::max(endJuncCandAngleDists[0], endJuncCandAngleDists[1]);
                        while (true) {
                            auto & nextBndCands = segtopo.junc2bnds[lastJunc];
                            bool hasNext = false;
                            for (int nextBndCand : nextBndCands) {
                                if (nextBndCand == lastBnd) {
                                    continue;
                                }
                                int nextJunc = segtopo.bnd2juncs[nextBndCand].first;
                                if (nextJunc == lastJunc) {
                                    nextJunc = segtopo.bnd2juncs[nextBndCand].second;
                                }
                                double nextJuncAngleDist = AngleBetweenDirections(jcenterdir, cam.toSpace(segtopo.juncpositions[nextJunc]));
                                if (nextJuncAngleDist <= lastJuncAngleDist) {
                                    continue;
                                }
                                auto & bndCandSamples = bndsamples[nextBndCand];
                                if (AllAlong(std::vector<std::vector<Vec3>>{branch.samples, bndCandSamples}, 
                                    jcenterdir, vps[vpid], angleThres)) {
                                    branch.bnds.push_back(nextBndCand);
                                    branch.samples.insert(branch.samples.end(), bndCandSamples.begin(), bndCandSamples.end());
                                    branch.endJuncId = nextJunc;
                                    
                                    lastBnd = nextBndCand;
                                    lastJunc = nextJunc;
                                    lastJuncAngleDist = nextJuncAngleDist;
                                    
                                    hasNext = true;
                                    break;
                                }
                            }
                            if (!hasNext) {
                                break;
                            }
                        }
                    }
                }


                std::vector<std::pair<int, Branch>> allBranches;
                for (auto & bs : goodBranches) {
                    for (auto & branch : bs.second) {
                        allBranches.emplace_back(bs.first, std::move(branch));
                    }
                }
                if (allBranches.size() < 3) {
                    continue;
                }

                std::vector<int> allBranchOrderedIds(allBranches.size());
                std::iota(allBranchOrderedIds.begin(), allBranchOrderedIds.end(), 0);
                std::sort(allBranchOrderedIds.begin(), allBranchOrderedIds.end(),
                    [&jcenterdir, &segtopo, &cam, &allBranches](int a, int b) {
                    return AngleBetweenDirections(jcenterdir, cam.toSpace(segtopo.juncpositions[allBranches[a].second.endJuncId]))
                        > AngleBetweenDirections(jcenterdir, cam.toSpace(segtopo.juncpositions[allBranches[b].second.endJuncId]));
                });
                
                std::vector<int> selectedBranchIds = { allBranchOrderedIds[0], allBranchOrderedIds[1] };
                if (allBranches[allBranchOrderedIds[0]].first == allBranches[allBranchOrderedIds[1]].first) {
                    for (int i = 2; i < allBranchOrderedIds.size(); i++) {
                        if (allBranches[allBranchOrderedIds[i]].first != allBranches[allBranchOrderedIds[0]].first) {
                            selectedBranchIds.push_back(allBranchOrderedIds[i]);
                            break;
                        }
                    }
                    std::swap(selectedBranchIds[0], selectedBranchIds[2]); // short long long
                } else {
                    for (int i = 2; i < allBranchOrderedIds.size(); i++) {
                        if (allBranches[allBranchOrderedIds[i]].first == allBranches[allBranchOrderedIds[0]].first) {
                            selectedBranchIds.push_back(allBranchOrderedIds[i]);
                            std::swap(selectedBranchIds[0], selectedBranchIds[1]);
                            break;
                        } else if (allBranches[allBranchOrderedIds[i]].first == allBranches[allBranchOrderedIds[1]].first) {
                            selectedBranchIds.push_back(allBranchOrderedIds[i]);
                            break;
                        }
                    }
                }
                assert(selectedBranchIds.size() == 3);
                std::set<int> vpidsset;
                for (auto & i : selectedBranchIds) {
                    vpidsset.insert(allBranches[i].first);
                }
                assert(vpidsset.size() == 2);


                TStructure ts;
                ts.centerJunctionId = jid;
                ts.shortBndIds = std::move(allBranches[selectedBranchIds[0]].second.bnds);
                ts.longBndIds[0] = std::move(allBranches[selectedBranchIds[1]].second.bnds);
                ts.longBndIds[1] = std::move(allBranches[selectedBranchIds[2]].second.bnds);
                ts.shortVPId = allBranches[selectedBranchIds[0]].first;
                ts.longVPId = allBranches[selectedBranchIds[1]].first;
                assert(allBranches[selectedBranchIds[1]].first == allBranches[selectedBranchIds[2]].first);

                ts.shortEndJunctionId = allBranches[selectedBranchIds[0]].second.endJuncId;
                ts.longEndJunctionIds[0] = allBranches[selectedBranchIds[1]].second.endJuncId;
                ts.longEndJunctionIds[1] = allBranches[selectedBranchIds[2]].second.endJuncId;

                if (std::set<int>({ ts.shortEndJunctionId, ts.longEndJunctionIds[0], ts.longEndJunctionIds[1], ts.centerJunctionId }).size() < 4) {
                    continue;
                }

                {
                    Vec3 center = normalize(cam.toSpace(segtopo.juncpositions[ts.centerJunctionId]));
                    Vec3 longEnds[] = {
                        normalize(cam.toSpace(segtopo.juncpositions[ts.longEndJunctionIds[0]])),
                        normalize(cam.toSpace(segtopo.juncpositions[ts.longEndJunctionIds[1]]))
                    };
                    Vec3 shortEnd = normalize(cam.toSpace(segtopo.juncpositions[ts.shortEndJunctionId]));

                    double angleSpan = std::max(AngleBetweenDirections(shortEnd, longEnds[0]), AngleBetweenDirections(shortEnd, longEnds[1]));
                    if (angleSpan < minSpanAngle) // too small!
                        continue;
                }

                tstructs.push_back(std::move(ts));

            }

            return tstructs;

        }




        namespace {
            template <class Vec3ArrayT>
            void MaskViewFromTriangle(const Vec3ArrayT & triangle, View<PartialPanoramicCamera, Imageub> & maskView, double focal) {
                Vec3 triangleCenter = normalize(std::accumulate(std::begin(triangle), std::end(triangle), Vec3()));
                double radiusAngle = 0;
                for (auto & c : triangle) {
                    radiusAngle = std::max(radiusAngle, AngleBetweenDirections(triangleCenter, c));
                }

                int ppcSize = std::ceil(2 * radiusAngle * focal);
                Vec3 x;
                std::tie(x, std::ignore) = ProposeXYDirectionsFromZDirection(triangleCenter);
                PartialPanoramicCamera ppc(ppcSize, ppcSize, focal, Point3(0, 0, 0), triangleCenter, x);

                Imageub mask = Imageub::zeros(ppc.screenSize());
                std::vector<Point2i> triangleProj(3);
                for (int k = 0; k < 3; k++) {
                    triangleProj[k] = ppc.toScreen(triangle[k]);
                }
                cv::fillPoly(mask, std::vector<std::vector<Point2i>>{triangleProj}, (uint8_t)1);

                maskView.camera = ppc;
                maskView.image = std::move(mask);
            }
        }


        void TStructureMaskView(const PanoramicCamera & cam,
            const SegmentationTopo & segtopo, const TStructure & ts,
            std::array<View<PartialPanoramicCamera, Imageub>, 2> & smallParts,
            View<PartialPanoramicCamera, Imageub> & largePart,
            double focal, double largePartWidthAngle) {

            Vec3 center = normalize(cam.toSpace(segtopo.juncpositions[ts.centerJunctionId]));
            Vec3 longEnds[] = {
                normalize(cam.toSpace(segtopo.juncpositions[ts.longEndJunctionIds[0]])),
                normalize(cam.toSpace(segtopo.juncpositions[ts.longEndJunctionIds[1]]))
            };
            Vec3 shortEnd = normalize(cam.toSpace(segtopo.juncpositions[ts.shortEndJunctionId]));

            // small parts
            for (int i = 0; i < 2; i++) {
                Vec3 triangle[] = { center, longEnds[i], shortEnd };
                MaskViewFromTriangle(triangle, smallParts[i], focal);
            }

            // large part
            Vec3 backBone = normalize(RotateDirection(center, shortEnd, -largePartWidthAngle));
            Vec3 triangle[] = { longEnds[0], longEnds[1], backBone };
            MaskViewFromTriangle(triangle, largePart, focal);
        }




        std::vector<int> GuessOccludedTStructures(const std::vector<Vec3> & vps,
            const SegmentationTopo & segtopo, const std::vector<TStructure> & tstructs, const PanoramicCamera & cam,
            const Image5d & gc) {
            
            std::vector<int> ots;
            
            int vertVPId = NearestDirectionId(vps);

            for (int i = 0; i < tstructs.size(); i++) {
                auto & ts = tstructs[i];

           /*     auto canvas = gui::MakeCanvas(Image3ub(gc.size(), Vec3ub()));
                ShowTStructures(canvas.color(gui::White).thickness(2), segtopo, { ts });
                auto im = canvas.image();*/

                std::array<View<PartialPanoramicCamera, Imageub>, 2> smallPartMaskViews;
                View<PartialPanoramicCamera, Imageub> largePartMaskView;
                TStructureMaskView(cam, segtopo, ts, smallPartMaskViews, largePartMaskView, 50, DegreesToRadians(5));

                auto gcMeanSmall1 = MeanInMask(MakeView(gc, cam), smallPartMaskViews[0]);
                auto gcMeanSmall2 = MeanInMask(MakeView(gc, cam), smallPartMaskViews[1]);
                auto gcMeanLarge = MeanInMask(MakeView(gc, cam), largePartMaskView);

                int shortVPId = ts.shortVPId;
                int longVPId = ts.longVPId;
                assert(shortVPId != longVPId);

                if (ts.longVPId == vertVPId) {
                    /*  GeometricContextIndex gciSmall1 = MaxGeometricIndex(gcMeanSmall1);
                      GeometricContextIndex gciSmall2 = MaxGeometricIndex(gcMeanSmall2);
                      GeometricContextIndex gciLarge = MaxGeometricIndex(gcMeanLarge);
                      if ((gciSmall1 == GeometricContextIndex::FloorOrGround || gciSmall2 == GeometricContextIndex::FloorOrGround) &&
                      gciLarge == GeometricContextIndex::Vertical) {*/
                        ots.push_back(i);
                    //}
                }
            }

            return ots;
        }


        void ApplyOccludedTStructure(const RLGraph & mg, RLGraphControls & controls, 
            const SegmentationTopo & segtopo,
            const std::vector<RegionBoundaryHandle> & bhs,
            const PanoramicCamera & cam, const std::vector<TStructure> & occtstructs) {

            int i = 0;
            for (auto & ts : occtstructs) {
                if (i == 14) {
                    std::cout << std::endl;
                }
                i++;

                Vec3 center = normalize(cam.toSpace(segtopo.juncpositions[ts.centerJunctionId]));
                Vec3 longEnds[] = {
                    normalize(cam.toSpace(segtopo.juncpositions[ts.longEndJunctionIds[0]])),
                    normalize(cam.toSpace(segtopo.juncpositions[ts.longEndJunctionIds[1]]))
                };
                Vec3 shortEnd = normalize(cam.toSpace(segtopo.juncpositions[ts.shortEndJunctionId]));
                Vec3 longEdgeNormal = normalize(longEnds[0].cross(longEnds[1]));


                for (auto & bndids : ts.longBndIds) {
                    for (int bndid : bndids) {
                        auto & bh = bhs[bndid];
                        if (bh.invalid()) {
                            continue;
                        }
                        controls[bh].used = false;

                        // the front side region
                        auto rh1 = mg.topo(bh).component<0>();
                        auto rh2 = mg.topo(bh).component<1>();
                        auto & rh1center = mg.data(rh1).normalizedCenter;
                        auto & rh2center = mg.data(rh2).normalizedCenter;

                        // which lies at the short end side from center?
                        double dot1 = rh1center.dot(longEdgeNormal);
                        double dot2 = rh2center.dot(longEdgeNormal);

                        double dotShortEnd = shortEnd.dot(longEdgeNormal);

                        // rh1 | rh2
                        bool rh1isFront = true;
                        if (dotShortEnd > 0) {
                            rh1isFront = dot1 < dot2;
                        } else {
                            rh1isFront = dot1 > dot2;
                        }

                        // related lh ?
                        for (auto & l : mg.components<LineData>()) {
                            bool hasRh1 = false, hasRh2 = false;
                            RegionLineConnectionHandle rlh1, rlh2;
                            for (RegionLineConnectionHandle rlh : l.topo.constraints<RegionLineConnectionData>()) {
                                if (mg.topo(rlh).component<0>() == rh1) {
                                    hasRh1 = true;
                                    rlh1 = rlh;
                                }
                                if (mg.topo(rlh).component<0>() == rh2) {
                                    hasRh2 = true;
                                    rlh2 = rlh;
                                }
                                if (hasRh1 && hasRh2) {
                                    break;
                                }
                            }
                            if (hasRh1 && hasRh2) {
                                if (rh1isFront) {
                                    controls[rlh2].used = false;
                                    controls[rlh1].used = true;
                                } else {
                                    controls[rlh1].used = false;
                                    controls[rlh2].used = true;
                                }

                                // spread to all rhs related to this lh
                                if (false) { ///// todo
                                    Vec3 lineNormal = normalize(l.data.line.first.cross(l.data.line.second));
                                    double dot1 = mg.data(rh1).normalizedCenter.dot(longEdgeNormal);
                                    double dot2 = mg.data(rh2).normalizedCenter.dot(longEdgeNormal);
                                    if (dot1 * dot2 >= 0) {
                                        continue; // confused case...
                                    }
                                    for (RegionLineConnectionHandle rlh : l.topo.constraints<RegionLineConnectionData>()) {
                                        double dotHere = mg.data(mg.topo(rlh).component<0>()).normalizedCenter.dot(longEdgeNormal);
                                        bool tendToRh1 = (dot1 * dotHere >= 0);
                                        if (tendToRh1 && !rh1isFront || !tendToRh1 && rh1isFront) {
                                            controls[rlh].used = false;
                                        }
                                    }
                                }
                            }
                        }
                    }
                }

            }

        }


        void ShowTStructures(gui::Canvas3ub & canvas,
            const SegmentationTopo & segtopo, const std::vector<TStructure> & occtstructs){

            auto ctable = gui::CreateRandomColorTableWithSize(occtstructs.size());
            int i = 0;
            for (auto & ts : occtstructs) {
                auto & c = ctable[i];
                canvas.color(c);
                for (auto & bid : ts.shortBndIds) {
                    canvas.add(segtopo.bndpixels[bid]);
                }
                for (auto & bids : ts.longBndIds) {
                    for (auto & bid : bids) {
                        canvas.add(segtopo.bndpixels[bid]);
                    }
                }
                canvas.add(NoteAs(segtopo.juncpositions[ts.centerJunctionId], std::to_string(i)));
                i++;
            }

        }

      

    }
}