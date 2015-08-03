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
            double testSpanAngle = minSpanAngle / 3.0;

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
                            if (samples[i].size() < allSamplePointCount * 0.7) {
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

            for (auto & ts : occtstructs) {

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
                auto & c = ctable[i++];
                canvas.color(c);
                for (auto & bid : ts.shortBndIds) {
                    canvas.add(segtopo.bndpixels[bid]);
                }
                for (auto & bids : ts.longBndIds) {
                    for (auto & bid : bids) {
                        canvas.add(segtopo.bndpixels[bid]);
                    }
                }
            }

        }

      

    }
}