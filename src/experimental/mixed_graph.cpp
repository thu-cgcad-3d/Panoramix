#include "../core/algorithms.hpp"
#include "../core/utilities.hpp"
#include "../core/containers.hpp"
#include "../ml/factor_graph.hpp"
#include "../gui/visualize2d.hpp"
#include "../gui/visualizers.hpp"
#include "mixed_graph.hpp"

namespace panoramix {
    namespace experimental {


        MixedGraph::MixedGraph(const Image & image,
            const Point2 & cameraCenterPosition,
            double cameraFocal) {

            core::LineSegmentExtractor lineExtractor;
            lineExtractor.params().algorithm = core::LineSegmentExtractor::LSD;
            lines = core::ClassifyEachAs(lineExtractor(image, 1), -1);

            view.camera = core::PerspectiveCamera(image.cols, image.rows, cameraCenterPosition, cameraFocal,
                core::Point3(0, 0, 0), core::Point3(1, 0, 0), core::Point3(0, 0, -1));
            view.image = image;

            // classify lines
            vps = core::EstimateVanishingPointsAndClassifyLines(view.camera, lines);
            // remove non-manhattan vps
            vps.erase(vps.begin() + 3, vps.end());
            for (auto & l : lines){
                if (l.claz >= 3){
                    l.claz = -1;
                }
            }

            // append lines
            core::AppendLines(regionLineGraph, lines, view.camera, vps);

            // append regions
            core::SegmentationExtractor segmenter;
            segmenter.params().algorithm = core::SegmentationExtractor::GraphCut;
            segmenter.params().sigma = 0.5;
            segmenter.params().c = 30.0;
            segmenter.params().superpixelSizeSuggestion = -1;
            segmenter.params().superpixelNumberSuggestion = 500;
            std::vector<core::Line2> pureLines(lines.size());
            for (int i = 0; i < lines.size(); i++)
                pureLines[i] = lines[i].component;
            int segmentsNum = 0;
            std::tie(regions, segmentsNum) = segmenter(image, pureLines, image.cols / 100.0);

            regionIds2Handles = core::AppendRegions(regionLineGraph, regions, view.camera, 0.001, 0.001, 3, 1);
            

            // boundary junctions
            // add region boundary handle
            auto bjunctions = ExtractBoundaryJunctions(regions);
            for (auto & junc : bjunctions){
                auto & regionIds = junc.first;
                auto & ps = junc.second;
                std::set<RegionHandle> rhs;
                for (int regionId : regionIds){
                    RegionHandle rh = regionIds2Handles[regionId];
                    if (rh.invalid())
                        continue;
                    rhs.insert(rh);
                }
                if (rhs.size() < 3)
                    continue;
                // locate boundary handles among rhs
                std::set<RegionBoundaryHandle> bhs;
                for (RegionHandle rh : rhs){
                    auto & relatedBhs = regionLineGraph.topo(rh).constraints<RegionBoundaryData>();
                    for (RegionBoundaryHandle bh : relatedBhs){
                        auto anotherRh = regionLineGraph.topo(bh).component<0>();
                        if (anotherRh == rh){
                            anotherRh = regionLineGraph.topo(bh).component<1>();
                        }
                        if (Contains(rhs, anotherRh)){
                            bhs.insert(bh);
                        }
                    }
                }
                std::vector<Vec3> positions(ps.size());
                for (int i = 0; i < ps.size(); i++){
                    positions[i] = normalize(view.camera.spatialDirection(ps[i]));
                }
                boundaryJunctions.push_back(BoundaryJunction{ std::move(bhs), std::move(positions) });
            }

        }

        void MixedGraph::installGCResponse(const Image7d & gc){
            assert(gc.size() == regions.size());
            gcResponse = regionLineGraph.createComponentTable<RegionData>(Vec7());
            auto regionAreas = regionLineGraph.createComponentTable<RegionData>(0.0);
            for (auto it = gc.begin(); it != gc.end(); ++it){
                auto p = it.pos();
                const Vec7 & v = gc(p);
                RegionHandle rh(regions(p));
                if (rh.invalid())
                    continue;
                gcResponse[rh] += v;
                regionAreas[rh] += 1.0;
            }
            for (int i = 0; i < gcResponse.data.size(); i++){
                gcResponse.data[i] /= regionAreas.data[i];
                double sum = std::accumulate(
                    std::begin(gcResponse.data[i].val), 
                    std::end(gcResponse.data[i].val), 
                    0.0);
                assert(IsFuzzyZero(sum - 1.0, 0.1));
            }
        }


        double Distance(const Point2 & p, const std::vector<PixelLoc> & edge){
            double mind = std::numeric_limits<double>::max();
            for (const auto & c : edge){
                auto cc = core::vec_cast<double>(c);
                double d = core::Distance(p, cc);
                if (d < mind){
                    mind = d;
                }
            }
            return mind;
        }


        void MixedGraph::installOcclusionResponce(const std::vector<std::vector<PixelLoc>> & edges,
            const std::vector<double> & scores){
            assert(edges.size() == scores.size());
            occlusionResponse = regionLineGraph.createConstraintTable<RegionBoundaryData>(0.0);
            auto getBB = [&edges](int eid)->Box2 {
                return BoundingBoxOfContainer(edges.at(eid));
            };
            std::vector<int> eids(edges.size());
            std::iota(eids.begin(), eids.end(), 0);
            RTreeWrapper<int, decltype(getBB)> rtree(eids.begin(), eids.end(), getBB);
            for (auto & bd : regionLineGraph.constraints<RegionBoundaryData>()){
                double & resp = occlusionResponse[bd.topo.hd];
                double allSampleScoreSum = 0;
                int allSampleNum = 0;
                for (auto & ps : bd.data.normalizedEdges){
                    for (auto & p : ps){
                        allSampleNum++;
                        auto sample = view.camera.screenProjection(p);
                        double sampleScore = 0.0;
                        int nearbyEdgeNum = 0;
                        rtree.search(Box2(sample, sample).expand(5.0), 
                            [&scores, &sampleScore, &nearbyEdgeNum, &edges, &sample](int eid) -> bool {
                            double dist = Distance(sample, edges[eid]);
                            if (dist <= 3.0){
                                sampleScore += scores[eid] > 0 ? 1.0 : (scores[eid] == 0.0 ? 0.0 : -1.0);
                                nearbyEdgeNum++;
                            }
                            return true;
                        });
                        sampleScore /= std::max(nearbyEdgeNum, 1);
                        allSampleScoreSum += sampleScore;
                    }
                }
                occlusionResponse[bd.topo.hd] = allSampleScoreSum / std::max(allSampleNum, 1);
            }
        }

        void MixedGraph::installOcclusionResponce(const std::vector<std::vector<int>> & edges,
            const std::vector<double> & scores){
            assert(edges.size() == scores.size());
            std::vector<std::vector<PixelLoc>> pixels(edges.size());
            for (int i = 0; i < edges.size(); i++){
                auto & ps = pixels[i];
                ps.reserve(edges[i].size());
                for (int id : edges[i]){
                    ps.push_back(DecodeIndexToSubscript(id, regions.size()));
                }
            }
            installOcclusionResponce(pixels, scores);
        }

        void MixedGraph::showGCResponse() const{
            
        }

        void MixedGraph::showOcclusionResponse() const{
            gui::Visualizer2D vis(view.image);
            for (auto & b : regionLineGraph.constraints<RegionBoundaryData>()){
                double score = occlusionResponse[b.topo.hd];
                if (score <= 0)
                    continue;
                vis.params.color = gui::Red;
                vis.params.thickness = 1;
                for (auto & edge : b.data.normalizedEdges){
                    for (int i = 0; i < edge.size() - 1; i++){
                        vis << Line2(view.camera.screenProjection(edge[i]), view.camera.screenProjection(edge[i + 1]));
                    }
                }
            }
            vis << gui::manip2d::Show();
        }

        void MixedGraph::showBoundaryJunctions() const{
            int regionsNum = MinMaxValOfImage(regions).second + 1;
            auto ctable = gui::CreateRandomColorTableWithSize(regionsNum);
            auto segs = ctable(regions);
            std::vector<Point2> ps;
            for (auto & j : boundaryJunctions){
                for (auto & d : j.positions){
                    ps.push_back(view.camera.screenProjection(d));
                }
            }
            gui::Visualizer2D(segs)
                << gui::manip2d::SetColor(gui::Blue)
                << ps
                << gui::manip2d::Show();
        }

        void MixedGraph::showDetachableRegionLineConnections() const{
            int regionsNum = MinMaxValOfImage(regions).second + 1;
            auto ctable = gui::CreateRandomColorTableWithSize(regionsNum);
            auto segs = ctable(regions);
            std::vector<Line2> lines;
            lines.reserve(regionLineGraph.internalComponents<LineData>().size());
            for (auto & l : regionLineGraph.components<LineData>()){
                Line2 line(view.camera.screenProjection(l.data.line.first), view.camera.screenProjection(l.data.line.second));
                lines.push_back(line);
            }
            std::vector<Line2> detachableConnections, undetachableConnections;
            for (auto & rl : regionLineGraph.constraints<RegionLineConnectionData>()){
                auto rh = regionLineGraph.topo(rl.topo.hd).component<0>();
                auto lh = regionLineGraph.topo(rl.topo.hd).component<1>();
                auto & rcenter = regionLineGraph.data(rh).normalizedCenter;
                Point2 rcenterProj = view.camera.screenProjection(rcenter);
                std::vector<Point2> anchorProjs;
                anchorProjs.reserve(rl.data.normalizedAnchors.size());
                for (auto & a : rl.data.normalizedAnchors){
                    anchorProjs.push_back(view.camera.screenProjection(a));
                }
                Point2 last = anchorProjs.front();
                for (int i = 1; i < anchorProjs.size(); i++){
                    double dist = core::Distance(last, anchorProjs[i]);
                    if (dist < 5) {
                        continue;
                    }
                    (rl.data.detachable ? detachableConnections : undetachableConnections).emplace_back(rcenterProj, anchorProjs[i]);
                    last = anchorProjs[i];
                }
            }
            gui::Visualizer2D(segs)               
                << gui::manip2d::SetThickness(1) << gui::manip2d::SetColor(gui::Black) << detachableConnections
                << gui::manip2d::SetColor(gui::White) << undetachableConnections
                << gui::manip2d::SetThickness(1) << gui::manip2d::SetColor(gui::Yellow) << lines
                << gui::manip2d::Show();
        }


        void MixedGraph::showSegmentations() const{
            int regionsNum = MinMaxValOfImage(regions).second + 1;
            auto ctable = gui::CreateRandomColorTableWithSize(regionsNum);
            auto segs = ctable(regions);
            gui::Visualizer2D(segs) << gui::manip2d::Show();
        }

        void MixedGraph::showVanishingPointsAndLines() const{
            
        }

        void MixedGraph::showRegionOrientationConstriants() const{

        }

        namespace {

            //int ToRegionOrientationFlag(int orientationClaz, int orientatinNotClaz, int vpnum, int vertVPId){
            //    assert(orientationClaz >= -1 && orientationClaz < vpnum);
            //    assert(orientatinNotClaz >= -1 && orientatinNotClaz < vpnum);
            //    if (orientationClaz != -1){
            //        return orientationClaz;
            //    }
            //    else if(orientatinNotClaz != -1){
            //        assert(orientatinNotClaz == vertVPid); 
            //        IMPROVABLE_HERE("donnot support constraints that are along other directions crrently");
            //        return vpnum;
            //    }
            //    else{
            //        return vpnum + 1;
            //    }
            //}

            //void FromRegionOrientationFlag(int flag, int & orientationClaz, int & orientatinNotClaz, int vpnum, int vertVPId) {
            //    assert(flag >= 0 && flag <= vpnum + 1);
            //    if (flag < vpnum){
            //        orientationClaz = flag;
            //        orientatinNotClaz = -1;
            //    }
            //    else if (flag == vpnum){
            //        orientationClaz = -1;
            //        orientatinNotClaz = vertVPId;
            //    }
            //    else if (flag == vpnum + 1){
            //        orientationClaz = orientatinNotClaz = -1;
            //    }
            //}

            //inline int ToLineOrientationFlag(int orientationClaz){
            //    return orientationClaz + 1;
            //}

            //inline void FromLineOrientationFlag(int flag, int & orientationClaz){
            //    orientationClaz = flag - 1;
            //}

        }

        HandledTable<RegionHandle, Plane3> MixedGraph::solve() const {           
            RLGraphPropertyTable props = core::MakeRLGraphPropertyTable(regionLineGraph, vps);

            size_t vpnum = vps.size();
            int vertVPId = GetVerticalDirectionId(vps);            

            // some precalculated data
            auto regionIntersectWithVPs = regionLineGraph.createComponentTable<RegionData, int>(-1);
            for (int i = 0; i < vpnum; i++){
                for (auto & p : { view.camera.screenProjection(vps[i]), view.camera.screenProjection(-vps[i]) }){
                    PixelLoc pixel = ToPixelLoc(p);
                    int s = 1;
                    for (int x = -s; x <= s; x++){
                        for (int y = -s; y <= s; y++){
                            PixelLoc pp(pixel.x + x, pixel.y + y);
                            if (Contains(regions, pixel)){
                                regionIntersectWithVPs[regionIds2Handles[regions(pp)]] = i;
                            }
                        }
                    }                   
                }
            }
            auto regionIntersectWithHorizon = regionLineGraph.createComponentTable<RegionData, int8_t>(false);
            for (auto & r : regionLineGraph.components<RegionData>()){
                auto h = r.topo.hd;
                auto & contours = r.data.normalizedContours;
                bool intersected = false;
                for (auto & cs : r.data.normalizedContours){
                    if (intersected)
                        break;
                    for (auto & c : cs){
                        double angle = M_PI_2 - AngleBetweenUndirectedVectors(c, vps[vertVPId]);
                        if (angle <= M_PI / 100.0){
                            intersected = true;
                            break;
                        }
                    }
                }
                if (intersected){
                    regionIntersectWithHorizon[h] = true;
                }
            }

            auto lineLeftDetachableRegionConnections = 
                regionLineGraph.createComponentTable<LineData, std::vector<RegionLineConnectionHandle>>();
            auto lineRightDetachableRegionConnections = 
                regionLineGraph.createComponentTable<LineData, std::vector<RegionLineConnectionHandle>>();
            for (auto & rl : regionLineGraph.constraints<RegionLineConnectionData>()){
                if (!rl.data.detachable)
                    continue;
                RegionHandle rh = rl.topo.component<0>();
                LineHandle lh = rl.topo.component<1>();
                const Line3 & lineProj = regionLineGraph.data(lh).line;
                Vec3 rightDir = normalize(lineProj.first.cross(lineProj.second));
                const Vec3 & regionCenter = regionLineGraph.data(rh).normalizedCenter;
                bool isOnRight = (regionCenter - lineProj.center()).dot(rightDir);
                (isOnRight ? lineRightDetachableRegionConnections : lineLeftDetachableRegionConnections)[lh].push_back(rl.topo.hd);
            }


            size_t regionNum = regionLineGraph.internalComponents<RegionData>().size();
            size_t lineNum = regionLineGraph.internalComponents<LineData>().size();
            size_t boundaryNum = regionLineGraph.internalConstraints<RegionBoundaryData>().size();
            size_t rlConNum = regionLineGraph.internalConstraints<RegionLineConnectionData>().size();
            size_t llConNum = regionLineGraph.internalConstraints<LineRelationData>().size();
            size_t bjuncNum = boundaryJunctions.size();

            // create factor graph
            ml::FactorGraph fg;
            fg.reserveVarCategories(5);
            fg.reserveVars(regionNum + lineNum + boundaryNum + rlConNum + llConNum);

            /// add vars

            // add region orientation constraint flags
            auto regionOrientationVc = fg.addVarCategory(ml::FactorGraph::VarCategory{ 3, 0.5 }); // 3: {no constraints, vertical, horizontal}
            auto regionVhs = regionLineGraph.createComponentTable<RegionData, ml::FactorGraph::VarHandle>();
            for (auto & r : regionLineGraph.components<RegionData>()){
                regionVhs[r.topo.hd] = fg.addVar(regionOrientationVc);
            }

            // add line orientation constraint flags
            auto lineOrientationVc = fg.addVarCategory(ml::FactorGraph::VarCategory{ 2, 0.5 }); // 2: {no constraints, oriented}
            auto lineVhs = regionLineGraph.createComponentTable<LineData, ml::FactorGraph::VarHandle>();
            for (auto & l : regionLineGraph.components<LineData>()){
                lineVhs[l.topo.hd] = fg.addVar(lineOrientationVc);
            }
            
            // add boundary flags
            auto boundaryOcclusionVc = fg.addVarCategory(ml::FactorGraph::VarCategory{ 2, 0.5 }); // 2: {not connected, connected}
            auto boundaryVhs = regionLineGraph.createConstraintTable<RegionBoundaryData, ml::FactorGraph::VarHandle>();
            for (auto & b : regionLineGraph.constraints<RegionBoundaryData>()){
                boundaryVhs[b.topo.hd] = fg.addVar(boundaryOcclusionVc);
            }

            // add rl connection flags
            auto rlConVc = fg.addVarCategory(ml::FactorGraph::VarCategory{ 2, 0.5 }); // 2: {not connected, connected}
            auto rlConVhs = regionLineGraph.createConstraintTable<RegionLineConnectionData, ml::FactorGraph::VarHandle>();
            for (auto & c : regionLineGraph.constraints<RegionLineConnectionData>()){
                rlConVhs[c.topo.hd] = fg.addVar(rlConVc);
            }

            // add ll connection flags
            auto llConVc = fg.addVarCategory(ml::FactorGraph::VarCategory{ 2, 0.5 }); // 2: {not connected, connected}
            auto llConVhs = regionLineGraph.createConstraintTable<LineRelationData, ml::FactorGraph::VarHandle>();
            for (auto & r : regionLineGraph.constraints<LineRelationData>()){
                llConVhs[r.topo.hd] = fg.addVar(llConVc);
            }

            /// add factors

            // add region orientation unary cost
            for (auto & r : regionLineGraph.components<RegionData>()){
                auto rh = r.topo.hd;
                bool intersectWithHorizon = regionIntersectWithHorizon[rh];
                bool intersectWithVertical = regionIntersectWithVPs[rh] == vertVPId;
                auto regionOrientationFc =
                    fg.addFactorCategory(
                    ml::FactorGraph::FactorCategory{ [this, rh, &props, vertVPId,
                    intersectWithHorizon, intersectWithVertical](const int * labels, size_t nvars) -> double {
                    assert(nvars == 1);
                    int rolabel = labels[0];
                    
                    // 0: front, 1: left, 2: right, 3: floor, 4: ceiling, 5: clutter, 6: unknown
                    const Vec7 & gcResp = gcResponse[rh]; 

                    double gcVertResp = std::max({ gcResp[0], gcResp[1], gcResp[2] });
                    double gcHoriResp = std::max({ gcResp[3], gcResp[4] });
                    assert(IsBetween(gcVertResp, 0.0, 1.1));
                    assert(IsBetween(gcHoriResp, 0.0, 1.1));
                    
                    if (rolabel == 0){ // no constraints
                        return 0.5;
                    }
                    else if (rolabel == 1){ // vertical
                        // consistency with gcVertResp
                        // intersect horizon ?
                        // consistency with current props
                        Plane3 plane = Instance(regionLineGraph, props, rh);
                        bool maybeVertical = M_PI_2 - AngleBetweenUndirectedVectors(plane.normal, vps[vertVPId]) < M_PI_4;
                        // todo
                        return (1.0 - gcVertResp) * (intersectWithHorizon ? 0.1 : 1.0) * (maybeVertical ? 1.0 : 2.0);
                    }
                    else if (rolabel == 2){ // horizontal
                        // consistency with gcHoriResp
                        // todo
                        Plane3 plane = Instance(regionLineGraph, props, rh);
                        bool maybeHorizontal = AngleBetweenUndirectedVectors(plane.normal, vps[vertVPId]) < M_PI_4;
                        // todo
                        return (1.0 - gcHoriResp) * (intersectWithVertical ? 0.1 : 1.0) * (maybeHorizontal ? 1.0 : 2.0);
                    }
                    assert(false);
                    return 1e5;
                }, 0.5 });
                fg.addFactor({ regionVhs[rh] }, regionOrientationFc);
            }

            // add line orientation unary cost
            for (auto & l : regionLineGraph.components<LineData>()){
                auto lh = l.topo.hd;
                auto lineOrientationFc = fg.addFactorCategory(
                    ml::FactorGraph::FactorCategory{ [this, lh, &props](const int * labels, size_t nvars) -> double{
                    assert(nvars == 1);
                    int lolabel = labels[0];
                    
                    auto & ld = regionLineGraph.data(lh);
                    Line3 line = Instance(regionLineGraph, props, lh);
                    // distortion
                    double angle = AngleBetweenUndirectedVectors(line.direction(), line.center());
                    double distortion = Gaussian(angle, M_PI_4 / 4.0);
                    if (lolabel == 0){ // no constraints
                        return 1.0 - distortion;
                    }
                    else if (lolabel == 1){ // abey the orientation class
                        return distortion;
                    }
                    assert(false);
                    return 1e5;
                }, 0.5 });
                fg.addFactor({ lineVhs[lh] }, lineOrientationFc);
            }

            // add boundary flag unary cost
            for (auto & b : regionLineGraph.constraints<RegionBoundaryData>()){
                auto bh = b.topo.hd;
                auto boundaryFlagFc = fg.addFactorCategory(
                    ml::FactorGraph::FactorCategory{ [this, bh](const int * labels, size_t nvars) -> double{
                    assert(nvars == 1);
                    const double & occResp = occlusionResponse[bh];
                    
                    int flag = labels[0];
                    if (flag == 0){ // occlusion
                        return occResp < 0 ? 1.0 : 0.0;
                    }
                    else { // not occlusion
                        return occResp > 0 ? 1.0 : 0.0;
                    }

                }, 0.5 });
                fg.addFactor({ boundaryVhs[bh] }, boundaryFlagFc);
            }


            // add rl connection flag unary cost
            for (auto & c : regionLineGraph.constraints<RegionLineConnectionData>()){
                auto ch = c.topo.hd;
                auto rlConFlagFc = fg.addFactorCategory(
                    ml::FactorGraph::FactorCategory{ [this, ch](const int * labels, size_t nvars) -> double{
                    assert(nvars == 1);
                    // todo
                    return 0.0;
                }, 0.5 });
                fg.addFactor({ rlConVhs[ch] }, rlConFlagFc);
            }

            // add line-line connection flag unary cost
            for (auto & r : regionLineGraph.constraints<LineRelationData>()){
                auto rh = r.topo.hd;
                auto linelineRelFlagFc = fg.addFactorCategory(
                    ml::FactorGraph::FactorCategory{ [this, rh](const int * labels, size_t nvars) -> double{
                    assert(nvars == 1);
                    const auto & rd = regionLineGraph.data(rh);
                    // junctionweight, type ...
                    // todo
                    int flag = labels[0];
                    if (flag == 0){ // not connected
                        return rd.junctionWeight / 10.0;
                    }
                    else{ // connected
                        return 1.0 - rd.junctionWeight / 10.0;
                    }
                }, 0.5 });
                fg.addFactor({ llConVhs[rh] }, linelineRelFlagFc);
            }


            // add regionOrientation-lineOrientation-regionLineConnectionFlag consistency cost
            auto roLoRLCConsistencyCostFc = fg.addFactorCategory(
                ml::FactorGraph::FactorCategory{ [](const int * labels, size_t nvars) -> double{
                assert(nvars == 3);
                int roLabel = labels[0];
                int loLabel = labels[1];
                int rlConFlag = labels[2];
                if (roLabel == 0 || loLabel == 0){
                    return 0.0;
                }
                // todo
                
                return 0.0;
            }, 0.5 });
            for (auto & c : regionLineGraph.constraints<RegionLineConnectionData>()){
                auto ch = c.topo.hd;
                RegionHandle rh = regionLineGraph.topo(ch).component<0>();
                LineHandle lh = regionLineGraph.topo(ch).component<1>();
                fg.addFactor({ regionVhs[rh], lineVhs[lh], rlConVhs[ch] }, roLoRLCConsistencyCostFc);
            }


            // add boundary junction consistency cost
            auto boundaryJunctionConsistencyCostFc = fg.addFactorCategory(
                ml::FactorGraph::FactorCategory{ [](const int * labels, size_t nvars) -> double{
                assert(nvars >= 3);
                // judge validness
                int occlusionNum = 0;
                for (int i = 0; i < nvars; i++){
                    if (labels[i] == 0){
                        occlusionNum++;
                    }
                }
                return occlusionNum <= 2 ? 0.0 : 10.0;
            }, 0.5 });
            for (auto & bj : boundaryJunctions){
                assert(bj.bhs.size() >= 3);
                std::vector<ml::FactorGraph::VarHandle> bdVhs;
                bdVhs.reserve(bj.bhs.size());
                for (auto bh : bj.bhs){
                    bdVhs.push_back(boundaryVhs[bh]);
                }
                fg.addFactor(bdVhs.begin(), bdVhs.end(), boundaryJunctionConsistencyCostFc);
            }

            // add regionLineConnectionFlag consistency test cost on each line
            for (auto & l : regionLineGraph.components<LineData>()){
                auto lh = l.topo.hd;
                auto rlConConsistencyCostFc = fg.addFactorCategory(
                    ml::FactorGraph::FactorCategory{ [this, lh](const int * labels, size_t nvars) -> double{
                    // todo
                    return 0.0;
                }, 0.5 });
                std::vector<ml::FactorGraph::VarHandle> rlCVhs;
                rlCVhs.reserve(l.topo.constraints<RegionLineConnectionData>().size());
                for (RegionLineConnectionHandle rlch : l.topo.constraints<RegionLineConnectionData>()){
                    rlCVhs.push_back(rlConVhs[rlch]);
                }
                fg.addFactor(rlCVhs.begin(), rlCVhs.end(), rlConConsistencyCostFc);
            }


            auto results = regionLineGraph.createComponentTable<RegionData, Plane3>();
            return results;

        }

    }
}