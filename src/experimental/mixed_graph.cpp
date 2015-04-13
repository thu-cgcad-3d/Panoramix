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

            LineSegmentExtractor lineExtractor;
            lineExtractor.params().algorithm = core::LineSegmentExtractor::LSD;
            lines = core::ClassifyEachAs(lineExtractor(image, 1), -1);

            view.camera = core::PerspectiveCamera(image.cols, image.rows, cameraCenterPosition, cameraFocal,
                core::Point3(0, 0, 0), core::Point3(1, 0, 0), core::Point3(0, 0, -1));
            view.image = image;

            // classify lines
            vps = EstimateVanishingPointsAndClassifyLines(view.camera, lines);
            // remove non-manhattan vps
            vps.erase(vps.begin() + 3, vps.end());
            for (auto & l : lines){
                if (l.claz >= 3){
                    l.claz = -1;
                }
            }

            // append lines
            AppendLines(regionLineGraph, lines, view.camera, vps);

            // append regions
            SegmentationExtractor segmenter;
            segmenter.params().algorithm = core::SegmentationExtractor::SLIC;
            segmenter.params().sigma = 0.5;
            segmenter.params().c = 30.0;
            segmenter.params().superpixelSizeSuggestion = -1;
            segmenter.params().superpixelNumberSuggestion = 1000;
            std::vector<core::Line2> pureLines(lines.size());
            for (int i = 0; i < lines.size(); i++)
                pureLines[i] = lines[i].component;
            int segmentsNum = 0;
            std::tie(regions, segmentsNum) = segmenter(image/*, pureLines, image.cols / 100.0*/);

            regionIds2Handles = AppendRegions(regionLineGraph, regions, view.camera, 0.001, 0.001, 3, 1);
            

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
                if (bhs.size() < 3)
                    continue;
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
                RegionHandle rh = regionIds2Handles[regions(p)];
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

            inline std::tuple<int, int, int> GetVerticalFrontSidedVPIds(const std::vector<Vec3> & vps, 
                const PerspectiveCamera & cam){
                Vec3 dirs[] = { 
                    cam.upward(),
                    cam.forward(),
                    cam.rightward() 
                };
                int vids[] = { -1, - 1, -1 };
                double minAngles[] = { M_PI, M_PI, M_PI };
                for (int i = 0; i < vps.size(); i++){
                    for (int k = 0; k < 3; k++){
                        double angle = AngleBetweenUndirectedVectors(vps[i], dirs[k]);
                        if (angle < minAngles[k]){
                            minAngles[k] = angle;
                            vids[k] = i;
                        }
                    }
                }
                return std::make_tuple(vids[0], vids[1], vids[2]);
            }


        }


        HandledTable<RegionHandle, Plane3> MixedGraph::solve() const {           
            RLGraphControls controls(regionLineGraph, vps);

            size_t vpnum = vps.size();
            int vertVPId = -1, frontVPId = -1, sideVPId = -1;
            std::tie(vertVPId, frontVPId, sideVPId) = GetVerticalFrontSidedVPIds(vps, view.camera);
            assert(!AllSameInContainer({ vertVPId, frontVPId, sideVPId }));

            // some precalculated data
            
            // region data
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
            
            // line data
            auto lineScoresForVPs =
                regionLineGraph.createComponentTable<LineData>(std::vector<double>(vpnum, 0.0));
            for (auto & l : regionLineGraph.components<LineData>()){
                auto & line = l.data.line;
                for (int i = 0; i < vpnum; i++){
                    auto normal = line.first.cross(line.second);
                    double angleToVP = M_PI_2 - AngleBetweenUndirectedVectors(normal, vps[i]);
                    double distanceToVP = AngleBetweenUndirectedVectors(line.center(), vps[i]);
                    lineScoresForVPs[l.topo.hd][i] =
                        Gaussian(angleToVP, M_PI_4 / 4.0) * Pitfall(distanceToVP, M_PI_4 / 4.0);
                }
            }


            auto lineLeftRegionConnections = 
                regionLineGraph.createComponentTable<LineData, std::set<RegionLineConnectionHandle>>();
            auto lineRightRegionConnections = 
                regionLineGraph.createComponentTable<LineData, std::set<RegionLineConnectionHandle>>();
            for (auto & rl : regionLineGraph.constraints<RegionLineConnectionData>()){
                RegionHandle rh = rl.topo.component<0>();
                LineHandle lh = rl.topo.component<1>();
                const Line3 & lineProj = regionLineGraph.data(lh).line;
                Vec3 rightDir = normalize(lineProj.first.cross(lineProj.second));
                const Vec3 & regionCenter = regionLineGraph.data(rh).normalizedCenter;
                bool isOnRight = (regionCenter - lineProj.center()).dot(rightDir) >= 0;
                (isOnRight ? lineRightRegionConnections : lineLeftRegionConnections)[lh].insert(rl.topo.hd);
            }

            float junctionWeightMax = 0.0f;
            for (auto & ll : regionLineGraph.constraints<LineRelationData>()){
                junctionWeightMax = std::max(ll.data.junctionWeight, junctionWeightMax);
            }
            assert(junctionWeightMax > 0.0);
            double regionAreaMax = 0.0;
            for (auto & r : regionLineGraph.components<RegionData>()){
                regionAreaMax = std::max(regionAreaMax, r.data.area);
            }

            size_t regionNum = regionLineGraph.internalComponents<RegionData>().size();
            size_t lineNum = regionLineGraph.internalComponents<LineData>().size();
            size_t boundaryNum = regionLineGraph.internalConstraints<RegionBoundaryData>().size();
            size_t rlConNum = regionLineGraph.internalConstraints<RegionLineConnectionData>().size();
            size_t llConNum = regionLineGraph.internalConstraints<LineRelationData>().size();
            size_t bjuncNum = boundaryJunctions.size();



            // data that MUST be updated after each inverse depth optimization !!!!
            struct ContinuousCaches {
                HandledTable<RegionHandle, Plane3> regionPlanes;
                HandledTable<LineHandle, Line3> spatialLines;
                HandledTable<RegionBoundaryHandle, double> rrDistances;
                HandledTable<RegionLineConnectionHandle, double> rlDistances;
                HandledTable<LineRelationHandle, double> llDistances;

                Plane3 & operator[](RegionHandle rh) { return regionPlanes[rh]; }
                Line3 & operator[](LineHandle lh) { return spatialLines[lh]; }
                double & operator[](RegionBoundaryHandle bh) { return rrDistances[bh]; }
                double & operator[](RegionLineConnectionHandle bh) { return rlDistances[bh]; }
                double & operator[](LineRelationHandle bh) { return llDistances[bh]; }

                void update(const RLGraph & mg, const RLGraphControls & controls, const RLGraphVars & vars){
                    // update continuous cache
                    for (auto & r : mg.components<RegionData>()){
                        (*this)[r.topo.hd] = Instance(mg, controls, vars, r.topo.hd);
                    }
                    for (auto & l : mg.components<LineData>()){
                        (*this)[l.topo.hd] = Instance(mg, controls, vars, l.topo.hd);
                    }

                    for (auto & c : mg.constraints<RegionBoundaryData>()){
                        double dist = 0.0;
                        int num = 0;
                        auto & samples = c.data.normalizedSampledPoints;
                        auto & inst1 = (*this)[c.topo.component<0>()];
                        auto & inst2 = (*this)[c.topo.component<1>()];
                        for (auto & ss : samples){
                            for (auto & s : ss){
                                double d1 = DepthAt(s, inst1);
                                double d2 = DepthAt(s, inst2);
                                dist += abs(d1 - d2);
                                num++;
                            }
                        }
                        dist /= num;
                        (*this)[c.topo.hd] = dist;
                    }
                    for (auto & c : mg.constraints<RegionLineConnectionData>()){
                        double dist = 0.0;
                        int num = 0;
                        auto & inst1 = (*this)[c.topo.component<0>()];
                        auto & inst2 = (*this)[c.topo.component<1>()];
                        for (auto & s : c.data.normalizedAnchors){
                            double d1 = DepthAt(s, inst1);
                            double d2 = DepthAt(s, inst2);
                            dist += abs(d1 - d2);
                            num++;
                        }
                        dist /= num;
                        (*this)[c.topo.hd] = dist;
                    }
                }
            };
            ContinuousCaches cache = {
                regionLineGraph.createComponentTable<RegionData, Plane3>(),
                regionLineGraph.createComponentTable<LineData, Line3>(),
                regionLineGraph.createConstraintTable<RegionBoundaryData>(0.0),
                regionLineGraph.createConstraintTable<RegionLineConnectionData>(0.0),
                regionLineGraph.createConstraintTable<LineRelationData>(0.0)
            };

            

            // create factor graph
            ml::FactorGraph fg;
            fg.reserveVarCategories(5);
            fg.reserveVars(regionNum + lineNum + boundaryNum + rlConNum + llConNum);

            /// add vars

            // add region orientation constraint flags
            // 3: {no constraints, vertical, tofront, ceilingfloor, side}
            auto regionOrientationVc = fg.addVarCategory(ml::FactorGraph::VarCategory{ 5, 0.5 });
            auto regionVhs = regionLineGraph.createComponentTable<RegionData, ml::FactorGraph::VarHandle>();
            for (auto & r : regionLineGraph.components<RegionData>()){
                regionVhs[r.topo.hd] = fg.addVar(regionOrientationVc);
            }

            // add line orientation constraint flags
            // 2: {no constraints, vp0, vp1, ...}
            auto lineOrientationVc = fg.addVarCategory(ml::FactorGraph::VarCategory{ vpnum + 1, 0.5 }); 
            auto lineVhs = regionLineGraph.createComponentTable<LineData, ml::FactorGraph::VarHandle>();
            for (auto & l : regionLineGraph.components<LineData>()){
                lineVhs[l.topo.hd] = fg.addVar(lineOrientationVc);
            }
            
            // add boundary flags
            // 2: {not connected, connected}
            auto boundaryOcclusionVc = fg.addVarCategory(ml::FactorGraph::VarCategory{ 2, 0.5 });
            auto boundaryVhs = regionLineGraph.createConstraintTable<RegionBoundaryData, ml::FactorGraph::VarHandle>();
            for (auto & b : regionLineGraph.constraints<RegionBoundaryData>()){
                boundaryVhs[b.topo.hd] = fg.addVar(boundaryOcclusionVc);
            }

            // add line connection flags
            // 3 : {both, left, right}
            auto lineConnectionSidesVc = fg.addVarCategory(ml::FactorGraph::VarCategory{ 3, 0.5 }); 
            auto lineConnectionSidesVhs = regionLineGraph.createComponentTable<LineData, ml::FactorGraph::VarHandle>();
            for (auto & l : regionLineGraph.components<LineData>()){
                lineConnectionSidesVhs[l.topo.hd] = fg.addVar(lineConnectionSidesVc);
            }

            // add ll connection flags
            // 2: {not connected, connected}
            auto llConVc = fg.addVarCategory(ml::FactorGraph::VarCategory{ 2, 0.5 });
            auto llConVhs = regionLineGraph.createConstraintTable<LineRelationData, ml::FactorGraph::VarHandle>();
            for (auto & r : regionLineGraph.constraints<LineRelationData>()){
                llConVhs[r.topo.hd] = fg.addVar(llConVc);
            }

            /// add factors

            // add region orientation unary cost (discrete and mixed)
            for (auto & r : regionLineGraph.components<RegionData>()){
                auto rh = r.topo.hd;
                bool intersectWithHorizon = regionIntersectWithHorizon[rh];
                int intersectWithVPs = regionIntersectWithVPs[rh];
                auto regionOrientationFc =
                    fg.addFactorCategory(
                    ml::FactorGraph::FactorCategory{ [this, rh, vertVPId, frontVPId, sideVPId, &cache,
                    intersectWithHorizon, intersectWithVPs](const int * labels, size_t nvars) -> double {
                    assert(nvars == 1);
                    assert(labels[0] >= 0 && labels[0] <= 4);
                    int rolabel = labels[0];
                    
                    // 0: front, 1: left, 2: right, 3: floor, 4: ceiling, 5: clutter, 6: unknown
                    const Vec7 & gcResp = gcResponse[rh]; 

                    double gcVertResp = std::max({ gcResp[0], gcResp[1], gcResp[2] });
                    double gcHoriResp = std::max({ gcResp[3], gcResp[4] });
                    assert(IsBetween(gcVertResp, 0.0, 1.1));
                    assert(IsBetween(gcHoriResp, 0.0, 1.1));
                    
                    if (rolabel == 0){ // no constraints
                        return 2;
                    }
                    else if (rolabel == 1){ // vertical
                        double violationWithGCVertResp = 1.0 - gcVertResp; // 0-1
                        double violationOfFarFromHorizon = intersectWithHorizon ? 0.1 : 1.0; // [0 1]

                        const Plane3 & plane = cache[rh];
                        double angleToVertical = M_PI_2 - AngleBetweenUndirectedVectors(plane.normal, vps[vertVPId]);
                        double violationThatCurrentPlaneTendNotToBeVertical = 1.0 - Gaussian(angleToVertical, M_PI_4);

                        return (violationWithGCVertResp * 2 + violationOfFarFromHorizon * 1 + violationThatCurrentPlaneTendNotToBeVertical);
                    }
                    else if (rolabel == 2){ // toward front
                        double violationWithGC = 1.0 - gcResp[0];
                        double violationOfFarFromVP = intersectWithVPs == frontVPId ? 0 : 1.0;
                        return violationWithGC * 3 + violationOfFarFromVP;
                    }
                    else if (rolabel == 3){ // floor/ceiling
                        double violationWithGC = 1.0 - std::max(gcResp[3], gcResp[4]);
                        double violationOfFarFromVP = intersectWithVPs == vertVPId ? 0 : 1.0;
                        return violationWithGC * 3 + violationOfFarFromVP;
                    }
                    else { // sided
                        double violationWithGC = 1.0 - std::max(gcResp[1], gcResp[2]);
                        double violationOfFarFromVP = intersectWithVPs == sideVPId ? 0 : 1.0;
                        return violationWithGC * 3 + violationOfFarFromVP;
                    }
                }, 0.5 });
                fg.addFactor({ regionVhs[rh] }, regionOrientationFc);
            }

            // add line orientation cost
            for (auto & l : regionLineGraph.components<LineData>()){
                auto lh = l.topo.hd;
                auto lineOrientationFc = fg.addFactorCategory(
                    ml::FactorGraph::FactorCategory{ [this, lh, &lineScoresForVPs, &cache](const int * labels, size_t nvars) -> double{
                    assert(nvars == 1);
                    int lolabel = labels[0];
                    assert(lolabel >= 0 && lolabel <= vps.size());

                    double maxDistanceToRegions = 0.0;
                    auto & rlConnections = regionLineGraph.topo(lh).constraints<RegionLineConnectionData>();
                    for (auto rlh : rlConnections){
                        maxDistanceToRegions = std::max(cache[rlh], maxDistanceToRegions);
                    }
                    assert(maxDistanceToRegions > 0);
                    
                    if (lolabel == 0){ // no constraints
                        return 5;
                    }
                    else {
                        int vpclaz = lolabel - 1;
                        assert(0 <= lineScoresForVPs[lh][vpclaz] && lineScoresForVPs[lh][vpclaz] <= 1.0);
                        return 10 * (1.0 - lineScoresForVPs[lh][vpclaz]) * 
                            maxDistanceToRegions * 5.0 / normalize(regionLineGraph.data(lh).line).length();
                    }
                }, 0.5 });
                fg.addFactor({ lineVhs[lh] }, lineOrientationFc);
            }

            // add boundary flag unary cost
            for (auto & b : regionLineGraph.constraints<RegionBoundaryData>()){
                auto bh = b.topo.hd;
                auto boundaryFlagFc = fg.addFactorCategory(
                    ml::FactorGraph::FactorCategory{ [this, bh, &cache](const int * labels, size_t nvars) -> double{
                    assert(nvars == 1);
                    //double distance = cache[bh];
                    const double & occResp = occlusionResponse[bh];                    
                    if (labels[0] == 0){ // occlusion
                        return occResp < 0 ? 10.0 : 0.0;
                    }
                    else { // not occlusion
                        return occResp > 0 ? 10.0 : 0.0;
                    }
                }, 0.5 });
                fg.addFactor({ boundaryVhs[bh] }, boundaryFlagFc);
            }

            // add line connection side cost
            for (auto & l : regionLineGraph.components<LineData>()){
                auto lh = l.topo.hd;
                auto & leftRegionCons = lineLeftRegionConnections[lh];
                auto & rightRegionCons = lineRightRegionConnections[lh];
                auto lineConnectionSideFlagFc = fg.addFactorCategory(
                    ml::FactorGraph::FactorCategory{ [this, lh, &cache, 
                    &leftRegionCons, &rightRegionCons](const int * labels, size_t nvars) -> double{
                    assert(nvars == 1);
                    assert(labels[0] == 0 || labels[0] == 1 || labels[0] == 2);

                    // check consistency of regions on both sides
                    double meanRegionDistances[] = { 0, 0 }; // left, right               
                    
                    for (RegionLineConnectionHandle ch : leftRegionCons){
                        assert(!IsInfOrNaN(cache[ch]));
                        meanRegionDistances[0] += cache[ch];
                    }
                    meanRegionDistances[0] /= std::max(leftRegionCons.size(), 1ull);
                    for (RegionLineConnectionHandle ch : rightRegionCons){
                        assert(!IsInfOrNaN(cache[ch]));
                        meanRegionDistances[1] += cache[ch];
                    }
                    meanRegionDistances[1] /= std::max(rightRegionCons.size(), 1ull);

                    double maxx = std::max(meanRegionDistances[0], meanRegionDistances[1]);
                    if (maxx > 0){
                        meanRegionDistances[0] /= maxx;
                        meanRegionDistances[1] /= maxx;
                    }
                    if (labels[0] == 0){ // both
                        return meanRegionDistances[0] + meanRegionDistances[1];
                    }
                    else if (labels[0] == 1){ // left
                        return meanRegionDistances[0] * 3;
                    }
                    else{ // right
                        return meanRegionDistances[1] * 3;
                    }
                }, 0.5 });
                fg.addFactor({ lineConnectionSidesVhs[lh] }, lineConnectionSideFlagFc);
            }


            // add line-line connection flag unary cost
            for (auto & r : regionLineGraph.constraints<LineRelationData>()){
                auto rh = r.topo.hd;
                double junctionWeightRatio = r.data.junctionWeight / junctionWeightMax;
                assert(!IsInfOrNaN(junctionWeightRatio));
                assert(junctionWeightRatio >= 0 && junctionWeightRatio <= 1.0);
                auto linelineRelFlagFc = fg.addFactorCategory(
                    ml::FactorGraph::FactorCategory{ [this, rh, junctionWeightRatio](const int * labels, size_t nvars) -> double{
                    assert(nvars == 1);
                    int flag = labels[0];
                    if (flag == 0){ // not connected
                        return junctionWeightRatio * 5.0;
                    }
                    else{ // connected
                        return (1.0 - junctionWeightRatio);
                    }
                }, 0.5 });
                fg.addFactor({ llConVhs[rh] }, linelineRelFlagFc);
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
                return occlusionNum <= 2 ? 0.0 : 100.0;
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


            

            // initialize
            LOG("initialize");
            AttachPrincipleDirectionConstraints(regionLineGraph, controls, M_PI / 50.0);
            AttachWallConstriants(regionLineGraph, controls, M_PI / 80.0);
            AttachAnchorToCenterOfLargestLineIfNoAnchorExists(regionLineGraph, controls);

            const auto defaultControls = controls;

            // solve continuous vars
            auto vars = SolveVariables(regionLineGraph, controls);
            NormalizeVariables(regionLineGraph, controls, vars);
            cache.update(regionLineGraph, controls, vars);

            //Visualize(view, regionLineGraph, controls, vars);

            for(int epoch = 0; epoch < 1; epoch ++){

                LOG("epoch " + std::to_string(epoch));          

                // solve discrete vars
                auto results = fg.solve(100, 3,
                    [](int epoch, double energy, double denergy, const ml::FactorGraph::ResultTable & ){
                    std::cout << "energy: " << energy << std::endl;
                    return denergy < 0;
                });

                // install discrete vars into props
                // region orientation
                for (auto & r : regionLineGraph.components<RegionData>()){
                    int rolabel = results[regionVhs[r.topo.hd]];
                    auto & prop = controls[r.topo.hd];
                    if (rolabel == 0){
                        prop.orientationClaz = prop.orientationNotClaz = -1;
                    }
                    else if (rolabel == 1){ // vertical
                        prop.orientationClaz = -1;
                        prop.orientationNotClaz = vertVPId;
                    }
                    else if (rolabel == 2){ // tofront
                        prop.orientationClaz = frontVPId;
                        prop.orientationNotClaz = -1;
                    }
                    else if (rolabel == 3){ // floorceiling
                        prop.orientationClaz = vertVPId;
                        prop.orientationNotClaz = -1;
                    }
                    else{
                        prop.orientationClaz = sideVPId;
                        prop.orientationNotClaz = -1;
                    }
                }
                // line orientation
                for (auto & l : regionLineGraph.components<LineData>()){
                    int lolabel = results[lineVhs[l.topo.hd]];
                    auto & prop = controls[l.topo.hd];
                    if (lolabel == 0){
                        prop.orientationClaz = prop.orientationNotClaz = -1;
                    }
                    else{
                        assert(lolabel >= 1 && lolabel <= vpnum);
                        prop.orientationClaz = lolabel - 1;
                        prop.orientationNotClaz = lolabel - 1;
                    }
                }

                // boundary occlusion
                for (auto & b : regionLineGraph.constraints<RegionBoundaryData>()){
                    controls[b.topo.hd].weight = results[boundaryVhs[b.topo.hd]] == 1 ? defaultControls[b.topo.hd].weight : 1e-3;
                }

                std::array<int, 3> lineConnectionSideLabels = { { 0, 0, 0 } };
                for (auto & l : regionLineGraph.components<LineData>()){
                    int label = results[lineConnectionSidesVhs[l.topo.hd]];
                    lineConnectionSideLabels[label] ++;
                    for (RegionLineConnectionHandle leftConH : lineLeftRegionConnections[l.topo.hd]){
                        auto & prop = controls[leftConH];
                        auto & defaultProp = defaultControls[leftConH];
                        /*if (!regionLineGraph.data(leftConH).detachable){
                            prop.weight = defaultProp.weight;
                            continue;
                        }*/
                        prop.weight = (label == 0 || label == 1) ? defaultProp.weight : 1e-3;
                    }
                    for (RegionLineConnectionHandle rightConH : lineRightRegionConnections[l.topo.hd]){
                        auto & prop = controls[rightConH];
                        auto & defaultProp = defaultControls[rightConH];
                        /*if (!regionLineGraph.data(rightConH).detachable){
                            prop.weight = defaultProp.weight;
                            continue;
                        }*/
                        prop.weight = (label == 0 || label == 2) ? defaultProp.weight : 1e-3;
                    }
                }

                for (auto & r : regionLineGraph.constraints<LineRelationData>()) {
                    controls[r.topo.hd].weight = results[llConVhs[r.topo.hd]] == 1 ? defaultControls[r.topo.hd].weight : 1e-3;
                }

                // solve continuous vars
                vars = SolveVariables(regionLineGraph, controls);
                NormalizeVariables(regionLineGraph, controls, vars);
                cache.update(regionLineGraph, controls, vars);

                //Visualize(view, regionLineGraph, controls, vars);
            }

            return std::move(cache.regionPlanes);

        }

    }
}