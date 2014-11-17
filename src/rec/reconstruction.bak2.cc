
            #pragma region function_declarations

            // vertex functions
            inline MixedGraphVertex CreateRegionCCVertex(int regionCCId, const RecContext & context);
            inline MixedGraphVertex CreateLineCCVertex(int lineCCId, const RecContext & context);
            inline Rational ComputeDeterminedAnchorsRatio(const MixedGraph & g,
                const MixedGraphVertHandle & vh);
            inline std::vector<Point3> CollectDeterminedAnchors(const MixedGraph & g,
                const MixedGraphVertHandle & vh);

            inline void BuildCandidates(const RecContext & context, MixedGraph & g, const MixedGraphVertHandle & vh);
            inline void RegisterChoices(const RecContext & context,
                MixedGraph & g, const MixedGraphVertHandle & vh,
                std::vector<Choice> & choices, std::vector<double> & probabilities,
                double baseProb, int maxChoiceNum);
            inline void PickChoice(const RecContext & context,
                MixedGraph & g, const MixedGraphVertHandle & vh,
                const Choice & choice);


            // edge functions
            inline MixedGraphEdge CreateRegionRegionEdge(const RegionIndex & ri1, const RegionIndex & ri2, const std::vector<Point3> & anchors);
            inline MixedGraphEdge CreateRegionLineEdge(const RegionIndex & ri, const LineIndex & li, const std::vector<Point3> & anchors);
            

            // holistic functions
            inline void UpdateNeighbors(const RecContext & context, MixedGraph & g, const MixedGraphVertHandle & vh);

            #pragma endregion function_declarations


            #pragma region function_implementations
            // implementations
            // vertex functions

           
            // region cc vertex method implementations
            void RegionCCVertexData::buildCandidates(const RecContext & context,
                const MixedGraph & g, const MixedGraphVertHandle & selfHandle) {

                // make candidate planes
                candidatePlanesByRoot.clear();

                double scale = context.initialBoundingBox.outerSphere().radius;
                auto determinedAnchors = CollectDeterminedAnchors(g, selfHandle);

                //// merge near anchors
                //std::vector<decltype(determinedAnchors.begin())> mergedAnchorIters;
                //core::MergeNearRTree(determinedAnchors.begin(), determinedAnchors.end(), std::back_inserter(mergedAnchorIters), 
                //    core::no(), scale / 100.0);

                // iterate over merged anchors to collect plane candidates
                for (auto & anchor : determinedAnchors/*mergedAnchorIters*/){
                    //const Point3 & anchor = *i;
                    for (auto & vp : context.vanishingPoints){
                        Plane3 plane(anchor, vp);
                        if (OPT_IgnoreTooSkewedPlanes){
                            if (norm(plane.root()) <= scale / 5.0)
                                continue;
                        }
                        if (OPT_IgnoreTooFarAwayPlanes){
                            bool valid = true;
                            for (auto & ri : regionIndices){
                                if (!valid)
                                    break;
                                auto & rd = GetData(ri, context.regionsNets);
                                if (rd.contours.back().size() < 3)
                                    continue;
                                auto & cam = context.views[ri.viewId].camera;
                                for (int i = 0; i < rd.contours.back().size(); i++){
                                    auto dir = cam.spatialDirection(ToPoint2(rd.contours.back()[i]));
                                    auto intersectionOnPlane = IntersectionOfLineAndPlane(InfiniteLine3(Point3(0, 0, 0), dir), plane).position;
                                    if (norm(intersectionOnPlane) > scale * 5.0){
                                        valid = false;
                                        break;
                                    }
                                }
                            }
                            if (!valid)
                                continue;
                        }

                        static const double distFromPointToPlaneThres = scale / 12.0;

                        // insert new root data
                        auto & pcd = candidatePlanesByRoot[plane.root()];
                        pcd.plane = plane;

                        // collect distance votes
                        double distVotes = 0.0;
                        std::vector<Vec3> nearbyAnchors;
                        for (int i = 0; i < /*mergedAnchorIters*/determinedAnchors.size(); i++){
                            double distanceToPlane = plane.distanceTo(/**mergedAnchorIters[i]*/determinedAnchors[i]);
                            if (distanceToPlane > distFromPointToPlaneThres)
                                continue;
                            distVotes += Gaussian(distanceToPlane, distFromPointToPlaneThres);
                            pcd.inlierAnchors.push_back(i);
                            nearbyAnchors.push_back(/**mergedAnchorIters[i]*/determinedAnchors[i]);
                        }
                        pcd.regionInlierAnchorsDistanceVotesSum = distVotes;
                        pcd.regionInlierAnchorsConvexContourVisualArea =
                            ComputeVisualAreaOfDirections(tangentialPlane,
                                xOnTangentialPlane, yOnTangentialPlane, nearbyAnchors, true);
                    }
                }
               
            }

            void RegionCCVertexData::registerChoices(const RecContext & context,
                const MixedGraph & g, const MixedGraphVertHandle & selfHandle,
                std::vector<Choice> & choices, std::vector<double> & probabilities,
                double baseProb, int maxChoiceNum) const{

                std::vector<Scored<Choice>> newChoices;

                int planeId = 0;
                double maxMeanVote = 0.0;
                for (auto & c : candidatePlanesByRoot){
                    maxMeanVote = std::max(maxMeanVote, c.second.regionInlierAnchorsDistanceVotesSum / c.second.inlierAnchors.size());
                }

                double fullCompleteness = ComputeDeterminedAnchorsRatio(g, selfHandle).value(0.0);

                // collect all choices and probabilities
                for (auto & c : candidatePlanesByRoot){
                    auto & candidatePlaneData = c.second;
                    Choice choice = { selfHandle, planeId++ };
                    auto inlierOccupationRatio =
                        candidatePlaneData.regionInlierAnchorsConvexContourVisualArea / regionConvexContourVisualArea;
                    double probability =
                        (fullCompleteness *
                        (inlierOccupationRatio > 0.4 ? 1.0 : 1e-4)) *
                        c.second.regionInlierAnchorsDistanceVotesSum / c.second.inlierAnchors.size();

                    assert(c.second.regionInlierAnchorsDistanceVotesSum / c.second.inlierAnchors.size() <= 1.0);
                    newChoices.push_back(ScoreAs(choice, probability + baseProb));
                }

                // select best [maxChoiceNum] choices
                std::sort(newChoices.begin(), newChoices.end(), std::greater<void>());
                for (int i = 0; i < std::min<int>(maxChoiceNum, newChoices.size()); i++){
                    choices.push_back(newChoices[i].component);
                    probabilities.push_back(newChoices[i].score);
                }
            }

            void RegionCCVertexData::pickChoice(const RecContext & context,
                MixedGraph & g,
                const Choice & choice) {
                assert(&(g.data(choice.vertHandle).regionCCVD()) == this);
                chosenPlane = (candidatePlanesByRoot.begin() + choice.choiceId)->second.plane;
                //// update related edge anchors and notify the vertex on the other side
                //for (const MixedGraphEdgeHandle & eh : g.topo(choice.vertHandle).uppers){
                //    auto & ed = g.data(eh);
                //    assert(
                //        (ed.connectsRegionAndLine() && context.regionConnectedComponentIds.at(ed.rili().first) == ccId) ||
                //        (ed.connectsRegionAndRegion() && 
                //            (context.regionConnectedComponentIds.at(ed.riri().first) == ccId || context.regionConnectedComponentIds.at(ed.riri().second) == ccId)
                //        ));
                //    for (auto & anchor : ed.anchors()){
                //        anchor = IntersectionOfLineAndPlane(InfiniteLine3(Point3(0, 0, 0), anchor), chosenPlane).position;
                //    }
                //    ed.determined = true;
                //}
            }

            // line cc vertex method implementations
            void LineCCVertexData::buildCandidates(const RecContext & context,
                const MixedGraph & g, const MixedGraphVertHandle & selfHandle) {

                candidateDepthFactors.clear();
                candidateDepthFactors[1.0] = 0.1;
                std::vector<double> depthFactors;
                for (const MixedGraphEdgeHandle & eh : g.topo(selfHandle).uppers){
                    auto & ed = g.data(eh);
                    if (ed.determined()){
                        assert(ed.connectsRegionAndLine());
                        assert(context.lineConnectedComponentIds.at(ed.rili().second) == ccId);
                        auto & li = ed.rili().second;
                        auto & line = context.reconstructedLines.at(li);
                        for (auto & anchor : ed.anchors()){
                            double depthVarOnLine = norm(DistanceBetweenTwoLines(line.infinieLine(), InfiniteLine3(Point3(0, 0, 0), anchor))
                                .second.second);
                            double depthValueOnRegion = norm(anchor);
                            if (!IsInfOrNaN(depthVarOnLine) && !IsInfOrNaN(depthValueOnRegion))
                                candidateDepthFactors[depthValueOnRegion / depthVarOnLine] += 1.0;
                        }
                    }
                }

            }

            void LineCCVertexData::registerChoices(const RecContext & context,
                const MixedGraph & g, const MixedGraphVertHandle & selfHandle,
                std::vector<Choice> & choices, std::vector<double> & probabilities,
                double baseProb, int maxChoiceNum) const {

                std::vector<Scored<Choice>> newChoices;
                int depthId = 0;
                double maxVote = 0.0;
                for (auto & c : candidateDepthFactors)
                    maxVote = maxVote < c.second ? c.second : maxVote;

                double fullCompleteness = ComputeDeterminedAnchorsRatio(g, selfHandle).value(0.0);

                for (auto & c : candidateDepthFactors){
                    auto & candidateDepthVote = c.second;
                    Choice choice = { selfHandle, depthId++ };
                    double probability = (fullCompleteness * 0.9
                        + double(lineIndices.size()) / context.reconstructedLines.size() * 0.1)
                        * candidateDepthVote / maxVote;
                    newChoices.push_back(ScoreAs(choice, probability + baseProb));
                }

                std::sort(newChoices.begin(), newChoices.end(), std::greater<void>());
                for (int i = 0; i < std::min<int>(maxChoiceNum, newChoices.size()); i++){
                    choices.push_back(newChoices[i].component);
                    probabilities.push_back(newChoices[i].score);
                }
            }

            void LineCCVertexData::pickChoice(const RecContext & context,
                MixedGraph & g,
                const Choice & choice) {
                assert(&(g.data(choice.vertHandle).lineCCVD()) == this);
                chosenDepthFactor = (candidateDepthFactors.begin() + choice.choiceId)->first[0];
                // update related edge anchors
                /* for (const MixedGraphEdgeHandle & eh : g.topo(selfHandle).uppers){
                auto & ed = g.data(eh);
                assert(ed.connectsRegionAndLine() &&
                context.lineConnectedComponentIds.at(ed.rili.second) == ccId);
                const auto & line = context.reconstructedLines.at(ed.rili.second);
                for (auto & anchor : ed.anchors){
                auto pOnLine = DistanceBetweenTwoLines(line.infinieLine(), InfiniteLine3(Point3(0, 0, 0), anchor))
                .second.second;
                anchor = pOnLine * depthFactor;
                }
                ed.determined = true;
                }*/
            }



            // unified methods for vertex
            inline Rational ComputeDeterminedAnchorsRatio(const MixedGraph & g,
                const MixedGraphVertHandle & selfHandle) {
                Rational r(0.0, 0.0);
                for (const MixedGraphEdgeHandle & eh : g.topo(selfHandle).uppers){
                    const auto & ed = g.data(eh);
                    r.denominator += ed.anchors.size();
                    r.numerator += ed.determined() ? ed.anchors.size() : 0.0;
                }
                return r;
            }

            inline std::vector<Point3> CollectDeterminedAnchors(const MixedGraph & g,
                const MixedGraphVertHandle & selfHandle)  {
                std::vector<Point3> ps;
                for (const MixedGraphEdgeHandle & eh : g.topo(selfHandle).uppers){
                    const auto & ed = g.data(eh);
                    if (ed.determined())
                        ps.insert(ps.end(), ed.anchors.begin(), ed.anchors.end());
                }
                return ps;
            }


            inline void BuildCandidates(const RecContext & context, MixedGraph & g, const MixedGraphVertHandle & vh){
                auto & vd = g.data(vh);
                if (vd.isRegionCC())
                    vd.regionCCVD().buildCandidates(context, g, vh);
                else if (vd.isLineCC())
                    vd.lineCCVD().buildCandidates(context, g, vh);
            }

            inline void RegisterChoices(const RecContext & context,
                MixedGraph & g, const MixedGraphVertHandle & vh,
                std::vector<Choice> & choices, std::vector<double> & probabilities,
                double baseProb, int maxChoiceNum){
                auto & vd = g.data(vh);
                if (vd.isRegionCC())
                    vd.regionCCVD().registerChoices(context, g, vh, choices, probabilities, baseProb, maxChoiceNum);
                else if (vd.isLineCC())
                    vd.lineCCVD().registerChoices(context, g, vh, choices, probabilities, baseProb, maxChoiceNum);
            }

            inline void PickChoice(const RecContext & context,
                MixedGraph & g, const MixedGraphVertHandle & vh,
                const Choice & choice){
                auto & vd = g.data(vh);
                if (vd.isRegionCC())
                    vd.regionCCVD().pickChoice(context, g, choice);
                else if (vd.isLineCC())
                    vd.lineCCVD().pickChoice(context, g, choice);
            }





            inline void UpdateNeighbors(const RecContext & context, MixedGraph & g, const MixedGraphVertHandle & vh){
                auto & vd = g.data(vh);
                auto & ehs = g.topo(vh).uppers;
                for (auto & eh : ehs){
                    auto anotherVh = g.topo(eh).lowers[0];
                    if (anotherVh == vh)
                        anotherVh = g.topo(eh).lowers[1];
                    auto & ed = g.data(eh);
                    if (vd.isRegionCC()){
                        auto & plane = vd.regionCCVD().chosenPlane;

                    }
                    else if (vd.isLineCC()){
                        auto depthFactor = vd.lineCCVD().chosenDepthFactor;
                        
                    }
                    ed._determined = true;
                }
            }
            
          
            #pragma endregion function_implementations





            
            void InitializSpatialRegionPlanes(const RecContext & context,
                MixedGraph & graph, 
                const std::vector<MixedGraphVertHandle> & regionCCIdToVHandles,
                const std::vector<MixedGraphVertHandle> & lineCCidToVHandles,
                std::vector<Plane3> & resultRegionConnectedComponentPlanes,
                std::vector<double> & resultLineConnectedComponentDepthFactors,
                int trialNum,
                bool useWeightedRandomSelection) {

                double scale = context.initialBoundingBox.outerSphere().radius;

                // initial cc ids status
                std::vector<int> regionCCIds(context.regionConnectedComponentsNum);
                std::iota(regionCCIds.begin(), regionCCIds.end(), 0);
                const std::set<int> initialRegionCCIdsNotDeterminedYet(regionCCIds.begin(), regionCCIds.end());

                std::vector<int> lineCCIds(context.lineConnectedComponentsNum);
                std::iota(lineCCIds.begin(), lineCCIds.end(), 0);
                const std::set<int> initialLineCCIdsNotDeterminedYet(lineCCIds.begin(), lineCCIds.end());


                // initialize region planes
                std::vector<Plane3> initialRegionConnectedComponentPlanes(context.regionConnectedComponentsNum);
                for (auto & r : context.regionConnectedComponentIds){
                    auto & ri = r.first;
                    auto & rd = GetData(ri, context.regionsNets);
                    Vec3 centerDir = context.views[ri.viewId].camera.spatialDirection(rd.center);
                    int regionCCId = context.regionConnectedComponentIds.at(ri);
                    initialRegionConnectedComponentPlanes[regionCCId].anchor = normalize(centerDir) * scale;
                    initialRegionConnectedComponentPlanes[regionCCId].normal = normalize(centerDir);
                }
                // initialize line depth factors
                std::vector<double> initialLineConnectedComponentDepthFactors(context.lineConnectedComponentsNum, 1.0);


                // start MC reasoning
                std::vector<Scored<std::pair<std::vector<Plane3>, std::vector<double>>>>
                    candidates(trialNum,
                    ScoreAs(std::make_pair(initialRegionConnectedComponentPlanes, initialLineConnectedComponentDepthFactors),
                    0.0));

                const auto & constGraph = graph;
                auto task = [&candidates, &context, 
                    &initialRegionCCIdsNotDeterminedYet, &initialLineCCIdsNotDeterminedYet, 
                    &useWeightedRandomSelection, &constGraph, &regionCCIdToVHandles, &lineCCidToVHandles](int t){
                    std::cout << "task: " << t << std::endl;
                    auto & candidate = candidates[t];

                    // random engine initialized
                    std::random_device rd;
                    std::default_random_engine gen(rd());

                    // copy graph
                    MixedGraph g = constGraph;

                    // undetermined checkers
                    std::set<int> regionCCIdsNotDeterminedYet = initialRegionCCIdsNotDeterminedYet;
                    std::set<int> lineCCIdsNotDeterminedYet = initialLineCCIdsNotDeterminedYet;

                    // choices and probabilities
                    std::vector<Choice> choices;
                    std::vector<double> choiceProbabilities;
                    choices.reserve(regionCCIdsNotDeterminedYet.size() + lineCCIdsNotDeterminedYet.size());
                    choiceProbabilities.reserve(regionCCIdsNotDeterminedYet.size() + lineCCIdsNotDeterminedYet.size());

                    // start expansion
                    std::cout << "start expansion" << std::endl;
                    while ((regionCCIdsNotDeterminedYet.size() + lineCCIdsNotDeterminedYet.size()) > 0){

                        // collect choices and probabilities
                        choices.clear();
                        choiceProbabilities.clear();
                        for (int regionCCId : regionCCIdsNotDeterminedYet){
                            auto & vh = regionCCIdToVHandles[regionCCId];
                            auto & vd = g.data(vh);
                            if (vd.determined)
                                continue;
                            vd.regionCCVD().buildCandidates(context, g, vh);
                            vd.regionCCVD().registerChoices(context, g, vh, choices, choiceProbabilities, 1e-5, 1);
                        }
                        for (int lineCCId : lineCCIdsNotDeterminedYet){
                            auto & vh = lineCCidToVHandles[lineCCId];
                            auto & vd = g.data(vh);
                            if (vd.determined)
                                continue;
                            vd.lineCCVD().buildCandidates(context, g, vh);
                            vd.lineCCVD().registerChoices(context, g, vh, choices, choiceProbabilities, 1e-5, 1);
                        }

                        assert(choices.size() == choiceProbabilities.size());

                        if (std::accumulate(choiceProbabilities.begin(), choiceProbabilities.end(), 0.0) == 0.0){
                            std::cerr << "all zero probabilities!" << std::endl;
                            break;
                        }

                        int selected = -1;

                        if (useWeightedRandomSelection){
                            // a workaround constructor since VS2013 lacks the range iterator constructor for std::discrete_distribution
                            int ord = 0;
                            std::discrete_distribution<int> distribution(choiceProbabilities.size(),
                                0.0, 1000.0,
                                [&choiceProbabilities, &ord](double){
                                return choiceProbabilities[ord++];
                            });
                            selected = distribution(gen);
                        }
                        else{
                            selected = std::distance(choiceProbabilities.begin(),
                                std::max_element(choiceProbabilities.begin(), choiceProbabilities.end()));
                        }

                        // made choice
                        const Choice & choice = choices[selected];
                        auto & exeVD = g.data(choice.vertHandle);
                        if (exeVD.isRegionCC()){
                            if (OPT_DisplayMessages)
                                std::cout << "chosen unit - region cc: " << exeVD.regionCCVD().ccId << std::endl;
                            exeVD.regionCCVD().pickChoice(context, g, choice.vertHandle, choice, 
                                candidate.component.first[exeVD.regionCCVD().ccId]);
                            regionCCIdsNotDeterminedYet.erase(exeVD.regionCCVD().ccId);
                        }
                        else if(exeVD.isLineCC()){
                            if (OPT_DisplayMessages)
                                std::cout << "chosen unit - line cc: " << exeVD.lineCCVD().ccId << std::endl;
                            exeVD.lineCCVD().pickChoice(context, g, choice.vertHandle, choice,
                                candidate.component.second[exeVD.lineCCVD().ccId]);
                            lineCCIdsNotDeterminedYet.erase(exeVD.lineCCVD().ccId);
                        }
                        exeVD.determined = true;


                    } // while
                    std::cout << "expansion done" << std::endl;


                    // score this candidate
                    double distanceSumOfRegionRegionConnections = 0.0;
                    double distanceSumOfRegionLineConnections = 0.0;

                    auto & regionConnectedComponentPlanes = candidate.component.first;
                    auto & lineConnectedComponentDepthFactors = candidate.component.second;

                    // region region
                    for (int i = 0; i < context.views.size(); i++){
                        auto & cam = context.views[i].camera;
                        for (auto & b : context.regionsNets[i].regions().elements<1>()){
                            auto ri1 = RegionIndex{ i, b.topo.lowers[0] };
                            auto ri2 = RegionIndex{ i, b.topo.lowers[1] };
                            int regionCCId1 = context.regionConnectedComponentIds.at(ri1);
                            int regionCCId2 = context.regionConnectedComponentIds.at(ri2);
                            for (auto & pts : b.data.sampledPoints){
                                for (auto & p : pts){
                                    auto dir = cam.spatialDirection(p);
                                    auto anchor1 = IntersectionOfLineAndPlane(InfiniteLine3(Point3(0, 0, 0), dir),
                                        regionConnectedComponentPlanes[regionCCId1]).position;
                                    auto anchor2 = IntersectionOfLineAndPlane(InfiniteLine3(Point3(0, 0, 0), dir),
                                        regionConnectedComponentPlanes[regionCCId2]).position;
                                    double distance = Distance(anchor1, anchor2);
                                    distanceSumOfRegionRegionConnections += distance;
                                }
                            }
                        }
                    }
                    // region line
                    for (auto & pp : context.regionLineConnections){
                        LineIndex li = pp.first.second;
                        RegionIndex ri = pp.first.first;
                        int lineCCId = context.lineConnectedComponentIds.at(li);
                        int regionCCId = context.regionConnectedComponentIds.at(ri);
                        auto & samplePoints = pp.second;
                        auto line = context.reconstructedLines.at(li);
                        line.first *= lineConnectedComponentDepthFactors[lineCCId];
                        line.second *= lineConnectedComponentDepthFactors[lineCCId];
                        for (auto & p : samplePoints){ // insert anchors into the rec info of the related region
                            auto pOnLine = DistanceBetweenTwoLines(line.infinieLine(), InfiniteLine3(Point3(0, 0, 0), p))
                                .second.second;
                            auto pOnRegion = IntersectionOfLineAndPlane(InfiniteLine3(Point3(0, 0, 0), p),
                                regionConnectedComponentPlanes[regionCCId]).position;
                            double distance = Distance(pOnLine, pOnRegion);
                            distanceSumOfRegionLineConnections += distance;
                        }
                    }

                    std::cout << "distance sum of region-region connections: " << distanceSumOfRegionRegionConnections << std::endl;
                    std::cout << "distance sum of region-line connections: " << distanceSumOfRegionLineConnections << std::endl;

                    candidate.score = -(distanceSumOfRegionRegionConnections + distanceSumOfRegionLineConnections);

                    if (OPT_DisplayOnEachTrial){
                        IF_DEBUG_USING_VISUALIZERS{
                        DisplayReconstruction(-1, -1, {}, {},
                        candidate.component.first, candidate.component.second, context);
                    }
                    }

                }; // task

                // run tasks
                int threadsNum = std::min<int>(std::max(1u, std::thread::hardware_concurrency() - 1), trialNum);
                std::cout << "threads num: " << threadsNum << std::endl;
                if (threadsNum == 1){
                    task(0);
                }else{
                    for (int i = 0; i < trialNum; i += threadsNum){
                        std::vector<std::thread> parallelThreads(threadsNum);
                        for (int t = i; t < std::min(trialNum, i + threadsNum); t++){
                            parallelThreads[t - i] = std::thread(task, t);
                        }
                        for (int t = i; t < std::min(trialNum, i + threadsNum); t++){
                            parallelThreads[t - i].join();
                        }
                    }
                }


                // select best candidate
                const auto & result = std::max_element(candidates.begin(), candidates.end())->component;
                resultRegionConnectedComponentPlanes = result.first;
                resultLineConnectedComponentDepthFactors = result.second;

                // visualize result of this task
                if (OPT_DisplayAtLast){
                    IF_DEBUG_USING_VISUALIZERS{
                    DisplayReconstruction(-1, -1, {}, {},
                    resultRegionConnectedComponentPlanes, resultLineConnectedComponentDepthFactors, context);
                }
                }

            }


            void OptimizeSpatialRegionPlanes(const RecContext & context,
                std::vector<Plane3> & resultRegionConnectedComponentPlanes,
                std::vector<double> & resultLineConnectedComponentDepthFactors) {

                // Simulated Annealing



            }