
#include <Eigen/Dense>
#include <Eigen/Sparse>

#include "../core/algorithms.hpp"
#include "../core/utilities.hpp"
#include "../vis/visualize3d.hpp"
#include "reconstruction.hpp"

#include "../core/debug.hpp"

namespace panoramix {
    namespace rec {


        namespace {

            template <class T, int N>
            inline Line<T, N> NormalizeLine(const Line<T, N> & l) {
                return Line<T, N>(normalize(l.first), normalize(l.second));
            }

            void CollectIndices(const std::vector<View<PerspectiveCamera>> & views,
                const std::vector<RegionsGraph> & regionsNets,
                std::vector<RegionIndex> & regionIndices,
                ComponentIndexHashMap<RegionIndex, int> & regionIndexToId){
                regionIndices.clear();
                regionIndexToId.clear();
                for (int i = 0; i < views.size(); i++){
                    RegionIndex ri;
                    ri.viewId = i;
                    for (auto & rd : regionsNets[i].elements<0>()){
                        ri.handle = rd.topo.hd;
                        regionIndices.push_back(ri);
                        regionIndexToId[ri] = regionIndices.size() - 1;
                    }
                }
            }

            void CollectIndices(const std::vector<View<PerspectiveCamera>> & views,
                const std::vector<LinesGraph> & linesNets,
                std::vector<LineIndex> & lineIndices,
                ComponentIndexHashMap<LineIndex, int> & lineIndexToIds){
                lineIndices.clear();
                lineIndexToIds.clear();
                for (int i = 0; i < views.size(); i++) {
                    LineIndex li;
                    li.viewId = i;
                    for (auto & ld : linesNets[i].elements<0>()) {
                        li.handle = ld.topo.hd;
                        lineIndices.push_back(li);
                        lineIndexToIds[li] = lineIndices.size() - 1;
                    }
                }
            }

            // line depth ratio
            double ComputeDepthRatioOfPointOnSpatialLine(Vec3 lineFirstPointDir,
                Vec3 p, Vec3 vp) {
                // firstp -> p vp
                //  \      /
                //   \    /
                //    center
                lineFirstPointDir /= norm(lineFirstPointDir);
                p /= norm(p);
                vp /= norm(vp);

                if ((p - lineFirstPointDir).dot(vp) < 0)
                    vp = -vp;
                double angleCenter = AngleBetweenDirections(lineFirstPointDir, p);
                double angleFirstP = AngleBetweenDirections(-lineFirstPointDir, vp);
                double angleP = AngleBetweenDirections(-p, -vp);
                //assert(FuzzyEquals(angleCenter + angleFirstP + angleP, M_PI, 0.1));
                return sin(angleFirstP) / sin(angleP);
            }



            static const double MinimumJunctionWeght = 1e-5;

            void EstimateSpatialLineDepthsOnce(const std::vector<View<PerspectiveCamera>> & views,
                const std::vector<LinesGraph> & linesNets,
                const std::vector<Vec3> & vanishingPoints,
                const std::vector<LineIndex> & lineIndices,
                const std::vector<LineRelationIndex> & lineRelationIndices,
                const ComponentIndexHashMap<std::pair<LineIndex, LineIndex>, Vec3> & interViewLineIncidences,
                int lineConnectedComponentsNum, const ComponentIndexHashMap<LineIndex, int> & lineConnectedComponentIds,
                ComponentIndexHashMap<LineIndex, Line3> & reconstructedLines,
                double constantEtaForFirstLineInEachConnectedComponent,
                bool useWeights){

                ComponentIndexHashMap<LineIndex, int> lineIndexToIds;
                for (int i = 0; i < lineIndices.size(); i++)
                    lineIndexToIds[lineIndices[i]] = i;

                using namespace Eigen;
                Eigen::SparseMatrix<double> A, W;
                VectorXd B;

                // try minimizing ||W(AX-B)||^2

                // pick the first line id in each connected component
                ComponentIndexHashSet<LineIndex> firstLineIndexInConnectedComponents;
                std::set<int> ccIdsRecorded;
                for (auto & lineIndexAndItsCCId : lineConnectedComponentIds) {
                    int ccid = lineIndexAndItsCCId.second;
                    if (ccIdsRecorded.find(ccid) == ccIdsRecorded.end()) { // not recorded yet
                        firstLineIndexInConnectedComponents.insert(lineIndexAndItsCCId.first);
                        ccIdsRecorded.insert(ccid);
                    }
                }

                std::cout << "anchor size: " << firstLineIndexInConnectedComponents.size() << std::endl;
                for (auto & ccId : ccIdsRecorded) {
                    std::cout << "ccid: " << ccId << std::endl;
                }


                // setup matrices
                int n = lineIndices.size(); // var num
                int m = lineRelationIndices.size() + interViewLineIncidences.size();  // cons num

                A.resize(m, n);
                W.resize(m, m);
                B.resize(m);

                // write equations
                int curEquationNum = 0;

                // write intersection/incidence constraint equations in same view
                for (const LineRelationIndex & lri : lineRelationIndices) {
                    auto & lrd = GetData(lri, linesNets);
                    auto & relationCenter = lrd.relationCenter;
                    //auto & weightDistribution = _views.data(lri.viewHandle).lineNet->lineVotingDistribution();

                    auto & topo = linesNets[lri.viewId].topo(lri.handle);
                    auto & camera = views[lri.viewId].camera;
                    LineIndex li1 = { lri.viewId, topo.lowers[0] };
                    LineIndex li2 = { lri.viewId, topo.lowers[1] };

                    int lineId1 = lineIndexToIds[li1];
                    int lineId2 = lineIndexToIds[li2];

                    auto & line1 = GetData(li1, linesNets).line;
                    auto & line2 = GetData(li2, linesNets).line;

                    auto & vp1 = vanishingPoints[line1.claz];
                    auto & vp2 = vanishingPoints[line2.claz];

                    double ratio1 = ComputeDepthRatioOfPointOnSpatialLine(
                        camera.spatialDirection(line1.component.first),
                        camera.spatialDirection(relationCenter), vp1);
                    double ratio2 = ComputeDepthRatioOfPointOnSpatialLine(
                        camera.spatialDirection(line2.component.first),
                        camera.spatialDirection(relationCenter), vp2);

                    if (!core::Contains(firstLineIndexInConnectedComponents, li1) &&
                        !core::Contains(firstLineIndexInConnectedComponents, li2)) {
                        // eta1 * ratio1 - eta2 * ratio2 = 0
                        A.insert(curEquationNum, lineId1) = ratio1;
                        A.insert(curEquationNum, lineId2) = -ratio2;
                        B(curEquationNum) = 0;
                    }
                    else if (core::Contains(firstLineIndexInConnectedComponents, li1)) {
                        // const[eta1] * ratio1 - eta2 * ratio2 = 0 -> 
                        // eta2 * ratio2 = const[eta1] * ratio1
                        A.insert(curEquationNum, lineId2) = ratio2;
                        B(curEquationNum) = constantEtaForFirstLineInEachConnectedComponent * ratio1;
                    }
                    else if (core::Contains(firstLineIndexInConnectedComponents, li2)) {
                        // eta1 * ratio1 - const[eta2] * ratio2 = 0 -> 
                        // eta1 * ratio1 = const[eta2] * ratio2
                        A.insert(curEquationNum, lineId1) = ratio1;
                        B(curEquationNum) = constantEtaForFirstLineInEachConnectedComponent * ratio2;
                    }

                    // set junction weights
                    W.insert(curEquationNum, curEquationNum) = lrd.junctionWeight < MinimumJunctionWeght ? 0.0 : lrd.junctionWeight;

                    curEquationNum++;
                }

                // write inter-view incidence constraints
                for (auto & lineIncidenceAcrossView : interViewLineIncidences) {
                    auto & li1 = lineIncidenceAcrossView.first.first;
                    auto & li2 = lineIncidenceAcrossView.first.second;
                    auto & relationCenter = lineIncidenceAcrossView.second;

                    auto & camera1 = views[li1.viewId].camera;
                    auto & camera2 = views[li2.viewId].camera;

                    int lineId1 = lineIndexToIds[li1];
                    int lineId2 = lineIndexToIds[li2];

                    auto & line1 = GetData(li1, linesNets).line;
                    auto & line2 = GetData(li2, linesNets).line;

                    auto & vp1 = vanishingPoints[line1.claz];
                    auto & vp2 = vanishingPoints[line2.claz];

                    double ratio1 = ComputeDepthRatioOfPointOnSpatialLine(
                        normalize(camera1.spatialDirection(line1.component.first)),
                        normalize(relationCenter), vp1);
                    double ratio2 = ComputeDepthRatioOfPointOnSpatialLine(
                        normalize(camera2.spatialDirection(line2.component.first)),
                        normalize(relationCenter), vp2);

                    if (ratio1 == 0.0 || ratio2 == 0.0) {
                        std::cout << "!!!!!!!ratio is zero!!!!!!!!" << std::endl;
                    }

                    if (!core::Contains(firstLineIndexInConnectedComponents, li1) &&
                        !core::Contains(firstLineIndexInConnectedComponents, li2)) {
                        // eta1 * ratio1 - eta2 * ratio2 = 0
                        A.insert(curEquationNum, lineId1) = ratio1;
                        A.insert(curEquationNum, lineId2) = -ratio2;
                        B(curEquationNum) = 0;
                    }
                    else if (core::Contains(firstLineIndexInConnectedComponents, li1)) {
                        // const[eta1] * ratio1 - eta2 * ratio2 = 0 -> 
                        // eta2 * ratio2 = const[eta1] * ratio1
                        A.insert(curEquationNum, lineId2) = ratio2;
                        B(curEquationNum) = constantEtaForFirstLineInEachConnectedComponent * ratio1;
                    }
                    else if (core::Contains(firstLineIndexInConnectedComponents, li2)) {
                        // eta1 * ratio1 - const[eta2] * ratio2 = 0 -> 
                        // eta1 * ratio1 = const[eta2] * ratio2
                        A.insert(curEquationNum, lineId1) = ratio1;
                        B(curEquationNum) = constantEtaForFirstLineInEachConnectedComponent * ratio2;
                    }

                    double junctionWeight = 5.0;
                    W.insert(curEquationNum, curEquationNum) = junctionWeight;

                    curEquationNum++;
                }

                // solve the equation system
                VectorXd X;
                SparseQR<SparseMatrix<double>, COLAMDOrdering<int>> solver;
                static_assert(!(Eigen::SparseMatrix<double>::IsRowMajor), "COLAMDOrdering only supports column major");
                Eigen::SparseMatrix<double> WA = W * A;
                A.makeCompressed();
                WA.makeCompressed();
                solver.compute(useWeights ? WA : A);
                if (solver.info() != Success) {
                    assert(0);
                    std::cout << "computation error" << std::endl;
                    return;
                }
                VectorXd WB = W * B;
                X = solver.solve(useWeights ? WB : B);
                if (solver.info() != Success) {
                    assert(0);
                    std::cout << "solving error" << std::endl;
                    return;
                }

                // fill back all etas
                int k = 0;
                for (int i = 0; i < lineIndices.size(); i++) {
                    auto & li = lineIndices[i];
                    double eta = X(i);
                    if (firstLineIndexInConnectedComponents.find(li) != firstLineIndexInConnectedComponents.end()) { // is first of a cc
                        eta = constantEtaForFirstLineInEachConnectedComponent;
                        std::cout << "is the " << (++k) << "-th anchor!" << std::endl;
                    }
                    auto & line2 = linesNets[li.viewId].data(li.handle).line;
                    auto & camera = views[li.viewId].camera;
                    Line3 line3 = {
                        normalize(camera.spatialDirection(line2.component.first)),
                        normalize(camera.spatialDirection(line2.component.second))
                    };

                    //std::cout << "eta: " << eta << " --- " << "ccid: " << lineConnectedComponentIds.at(li) << std::endl;

                    double resizeScale = eta / norm(line3.first);
                    line3.first *= resizeScale;
                    line3.second *= (resizeScale *
                        ComputeDepthRatioOfPointOnSpatialLine(line3.first, line3.second, vanishingPoints[line2.claz]));

                    reconstructedLines[li] = line3;
                }


            }


        }


        void EstimateSpatialLineDepths(const Context & context,
            ComponentIndexHashMap<LineIndex, Line3> & reconstructedLines,
            double constantEtaForFirstLineInEachConnectedComponent,
            bool twiceEstimation){

            assert(context.views.size() == context.linesGraphs.size());

            // collect all lines
            std::vector<LineIndex> lineIndices;
            ComponentIndexHashMap<LineIndex, int> lineIndexToIds;
            CollectIndices(context.views, context.linesGraphs, lineIndices, lineIndexToIds);

            // collect all same view constraints
            std::vector<LineRelationIndex> lineRelationIndices; // constraint indices in same views
            for (int i = 0; i < context.views.size(); i++) {
                LineRelationIndex lri;
                lri.viewId = i;
                for (auto & ld : context.linesGraphs[i].elements<1>()) {
                    lri.handle = ld.topo.hd;
                    lineRelationIndices.push_back(lri);
                }
            }

            // reconstruct
            ComponentIndexHashMap<LineIndex, Line3> reconstructedLinesOriginal;
            EstimateSpatialLineDepthsOnce(context.views, context.linesGraphs, context.vanishingPoints, lineIndices, lineRelationIndices,
                context.lineIndencesAcrossViews, context.lineConnectedComponentsNum, context.lineConnectedComponentIds,
                reconstructedLinesOriginal, constantEtaForFirstLineInEachConnectedComponent, true);

            if (!twiceEstimation){
                reconstructedLines = std::move(reconstructedLinesOriginal);
                return;
            }

            // store all line constraints homogeneously
            struct ConstraintBetweenLines{
                enum { InnerView, InterView } type;
                LineRelationIndex lineRelationIndex;
                std::pair<LineIndex, LineIndex> linePairIndex;
                double distance;
            };
            std::vector<ConstraintBetweenLines> homogeneousConstaints;
            homogeneousConstaints.reserve(lineRelationIndices.size() + context.lineIndencesAcrossViews.size());
            for (auto & lri : lineRelationIndices){
                int viewId = lri.viewId;
                // ignore too light constraints, same in ComputeConnectedComponentsUsingRegionLineConstraints
                if (GetData(lri, context.linesGraphs).junctionWeight < MinimumJunctionWeght)
                    continue;
                auto lineHandles = context.linesGraphs.at(viewId).topo(lri.handle).lowers;
                auto & line1 = reconstructedLinesOriginal[LineIndex{ viewId, lineHandles[0] }];
                auto & line2 = reconstructedLinesOriginal[LineIndex{ viewId, lineHandles[1] }];
                auto nearestPoints = DistanceBetweenTwoLines(line1.infiniteLine(), line2.infiniteLine()).second;
                auto c = (nearestPoints.first + nearestPoints.second) / 2.0;
                double distance = abs((nearestPoints.first - nearestPoints.second).dot(normalize(c)))
                    / constantEtaForFirstLineInEachConnectedComponent;
                homogeneousConstaints.push_back(ConstraintBetweenLines{
                    ConstraintBetweenLines::InnerView, lri, std::pair<LineIndex, LineIndex>(), distance
                });
            }
            for (auto & ivl : context.lineIndencesAcrossViews){
                auto & line1 = reconstructedLinesOriginal[ivl.first.first];
                auto & line2 = reconstructedLinesOriginal[ivl.first.second];
                auto nearestPoints = DistanceBetweenTwoLines(line1.infiniteLine(), line2.infiniteLine()).second;
                auto c = (nearestPoints.first + nearestPoints.second) / 2.0;
                double distance = abs((nearestPoints.first - nearestPoints.second).dot(normalize(c)))
                    / constantEtaForFirstLineInEachConnectedComponent;
                homogeneousConstaints.push_back(ConstraintBetweenLines{
                    ConstraintBetweenLines::InterView, LineRelationIndex(), ivl.first, distance
                });
            }

            std::cout << "original line constraints num = " << homogeneousConstaints.size() << std::endl;
            std::vector<int> constraintIds(homogeneousConstaints.size());
            std::iota(constraintIds.begin(), constraintIds.end(), 0);

            // minimum spanning tree
            auto edgeVertsGetter = [&homogeneousConstaints, &context](int cid)->std::pair<LineIndex, LineIndex> {
                std::pair<LineIndex, LineIndex> verts;
                auto & c = homogeneousConstaints[cid];
                if (c.type == ConstraintBetweenLines::InnerView){
                    int viewId = c.lineRelationIndex.viewId;
                    auto lineHandles = context.linesGraphs.at(viewId)
                        .topo(c.lineRelationIndex.handle).lowers;
                    verts = std::make_pair(LineIndex{ viewId, lineHandles[0] }, LineIndex{ viewId, lineHandles[1] });
                }
                else if (c.type == ConstraintBetweenLines::InterView){
                    verts = c.linePairIndex;
                }
                return verts;
            };

            std::vector<int> reservedHomogeneousConstaintsIds;
            reservedHomogeneousConstaintsIds.reserve(homogeneousConstaints.size() / 2);
            core::MinimumSpanningTree(lineIndices.begin(), lineIndices.end(),
                constraintIds.begin(), constraintIds.end(),
                std::back_inserter(reservedHomogeneousConstaintsIds), edgeVertsGetter,
                [&homogeneousConstaints](int cid1, int cid2)->bool {
                return homogeneousConstaints[cid1].distance < homogeneousConstaints[cid2].distance;
            });

            std::cout << "line constraints num after MST = " << reservedHomogeneousConstaintsIds.size() << std::endl;


            // build trimmed line relation indices and inter-view-incidences
            std::vector<LineRelationIndex> trimmedLineRelationIndices;
            trimmedLineRelationIndices.reserve(reservedHomogeneousConstaintsIds.size() / 2);
            ComponentIndexHashMap<std::pair<LineIndex, LineIndex>, Vec3> trimmedInterViewLineIncidences;
            for (int i : reservedHomogeneousConstaintsIds){
                auto & c = homogeneousConstaints[i];
                if (c.type == ConstraintBetweenLines::InnerView){
                    trimmedLineRelationIndices.push_back(c.lineRelationIndex);
                }
                else if (c.type == ConstraintBetweenLines::InterView){
                    trimmedInterViewLineIncidences.emplace(c.linePairIndex, context.lineIndencesAcrossViews.at(c.linePairIndex));
                }
            }

            // reconstruct again
            EstimateSpatialLineDepthsOnce(context.views, context.linesGraphs, context.vanishingPoints, lineIndices, trimmedLineRelationIndices,
                trimmedInterViewLineIncidences, context.lineConnectedComponentsNum, context.lineConnectedComponentIds,
                reconstructedLines, constantEtaForFirstLineInEachConnectedComponent, false);


            // visualize ccids
            // display reconstructed lines
            IF_DEBUG_USING_VISUALIZERS {
                vis::Visualizer3D viz;
                viz << vis::manip3d::SetBackgroundColor(vis::ColorTag::White)
                    << vis::manip3d::SetDefaultColorTable(vis::CreateRandomColorTableWithSize(context.lineConnectedComponentsNum))
                    << vis::manip3d::SetDefaultLineWidth(2.0);
                for (auto & l : reconstructedLines) {
                    viz << core::ClassifyAs(NormalizeLine(l.second), context.lineConnectedComponentIds.at(l.first));
                }
                viz << vis::manip3d::SetDefaultLineWidth(4.0);
                for (auto & c : context.lineIndencesAcrossViews) {
                    auto & line1 = reconstructedLines[c.first.first];
                    auto & line2 = reconstructedLines[c.first.second];
                    auto nearest = DistanceBetweenTwoLines(NormalizeLine(line1), NormalizeLine(line2));
                    viz << vis::manip3d::SetDefaultForegroundColor(vis::ColorTag::Black)
                        << Line3(nearest.second.first.position, nearest.second.second.position);
                }
                viz << vis::manip3d::SetWindowName("not-yet-reconstructed lines with ccids");
                viz << vis::manip3d::Show(false, true);
            }
            IF_DEBUG_USING_VISUALIZERS{
                vis::Visualizer3D viz;
                viz << vis::manip3d::SetBackgroundColor(vis::ColorTag::White)
                    << vis::manip3d::SetDefaultColorTable(vis::CreateRandomColorTableWithSize(context.lineConnectedComponentsNum))
                    << vis::manip3d::SetDefaultLineWidth(4.0);
                for (auto & l : reconstructedLinesOriginal) {
                    viz << core::ClassifyAs(l.second, context.lineConnectedComponentIds.at(l.first));
                }
                viz << vis::manip3d::SetWindowName("reconstructed lines with ccids, 1st time");
                viz << vis::manip3d::Show(false, true);
            }
            IF_DEBUG_USING_VISUALIZERS{
                vis::Visualizer3D viz;
                viz << vis::manip3d::SetBackgroundColor(vis::ColorTag::White)
                    << vis::manip3d::SetDefaultColorTable(vis::CreateRandomColorTableWithSize(context.lineConnectedComponentsNum))
                    << vis::manip3d::SetDefaultLineWidth(4.0);
                for (auto & l : reconstructedLines) {
                    viz << core::ClassifyAs(l.second, context.lineConnectedComponentIds.at(l.first));
                }
                viz << vis::manip3d::SetWindowName("reconstructed lines with ccids, 2nd time");
                viz << vis::manip3d::Show(false, true);
            }

            IF_DEBUG_USING_VISUALIZERS{ // show interview constraints
                vis::Visualizer3D viz;
                viz << vis::manip3d::SetBackgroundColor(vis::ColorTag::White)
                    << vis::manip3d::SetDefaultColorTable(vis::CreateRandomColorTableWithSize(context.lineConnectedComponentsNum))
                    << vis::manip3d::SetDefaultLineWidth(2.0);
                for (auto & l : reconstructedLines) {
                    viz << core::ClassifyAs(l.second, context.lineConnectedComponentIds.at(l.first));
                }
                viz << vis::manip3d::SetDefaultLineWidth(4.0);
                for (auto & c : context.lineIndencesAcrossViews) {
                    auto & line1 = reconstructedLines[c.first.first];
                    auto & line2 = reconstructedLines[c.first.second];
                    auto nearest = DistanceBetweenTwoLines(line1, line2);
                    viz << vis::manip3d::SetDefaultForegroundColor(vis::ColorTag::Black)
                        << Line3(nearest.second.first.position, nearest.second.second.position);
                }
                viz << vis::manip3d::SetWindowName("reconstructed lines with interview constraints");
                viz << vis::manip3d::Show(true, true);
            }

        }

    }
}