
//extern "C" {
//    #include <mosek.h>
//}
#include <Eigen/Dense>
#include <Eigen/Sparse>
//
#include "algorithms.hpp"
#include "containers.hpp"
#include "matlab.hpp"
#include "utilities.hpp"
#include "mixed_graph.hpp"

#include "../vis/visualizers.hpp"
//#include "matlab.hpp"


namespace panoramix {
    namespace core{

        template <class FunT>
        inline void ForeachMixedGraphComponentHandle(const MixedGraph & mg, FunT && fun){
            for (auto & c : mg.components<LineData>()){
                fun(c.topo.hd);
            }
            for (auto & c : mg.components<RegionData>()){
                fun(c.topo.hd);
            }
        }

        template <class FunT>
        inline void ForeachMixedGraphConstraintHandle(const MixedGraph & mg, FunT && fun){
            for (auto & c : mg.constraints<RegionBoundaryData>()){
                fun(c.topo.hd);
            }
            for (auto & c : mg.components<LineRelationData>()){
                fun(c.topo.hd);
            }
            for (auto & c : mg.components<RegionLineConnectionData>()){
                fun(c.topo.hd);
            }
        }

        namespace {


            template <class T>
            inline size_t ElementsNum(const std::vector<T> & v){
                return v.size();
            }

            template <class T>
            inline size_t ElementsNum(const std::vector<std::vector<T>> & v){
                size_t n = 0;
                for (auto & vv : v){
                    n += vv.size();
                }
                return n;
            }

            struct ComparePixelLoc {
                inline bool operator ()(const PixelLoc & a, const PixelLoc & b) const {
                    if (a.x != b.x)
                        return a.x < b.x;
                    return a.y < b.y;
                }
            };


            inline Point2 ToPoint2(const PixelLoc & p) {
                return Point2(p.x, p.y);
            }

            template <class T>
            inline Point2 ToPoint2(const Point<T, 2> & p) {
                return Point2(static_cast<double>(p[0]), static_cast<double>(p[1]));
            }

            template <class T>
            inline Point<float, 2> ToPoint2f(const Point<T, 2> & p) {
                return Point<float, 2>(static_cast<float>(p[0]), static_cast<float>(p[1]));
            }

            inline PixelLoc ToPixelLoc(const Point2 & p) {
                return PixelLoc(static_cast<int>(p[0]), static_cast<int>(p[1]));
            }

            std::pair<double, double> ComputeSpanningArea(const Point2 & a, const Point2 & b, const InfiniteLine2 & line) {
                auto ad = SignedDistanceFromPointToLine(a, line);
                auto bd = SignedDistanceFromPointToLine(b, line);
                auto ap = DistanceFromPointToLine(a, line).second;
                auto bp = DistanceFromPointToLine(b, line).second;
                auto len = norm(ap - bp);
                if (ad * bd >= 0) {
                    return std::make_pair(len * std::abs(ad + bd) / 2.0, len);
                }
                ad = abs(ad);
                bd = abs(bd);
                return std::make_pair((ad * ad + bd * bd) * len / (ad + bd) / 2.0, len);
            }

        }







        float IncidenceJunctionWeight(bool acrossViews){
            return acrossViews ? 10.0f : 7.0f;
        }
        float OutsiderIntersectionJunctionWeight(){
            return 2.0f;
        }
        float ComputeIntersectionJunctionWeightWithLinesVotes(const Mat<float, 3, 2> & v){
            double junctionWeight = 0.0;
            // Y
            double Y = 0.0;
            for (int s = 0; s < 2; s++) {
                Y += v(0, s) * v(1, s) * v(2, s) * DiracDelta(v(0, 1 - s) + v(1, 1 - s) + v(2, 1 - s));
            }

            // W
            double W = 0.0;
            for (int i = 0; i < 3; i++) {
                for (int j = 0; j < 3; j++) {
                    if (i == j)
                        continue;
                    int k = 3 - i - j;
                    for (int s = 0; s < 2; s++) {
                        W += v(i, s) * v(j, 1 - s) * v(k, 1 - s) * DiracDelta(v(i, 1 - s) + v(j, s) + v(k, s));
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
                    K += v(i, 0) * v(i, 1) * v(j, 0) * v(k, 1) * DiracDelta(v(j, 1) + v(k, 0));
                    K += v(i, 0) * v(i, 1) * v(j, 1) * v(k, 0) * DiracDelta(v(j, 0) + v(k, 1));
                }
            }

            // compute X junction
            double X = 0.0;
            for (int i = 0; i < 3; i++) {
                for (int j = 0; j < 3; j++) {
                    if (i == j)
                        continue;
                    int k = 3 - i - j;
                    X += v(i, 0) * v(i, 1) * v(j, 0) * v(j, 1) * DiracDelta(v(k, 0) + v(k, 1));
                }
            }

            // compute T junction
            double T = 0.0;
            for (int i = 0; i < 3; i++) {
                for (int j = 0; j < 3; j++) {
                    if (i == j)
                        continue;
                    int k = 3 - i - j;
                    T += v(i, 0) * v(i, 1) * v(j, 0) * DiracDelta(v(j, 1) + v(k, 0) + v(k, 1));
                    T += v(i, 0) * v(i, 1) * v(j, 1) * DiracDelta(v(j, 0) + v(k, 0) + v(k, 1));
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
                            L += v(i, a) * v(j, b) * DiracDelta(v(i, nota) + v(j, notb) + v(k, 0) + v(k, 1));
                        }
                    }
                }
            }

            //std::cout << " Y-" << Y << " W-" << W << " K-" << K << 
            //    " X-" << X << " T-" << T << " L-" << L << std::endl; 
            static const double threshold = 1e-4;
            if (Y > threshold) {
                junctionWeight += 5.0;
            }
            else if (W > threshold) {
                junctionWeight += 5.0;
            }
            else if (L > threshold) {
                junctionWeight += 4.0;
            }
            else if (K > threshold) {
                junctionWeight += 3.0;
            }
            else if (X > threshold) {
                junctionWeight += 5.0;
            }
            else if (T > threshold) {
                junctionWeight += 0.0;
            }

            return junctionWeight;
        }


        std::vector<Vec3> EstimateVanishingPointsAndClassifyLines(const std::vector<PerspectiveCamera> & cams,
            std::vector<std::vector<Classified<Line2>>> & lineSegments){

            assert(cams.size() == lineSegments.size());
            std::vector<Vec3> lineIntersections;

            int linesNum = 0;
            for (int i = 0; i < cams.size(); i++){
                std::vector<Line2> pureLines(lineSegments[i].size());
                linesNum += lineSegments[i].size();
                for (int k = 0; k < pureLines.size(); k++){
                    pureLines[k] = lineSegments[i][k].component;
                }
                auto inters = ComputeLineIntersections(pureLines, nullptr, true, std::numeric_limits<double>::max());
                // insert line intersections
                for (auto & p : inters){
                    lineIntersections.push_back(normalize(cams[i].spatialDirection(p.value())));
                }
            }

            auto vanishingPoints = FindThreeOrthogonalPrinicipleDirections(lineIntersections);

            // project lines to space
            std::vector<Classified<Line3>> spatialLineSegments;
            spatialLineSegments.reserve(linesNum);
            for (int i = 0; i < cams.size(); i++){
                for (const auto & line : lineSegments[i]) {
                    auto & p1 = line.component.first;
                    auto & p2 = line.component.second;
                    auto pp1 = cams[i].spatialDirection(p1);
                    auto pp2 = cams[i].spatialDirection(p2);
                    Classified<Line3> cline3;
                    cline3.claz = -1;
                    cline3.component = Line3{ pp1, pp2 };
                    spatialLineSegments.push_back(cline3);
                }
            }

            // classify lines
            ClassifyLines(spatialLineSegments, vanishingPoints, M_PI / 3.0, 0.1, 0.8);

            int ii = 0;
            for (int i = 0; i < lineSegments.size(); i++){
                for (int j = 0; j < lineSegments[i].size(); j++){
                    lineSegments[i][j].claz = spatialLineSegments[ii].claz;
                    ii++;
                }
            }

            return vanishingPoints;

        }



        namespace {

            // lines graph in 2d
            struct LineData2D {
                Line2 line;
                int claz;
                template <class Archiver>
                void serialize(Archiver & ar) {
                    ar(line, claz);
                }
            };
            struct LineRelationData2D {
                Point2 relationCenter;
                float junctionWeight;
                LineRelationData::Type type;
                template <class Archiver>
                void serialize(Archiver & ar) {
                    ar(relationCenter, junctionWeight, type);
                }
            };
            using LinesGraph2D = HomogeneousGraph02<LineData2D, LineRelationData2D>;
            using LineHandle2D = HandleOfTypeAtLevel<LinesGraph2D, 0>;
            using LineRelationHandle2D = HandleOfTypeAtLevel<LinesGraph2D, 1>;

            LinesGraph2D CreateLinesGraph2D(const std::vector<Classified<Line2>> & lines,
                const std::vector<HPoint2> & vps,
                double intersectionDistanceThreshold,
                double incidenceDistanceAlongDirectionThreshold,
                double incidenceDistanceVerticalDirectionThreshold){

                /// TODO include all line relations is NOT STABLE!

                LinesGraph2D graph;

                // insert lines
                std::vector<LineHandle2D> handles;
                handles.reserve(lines.size());
                graph.internalElements<0>().reserve(lines.size());

                for (auto & line : lines){
                    LineData2D ld;
                    ld.line = line.component;
                    ld.claz = line.claz;
                    handles.push_back(graph.add(std::move(ld)));
                }

                // construct incidence/intersection relations
                auto & linesData = graph.internalElements<0>();
                for (int i = 0; i < linesData.size(); i++){
                    auto & linei = linesData[i].data.line;
                    int clazi = linesData[i].data.claz;

                    for (int j = i + 1; j < linesData.size(); j++){
                        auto & linej = linesData[j].data.line;
                        int clazj = linesData[j].data.claz;

                        auto nearest = DistanceBetweenTwoLines(linei, linej);
                        double d = nearest.first;

                        if (clazi == clazj && clazi >= 0 ){ // incidences for classified lines
                            auto conCenter = (nearest.second.first.position + nearest.second.second.position) / 2.0;
                            auto conDir = (nearest.second.first.position - nearest.second.second.position);

                            auto & vp = vps[clazi];

                            if (Distance(vp.value(), conCenter) < intersectionDistanceThreshold)
                                continue;

                            auto dir = normalize((vp - HPoint2(conCenter)).numerator);
                            double dAlong = abs(conDir.dot(dir));
                            double dVert = sqrt(Square(norm(conDir)) - dAlong*dAlong);

                            if (dAlong < incidenceDistanceAlongDirectionThreshold &&
                                dVert < incidenceDistanceVerticalDirectionThreshold){
                                LineRelationData2D lrd;
                                lrd.type = LineRelationData::Type::Incidence;
                                lrd.relationCenter = conCenter;

                                if (HasValue(lrd.relationCenter, IsInfOrNaN<double>))
                                    continue;

                                graph.add<1>({ handles[i], handles[j] }, std::move(lrd));
                            }
                        }
                        else if(clazi != clazj && clazi >= 0 && clazj >= 0) { // intersections for classified lines
                            if (d < intersectionDistanceThreshold){
                                auto conCenter = HPointFromVector(GetCoeffs(linei.infiniteLine())
                                    .cross(GetCoeffs(linej.infiniteLine()))).value();

                                assert(conCenter != Point2(0.0, 0.0));

                                if (DistanceFromPointToLine(conCenter, linei).first > intersectionDistanceThreshold * 4 ||
                                    DistanceFromPointToLine(conCenter, linej).first > intersectionDistanceThreshold * 4)
                                    continue;

                                LineRelationData2D lrd;
                                lrd.type = LineRelationData::Type::Intersection;
                                lrd.relationCenter = conCenter;

                                if (HasValue(lrd.relationCenter, IsInfOrNaN<double>))
                                    continue;

                                graph.add<1>({ handles[i], handles[j] }, std::move(lrd));
                            }
                        }
                        else if (clazi == clazj && clazi == -1){ // incidence for unclassified lines
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
                            lrd.type = LineRelationData::Type::Incidence;
                            lrd.relationCenter = (nearest.second.first.position + nearest.second.second.position) / 2.0;

                            if (HasValue(lrd.relationCenter, IsInfOrNaN<double>))
                                continue;

                            graph.add<1>({ handles[i], handles[j] }, std::move(lrd));
                        }

                    }
                }

                static const double angleThreshold = M_PI / 32;
                static const double sigma = 0.1;

                enum LineVotingDirection : int {
                    TowardsVanishingPoint = 0,
                    TowardsOppositeOfVanishingPoint = 1
                };
                enum class JunctionType : int {
                    L, T, Y, W, X
                };

                Box2 linesBBox = BoundingBoxOfContainer(lines);
                linesBBox.expand(linesBBox.outerSphere().radius / 2.0);

                // compute junction weights
                for (auto & lr : graph.elements<1>()){
                    auto & lrd = lr.data;
                    if (lrd.type == LineRelationData::Type::Incidence){
                        lrd.junctionWeight = IncidenceJunctionWeight(false);
                    }
                    else if (lrd.type == LineRelationData::Type::Intersection){
                        Mat<float, 3, 2> votingData;
                        std::fill(std::begin(votingData.val), std::end(votingData.val), 0);

                        for (auto & ld : graph.elements<0>()){
                            auto & line = ld.data.line;
                            int claz = ld.data.claz;
                            if (claz == -1)
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

                            if (AngleBetweenDirections(center2vp, line.direction()) < M_PI_2){ // first-second-vp
                                votingData(claz, TowardsVanishingPoint) += angleScore * line.length() * (1 - projRatio);
                                votingData(claz, TowardsOppositeOfVanishingPoint) += angleScore * line.length() * projRatio;
                            }
                            else{ // vp-first-second
                                votingData(claz, TowardsOppositeOfVanishingPoint) += angleScore * line.length() * (1 - projRatio);
                                votingData(claz, TowardsVanishingPoint) += angleScore * line.length() * projRatio;
                            }
                        }


                        auto p = ToPixelLoc(lrd.relationCenter);
                        if (core::IsBetween(lrd.relationCenter[0], 0, linesBBox.size()[0] - 1) &&
                            core::IsBetween(lrd.relationCenter[1], 0, linesBBox.size()[1] - 1)){
                            lrd.junctionWeight = ComputeIntersectionJunctionWeightWithLinesVotes(
                                votingData);
                        }
                        else{
                            lrd.junctionWeight = OutsiderIntersectionJunctionWeight();
                        }
                    }
                    else{
                        assert(false && "invalid branch!");
                    }
                }

                return graph;

            }


            template <class T, int N>
            inline Line<T, N> NormalizeLine(const Line<T, N> & l) {
                return Line<T, N>(normalize(l.first), normalize(l.second));
            }
        }


        void AppendLines(MixedGraph & mg, const std::vector<Classified<Line2>> & lineSegments,
            const PerspectiveCamera & cam,
            const std::vector<Vec3> & vps,
            double intersectionAngleThreshold /*= 0.04*/,
            double incidenceAngleAlongDirectionThreshold /*= 0.1*/,
            double incidenceAngleVerticalDirectionThreshold /*= 0.02*/,
            double interViewIncidenceAngleAlongDirectionThreshold /*= 0.15*/, // for new line-line incidence recognition
            double interViewIncidenceAngleVerticalDirectionThreshold /*= 0.03*/){

            if (lineSegments.empty())
                return;

            assert(mg.internalComponents<RegionData>().empty() && "Regions must be added AFTER all lines!!!");

            std::vector<HPoint2> vps2d(vps.size());
            for (int i = 0; i < vps.size(); i++){
                vps2d[i] = cam.screenProjectionInHPoint(vps[i]);
            }
            LinesGraph2D graph2d = CreateLinesGraph2D(lineSegments, vps2d, 
                cam.focal() * intersectionAngleThreshold,
                cam.focal() * incidenceAngleAlongDirectionThreshold,
                cam.focal() * incidenceAngleVerticalDirectionThreshold);

            mg.internalComponents<LineData>().reserve(mg.internalComponents<LineData>().size() + 
                graph2d.internalElements<0>().size());
            mg.internalConstraints<LineRelationData>().reserve(mg.internalConstraints<LineRelationData>().size() + 
                graph2d.internalElements<1>().size());

            assert(graph2d.isDense());

            std::unordered_set<LineHandle> newLineHandles;
            std::unordered_map<LineHandle2D, LineHandle> lh2dToLh;

            for (auto & l2d : graph2d.elements<0>()){
                auto & ld2d = l2d.data;
                LineData ld;
                ld.initialClaz = ld2d.claz;
                ld.line.first = normalize(cam.spatialDirection(ld2d.line.first));
                ld.line.second = normalize(cam.spatialDirection(ld2d.line.second));
                lh2dToLh[l2d.topo.hd] = mg.addComponent(std::move(ld));
                newLineHandles.insert(lh2dToLh[l2d.topo.hd]);
            }

            for (auto & r2d : graph2d.elements<1>()){
                auto & rd2d = r2d.data;
                LineRelationData rd;
                rd.type = rd2d.type;
                rd.junctionWeight = rd2d.junctionWeight;
                rd.normalizedRelationCenter = normalize(cam.spatialDirection(rd2d.relationCenter));
                mg.addConstraint(std::move(rd), lh2dToLh.at(r2d.topo.lowers.front()), lh2dToLh.at(r2d.topo.lowers.back()));
            }

            // line incidence relations between old lines and new lines
            std::map<std::pair<LineHandle, LineHandle>, Vec3> interViewLineIncidences;

            // build rtree for lines
            auto lookupLineNormal = [&mg](const LineHandle & li) -> Box3 {
                auto normal = mg.data(li).line.first.cross(mg.data(li).line.second);
                Box3 b = BoundingBox(normalize(normal));
                static const double s = 0.2;
                b.minCorner = b.minCorner - Vec3(s, s, s);
                b.maxCorner = b.maxCorner + Vec3(s, s, s);
                return b;
            };

            RTreeWrapper<LineHandle, decltype(lookupLineNormal)> linesRTree(lookupLineNormal);
            for (auto & l : mg.components<LineData>()) {
                linesRTree.insert(l.topo.hd);
            }

            // recognize incidence constraints between lines of different views
            for (auto & l : mg.components<LineData>()) {
                auto lh = l.topo.hd;
                if (Contains(newLineHandles, lh))
                    continue; // stick old lh, find new lh
                auto & lineData = l.data;
                linesRTree.search(lookupLineNormal(lh),
                    [interViewIncidenceAngleAlongDirectionThreshold, interViewIncidenceAngleVerticalDirectionThreshold,
                    &l, &mg, &interViewLineIncidences, &newLineHandles](const LineHandle & relatedLh) -> bool {
                    if (!Contains(newLineHandles, relatedLh))
                        return true;

                    assert(l.topo.hd < relatedLh);

                    auto & line1 = l.data.line;
                    auto & line2 = mg.data(relatedLh).line;
                    
                    if (l.data.initialClaz != mg.data(relatedLh).initialClaz) // only incidence relations are recognized here
                        return true;

                    auto normal1 = normalize(line1.first.cross(line1.second));
                    auto normal2 = normalize(line2.first.cross(line2.second));

                    if (std::min(std::abs(AngleBetweenDirections(normal1, normal2)),
                        std::abs(AngleBetweenDirections(normal1, -normal2))) <
                        interViewIncidenceAngleVerticalDirectionThreshold) {

                        auto nearest = DistanceBetweenTwoLines(NormalizeLine(line1), NormalizeLine(line2));
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
            for (auto & ivli : interViewLineIncidences){
                LineRelationData lrd;
                lrd.type = LineRelationData::Incidence;
                lrd.junctionWeight = IncidenceJunctionWeight(true);
                lrd.normalizedRelationCenter = ivli.second;
                mg.addConstraint(std::move(lrd), ivli.first.first, ivli.first.second);
            }

        }


        namespace {


            template <class T, int N>
            std::vector<Point<T, N>> NormalizedSamplePointsOnLine(const Line<T, N> & line, const T & stepLen){
                T len = line.length();
                auto dir = normalize(line.direction());
                std::vector<Point<T, N>> pts;
                pts.reserve(len / stepLen);
                for (T l = 0; l <= len; l += stepLen){
                    pts.push_back(normalize(line.first + dir * l));
                }
                return pts;
            }

            std::vector<RegionHandle> CollectRegionsFromSegmentation(MixedGraph & mg, const Imagei & segmentedRegions, const PerspectiveCamera & cam){
                double minVal, maxVal;
                std::tie(minVal, maxVal) = MinMaxValOfImage(segmentedRegions);
                assert(minVal == 0.0);

                int regionNum = static_cast<int>(maxVal)+1;
                mg.internalComponents<RegionData>().reserve(mg.internalComponents<RegionData>().size() + regionNum);
                std::vector<RegionHandle> regionHandles(regionNum);

                // calculate contours and tangential projected areas for each region data
                for (int i = 0; i < regionNum; i++){
                    Image regionMask = (segmentedRegions == i);

                    // find contour of the region
                    std::vector<std::vector<PixelLoc>> contours;

                    cv::findContours(regionMask, contours, CV_RETR_EXTERNAL, CV_CHAIN_APPROX_SIMPLE); // CV_RETR_EXTERNAL: get only the outer contours
                    std::sort(contours.begin(), contours.end(),
                        [](const std::vector<PixelLoc> & ca, const std::vector<PixelLoc> & cb){return ca.size() > cb.size(); });

                    if (contours.size() >= 2){
                        std::cout << "multiple contours for one region in perspective projection!";
                    }

                    auto iter = std::find_if(contours.begin(), contours.end(), [](const std::vector<PixelLoc> & c){return c.size() <= 2; });
                    contours.erase(iter, contours.end());

                    //assert(!contours.empty() && "no contour? impossible~");
                    if (contours.empty() || contours.front().size() <= 2){
                        continue;
                    }

                    for (auto & cs : contours){
                        assert(cs.size() > 2);
                    }

                    RegionData rd;
                    Vec3 center(0, 0, 0);
                    rd.normalizedContours.resize(contours.size());
                    for (int k = 0; k < contours.size(); k++){
                        rd.normalizedContours[k].reserve(contours[k].size());
                        for (auto & p : contours[k]){
                            rd.normalizedContours[k].push_back(normalize(cam.spatialDirection(p)));
                            center += rd.normalizedContours[k].back();
                        }
                    }
                    rd.normalizedCenter = normalize(center);

                    // project the contours onto tangential plane 
                    Vec3 x, y;
                    std::tie(x, y) = core::ProposeXYDirectionsFromZDirection(rd.normalizedCenter);
                    rd.area = 0.0;
                    for (int k = 0; k < contours.size(); k++){
                        std::vector<Point2f> tangentialContour;
                        tangentialContour.reserve(contours[k].size());
                        for (const Vec3 & d : rd.normalizedContours[k]){
                            tangentialContour.emplace_back((d - rd.normalizedCenter).dot(x), (d - rd.normalizedCenter).dot(y));
                        }
                        double a = cv::contourArea(tangentialContour);
                        assert(a >= 0.0);
                        rd.area += a;
                    }

                   /* if (rd.area < 1e-6)
                        continue;*/

                    regionHandles[i] = mg.addComponent(std::move(rd));
                }

                return regionHandles;
            }

            std::vector<RegionHandle> CollectRegionsFromSegmentation(MixedGraph & mg, const Imagei & segmentedRegions, const PanoramicCamera & cam){
                double minVal, maxVal;
                std::tie(minVal, maxVal) = MinMaxValOfImage(segmentedRegions);
                assert(minVal == 0.0);

                int regionNum = static_cast<int>(maxVal)+1;

                mg.internalComponents<RegionData>().reserve(mg.internalComponents<RegionData>().size() + regionNum);
                std::vector<RegionHandle> regionHandles(regionNum);

                // calculate contours and tangential projected areas for each region data
                for (int i = 0; i < regionNum; i++){

                    Image regionMask = (segmentedRegions == i);

                    // find contour of the region
                    std::vector<std::vector<PixelLoc>> contours;
                    cv::findContours(regionMask, contours, CV_RETR_EXTERNAL, CV_CHAIN_APPROX_SIMPLE); // CV_RETR_EXTERNAL: get only the outer contours
                    if (contours.empty()){
                        continue;
                    }

                    Vec3 centerDirection(0, 0, 0);
                    std::vector<Vec3> directions;
                    directions.reserve(ElementsNum(contours));
                    for (auto & cs : contours){
                        for (auto & c : cs){
                            directions.push_back(normalize(cam.spatialDirection(c)));
                            centerDirection += directions.back();
                        }
                    }
                    centerDirection /= norm(centerDirection);
                    // get max angle distance from center direction
                    double radiusAngle = 0.0;
                    for (auto & d : directions){
                        double a = AngleBetweenDirections(centerDirection, d);
                        if (radiusAngle < a){
                            radiusAngle = a;
                        }
                    }

                    // perform a more precise sample !
                    int newSampleSize = cam.focal() * radiusAngle * 2 + 2;
                    PartialPanoramicCamera sCam(newSampleSize, newSampleSize, cam.focal(), cam.eye(), centerDirection, 
                        ProposeXYDirectionsFromZDirection(centerDirection).second);
                    Imagei sampledSegmentedRegions = MakeCameraSampler(sCam, cam)(segmentedRegions);

                    // collect better contours
                    contours.clear();
                    regionMask = (sampledSegmentedRegions == i);
                    cv::findContours(regionMask, contours, CV_RETR_EXTERNAL, CV_CHAIN_APPROX_SIMPLE); // CV_RETR_EXTERNAL: get only the outer contours
                    std::sort(contours.begin(), contours.end(),
                        [](const std::vector<PixelLoc> & ca, const std::vector<PixelLoc> & cb){return ca.size() > cb.size(); });

                    if (contours.size() >= 2){
                        //std::cout << "multiple contours for one region in projection!";
                    }

                    auto iter = std::find_if(contours.begin(), contours.end(), [](const std::vector<PixelLoc> & c){return c.size() <= 2; });
                    contours.erase(iter, contours.end());

                    //assert(!contours.empty() && "no contour? impossible~");
                    if (contours.empty() || contours.front().size() <= 2){
                        continue;
                    }


                    RegionData rd;
                    Vec3 center(0, 0, 0);
                    rd.normalizedContours.resize(contours.size());
                    rd.area = 0.0;
                    for (int k = 0; k < contours.size(); k++){
                        rd.normalizedContours[k].reserve(contours[k].size());
                        for (auto & p : contours[k]){
                            rd.normalizedContours[k].push_back(normalize(sCam.spatialDirection(p)));
                            center += rd.normalizedContours[k].back();
                        }
                        std::vector<Point2f> contourf(contours[k].size());
                        for (int kk = 0; kk < contours[k].size(); kk++){
                            contourf[kk] = vec_cast<float>(contours[k][kk]);
                        }
                        rd.area += cv::contourArea(contourf);
                    }
                    rd.normalizedCenter = normalize(center);
                    //std::cout << "region area: " << rd.area << std::endl;
                                       
                    regionHandles[i] = mg.addComponent(std::move(rd));
                }

                return regionHandles;
            }


            template <class T>
            inline std::pair<T, T> MakeOrderedPair(const T & a, const T & b) {
                return a < b ? std::make_pair(a, b) : std::make_pair(b, a);
            }

            void FindContoursOfRegionsAndBoundaries(const Imagei & segRegions, int regionNum,
                std::map<std::pair<int, int>, std::vector<std::vector<PixelLoc>>> & boundaryEdges) {
                static const int detectionSize = 3;
                std::map<std::pair<int, int>, std::set<PixelLoc, ComparePixelLoc>> boundaryPixels;

                int width = segRegions.cols;
                int height = segRegions.rows;
                for (int y = 0; y < height - 1; y++) {
                    for (int x = 0; x < width - 1; x++) {
                        int originalRegionId = segRegions(PixelLoc(x, y));
                        for (int xx = std::max(x - detectionSize, 0); xx <= std::min(x + detectionSize, width - 1); xx++){
                            for (int yy = std::max(y - detectionSize, 0); yy <= std::min(y + detectionSize, height - 1); yy++){
                                int regionIdHere = segRegions(PixelLoc(xx, yy));
                                if (originalRegionId != regionIdHere){
                                    boundaryPixels[MakeOrderedPair(originalRegionId, regionIdHere)].insert(PixelLoc((x + xx) / 2, (y + yy) / 2));
                                }
                            }
                        }
                    }
                }

                for (auto & bpp : boundaryPixels) {
                    int rid1 = bpp.first.first;
                    int rid2 = bpp.first.second;
                    auto & pixels = bpp.second;

                    if (pixels.empty())
                        continue;

                    PixelLoc p = *pixels.begin();

                    static const int xdirs[] = { 1, 0, -1, 0, -1, 1, 1, -1, 0, 0, 2, -2 };
                    static const int ydirs[] = { 0, 1, 0, -1, 1, -1, 1, -1, 2, -2, 0, 0 };

                    std::vector<std::vector<PixelLoc>> edges;
                    edges.push_back({ p });
                    pixels.erase(p);

                    while (true) {
                        auto & curEdge = edges.back();
                        auto & curTail = curEdge.back();

                        bool foundMore = false;
                        for (int i = 0; i < std::distance(std::begin(xdirs), std::end(xdirs)); i++) {
                            PixelLoc next = curTail;
                            next.x += xdirs[i];
                            next.y += ydirs[i];
                            if (!IsBetween(next.x, 0, width - 1) || !IsBetween(next.y, 0, height - 1))
                                continue;
                            if (pixels.find(next) == pixels.end()) // not a boundary pixel or already recorded
                                continue;

                            curEdge.push_back(next);
                            pixels.erase(next);
                            foundMore = true;
                            break;
                        }

                        if (!foundMore) {
                            // simplify current edge
                            if (edges.back().size() <= 1) {
                                edges.pop_back();
                            }
                            else {
                                bool closed = Distance(edges.back().front(), edges.back().back()) <= 1.5;
                                cv::approxPolyDP(edges.back(), edges.back(), 2, closed);
                                if (edges.back().size() <= 1)
                                    edges.pop_back();
                            }

                            if (pixels.empty()) { // no more pixels
                                break;
                            }
                            else { // more pixels
                                PixelLoc p = *pixels.begin();
                                edges.push_back({ p });
                                pixels.erase(p);
                            }
                        }
                    }

                    if (!edges.empty()) {
                        boundaryEdges[MakeOrderedPair(rid1, rid2)] = edges;
                    }
                }
            }


            template <class CameraT>
            void AppendRegionsTemplate(MixedGraph & mg, const Imagei & segmentedRegions, const CameraT & cam,
                double samplingStepAngleOnBoundary, double samplingStepAngleOnLine){

                auto regionHandles = CollectRegionsFromSegmentation(mg, segmentedRegions, cam);
                int regionNum = regionHandles.size();

                // add region boundary constraints
                std::map<std::pair<int, int>, std::vector<std::vector<PixelLoc>>> boundaryEdges;
                FindContoursOfRegionsAndBoundaries(segmentedRegions, regionNum, boundaryEdges);

                for (auto & bep : boundaryEdges) {
                    auto & rids = bep.first;
                    auto & edges = bep.second;

                    if (regionHandles[rids.first].invalid() || regionHandles[rids.second].invalid())
                        continue;

                    RegionBoundaryData bd;
                    bd.normalizedEdges.resize(edges.size());
                    for (int k = 0; k < edges.size(); k++){
                        bd.normalizedEdges[k].reserve(edges[k].size());
                        for (auto & p : edges[k]){
                            bd.normalizedEdges[k].push_back(normalize(cam.spatialDirection(p)));
                        }
                    }

                    bd.length = 0;
                    for (auto & e : bd.normalizedEdges) {
                        assert(!e.empty() && "edges should never be empty!");
                        // get edge point projections
                        for (int i = 0; i < e.size() - 1; i++) {
                            bd.length += AngleBetweenDirections(e[i], e[i + 1]);
                        }
                    }

                    bd.normalizedSampledPoints.resize(bd.normalizedEdges.size());
                    for (int k = 0; k < bd.normalizedEdges.size(); k++){
                        ForeachCompatibleWithLastElement(bd.normalizedEdges[k].begin(), bd.normalizedEdges[k].end(),
                            std::back_inserter(bd.normalizedSampledPoints[k]),
                            [samplingStepAngleOnBoundary](const Vec3 & a, const Vec3 & b){
                            return AngleBetweenDirections(a, b) >= samplingStepAngleOnBoundary;
                        });
                    }

                    mg.addConstraint(std::move(bd), regionHandles[rids.first], regionHandles[rids.second]);
                }


                // add region-line connections
                std::map<std::pair<RegionHandle, LineHandle>, std::vector<Vec3>> regionLineConnections;

                for (auto & ld : mg.components<LineData>()){
                    auto & line = ld.data.line;
                    auto samples = NormalizedSamplePointsOnLine(line, samplingStepAngleOnLine);
                    for (int k = 0; k < samples.size(); k++){
                        if (!cam.isVisibleOnScreen(samples[k]))
                            continue;
                        PixelLoc originalP = ToPixelLoc(cam.screenProjection(samples[k]));
                        // collect neighbors
                        static const int sz = 3;
                        for (int x = originalP.x - sz; x <= originalP.x + sz; x++){
                            for (int y = originalP.y - sz; y <= originalP.y + sz; y++){
                                PixelLoc p(x, y);
                                if (p.x < 0 || p.x >= segmentedRegions.cols || p.y < 0 || p.y >= segmentedRegions.rows)
                                    continue;
                                int regionId = segmentedRegions(p);
                                auto rh = regionHandles[regionId];
                                if (rh.invalid())
                                    continue;
                                regionLineConnections[std::make_pair(rh, ld.topo.hd)].push_back(samples[k]);
                            }
                        }                      
                        
                    }
                }

                for (auto & rlc : regionLineConnections){
                    RegionLineConnectionData rlcd;
                    rlcd.normalizedAnchors = std::move(rlc.second);
                    mg.addConstraint(std::move(rlcd), rlc.first.first, rlc.first.second);
                }

            }


        }

        void AppendRegions(MixedGraph & mg, const Imagei & segmentedRegions, const PerspectiveCamera & cam,
            double samplingStepAngleOnBoundary, double samplingStepAngleOnLine){
            AppendRegionsTemplate(mg, segmentedRegions, cam, samplingStepAngleOnBoundary, samplingStepAngleOnLine);
        }
        void AppendRegions(MixedGraph & mg, const Imagei & segmentedRegions, const PanoramicCamera & cam,
            double samplingStepAngleOnBoundary, double samplingStepAngleOnLine){
            AppendRegionsTemplate(mg, segmentedRegions, cam, samplingStepAngleOnBoundary, samplingStepAngleOnLine);
        }











        MixedGraphPropertyTable MakeMixedGraphPropertyTable(const MixedGraph & mg, const std::vector<Vec3> & vps) {
            auto compProps = MakeHandledTableForComponents<MixedGraphComponentProperty>(mg);
            auto consProps = MakeHandledTableForConstraints<MixedGraphConstraintProperty>(mg);
            for (auto & l : mg.internalComponents<LineData>()){
                compProps[l.topo.hd].orientationClaz = l.data.initialClaz;
                compProps[l.topo.hd].orientationNotClaz = -1;
            }
            for (auto & r : mg.internalComponents<RegionData>()){
                compProps[r.topo.hd].orientationClaz = compProps[r.topo.hd].orientationNotClaz = -1;
            }
            for (auto & dtable : consProps.data){
                for (auto & d : dtable){
                    d.used = true;
                }
            }
            return MixedGraphPropertyTable{ vps, std::move(compProps), std::move(consProps) };
        }





        namespace {

            int SwappedComponent(const Vec3 & orientation){
                for (int i = 0; i < 2; i++){
                    if (abs(orientation[i]) >= 1e-8){
                        return i;
                    }
                }
                return 2;
            }

            struct InitializeVariablesForEachHandle {
                const MixedGraph & mg;
                MixedGraphPropertyTable & props;

                void operator()(const LineHandle & lh) const {
                    auto & lp = props.componentProperties[lh];
                    if (lp.orientationClaz == -1){
                        // (1/cornerDepth1, 1/cornerDepth2) for LineFree
                        lp.variables = { 1.0, 1.0 };
                    }
                    else{
                        // 1/centerDepth for LineOriented,
                        lp.variables = { 1.0 };
                    }
                }

                void operator()(const RegionHandle & rh) const {
                    auto & rd = mg.data(rh);
                    auto & rp = props.componentProperties[rh];
                    if (rp.orientationClaz == -1 && rp.orientationNotClaz == -1){
                        // (a, b, c) for RegionFree ax+by+c=1, 
                        rp.variables = { rd.normalizedCenter[0], rd.normalizedCenter[1], rd.normalizedCenter[2] };
                    }
                    else if (rp.orientationClaz == -1 && rp.orientationNotClaz >= 0){
                        // (a, b), {or (b, c) or (a, c)} for RegionAlongFixedAxis  ax+by+c=1,
                        Vec3 anotherAxis = rd.normalizedCenter.cross(normalize(props.vanishingPoints[rp.orientationNotClaz]));
                        Vec3 trueNormal = props.vanishingPoints[rp.orientationNotClaz].cross(anotherAxis);
                        Plane3 plane(rd.normalizedCenter, trueNormal);
                        auto eq = Plane3ToEquation(plane);
                        int c = SwappedComponent(normalize(props.vanishingPoints[rp.orientationNotClaz]));
                        std::swap(eq[c], eq[2]);
                        rp.variables = { eq[0], eq[1] };
                    }
                    else {
                        // 1/centerDepth for RegionWithFixedNormal
                        rp.variables = { 1.0 };
                    }
                }
            };

            
        }


        void InitializeVariables(const MixedGraph & mg, MixedGraphPropertyTable & props){
            ForeachMixedGraphComponentHandle(mg, InitializeVariablesForEachHandle{ mg, props });
        }


        Line3 Instance(const MixedGraph & mg, const MixedGraphPropertyTable & props, const LineHandle & lh){
            auto & ld = mg.data(lh);
            auto & lp = props.componentProperties.at(lh);
            if (lp.orientationClaz >= 0){
                assert(lp.variables.size() == 1);
                InfiniteLine3 infLine(normalize(ld.line.center()) / lp.variables[0], props.vanishingPoints[lp.orientationClaz]);
                return Line3(DistanceBetweenTwoLines(InfiniteLine3(Point3(0, 0, 0), normalize(ld.line.first)), infLine).second.second,
                    DistanceBetweenTwoLines(InfiniteLine3(Point3(0, 0, 0), normalize(ld.line.second)), infLine).second.second);
            }
            else /*if (line.type == MGUnary::LineFree)*/{
                assert(lp.variables.size() == 2);
                /*           | sin(theta) | | p | | q |
                len:---------------------------------   => 1/len: [(1/q)sin(phi) - (1/p)sin(phi-theta)] / sin(theta)
                | p sin(phi) - q sin(phi - theta) |
                */
                // variables[0] -> 1/p
                // variables[1] -> 1/q
                return Line3(normalize(ld.line.first) / lp.variables[0], normalize(ld.line.second) / lp.variables[1]);
            }
        }



        Plane3 Instance(const MixedGraph & mg, const MixedGraphPropertyTable & props, const RegionHandle & rh){
            auto & rd = mg.data(rh);
            auto & rp = props.componentProperties.at(rh);
            if (rp.orientationClaz >= 0){
                assert(rp.variables.size() == 1);
                return Plane3(rd.normalizedCenter / rp.variables[0], props.vanishingPoints[rp.orientationClaz]);
            }
            else if (rp.orientationClaz == -1 && rp.orientationNotClaz >= 0){
                assert(rp.variables.size() == 2);
                double vs[] = { rp.variables[0], rp.variables[1], 0.0 }; // fake vs
                // v1 * o1 + v2 * o2 + v3 * o3 = 0, since region.orientation is orthogonal to normal
                auto orientation = normalize(props.vanishingPoints[rp.orientationNotClaz]);
                int c = SwappedComponent(orientation);
                std::swap(orientation[c], orientation[2]); // now fake orientation
                vs[2] = (-vs[0] * orientation[0] - vs[1] * orientation[1])
                    / orientation[2];
                std::swap(vs[c], vs[2]); // now real vs
                return Plane3FromEquation(vs[0], vs[1], vs[2]);
            }
            else /*if (region.type == MGUnary::RegionWithFixedNormal)*/{
                assert(rp.variables.size() == 3);
                return Plane3FromEquation(rp.variables[0], rp.variables[1], rp.variables[2]);
            }
        }



        std::vector<double> VariableCoefficientsForInverseDepthAtDirection(const MixedGraph & mg,
            const MixedGraphPropertyTable & props, const Vec3 & direction, const LineHandle & lh){
            auto & ld = mg.data(lh);
            auto & lp = props.componentProperties.at(lh);
            if (lp.orientationClaz >= 0){
                assert(lp.variables.size() == 1);
                InfiniteLine3 infLine(normalize(ld.line.center()), props.vanishingPoints[lp.orientationClaz]);
                 // variable is 1.0/centerDepth
                 // corresponding coeff is 1.0/depthRatio
                 // so that 1.0/depth = 1.0/centerDepth * 1.0/depthRatio -> depth = centerDepth * depthRatio
                 double depthRatio = norm(DistanceBetweenTwoLines(InfiniteLine3(Point3(0, 0, 0), direction), infLine).second.first);
                 return std::vector<double>{1.0 / depthRatio};
             }
             else /*if(u.type == MGUnary::LineFree)*/{
                 assert(lp.variables.size() == 2);
                 double theta = AngleBetweenDirections(normalize(ld.line.first), normalize(ld.line.second));
                 double phi = AngleBetweenDirections(normalize(ld.line.first), direction);
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


        std::vector<double> VariableCoefficientsForInverseDepthAtDirection(const MixedGraph & mg,
            const MixedGraphPropertyTable & props, const Vec3 & direction, const RegionHandle & rh){
            auto & rd = mg.data(rh);
            auto & rp = props.componentProperties.at(rh);
            if (rp.orientationClaz >= 0){
                assert(rp.variables.size() == 1);
                Plane3 plane(rd.normalizedCenter, props.vanishingPoints[rp.orientationClaz]);
                // variable is 1.0/centerDepth
                // corresponding coeff is 1.0/depthRatio
                // so that 1.0/depth = 1.0/centerDepth * 1.0/depthRatio -> depth = centerDepth * depthRatio
                double depthRatio = norm(IntersectionOfLineAndPlane(InfiniteLine3(Point3(0, 0, 0), direction), plane).position);
                return std::vector<double>{1.0 / depthRatio};
            }
            else if (rp.orientationClaz == -1 && rp.orientationNotClaz >= 0){
                assert(rp.variables.size() == 2);
                // depth = 1.0 / (ax + by + cz) where (x, y, z) = direction, (a, b, c) = variables
                // -> 1.0/depth = ax + by + cz
                // -> 1.0/depth = ax + by + (- o1a - o2b)/o3 * z = a(x-o1*z/o3) + b(y-o2*z/o3)
                auto orientation = normalize(props.vanishingPoints[rp.orientationNotClaz]);
                int c = SwappedComponent(orientation);
                Vec3 forientation = orientation;
                std::swap(forientation[c], forientation[2]);
                Vec3 fdirection = direction;
                std::swap(fdirection[c], fdirection[2]);
                return std::vector<double>{
                    fdirection[0] - forientation[0] * fdirection[2] / forientation[2],
                        fdirection[1] - forientation[1] * fdirection[2] / forientation[2]
                };
            }
            else{
                assert(rp.variables.size() == 3);
                // depth = 1.0 / (ax + by + cz) where (x, y, z) = direction, (a, b, c) = variables
                // -> 1.0/depth = ax + by + cz
                return std::vector<double>{direction[0], direction[1], direction[2]};
            }
        }




        namespace {

            inline double DepthAt(const Vec3 & direction, const Plane3 & plane){
                InfiniteLine3 ray(Point3(0, 0, 0), direction);
                return norm(IntersectionOfLineAndPlane(ray, plane).position);
            }

            inline double DepthAt(const Vec3 & direction, const Line3 & line){
                InfiniteLine3 ray(Point3(0, 0, 0), direction);
                return norm(DistanceBetweenTwoLines(ray, line.infiniteLine()).second.first);
            }



            std::vector<Vec3> NecessaryAnchorsForBinary(const MixedGraph & mg, RegionBoundaryHandle bh, double & weightForEachAnchor){
                size_t n = ElementsNum(mg.data(bh).normalizedSampledPoints);

                assert(n > 0);
                const auto & points = mg.data(bh).normalizedSampledPoints;

                if (n == 1){
                    weightForEachAnchor = 1.0;
                    return{ points.front().front() };
                }

                auto & p1 = points.front().front();
                auto pp2 = &p1;
                for (auto & ps : points){
                    for (auto & p : ps){
                        if (AngleBetweenDirections(p, p1) > AngleBetweenDirections(p, *pp2)){
                            pp2 = &p;
                        }
                    }
                }
                auto & p2 = *pp2;

                if (n == 2){
                    weightForEachAnchor = 1.0;
                    return{ p1, p2 };
                }
                
                auto normal12 = p1.cross(p2);
                auto pp3 = &p1;
                for (auto & ps : points){
                    for (auto & p : ps){
                        if (abs(p.dot(normal12)) > abs(pp3->dot(normal12))){
                            pp3 = &p;
                        }
                    }
                }
                auto & p3 = *pp3;
                weightForEachAnchor = n / 3.0;
                IMPROVABLE_HERE(?);
                return{ p1, p2, p3 };
            }

            inline std::vector<Vec3> NecessaryAnchorsForBinary(const MixedGraph & mg, LineRelationHandle bh, double & weightForEachAnchor){
                weightForEachAnchor = 5.0;
                return{ mg.data(bh).normalizedRelationCenter };
            }

            inline std::vector<Vec3> NecessaryAnchorsForBinary(const MixedGraph & mg, RegionLineConnectionHandle bh, double & weightForEachAnchor){
                weightForEachAnchor = mg.data(bh).normalizedAnchors.size() / 2.0;
                return{ mg.data(bh).normalizedAnchors.front(), mg.data(bh).normalizedAnchors.back() };
            }


            template <class DataT>
            inline void RegisterConstraintEquations(int & eid, const MixedGraph & mg, MixedGraphPropertyTable & props,
                ComponentHandledTableFromConstraintGraph<int, MixedGraph>::type & uh2varStartPosition,
                ConstraintHandledTableFromConstraintGraph<std::vector<Vec3>, MixedGraph>::type & appliedBinaryAnchors,
                ConstraintHandledTableFromConstraintGraph<double, MixedGraph>::type & weightsForEachAppliedBinaryAnchor,
                std::vector<SparseMatElement<double>> & Atriplets,
                std::vector<SparseMatElement<double>> & Wtriplets,
                std::vector<double> & B) {

                for (auto & c : mg.constraints<DataT>()){
                    if (!props[c.topo.hd].used)
                        continue;
                    auto bh = c.topo.hd;
                    auto uh1 = mg.topo(bh).component<0>();
                    auto uh2 = mg.topo(bh).component<1>();
                    auto & u1 = mg.data(uh1);
                    auto & u2 = mg.data(uh2);

                    int u1VarStartPosition = uh2varStartPosition.at(uh1);
                    int u1VarNum = props[uh1].variables.size();

                    int u2VarStartPosition = uh2varStartPosition.at(uh2);
                    int u2VarNum = props[uh2].variables.size();

                    for (auto & a : appliedBinaryAnchors.at(bh)){
                        B[eid] = 0.0;
                        Wtriplets.emplace_back(eid, eid, weightsForEachAppliedBinaryAnchor.at(bh));
                        {
                            auto u1VarCoeffs = VariableCoefficientsForInverseDepthAtDirection(mg, props, a, uh1);
                            assert(u1VarCoeffs.size() == u1VarNum);
                            for (int i = 0; i < u1VarCoeffs.size(); i++){
                                //A.insert(eid, u1VarStartPosition + i) = u1VarCoeffs[i]; // pos
                                Atriplets.emplace_back(eid, u1VarStartPosition + i, u1VarCoeffs[i]);
                            }
                        }
                        {
                            auto u2VarCoeffs = VariableCoefficientsForInverseDepthAtDirection(mg, props, a, uh2);
                            assert(u2VarCoeffs.size() == u2VarNum);
                            for (int i = 0; i < u2VarCoeffs.size(); i++){
                                //A.insert(eid, u2VarStartPosition + i) = -u2VarCoeffs[i]; // neg
                                Atriplets.emplace_back(eid, u2VarStartPosition + i, -u2VarCoeffs[i]);
                            }
                        }
                        eid++;
                    }
                }
            }


            void FormulateComponentsAndConstraintsAsMatrices(const MixedGraph & mg, MixedGraphPropertyTable & props,
                ComponentHandledTableFromConstraintGraph<int, MixedGraph>::type & uh2varStartPosition,
                ConstraintHandledTableFromConstraintGraph<int, MixedGraph>::type & bh2consStartPosition,
                ConstraintHandledTableFromConstraintGraph<std::vector<Vec3>, MixedGraph>::type & appliedBinaryAnchors,
                ConstraintHandledTableFromConstraintGraph<double, MixedGraph>::type & weightsForEachAppliedBinaryAnchor,
                int & varNum, int & consNum,
                std::vector<SparseMatElement<double>> & Atriplets,
                std::vector<SparseMatElement<double>> & Wtriplets,
                std::vector<double> & X,
                std::vector<double> & B,
                bool addAnchor = true){

                IMPROVABLE_HERE("make sure mg is single connected");

                varNum = 0;
                for (auto & c : mg.components<LineData>()){
                    uh2varStartPosition[c.topo.hd] = varNum;
                    varNum += props[c.topo.hd].variables.size();
                }
                for (auto & c : mg.components<RegionData>()){
                    uh2varStartPosition[c.topo.hd] = varNum;
                    varNum += props[c.topo.hd].variables.size();
                }

                X.resize(varNum);
                for (auto & c : mg.components<LineData>()){
                    for (int i = 0; i < props[c.topo.hd].variables.size(); i++){
                        X[uh2varStartPosition.at(c.topo.hd) + i] = props[c.topo.hd].variables.at(i);
                    }
                }
                for (auto & c : mg.components<RegionData>()){
                    for (int i = 0; i < props[c.topo.hd].variables.size(); i++){
                        X[uh2varStartPosition.at(c.topo.hd) + i] = props[c.topo.hd].variables.at(i);
                    }
                }


                consNum = 0;

                if (addAnchor){
                    consNum++;
                }
                for (auto & c : mg.constraints<RegionBoundaryData>()){
                    if (!props[c.topo.hd].used)
                        continue;
                    bh2consStartPosition[c.topo.hd] = consNum;
                    appliedBinaryAnchors[c.topo.hd] = NecessaryAnchorsForBinary(mg, c.topo.hd, 
                        weightsForEachAppliedBinaryAnchor[c.topo.hd]);
                    consNum += appliedBinaryAnchors[c.topo.hd].size();
                }
                for (auto & c : mg.constraints<LineRelationData>()){
                    if (!props[c.topo.hd].used)
                        continue;
                    bh2consStartPosition[c.topo.hd] = consNum;
                    appliedBinaryAnchors[c.topo.hd] = NecessaryAnchorsForBinary(mg, c.topo.hd,
                        weightsForEachAppliedBinaryAnchor[c.topo.hd]);
                    consNum += appliedBinaryAnchors[c.topo.hd].size();
                }
                for (auto & c : mg.constraints<RegionLineConnectionData>()){
                    if (!props[c.topo.hd].used)
                        continue;
                    bh2consStartPosition[c.topo.hd] = consNum;
                    appliedBinaryAnchors[c.topo.hd] = NecessaryAnchorsForBinary(mg, c.topo.hd,
                        weightsForEachAppliedBinaryAnchor[c.topo.hd]);
                    consNum += appliedBinaryAnchors[c.topo.hd].size();
                }

                B.resize(consNum);

                Atriplets.reserve(consNum * 6);
                Wtriplets.reserve(consNum);

                // write equations
                int eid = 0;
                if (addAnchor){ // the anchor constraint
                    IMPROVABLE_HERE("find a most connected component to set the anchor");

                    // largest plane
                    double maxArea = 0.0;
                    RegionHandle rhAnchored;
                    for (auto & c : mg.components<RegionData>()){
                        if (c.data.area > maxArea){
                            rhAnchored = c.topo.hd;
                            maxArea = c.data.area;
                        }
                    }

                    int uhVarNum = props[rhAnchored].variables.size();
                    Vec3 uhCenter = mg.data(rhAnchored).normalizedCenter;
                    auto uhVarCoeffsAtCenter = VariableCoefficientsForInverseDepthAtDirection(mg, props, uhCenter, rhAnchored);
                    assert(uhVarCoeffsAtCenter.size() == uhVarNum);
                    int uhVarStartPosition = uh2varStartPosition.at(rhAnchored);
                    for (int i = 0; i < uhVarCoeffsAtCenter.size(); i++){
                        //A.insert(eid, uhVarStartPosition + i) = uhVarCoeffsAtCenter[i];
                        Atriplets.emplace_back(eid, uhVarStartPosition + i, uhVarCoeffsAtCenter[i]);
                    }
                    B[eid] = 1.0;
                    //W.insert(eid, eid) = 1.0;
                    Wtriplets.emplace_back(eid, eid, 10.0);
                    eid++;
                }


                RegisterConstraintEquations<RegionBoundaryData>(eid, mg, props, uh2varStartPosition,
                    appliedBinaryAnchors, weightsForEachAppliedBinaryAnchor, Atriplets, Wtriplets, B);
                RegisterConstraintEquations<LineRelationData>(eid, mg, props, uh2varStartPosition,
                    appliedBinaryAnchors, weightsForEachAppliedBinaryAnchor, Atriplets, Wtriplets, B);
                RegisterConstraintEquations<RegionLineConnectionData>(eid, mg, props, uh2varStartPosition,
                    appliedBinaryAnchors, weightsForEachAppliedBinaryAnchor, Atriplets, Wtriplets, B);

                assert(eid == consNum);

            }

        }



        void SolveVariables(const MixedGraph & mg, MixedGraphPropertyTable & props){            

            auto uh2varStartPosition = MakeHandledTableForComponents<int>(mg);
            auto bh2consStartPosition = MakeHandledTableForConstraints<int>(mg);
            auto appliedBinaryAnchors = MakeHandledTableForConstraints<std::vector<Vec3>>(mg);
            auto weightsForEachAppliedBinaryAnchor = MakeHandledTableForConstraints<double>(mg);

            int varNum = 0;
            int consNum = 0;
            std::vector<SparseMatElement<double>> Atriplets;
            std::vector<SparseMatElement<double>> Wtriplets;
            std::vector<double> Xdata;
            std::vector<double> Bdata;

            FormulateComponentsAndConstraintsAsMatrices(mg, props, uh2varStartPosition, bh2consStartPosition, 
                appliedBinaryAnchors, weightsForEachAppliedBinaryAnchor, 
                varNum, consNum, Atriplets, Wtriplets, Xdata, Bdata);
            
            Eigen::SparseMatrix<double> A, W;
            A.resize(consNum, varNum);
            W.resize(consNum, consNum);
            A.reserve(Atriplets.size());
            W.reserve(Wtriplets.size());

            for (auto & t : Atriplets){
                A.insert(t.row, t.col) = t.value;
            }
            for (auto & t : Wtriplets){
                W.insert(t.row, t.col) = t.value;
            }

            Eigen::Map<Eigen::VectorXd> B(Bdata.data(), Bdata.size());

            Eigen::SparseQR<Eigen::SparseMatrix<double>, Eigen::COLAMDOrdering<int>> solver;
            static_assert(!(Eigen::SparseMatrix<double>::IsRowMajor), "COLAMDOrdering only supports column major");
            Eigen::SparseMatrix<double> WA = W * A;
            A.makeCompressed();
            WA.makeCompressed();

            static const bool useWeights = false;
            solver.compute(useWeights ? WA : A);

            if (solver.info() != Eigen::Success) {
                assert(0);
                std::cout << "computation error" << std::endl;
                return;
            }

            Eigen::VectorXd WB = W * B;
            Eigen::VectorXd X = solver.solve(useWeights ? WB : B);
            if (solver.info() != Eigen::Success) {
                assert(0);
                std::cout << "solving error" << std::endl;
                return;
            }


            for (auto & c : mg.components<RegionData>()){
                int uhStartPosition = uh2varStartPosition.at(c.topo.hd);
                for (int i = 0; i < props[c.topo.hd].variables.size(); i++){
                    props[c.topo.hd].variables[i] = X(uhStartPosition + i);
                }
            }
            for (auto & c : mg.components<LineData>()){
                int uhStartPosition = uh2varStartPosition.at(c.topo.hd);
                for (int i = 0; i < props[c.topo.hd].variables.size(); i++){
                    props[c.topo.hd].variables[i] = X(uhStartPosition + i);
                }
            }
        }



        void NormalizeVariables(const MixedGraph & mg, MixedGraphPropertyTable & props){
            for (auto & c : mg.components<RegionData>()){
                for (auto & v : props[c.topo.hd].variables){
                    assert(!IsInfOrNaN(v));
                }
            }
            for (auto & c : mg.components<LineData>()){
                for (auto & v : props[c.topo.hd].variables){
                    assert(!IsInfOrNaN(v));
                }
            }

            double meanCenterDepth = 0.0;
            int componentsNum = 0;
            for (auto & c : mg.components<RegionData>()){
                meanCenterDepth += DepthAt(c.data.normalizedCenter, Instance(mg, props, c.topo.hd));
                assert(!IsInfOrNaN(meanCenterDepth));
                componentsNum++;
            }
            for (auto & c : mg.components<LineData>()){
                meanCenterDepth += DepthAt(normalize(c.data.line.center()), Instance(mg, props, c.topo.hd));
                assert(!IsInfOrNaN(meanCenterDepth));
                componentsNum++;
            }
            meanCenterDepth /= componentsNum;
            assert(!IsInfOrNaN(meanCenterDepth));

            // normalize variables
            for (auto & c : mg.components<RegionData>()){
                for (int i = 0; i < props[c.topo.hd].variables.size(); i++){
                    props[c.topo.hd].variables[i] /= meanCenterDepth;
                }
            }
            for (auto & c : mg.components<LineData>()){
                for (int i = 0; i < props[c.topo.hd].variables.size(); i++){
                    props[c.topo.hd].variables[i] /= meanCenterDepth;
                }
            }
        }


        double ComputeScore(const MixedGraph & mg, const MixedGraphPropertyTable & props){

            // manhattan fitness
            double sumOfComponentWeightedFitness = 0.0;
            double sumOfComponentWeights = 0.0;

            // lines
            double maxLineSpanAngle = 0.0;
            for (auto & l : mg.components<LineData>()){
                double lineSpanAngle = AngleBetweenDirections(l.data.line.first, l.data.line.second);
                if (lineSpanAngle > maxLineSpanAngle)
                    maxLineSpanAngle = lineSpanAngle;
            }

            HandledTable<LineHandle, Line3> lines(mg.internalComponents<LineData>().size());
            HandledTable<RegionHandle, Plane3> planes(mg.internalComponents<RegionData>().size());

            for (auto & l : mg.components<LineData>()){
                lines[l.topo.hd] = Instance(mg, props, l.topo.hd);

                double fitness = 0.0;
                double weight = 0.0;
                double lineSpanAngle = AngleBetweenDirections(l.data.line.first, l.data.line.second);
                auto & prop = props.componentProperties[l.topo.hd];
                if (prop.orientationClaz >= 0){
                    fitness = 1.0;
                    weight = lineSpanAngle / maxLineSpanAngle;
                }
                else{
                    static const double angleThreshold = DegreesToRadians(10);
                    double minAngle = angleThreshold;
                    int bestVPId = -1;
                    for (int i = 0; i < props.vanishingPoints.size(); i++){
                        double angle = AngleBetweenUndirectedVectors(props.vanishingPoints[i], lines[l.topo.hd].direction());
                        if (angle < minAngle){
                            minAngle = angle;
                            bestVPId = i;
                        }
                    }
                    fitness = Gaussian(minAngle / angleThreshold, 0.1);
                    weight = lineSpanAngle / maxLineSpanAngle;
                }
                sumOfComponentWeightedFitness += (fitness * weight);
                sumOfComponentWeights += weight;
            }

            // regions
            double maxRegionArea = 0.0;
            for (auto & r : mg.components<RegionData>()){
                if (r.data.area > maxRegionArea){
                    maxRegionArea = r.data.area;
                }
            }

            for (auto & r : mg.components<RegionData>()){
                planes[r.topo.hd] = Instance(mg, props, r.topo.hd);

                double fitness = 0.0;
                double weight = 1.0;
                auto & prop = props.componentProperties[r.topo.hd];
                if (prop.orientationClaz >= 0){
                    fitness = 1.0;
                    weight = r.data.area / maxRegionArea;
                }
                else{
                    static const double angleThreshold = DegreesToRadians(10);
                    double minAngle = angleThreshold;
                    int bestVPId = -1;
                    for (int i = 0; i < props.vanishingPoints.size(); i++){
                        double angle = AngleBetweenUndirectedVectors(props.vanishingPoints[i], planes[r.topo.hd].normal);
                        if (angle < minAngle){
                            minAngle = angle;
                            bestVPId = i;
                        }
                    }
                    fitness = Gaussian(minAngle / angleThreshold, 0.1);
                    weight = r.data.area / maxRegionArea;
                }
                sumOfComponentWeightedFitness += (fitness * weight);
                sumOfComponentWeights += weight;
            }


            // constraint fitness
            double sumOfConstraintWeightedFitness = 0.0;
            double sumOfConstraintWeights = 0.0;
            double sumOfNotUsedConstraintWeights = 0.0;
            for (auto & c : mg.constraints<RegionBoundaryData>()){
                static const double typeWeight = 1.0;
                if (!props.constraintProperties[c.topo.hd].used){
                    sumOfNotUsedConstraintWeights += ElementsNum(c.data.normalizedSampledPoints) * typeWeight;
                    continue;
                }

                auto & plane1 = planes[c.topo.component<0>()];
                auto & plane2 = planes[c.topo.component<1>()];
                for (auto & ss : c.data.normalizedSampledPoints){
                    for (auto & s : ss){
                        double d1 = DepthAt(s, plane1);
                        double d2 = DepthAt(s, plane2);
                        assert(d1 > 0 && d2 > 0);
                        
                        double fitness = Gaussian(abs(d1 - d2) / (d1 + d2), 0.1);
                        double weight = typeWeight;

                        sumOfConstraintWeightedFitness += (fitness * weight);
                        sumOfConstraintWeights += weight;
                    }
                }
            }
            double maxJunctionWeight = 0.0;
            for (auto & c : mg.constraints<LineRelationData>()){
                if (c.data.junctionWeight > maxJunctionWeight){
                    maxJunctionWeight = c.data.junctionWeight;
                }
            }
            for (auto & c : mg.constraints<LineRelationData>()){
                static const double typeWeight = 8.0;
                if (!props.constraintProperties[c.topo.hd].used){
                    sumOfNotUsedConstraintWeights += c.data.junctionWeight / maxJunctionWeight * typeWeight;
                    continue;
                }

                auto & line1 = lines[c.topo.component<0>()];
                auto & line2 = lines[c.topo.component<1>()];
                double d1 = DepthAt(c.data.normalizedRelationCenter, line1);
                double d2 = DepthAt(c.data.normalizedRelationCenter, line2);

                double fitness = Gaussian(abs(d1 - d2) / (d1 + d2), 0.1);
                double weight = typeWeight;

                sumOfConstraintWeightedFitness += (fitness * weight);
                sumOfConstraintWeights += weight;
            }
            for (auto & c : mg.constraints<RegionLineConnectionData>()){
                static const double typeWeight = 1.0;
                if (!props.constraintProperties[c.topo.hd].used){
                    sumOfNotUsedConstraintWeights += ElementsNum(c.data.normalizedAnchors) * typeWeight;
                    continue;
                }

                auto & plane = planes[c.topo.component<0>()];
                auto & line = lines[c.topo.component<1>()];
                for (auto & a : c.data.normalizedAnchors){
                    double d1 = DepthAt(a, plane);
                    double d2 = DepthAt(a, line);
                    assert(d1 > 0 && d2 > 0);

                    double fitness = Gaussian(abs(d1 - d2) / (d1 + d2), 0.1);
                    double weight = typeWeight;

                    sumOfConstraintWeightedFitness += (fitness * weight);
                    sumOfConstraintWeights += weight;
                }
            }

            double score = sumOfComponentWeightedFitness / sumOfComponentWeights * 50
                + sumOfConstraintWeightedFitness / sumOfConstraintWeights * 1000
                - sumOfNotUsedConstraintWeights / (sumOfConstraintWeights + sumOfNotUsedConstraintWeights) * 5;

            return score;
        }




        namespace {

            template <class UhColorizerFunT = core::ConstantFunctor<vis::Color>>
            void ManuallyOptimizeMixedGraph(const core::Image & panorama,
                const MixedGraph & mg,
                MixedGraphPropertyTable & props,
                UhColorizerFunT uhColorizer = UhColorizerFunT(vis::ColorTag::White),
                bool optimizeInEachIteration = false) {

                bool modified = true;

                struct ComponentID {
                    int handleID;
                    bool isRegion;
                };

                auto sppCallbackFun = [&props, &mg, &modified](vis::InteractionID iid,
                    const std::pair<ComponentID, vis::Colored<vis::SpatialProjectedPolygon>> & spp) {
                    std::cout << (spp.first.isRegion ? "Region" : "Line") << spp.first.handleID << std::endl;
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
                        SolveVariables(mg, props);
                    }

                    vis::Visualizer viz("mixed graph optimizable");
                    viz.renderOptions.bwColor = 1.0;
                    viz.renderOptions.bwTexColor = 0.0;
                    viz.installingOptions.discretizeOptions.colorTable = vis::ColorTableDescriptor::RGB;
                    std::vector<std::pair<ComponentID, vis::Colored<vis::SpatialProjectedPolygon>>> spps;
                    std::vector<vis::Colored<core::Line3>> lines;

                    for (auto & c : mg.components<RegionData>()){
                        auto uh = c.topo.hd;
                        auto & region = c.data;
                        vis::SpatialProjectedPolygon spp;
                        // filter corners
                        core::ForeachCompatibleWithLastElement(c.data.normalizedContours.front().begin(), c.data.normalizedContours.front().end(),
                            std::back_inserter(spp.corners),
                            [](const core::Vec3 & a, const core::Vec3 & b) -> bool {
                            return core::AngleBetweenDirections(a, b) > M_PI / 100.0;
                        });
                        if (spp.corners.size() < 3)
                            continue;

                        spp.projectionCenter = core::Point3(0, 0, 0);
                        spp.plane = Instance(mg, props, uh);
                        spps.emplace_back(ComponentID{ uh.id, true }, std::move(vis::ColorAs(spp, uhColorizer(uh))));
                    }

                    for (auto & c : mg.components<LineData>()){
                        auto uh = c.topo.hd;
                        auto & line = c.data;
                        lines.push_back(vis::ColorAs(Instance(mg, props, uh), uhColorizer(uh)));
                    }

                    viz.begin(spps, sppCallbackFun).shaderSource(vis::OpenGLShaderSourceDescriptor::XPanorama).resource("texture").end();
                    viz.installingOptions.lineWidth = 4.0;
                    viz.add(lines);

                    std::vector<core::Line3> connectionLines;
                    for (auto & c : mg.constraints<RegionBoundaryData>()){
                        auto & samples = c.data.normalizedSampledPoints;
                        auto inst1 = Instance(mg, props, c.topo.component<0>());
                        auto inst2 = Instance(mg, props, c.topo.component<1>());
                        for (auto & ss : samples){
                            for (auto & s : ss){
                                double d1 = DepthAt(s, inst1);
                                double d2 = DepthAt(s, inst2);
                                connectionLines.emplace_back(normalize(s) * d1, normalize(s) * d2);
                            }
                        }
                    }
                    for (auto & c : mg.constraints<RegionLineConnectionData>()){
                        auto inst1 = Instance(mg, props, c.topo.component<0>());
                        auto inst2 = Instance(mg, props, c.topo.component<1>());
                        for (auto & s : c.data.normalizedAnchors){
                            double d1 = DepthAt(s, inst1);
                            double d2 = DepthAt(s, inst2);
                            connectionLines.emplace_back(normalize(s) * d1, normalize(s) * d2);
                        }
                    }

                    viz.installingOptions.discretizeOptions.color = vis::ColorTag::DarkGray;
                    viz.installingOptions.lineWidth = 2.0;
                    viz.renderOptions.renderMode = vis::RenderModeFlag::Triangles | vis::RenderModeFlag::Lines;
                    //viz.add(connectionLines);
                    viz.camera(core::PerspectiveCamera(800, 800, 500, { 1.0, 1.0, -1.0 }, { 0.0, 0.0, 0.0 }, { 0.0, 0.0, -1.0 }));
                    viz.renderOptions.backgroundColor = vis::ColorTag::White;
                    viz.renderOptions.bwColor = 0.0;
                    viz.renderOptions.bwTexColor = 1.0;
                    viz.camera(core::PerspectiveCamera(600, 600, 300, Point3(5, 5, 5), Point3(0, 0, 0)));
                    viz.show(true, false);

                    vis::ResourceStore::clear();
                }

            }


        }




        void Visualize(const View<PanoramicCamera> & texture, const PanoramicCamera & pcam,
            const MixedGraph & mg, MixedGraphPropertyTable & props){
            ManuallyOptimizeMixedGraph(texture.image, mg, props);
        }
        
        void Visualize(const View<PerspectiveCamera> & texture, const PanoramicCamera & pcam,
            const MixedGraph & mg, MixedGraphPropertyTable & props) {
            ManuallyOptimizeMixedGraph(texture.sampled(PanoramicCamera(250, Point3(0, 0, 0), Point3(1, 0, 0), Vec3(0, 0, -1))).image, mg, props);
        }

    }
}


#if 0
namespace panoramix {
    namespace core {

        namespace {
            int SwappedComponent(const MGUnary & u){
                for (int i = 0; i < 2; i++){
                    if (abs(u.normalizedOrientation[i]) >= 1e-8){
                        return i;
                    }
                }
                return 2;
            }
        }

        MGUnaryVariable::MGUnaryVariable(const MGUnary & u, bool f) : fixed(f) {
            if (u.type == MGUnary::RegionFree){
                // (a, b, c) for RegionFree ax+by+c=1, 
                variables = { u.normalizedCenter[0], u.normalizedCenter[1], u.normalizedCenter[2] };
            }
            else if (u.type == MGUnary::RegionAlongFixedAxis){
                // (a, b), {or (b, c) or (a, c)} for RegionAlongFixedAxis  ax+by+c=1,
                Vec3 anotherAxis = u.normalizedCenter.cross(u.normalizedOrientation);
                Vec3 trueNormal = u.normalizedOrientation.cross(anotherAxis);
                Plane3 plane(u.normalizedCenter, trueNormal);
                auto eq = Plane3ToEquation(plane);
                int c = SwappedComponent(u);
                std::swap(eq[c], eq[2]);
                variables = { eq[0], eq[1] };
            }
            else if (u.type == MGUnary::RegionWithFixedNormal){
                // 1/centerDepth for RegionWithFixedNormal
                variables = { 1.0 };
            }
            else if (u.type == MGUnary::LineFree){
                // (1/cornerDepth1, 1/cornerDepth2) for LineFree
                variables = { 1.0, 1.0 };
            }
            else if (u.type == MGUnary::LineOriented){
                // 1/centerDepth for LineOriented,
                variables = { 1.0 };
            }
        }


        double MGUnaryVariable::rawDepth() const{
            double sqSum = 0.0;
            for (auto v : variables){
                sqSum += v * v;
            }
            return 1.0 / sqrt(sqSum);
        }

        Plane3 MGUnaryVariable::interpretAsPlane(const MGUnary & region) const {
            assert(region.isRegion());
            if (region.type == MGUnary::RegionWithFixedNormal){
                assert(variables.size() == 1);
                return Plane3(region.normalizedCenter / variables[0], region.normalizedOrientation);
            }
            else if (region.type == MGUnary::RegionAlongFixedAxis){
                assert(variables.size() == 2);
                double vs[] = { variables[0], variables[1], 0.0 }; // fake vs
                // v1 * o1 + v2 * o2 + v3 * o3 = 0, since region.orientation is orthogonal to normal
                int c = SwappedComponent(region);
                Vec3 orientation = region.normalizedOrientation;
                std::swap(orientation[c], orientation[2]); // now fake orientation
                vs[2] = (-vs[0] * orientation[0] - vs[1] * orientation[1])
                    / orientation[2];
                std::swap(vs[c], vs[2]); // now real vs
                return Plane3FromEquation(vs[0], vs[1], vs[2]);
            }
            else /*if (region.type == MGUnary::RegionWithFixedNormal)*/{
                assert(variables.size() == 3);
                return Plane3FromEquation(variables[0], variables[1], variables[2]);
            }
        }

        Line3 MGUnaryVariable::interpretAsLine(const MGUnary & line) const {
            assert(line.isLine());
            if (line.type == MGUnary::LineOriented){
                assert(variables.size() == 1);
                InfiniteLine3 infLine(line.normalizedCenter / variables[0], line.normalizedOrientation);
                return Line3(DistanceBetweenTwoLines(InfiniteLine3(Point3(0, 0, 0), line.normalizedCorners.front()), infLine).second.second,
                    DistanceBetweenTwoLines(InfiniteLine3(Point3(0, 0, 0), line.normalizedCorners.back()), infLine).second.second);
            }
            else /*if (line.type == MGUnary::LineFree)*/{
                assert(variables.size() == 2);
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
            const MGUnary & u) const{
            if (u.type == MGUnary::RegionWithFixedNormal){
                assert(variables.size() == 1);
                Plane3 plane(u.normalizedCenter, u.normalizedOrientation);
                // variable is 1.0/centerDepth
                // corresponding coeff is 1.0/depthRatio
                // so that 1.0/depth = 1.0/centerDepth * 1.0/depthRatio -> depth = centerDepth * depthRatio
                double depthRatio = norm(IntersectionOfLineAndPlane(InfiniteLine3(Point3(0, 0, 0), direction), plane).position);
                return std::vector<double>{1.0 / depthRatio};
            }
            else if (u.type == MGUnary::RegionAlongFixedAxis){
                assert(variables.size() == 2);
                // depth = 1.0 / (ax + by + cz) where (x, y, z) = direction, (a, b, c) = variables
                // -> 1.0/depth = ax + by + cz
                // -> 1.0/depth = ax + by + (- o1a - o2b)/o3 * z = a(x-o1*z/o3) + b(y-o2*z/o3)
                int c = SwappedComponent(u);
                Vec3 forientation = u.normalizedOrientation;
                std::swap(forientation[c], forientation[2]);
                Vec3 fdirection = direction;
                std::swap(fdirection[c], fdirection[2]);
                return std::vector<double>{
                    fdirection[0] - forientation[0] * fdirection[2] / forientation[2],
                        fdirection[1] - forientation[1] * fdirection[2] / forientation[2]
                };
            }
            else if (u.type == MGUnary::RegionFree){
                assert(variables.size() == 3);
                // depth = 1.0 / (ax + by + cz) where (x, y, z) = direction, (a, b, c) = variables
                // -> 1.0/depth = ax + by + cz
                return std::vector<double>{direction[0], direction[1], direction[2]};
            }
            else if (u.type == MGUnary::LineOriented){
                assert(variables.size() == 1);
                InfiniteLine3 infLine(u.normalizedCenter, u.normalizedOrientation);
                // variable is 1.0/centerDepth
                // corresponding coeff is 1.0/depthRatio
                // so that 1.0/depth = 1.0/centerDepth * 1.0/depthRatio -> depth = centerDepth * depthRatio
                double depthRatio = norm(DistanceBetweenTwoLines(InfiniteLine3(Point3(0, 0, 0), direction), infLine).second.first);
                return std::vector<double>{1.0 / depthRatio};
            }
            else /*if(u.type == MGUnary::LineFree)*/{
                assert(variables.size() == 2);
                const auto & line = u;
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

        double MGUnaryVariable::inverseDepthAtDirection(const Vec3 & direction, const MGUnary & u) const {
            if (u.type == MGUnary::RegionWithFixedNormal){
                assert(variables.size() == 1);
                Plane3 plane(u.normalizedCenter, u.normalizedOrientation);
                // variable is 1.0/centerDepth
                // corresponding coeff is 1.0/depthRatio
                // so that 1.0/depth = 1.0/centerDepth * 1.0/depthRatio -> depth = centerDepth * depthRatio
                double depthRatio = norm(IntersectionOfLineAndPlane(InfiniteLine3(Point3(0, 0, 0), direction), plane).position);
                return variables[0] / depthRatio;
            }
            else if (u.type == MGUnary::RegionAlongFixedAxis){
                assert(variables.size() == 2);
                // depth = 1.0 / (ax + by + cz) where (x, y, z) = direction, (a, b, c) = variables
                // -> 1.0/depth = ax + by + cz
                // -> 1.0/depth = ax + by + (- o1a - o2b)/o3 * z = a(x-o1*z/o3) + b(y-o2*z/o3)
                int c = SwappedComponent(u);
                Vec3 forientation = u.normalizedOrientation;
                std::swap(forientation[c], forientation[2]);
                Vec3 fdirection = direction;
                std::swap(fdirection[c], fdirection[2]);
                double xx = fdirection[0] - forientation[0] * fdirection[2] / forientation[2];
                double yy = fdirection[1] - forientation[1] * fdirection[2] / forientation[2];
                return variables[0] * xx + variables[1] * yy;
            }
            else if (u.type == MGUnary::RegionFree){
                assert(variables.size() == 3);
                // depth = 1.0 / (ax + by + cz) where (x, y, z) = direction, (a, b, c) = variables
                // -> 1.0/depth = ax + by + cz
                return variables[0] * direction[0] + variables[1] * direction[1] + variables[2] * direction[2];
            }
            else if (u.type == MGUnary::LineOriented){
                assert(variables.size() == 1);
                const auto & line = u;
                InfiniteLine3 infLine(line.normalizedCenter, u.normalizedOrientation);
                // variable is 1.0/centerDepth
                // corresponding coeff is 1.0/depthRatio
                // so that 1.0/depth = 1.0/centerDepth * 1.0/depthRatio -> depth = centerDepth * depthRatio
                double depthRatio = norm(DistanceBetweenTwoLines(InfiniteLine3(Point3(0, 0, 0), direction), infLine).second.first);
                return variables[0] / depthRatio;
            }
            else  /*if(u.type == MGUnary::LineFree)*/{
                assert(variables.size() == 2);
                const auto & line = u;
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

        double MGUnaryVariable::depthAtCenter(const MGUnary & unary) const {
            return 1.0 / inverseDepthAtDirection(unary.normalizedCenter, unary);
        }

        void MGUnaryVariable::fitToClosestOrientation(MGUnary & u, const std::vector<Vec3> & vps, double angleThreshold){
            if (u.isRegion()){
                if (!u.hasOrientationConstraints()){
                    Plane3 plane = interpretAsPlane(u);
                    int claz = -1;
                    double minAngle = angleThreshold;
                    for (int i = 0; i < vps.size(); i++){
                        double angle = AngleBetweenUndirectedVectors(plane.normal, vps[i]);
                        assert(angle >= 0);
                        if (angle < minAngle){
                            claz = i;
                            minAngle = angle;
                        }
                    }
                    if (claz >= 0){
                        u.type = MGUnary::RegionWithFixedNormal;
                        u.normalizedOrientation = normalize(vps[claz]);
                        variables.resize(1);
                    }
                }
            }
            else /*if (unary.type == MGUnary::Line)*/{
                if (!u.hasOrientationConstraints()){
                    assert(variables.size() == 2);
                    Line3 line(u.normalizedCorners.front() / variables[0], u.normalizedCorners.back() / variables[1]);
                    Vec3 normal = normalize(line.first.cross(line.second));
                    Vec3 direction = line.direction();
                    int claz = -1;
                    double minAngle = angleThreshold;
                    for (int i = 0; i < vps.size(); i++){
                        if (abs(normal.dot(normalize(vps[i]))) > 0.2) // projection does not match!!!
                            continue;
                        double angle = AngleBetweenUndirectedVectors(direction, vps[i]);
                        assert(angle >= 0);
                        if (angle < minAngle){
                            claz = i;
                            minAngle = angle;
                        }
                    }
                    if (claz >= 0){
                        u.type = MGUnary::LineOriented;
                        u.normalizedOrientation = vps[claz];
                        variables.resize(1);
                    }
                }
            }
        }





        bool MGUnary::hasOrientationConstraints() const { return type == RegionAlongFixedAxis || type == RegionWithFixedNormal || type == LineOriented; }
        bool MGUnary::isRegion() const { return type == RegionFree || type == RegionAlongFixedAxis || type == RegionWithFixedNormal; }
        bool MGUnary::isLine() const { return type == LineFree || type == LineOriented; }





        template <class CameraT>
        MixedGraph BuildMixedGraphTemplate(const std::vector<View<PerspectiveCamera>> & views,
            const std::vector<RegionsGraph> & regionsGraphs, const std::vector<LinesGraph> & linesGraphs,
            const ComponentIndexHashMap<std::pair<RegionIndex, RegionIndex>, double> & regionOverlappingsAcrossViews,
            const ComponentIndexHashMap<std::pair<LineIndex, LineIndex>, Vec3> & lineIncidencesAcrossViews,
            const std::vector<std::map<std::pair<RegionHandle, LineHandle>, std::vector<Point2>>> & regionLineConnections,
            const std::vector<Vec3> & vps,
            const SurfaceLabels<CameraT> & surfaceLabels,
            double initialDepth){

            MixedGraph mg;
            ComponentIndexHashMap<RegionIndex, MGUnaryHandle> ri2mgh;
            ComponentIndexHashMap<LineIndex, MGUnaryHandle> li2mgh;
            std::unordered_map<MGUnaryHandle, RegionIndex> mgh2ri;
            std::unordered_map<MGUnaryHandle, LineIndex> mgh2li;

            // get the vertical vp
            int vertVPId = std::distance(vps.begin(), std::min_element(vps.begin(), vps.end(), [](const Vec3 & a, const Vec3 & b){
                return AngleBetweenUndirectedVectors(a, Vec3(0, 0, 1)) < AngleBetweenUndirectedVectors(b, Vec3(0, 0, 1));
            }));

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
                    SurfaceLabelNames slabel = surfaceLabels.mostLikelyLabelAt(normalizedContour);
                    if (slabel == SurfaceLabelNames::None) {
                        continue;
                    }

                    // initialization
                    MGUnary u;
                    u.type = MGUnary::Undefined;
                    u.normalizedCorners = std::move(normalizedContour);
                    u.normalizedCenter = normalize(cam.spatialDirection(rd.data.center));
                    u.normalizedOrientation = u.normalizedCenter;
                    u.material = MGUnary::Unknown;

                    // set members
                    if (slabel == SurfaceLabelNames::Horizontal){
                        u.type = MGUnary::RegionWithFixedNormal;
                        u.normalizedOrientation = normalize(vps.at(vertVPId));
                        u.material = MGUnary::Planar;
                    }
                    else if (slabel == SurfaceLabelNames::Vertical){
                        u.type = MGUnary::RegionAlongFixedAxis;
                        u.normalizedOrientation = normalize(vps.at(vertVPId));
                        u.material = MGUnary::Planar;
                    }
                    else if (slabel == SurfaceLabelNames::Planar){
                        u.type = MGUnary::RegionFree;
                        u.material = MGUnary::Planar;
                    }
                    else if (slabel == SurfaceLabelNames::NonPlanar){
                        u.type = MGUnary::RegionAlongFixedAxis;
                        u.normalizedOrientation = normalize(vps.at(vertVPId));
                        u.material = MGUnary::NonPlanar;
                    }
                    else if (slabel == SurfaceLabelNames::Void){
                        u.type = MGUnary::RegionWithFixedNormal;
                        u.normalizedOrientation = normalize(vps.at(vertVPId));
                        u.material = MGUnary::Void;
                    }
                    else{
                        SHOULD_NEVER_BE_CALLED();
                    }

                    ri2mgh[ri] = mg.add(std::move(u));
                    mgh2ri[ri2mgh[ri]] = ri;                    
                }
                // lines
                for (auto & ld : linesGraphs[i].elements<0>()){
                    auto li = LineIndex{ i, ld.topo.hd };

                    MGUnary u;
                    u.type = MGUnary::Undefined;
                    u.normalizedCorners = {
                        normalize(cam.spatialDirection(ld.data.line.component.first)),
                        normalize(cam.spatialDirection(ld.data.line.component.second))
                    };
                    u.normalizedCenter = normalize(cam.spatialDirection(ld.data.line.component.center()));
                    u.normalizedOrientation = u.normalizedCenter;
                    u.material = MGUnary::Unknown;

                    SurfaceLabelNames slabel = surfaceLabels.mostLikelyLabelAt(normalizedContour);
                    if (slabel == SurfaceLabelNames::None) {
                        continue;
                    }

                    // set members
                    if (slabel == SurfaceLabelNames::Horizontal){
                        u.type = MGUnary::RegionWithFixedNormal;
                        u.normalizedOrientation = normalize(vps.at(vertVPId));
                        u.material = MGUnary::Planar;
                    }
                    else if (slabel == SurfaceLabelNames::Vertical){
                        u.type = MGUnary::RegionAlongFixedAxis;
                        u.normalizedOrientation = normalize(vps.at(vertVPId));
                        u.material = MGUnary::Planar;
                    }
                    else if (slabel == SurfaceLabelNames::Planar){
                        u.type = MGUnary::RegionFree;
                        u.material = MGUnary::Planar;
                    }
                    else if (slabel == SurfaceLabelNames::NonPlanar){
                        u.type = MGUnary::RegionAlongFixedAxis;
                        u.normalizedOrientation = normalize(vps.at(vertVPId));
                        u.material = MGUnary::NonPlanar;
                    }
                    else if (slabel == SurfaceLabelNames::Void){
                        u.type = MGUnary::RegionWithFixedNormal;
                        u.normalizedOrientation = normalize(vps.at(vertVPId));
                        u.material = MGUnary::Void;
                    }
                    else{
                        SHOULD_NEVER_BE_CALLED();
                    }


                    li2mgh[li] = mg.add(MGUnary{
                        MGUnary::Line,
                        ,
                        ,
                        ld.data.line.claz
                    });
                    mgh2li[li2mgh[li]] = li;
                    /*unaryVars[li2mgh[li]].variables = ld.data.line.claz == -1 ?
                        std::vector<double>{ 1.0, 1.0 } : std::vector<double>{ 1.0 };
                    unaryVars[li2mgh[li]].fixed = false;
                    unaryVars[li2mgh[li]].orientationClaz = mg.data(li2mgh.at(li)).primaryOrientationClaz;*/
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
            /*for (auto & b : mg.internalElements<1>()){
                binaryVars[b.topo.hd].enabled = true;
            }*/

            return mg;
        }







        MixedGraph BuildMixedGraph(const std::vector<View<PerspectiveCamera>> & views,
            std::vector<Vec3> & vps,

            double initialDepth,
            const core::LineSegmentExtractor & lineseger,
            double intersectionDistanceThreshold,
            double incidenceDistanceAlongDirectionThreshold,
            double incidenceDistanceVerticalDirectionThreshold,

            const core::SegmentationExtractor & segmenter,
            double samplingStepLengthOnBoundary,
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
                initialDepth);

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


        MGPatch MakePatchOnBinary(const MixedGraph & mg, const MGBinaryHandle & bh, const std::vector<Vec3> & vps){
            MGPatch patch;
            patch.bhs[bh] = MGBinaryVariable(mg.data(bh));
            auto uhs = mg.topo(bh).lowers;
            patch.uhs[uhs[0]] = MGUnaryVariable(mg.data(uhs[0]), vps);
            patch.uhs[uhs[1]] = MGUnaryVariable(mg.data(uhs[1]), vps);
            assert(BinaryHandlesAreValidInPatch(mg, patch));
            assert(UnariesAreConnectedInPatch(mg, patch));
            return patch;
        }


        MGPatch MakeStarPatchAroundUnary(const MixedGraph & mg, const MGUnaryHandle & uh, const std::vector<Vec3> & vps){
            MGPatch patch;
            patch.uhs[uh] = MGUnaryVariable(mg.data(uh), vps);
            for (auto & bh : mg.topo(uh).uppers){
                auto anotherUh = mg.topo(bh).lowers.front();
                if (anotherUh == uh)
                    anotherUh = mg.topo(bh).lowers.back();
                patch.bhs[bh] = MGBinaryVariable(mg.data(bh));
                patch.uhs[anotherUh] = MGUnaryVariable(mg.data(uh), vps);
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



        std::vector<MGPatch> SplitMixedGraphIntoPatches(const MixedGraph & mg, const std::vector<Vec3> & vps){

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
                patches[ccid].uhs[uh] = MGUnaryVariable(mg.data(uh), vps);
            }
            for (auto & b : mg.elements<1>()){
                auto & uhs = mg.topo(b.topo.hd).lowers;
                assert(ccids[uhs.front()] == ccids[uhs.back()]);
                patches[ccids[uhs.front()]].bhs[b.topo.hd] = MGBinaryVariable(b.data);
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



            template <class SparseMatElementT, class T>
            void FormulateComponentsAndConstraintsAsMatrices(const MixedGraph & mg, MGPatch & patch,
                const std::vector<Vec3> & vanishingPoints,
                std::unordered_map<MGUnaryHandle, int> & uh2varStartPosition,
                std::unordered_map<MGBinaryHandle, int> & bh2consStartPosition,
                std::unordered_map<MGBinaryHandle, std::vector<Vec3>> & appliedBinaryAnchors,
                int & varNum, int & consNum,
                std::vector<SparseMatElementT> & Atriplets,
                std::vector<SparseMatElementT> & Wtriplets,
                std::vector<T> & X,
                std::vector<T> & B,
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

                X.resize(varNum);
                for (auto & uhv : patch.uhs){
                    if (uhv.second.fixed){
                        continue;
                    }
                    for (int i = 0; i < uhv.second.variables.size(); i++){
                        X[uh2varStartPosition.at(uhv.first) + i] = uhv.second.variables.at(i);
                    }
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
                        B[eid] = 1.0;
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

                        B[eid] = 0.0;
                        assert(mg.data(bh).weight >= 0.0);
                        //W.insert(eid, eid) = mg.data(bh).weight;
                        Wtriplets.emplace_back(eid, eid, mg.data(bh).weight);

                        if (u1IsFixed){
                            double inverseDepthAtA = u1Var.inverseDepthAtDirection(a, u1, vanishingPoints);
                            B[eid] = -inverseDepthAtA;
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
                            B[eid] = inverseDepthAtA;
                        }
                        else{
                            auto u2VarCoeffs = u2Var.variableCoeffsForInverseDepthAtDirection(a, u2, vanishingPoints);
                            assert(u2VarCoeffs.size() == u2VarNum);
                            for (int i = 0; i < u2VarCoeffs.size(); i++){
                                //A.insert(eid, u2VarStartPosition + i) = -u2VarCoeffs[i]; // neg
                                Atriplets.emplace_back(eid, u2VarStartPosition + i, -u2VarCoeffs[i]);
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
                std::vector<double> X, B;
                bool useWeights;

                std::unordered_map<MGUnaryHandle, int> uh2varStartPosition;
                std::unordered_map<MGBinaryHandle, int> bh2consStartPosition;

                std::unordered_map<MGBinaryHandle, std::vector<Vec3>> appliedBinaryAnchors;

                virtual void initialize(const MixedGraph & mg, MGPatch & patch,
                    const std::vector<Vec3> & vanishingPoints, bool useWeights) override {

                    int varNum, consNum;
                    std::vector<SparseMatElement<double>> Atriplets, Wtriplets;
                    FormulateComponentsAndConstraintsAsMatrices(mg, patch, vanishingPoints,
                        uh2varStartPosition, bh2consStartPosition, appliedBinaryAnchors, varNum, consNum, Atriplets, Wtriplets, X, B, true);
                    A = MakeSparseMatFromElements(consNum, varNum, Atriplets.begin(), Atriplets.end());
                    W = MakeSparseMatFromElements(consNum, consNum, Wtriplets.begin(), Wtriplets.end());
                }

                virtual bool optimize(const MixedGraph & mg, MGPatch & patch,
                    const std::vector<Vec3> & vanishingPoints)  override {

                    core::Matlab::RunScript("clear");
                    
                    core::Matlab::PutVariable("A", A);
                    core::Matlab::PutVariable("W", W);
                    core::Matlab::PutVariable("B", B);
                    core::Matlab::PutVariable("X0", X);
                    core::Matlab matlab;

                    matlab << "B = B(:); X0 = X0(:);";

                    matlab
                        << "cvx_setup;"
                        << "cvx_begin"
                        << "   variable X1(size(A, 2))"
                        << "   minimize(norm((A * X1 - B)))"
                        << "cvx_end";

                    //matlab
                    //    << "m = median(X1);"
                    //    << "outs = (X1 < 0.1 * m) | (X1 > 10 * m);";
                    //
                    //matlab
                    //    << "A2 = A(:, outs);"
                    //    << "B2 = B - A(:, ~outs) * X1(~outs);";

                    //matlab
                    //    << "cvx_setup;"
                    //    << "cvx_begin"
                    //    << "   variable X2(size(A2, 2))"
                    //    << "   minimize(norm((A2 * X2 - B2)) * 100 + norm(X2))"
                    //    //<< "     subject to"
                    //    //<< "        abs(X2) <= abs(X1(outs(:))) / 2.0"
                    //    << "cvx_end";

                    /*matlab
                        << "X1(find(outs)) = X2;";*/

                    matlab << "X1 = X1';";

                    //matlab << "save temp;";

                    
                    std::vector<double> X;
                    core::Matlab::GetVariable("X1", X, false);

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

            if (_at == Algorithm::EigenSparseQR){
                _internal = new MGPatchDepthsOptimizerInternalEigen;
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


        void FitUnariesToClosestOrientations(const MixedGraph & mg, MGPatch & patch, const std::vector<Vec3> & vps, double angleThreshold){
            for (auto & uhv : patch.uhs){
                uhv.second.fitToClosestOrientation(mg.data(uhv.first), vps, angleThreshold);
            }
        }


   	}
}

#endif