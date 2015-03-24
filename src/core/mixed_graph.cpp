
extern "C" {
    #include <gpc.h>
//    #include <mosek.h>
}
#include <Eigen/Dense>
#include <Eigen/Sparse>

#include <Eigen/SPQRSupport>

#include <unsupported/Eigen/NonLinearOptimization>
#include <unsupported/Eigen/NumericalDiff>

//
#include <GCoptimization.h>


#include "../misc/matlab.hpp"

//
#include "algorithms.hpp"
#include "containers.hpp"
#include "utilities.hpp"
#include "clock.hpp"
#include "mixed_graph.hpp"

#include "../gui/visualizers.hpp"
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
            for (auto & c : mg.constraints<LineRelationData>()){
                fun(c.topo.hd);
            }
            for (auto & c : mg.constraints<RegionLineConnectionData>()){
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

            std::pair<double, double> ComputeSpanningArea(const Point2 & a, const Point2 & b, const Ray2 & line) {
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

            auto vanishingPoints = FindOrthogonalPrinicipleDirections(lineIntersections, 1000, 500, true).unwrap();

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
            ClassifyLines(spatialLineSegments, vanishingPoints, M_PI / 3.0, 0.1, 0.8, M_PI / 18.0);

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
                double incidenceDistanceVerticalDirectionThreshold,
                bool includeIncidencesBetweenUnclassifiedLines = false,
                bool computeJunctionWeights = false){

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
                        else if (includeIncidencesBetweenUnclassifiedLines && clazi == clazj && clazi == -1){ // incidence for unclassified lines
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

                if (computeJunctionWeights) {
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
                }
                else {
                    for (auto & lr : graph.elements<1>()){
                        lr.data.junctionWeight = 1.0;
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


            //template <class T, int N>
            //std::vector<Point<T, N>> NormalizedSamplePointsOnLine(const Line<T, N> & line, const T & stepLen){
            //    T len = line.length();
            //    auto dir = normalize(line.direction());
            //    std::vector<Point<T, N>> pts;
            //    pts.reserve(len / stepLen);
            //    for (T l = 0; l <= len; l += stepLen){
            //        pts.push_back(normalize(line.first + dir * l));
            //    }
            //    return pts;
            //}

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
                        //std::cout << "multiple contours for one region in perspective projection!";
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
                std::map<std::pair<int, int>, std::vector<std::vector<PixelLoc>>> & boundaryEdges,
                int connectionExtendSize) {

                std::map<std::pair<int, int>, std::set<PixelLoc, ComparePixelLoc>> boundaryPixels;

                int width = segRegions.cols;
                int height = segRegions.rows;
                for (int y = 0; y < height - 1; y++) {
                    for (int x = 0; x < width - 1; x++) {
                        int originalRegionId = segRegions(PixelLoc(x, y));
                        for (int xx = std::max(x - connectionExtendSize, 0); xx <= std::min(x + connectionExtendSize, width - 1); xx++){
                            for (int yy = std::max(y - connectionExtendSize, 0); yy <= std::min(y + connectionExtendSize, height - 1); yy++){
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

                    IMPROVABLE_HERE("what if connectionExtendSize is too large? will it cause bugs here searching edges?");

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
                double samplingStepAngleOnBoundary, double samplingStepAngleOnLine, 
                int samplerSizeOnBoundary, int samplerSizeOnLine){

                auto regionHandles = CollectRegionsFromSegmentation(mg, segmentedRegions, cam);
                int regionNum = regionHandles.size();

                // add region boundary constraints
                std::map<std::pair<int, int>, std::vector<std::vector<PixelLoc>>> boundaryEdges;
                FindContoursOfRegionsAndBoundaries(segmentedRegions, regionNum, boundaryEdges, samplerSizeOnBoundary);

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
                    double spanAngle = AngleBetweenDirections(line.first, line.second);
                    int stepNum = static_cast<int>(std::ceil(spanAngle / samplingStepAngleOnLine));
                    assert(stepNum >= 1);
                    for (int step = 0; step <= stepNum; step ++){
                        double angle = step * samplingStepAngleOnLine;
                        Vec3 sample = RotateDirection(line.first, line.second, angle);
                        if (!cam.isVisibleOnScreen(sample))
                            continue;
                        PixelLoc originalP = ToPixelLoc(cam.screenProjection(sample));
                        // collect neighbors
                        std::set<int> connectedRegionIds;
                        for (int x = originalP.x - samplerSizeOnLine; x <= originalP.x + samplerSizeOnLine; x++){
                            for (int y = originalP.y - samplerSizeOnLine; y <= originalP.y + samplerSizeOnLine; y++){
                                PixelLoc p(x, y);
                                if (p.x < 0 || p.x >= segmentedRegions.cols || p.y < 0 || p.y >= segmentedRegions.rows)
                                    continue;
                                int regionId = segmentedRegions(p);
                                auto rh = regionHandles[regionId];
                                if (rh.invalid())
                                    continue;
                                connectedRegionIds.insert(regionId);
                            }
                        }                      
                        for (int regionId : connectedRegionIds){
                            regionLineConnections[std::make_pair(regionHandles[regionId], ld.topo.hd)].push_back(normalize(sample));
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
            double samplingStepAngleOnBoundary, double samplingStepAngleOnLine, 
            int samplerSizeOnBoundary, int samplerSizeOnLine){
            AppendRegionsTemplate(mg, segmentedRegions, cam, samplingStepAngleOnBoundary, samplingStepAngleOnLine, 
                samplerSizeOnBoundary, samplerSizeOnLine);
        }
        void AppendRegions(MixedGraph & mg, const Imagei & segmentedRegions, const PanoramicCamera & cam,
            double samplingStepAngleOnBoundary, double samplingStepAngleOnLine,
            int samplerSizeOnBoundary, int samplerSizeOnLine){
            AppendRegionsTemplate(mg, segmentedRegions, cam, samplingStepAngleOnBoundary, samplingStepAngleOnLine, 
                samplerSizeOnBoundary, samplerSizeOnLine);
        }



        namespace {

            template <class FunctorT>
            struct DataCostFunctorWrapper : GCoptimization::DataCostFunctor{
                inline DataCostFunctorWrapper(FunctorT && f) : fun(std::forward<FunctorT>(f)) {}
                virtual GCoptimization::EnergyTermType compute(GCoptimization::SiteID s, GCoptimization::LabelID l) override {
                    return fun(s, l);
                }
                FunctorT fun;
            };
            template <class FunctorT>
            inline DataCostFunctorWrapper<FunctorT> * AllocDataCostFunctor(FunctorT && f) {
                return new DataCostFunctorWrapper<FunctorT>(std::forward<FunctorT>(f));
            }

            template <class FunctorT>
            struct SmoothCostFunctorWrapper : GCoptimization::SmoothCostFunctor {
                inline SmoothCostFunctorWrapper(FunctorT && f) : fun(std::forward<FunctorT>(f)){}
                virtual GCoptimization::EnergyTermType compute(
                    GCoptimization::SiteID s1, GCoptimization::SiteID s2,
                    GCoptimization::LabelID l1, GCoptimization::LabelID l2) override {
                    return fun(s1, s2, l1, l2);
                }
                FunctorT fun;
            };
            template <class FunctorT>
            inline SmoothCostFunctorWrapper<FunctorT> * AllocSmoothCostFunctor(FunctorT && f) {
                return new SmoothCostFunctorWrapper<FunctorT>(std::forward<FunctorT>(f));
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
                        if (!lp.used){
                            lp.variables = {};
                            return;
                        }
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
                        if (!rp.used){
                            rp.variables = {};
                            return;
                        }
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

        }





        MixedGraphPropertyTable MakeMixedGraphPropertyTable(const MixedGraph & mg, const std::vector<Vec3> & vps) {
            auto compProps = MakeHandledTableForAllComponents<MixedGraphComponentProperty>(mg);
            auto consProps = MakeHandledTableForAllConstraints<MixedGraphConstraintProperty>(mg);
            for (auto & l : mg.internalComponents<LineData>()){
                compProps[l.topo.hd].orientationClaz = l.data.initialClaz;
                compProps[l.topo.hd].orientationNotClaz = -1;
            }
            for (auto & r : mg.internalComponents<RegionData>()){
                compProps[r.topo.hd].orientationClaz = compProps[r.topo.hd].orientationNotClaz = -1;
            }
            for (auto & dtable : compProps.data){
                for (auto & d : dtable){
                    d.used = true;
                }
            }
            for (auto & dtable : consProps.data){
                for (auto & d : dtable){
                    d.used = true;
                }
            }
            MixedGraphPropertyTable props{ vps, std::move(compProps), std::move(consProps) };
            ForeachMixedGraphComponentHandle(mg, InitializeVariablesForEachHandle{ mg, props });
            return props;
        }



        namespace {

            class GPCPolygon {
            public:
                GPCPolygon() : _p(nullptr) {}
                explicit GPCPolygon(gpc_polygon * p) : _p(p) {}
                template <class T>
                GPCPolygon(const std::vector<std::vector<Point<T, 2>>> & pts) {
                    _p = new gpc_polygon;
                    auto & poly = *_p;
                    poly.num_contours = pts.size();
                    poly.contour = new gpc_vertex_list[pts.size()];
                    for (int k = 0; k < pts.size(); k++){
                        poly.contour[k].num_vertices = pts[k].size();
                        poly.contour[k].vertex = new gpc_vertex[pts[k].size()];
                        for (int i = 0; i < pts[k].size(); i++){
                            poly.contour[k].vertex[i].x = pts[k][i][0];
                            poly.contour[k].vertex[i].y = pts[k][i][1];
                        }
                    }
                    poly.hole = new int[pts.size()];
                    std::fill(poly.hole, poly.hole + pts.size(), 0);
                }
                GPCPolygon(const GPCPolygon & p) = delete;
                GPCPolygon(GPCPolygon && p) {
                    _p = p._p;
                    p._p = nullptr;
                }
                GPCPolygon & operator = (const GPCPolygon & p) = delete;
                GPCPolygon & operator = (GPCPolygon && p) {
                    std::swap(_p, p._p);
                    return *this;
                }

                ~GPCPolygon(){
                    if (!_p)
                        return;
                    gpc_free_polygon(_p);
                    delete _p;
                    _p = nullptr;
                }
                template <class T>
                operator std::vector<std::vector<Point<T, 2>>>() const {
                    std::vector<std::vector<Point<T, 2>>> pts;
                    auto & poly = *_p;
                    pts.reserve(poly.num_contours);
                    for (int k = 0; k < poly.num_contours; k++){
                        if (poly.hole[k])
                            continue;

                        std::vector<Point<T, 2>> ps(poly.contour[k].num_vertices);
                        for (int i = 0; i < ps.size(); i++){
                            ps[i][0] = static_cast<T>(poly.contour[k].vertex[i].x);
                            ps[i][1] = static_cast<T>(poly.contour[k].vertex[i].y);
                        }
                        pts.push_back(std::move(ps));
                    }
                    return pts;
                }
            public:
                GPCPolygon clipWith(gpc_op op, const GPCPolygon & clipPoly) const {
                    GPCPolygon result(new gpc_polygon);
                    gpc_polygon_clip(op, _p, clipPoly._p, result._p);
                    return result;
                }

                inline GPCPolygon operator - (const GPCPolygon & clipPoly) const {
                    return clipWith(GPC_DIFF, clipPoly);
                }
                inline GPCPolygon operator & (const GPCPolygon & clipPoly) const {
                    return clipWith(GPC_INT, clipPoly);
                }
                inline GPCPolygon operator ^ (const GPCPolygon & clipPoly) const {
                    return clipWith(GPC_XOR, clipPoly);
                }
                inline GPCPolygon operator | (const GPCPolygon & clipPoly) const {
                    return clipWith(GPC_UNION, clipPoly);
                }

            private:
                gpc_polygon * _p;
            };



            enum class OrientationHint {
                Void,
                Horizontal,
                Vertical,
                OtherPlanar,
                NonPlanar,
                Count
            };

            inline OrientationHint ToOrientationHint(GeometricContextLabel label, bool leftFrontRightAsVertical){
                //THERE_ARE_BUGS_HERE("Only OtherPlanar is returned!");
                //return OrientationHint::OtherPlanar;
                switch (label){
                case GeometricContextLabel::Ceiling: return OrientationHint::Horizontal;
                case GeometricContextLabel::Floor: return OrientationHint::Horizontal;
                case GeometricContextLabel::Front:
                case GeometricContextLabel::Left:
                case GeometricContextLabel::Right:
                    return leftFrontRightAsVertical ? OrientationHint::Vertical : OrientationHint::OtherPlanar;
                case GeometricContextLabel::Furniture: return OrientationHint::Void;
                case GeometricContextLabel::Ground: return OrientationHint::Horizontal;
                case GeometricContextLabel::Sky: return OrientationHint::Void;
                case GeometricContextLabel::Vertical: return OrientationHint::OtherPlanar;
                case GeometricContextLabel::NotPlanar: return OrientationHint::NonPlanar;
                default:
                    return OrientationHint::Void;
                }
            }

            struct ConstraintUsableFlagUpdater {
                const MixedGraph & mg;
                MixedGraphPropertyTable & props;
                template <class T>
                inline void operator()(ConstraintHandle<T> h) const {
                    props[h].used = props[mg.topo(h).component<0>()].used && props[mg.topo(h).component<1>()].used;
                }
            };

            struct DanglingConstraintsDisabler {
                const MixedGraph & mg;
                MixedGraphPropertyTable & props;
                template <class T>
                inline void operator()(ConstraintHandle<T> h) const {
                    auto allEndsUsable = props[mg.topo(h).component<0>()].used && props[mg.topo(h).component<1>()].used;
                    if (!allEndsUsable)
                        props[h].used = false;
                }
            };

        }


        void ResetVariables(const MixedGraph & mg, MixedGraphPropertyTable & props){
            ForeachMixedGraphComponentHandle(mg, InitializeVariablesForEachHandle{ mg, props });
        }

        void UpdateConstraintUsabilities(const MixedGraph & mg, MixedGraphPropertyTable & props, bool enableDisabledConstraints){
            if (!enableDisabledConstraints){
                ForeachMixedGraphConstraintHandle(mg, DanglingConstraintsDisabler{ mg, props });
            }
            else{
                ForeachMixedGraphConstraintHandle(mg, ConstraintUsableFlagUpdater{ mg, props });
            }            
        }


        namespace {

            inline double NonZeroize(double d){
                return d == 0.0 ? 1e-6 : d;
            }

        }

        Line3 Instance(const MixedGraph & mg, const MixedGraphPropertyTable & props, const LineHandle & lh){
            auto & ld = mg.data(lh);
            auto & lp = props[lh];
            //if (!lp.used)
            //    return Line3();
            if (lp.orientationClaz >= 0){
                assert(lp.variables.size() == 1);
                Ray3 infLine(normalize(ld.line.center()) / NonZeroize(lp.variables[0]), props.vanishingPoints[lp.orientationClaz]);
                return Line3(DistanceBetweenTwoLines(Ray3(Point3(0, 0, 0), normalize(ld.line.first)), infLine).second.second,
                    DistanceBetweenTwoLines(Ray3(Point3(0, 0, 0), normalize(ld.line.second)), infLine).second.second);
            }
            else /*if (line.type == MGUnary::LineFree)*/{
                assert(lp.variables.size() == 2);
                /*           | sin(theta) | | p | | q |
                len:---------------------------------   => 1/len: [(1/q)sin(phi) - (1/p)sin(phi-theta)] / sin(theta)
                | p sin(phi) - q sin(phi - theta) |
                */
                // variables[0] -> 1/p
                // variables[1] -> 1/q
                return Line3(normalize(ld.line.first) / NonZeroize(lp.variables[0]), normalize(ld.line.second) / NonZeroize(lp.variables[1]));
            }
        }



        Plane3 Instance(const MixedGraph & mg, const MixedGraphPropertyTable & props, const RegionHandle & rh){
            auto & rd = mg.data(rh);
            auto & rp = props[rh];
            /*if (!rp.used)
                return Plane3();*/
            if (rp.orientationClaz >= 0){
                assert(rp.variables.size() == 1);
                return Plane3(rd.normalizedCenter / NonZeroize(rp.variables[0]), props.vanishingPoints[rp.orientationClaz]);
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
                Ray3 infLine(normalize(ld.line.center()), props.vanishingPoints[lp.orientationClaz]);
                 // variable is 1.0/centerDepth
                 // corresponding coeff is 1.0/depthRatio
                 // so that 1.0/depth = 1.0/centerDepth * 1.0/depthRatio -> depth = centerDepth * depthRatio
                 double depthRatio = norm(DistanceBetweenTwoLines(Ray3(Point3(0, 0, 0), direction), infLine).second.first);
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
                double depthRatio = norm(IntersectionOfLineAndPlane(Ray3(Point3(0, 0, 0), direction), plane).position);
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


        double DepthAtDirectionGivenVariables(const MixedGraph & mg, const double * variables,
            const MixedGraphPropertyTable & props, const Vec3 & direction, const LineHandle & lh){
            auto & ld = mg.data(lh);
            auto & lp = props.componentProperties.at(lh);
            if (lp.orientationClaz >= 0){
                //assert(lp.variables.size() == 1);
                Ray3 infLine(normalize(ld.line.center()), props.vanishingPoints[lp.orientationClaz]);
                // variable is 1.0/centerDepth
                // depths = depthRatio * centerDepth
                double depthRatio = norm(DistanceBetweenTwoLines(Ray3(Point3(0, 0, 0), direction), infLine).second.first);
                double inversedCenterDepth = variables[0];
                return depthRatio / inversedCenterDepth;
            }
            else /*if(u.type == MGUnary::LineFree)*/{
                //assert(lp.variables.size() == 2);
                double theta = AngleBetweenDirections(normalize(ld.line.first), normalize(ld.line.second));
                double phi = AngleBetweenDirections(normalize(ld.line.first), direction);
                /*           | sin(theta) | | p | | q |
                        len:---------------------------------   => 1/len: [(1/q)sin(phi) - (1/p)sin(phi-theta)] / sin(theta)
                            | p sin(phi) - q sin(phi - theta) |
                */
                // variables[0] -> 1/p
                // variables[1] -> 1/q
                return sin(theta) / (variables[1] * sin(phi) - variables[0] * sin(phi - theta));
            }
        }


        double DepthAtDirectionGivenVariables(const MixedGraph & mg, const double * variables,
            const MixedGraphPropertyTable & props, const Vec3 & direction, const RegionHandle & rh){
            auto & rd = mg.data(rh);
            auto & rp = props.componentProperties.at(rh);
            if (rp.orientationClaz >= 0){
                //assert(rp.variables.size() == 1);
                Plane3 plane(rd.normalizedCenter, props.vanishingPoints[rp.orientationClaz]);
                // variable is 1.0/centerDepth
                // depths = depthRatio * centerDepth
                double depthRatio = norm(IntersectionOfLineAndPlane(Ray3(Point3(0, 0, 0), direction), plane).position);
                double inversedCenterDepth = variables[0];
                return depthRatio / inversedCenterDepth;
            }
            else if (rp.orientationClaz == -1 && rp.orientationNotClaz >= 0){
                //assert(rp.variables.size() == 2);
                // depth = 1.0 / (ax + by + cz) where (x, y, z) = direction, (a, b, c) = variables
                // -> 1.0/depth = ax + by + cz
                // -> 1.0/depth = ax + by + (- o1a - o2b)/o3 * z = a(x-o1*z/o3) + b(y-o2*z/o3)
                auto orientation = normalize(props.vanishingPoints[rp.orientationNotClaz]);
                int c = SwappedComponent(orientation);
                Vec3 forientation = orientation;
                std::swap(forientation[c], forientation[2]);
                Vec3 fdirection = direction;
                std::swap(fdirection[c], fdirection[2]);
                double a = variables[0];
                double b = variables[1];

                NOT_TESTED_YET();
                
                return forientation[2] / (
                    a * (fdirection[0] * forientation[2] - forientation[0] * fdirection[2]) +
                    b * (fdirection[1] * forientation[2] - forientation[1] * fdirection[2]));
            }
            else{
                assert(rp.variables.size() == 3);
                // depth = 1.0 / (ax + by + cz) where (x, y, z) = direction, (a, b, c) = variables
                // -> 1.0/depth = ax + by + cz
                return 1.0 / (direction[0] * variables[0] + direction[1] * variables[1] + direction[2] * variables[2]);
            }
        }




        namespace {

            inline double DepthAt(const Vec3 & direction, const Plane3 & plane){
                Ray3 ray(Point3(0, 0, 0), direction);
                return norm(IntersectionOfLineAndPlane(ray, plane).position);
            }

            inline double DepthAt(const Vec3 & direction, const Line3 & line){
                Ray3 ray(Point3(0, 0, 0), direction);
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

                const Vec3 * pp1 = nullptr;
                const Vec3 * pp2 = nullptr;
                double maxAngle = 0.0;
                for (int i = 0; i < points.size(); i++){
                    for (int ii = 0; ii < points[i].size(); ii++){
                        for (int j = i; j < points.size(); j++){
                            for (int jj = 0; jj < points[j].size(); jj++){
                                if (std::tie(i, ii) >= std::tie(j, jj))
                                    continue;
                                double angle = AngleBetweenDirections(points[i][ii], points[j][jj]);
                                if (angle > maxAngle){
                                    maxAngle = angle;
                                    pp1 = &points[i][ii];
                                    pp2 = &points[j][jj];
                                }
                            }
                        }
                    }
                }
                auto & p1 = *pp1;
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
                weightForEachAnchor = mg.data(bh).junctionWeight * 2;
                return{ mg.data(bh).normalizedRelationCenter };
            }

            inline std::vector<Vec3> NecessaryAnchorsForBinary(const MixedGraph & mg, RegionLineConnectionHandle bh, double & weightForEachAnchor){
                weightForEachAnchor = mg.data(bh).normalizedAnchors.size() / 2.0;
                return{ mg.data(bh).normalizedAnchors.front(), mg.data(bh).normalizedAnchors.back() };
            }


            template <class DataT, class SparseMatElementT>
            inline void RegisterConstraintEquations(int & eid, const MixedGraph & mg, MixedGraphPropertyTable & props,
                ComponentHandledTableFromConstraintGraph<int, MixedGraph>::type & uh2varStartPosition,
                ConstraintHandledTableFromConstraintGraph<std::vector<Vec3>, MixedGraph>::type & appliedBinaryAnchors,
                ConstraintHandledTableFromConstraintGraph<double, MixedGraph>::type & weightsForEachAppliedBinaryAnchor,
                std::vector<SparseMatElementT> & Atriplets,
                std::vector<SparseMatElementT> & Wtriplets,
                std::vector<double> & B) {

                for (auto & c : mg.constraints<DataT>()){
                    if (!props[c.topo.hd].used)
                        continue;
                    auto bh = c.topo.hd;
                    auto uh1 = mg.topo(bh).component<0>();
                    auto uh2 = mg.topo(bh).component<1>();
                    assert(props[uh1].used || props[uh2].used);
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


            template <class SparseMatElementT>
            void FormulateComponentsAndConstraintsAsMatricesForInverseDepthSolution(const MixedGraph & mg, MixedGraphPropertyTable & props,
                ComponentHandledTableFromConstraintGraph<int, MixedGraph>::type & uh2varStartPosition,
                ConstraintHandledTableFromConstraintGraph<int, MixedGraph>::type & bh2consStartPosition,
                ConstraintHandledTableFromConstraintGraph<std::vector<Vec3>, MixedGraph>::type & appliedBinaryAnchors,
                ConstraintHandledTableFromConstraintGraph<double, MixedGraph>::type & weightsForEachAppliedBinaryAnchor,
                int & varNum, int & consNum,
                std::vector<SparseMatElementT> & Atriplets,
                std::vector<SparseMatElementT> & Wtriplets,
                std::vector<double> & X,
                std::vector<double> & B,
                bool addAnchor = true){

                SetClock();

                THERE_ARE_BUGS_HERE("make sure mg is single connected");

                varNum = 0;
                for (auto & c : mg.components<LineData>()){
                    if (!props[c.topo.hd].used)
                        continue;
                    uh2varStartPosition[c.topo.hd] = varNum;
                    varNum += props[c.topo.hd].variables.size();
                }
                for (auto & c : mg.components<RegionData>()){
                    if (!props[c.topo.hd].used)
                        continue;
                    uh2varStartPosition[c.topo.hd] = varNum;
                    varNum += props[c.topo.hd].variables.size();
                }

                X.resize(varNum);
                for (auto & c : mg.components<LineData>()){
                    if (!props[c.topo.hd].used)
                        continue;
                    for (int i = 0; i < props[c.topo.hd].variables.size(); i++){
                        X[uh2varStartPosition.at(c.topo.hd) + i] = props[c.topo.hd].variables.at(i);
                    }
                }
                for (auto & c : mg.components<RegionData>()){
                    if (!props[c.topo.hd].used)
                        continue;
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
                    assert(props[c.topo.component<0>()].used && props[c.topo.component<1>()].used);
                    bh2consStartPosition[c.topo.hd] = consNum;
                    appliedBinaryAnchors[c.topo.hd] = NecessaryAnchorsForBinary(mg, c.topo.hd, 
                        weightsForEachAppliedBinaryAnchor[c.topo.hd]);
                    consNum += appliedBinaryAnchors[c.topo.hd].size();
                }
                for (auto & c : mg.constraints<LineRelationData>()){
                    if (!props[c.topo.hd].used)
                        continue;
                    assert(props[c.topo.component<0>()].used && props[c.topo.component<1>()].used);
                    bh2consStartPosition[c.topo.hd] = consNum;
                    appliedBinaryAnchors[c.topo.hd] = NecessaryAnchorsForBinary(mg, c.topo.hd,
                        weightsForEachAppliedBinaryAnchor[c.topo.hd]);
                    consNum += appliedBinaryAnchors[c.topo.hd].size();
                }
                for (auto & c : mg.constraints<RegionLineConnectionData>()){
                    if (!props[c.topo.hd].used)
                        continue;
                    assert(props[c.topo.component<0>()].used && props[c.topo.component<1>()].used);
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
                        if (!props[c.topo.hd].used)
                            continue;
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



        void SolveVariablesUsingInversedDepths(const MixedGraph & mg, MixedGraphPropertyTable & props, bool useWeights){

            SetClock();

            auto uh2varStartPosition = MakeHandledTableForAllComponents<int>(mg);
            auto bh2consStartPosition = MakeHandledTableForAllConstraints<int>(mg);
            auto appliedBinaryAnchors = MakeHandledTableForAllConstraints<std::vector<Vec3>>(mg);
            auto weightsForEachAppliedBinaryAnchor = MakeHandledTableForAllConstraints<double>(mg);

            int varNum = 0;
            int consNum = 0;
            std::vector<Eigen::Triplet<double>> Atriplets;
            std::vector<Eigen::Triplet<double>> Wtriplets;
            std::vector<double> Xdata;
            std::vector<double> Bdata;

            FormulateComponentsAndConstraintsAsMatricesForInverseDepthSolution(mg, props, 
                uh2varStartPosition, bh2consStartPosition,
                appliedBinaryAnchors, weightsForEachAppliedBinaryAnchor, 
                varNum, consNum, Atriplets, Wtriplets, Xdata, Bdata);
            

            Eigen::SparseMatrix<double> A;
            {
                Clock clock("form matrix A");
                A.resize(consNum, varNum);
                A.setFromTriplets(Atriplets.begin(), Atriplets.end());
            }
            Eigen::Map<const Eigen::VectorXd> B(Bdata.data(), Bdata.size());


            Eigen::SparseMatrix<double> WA;
            Eigen::VectorXd WB;
            if(useWeights) {
                Clock clock("form matrix WA");
                Eigen::SparseMatrix<double> W;
                W.resize(consNum, consNum);
                W.setFromTriplets(Wtriplets.begin(), Wtriplets.end());
                WA = W * A;
                WB = W * B;
            }

            static const bool useSPQR = true;

            Eigen::VectorXd X;
            if (!useSPQR) {
                Clock clock("solve equations using Eigen::SparseQR");

                A.makeCompressed();
                WA.makeCompressed();

                Eigen::SparseQR<Eigen::SparseMatrix<double>, Eigen::COLAMDOrdering<int>> solver;
                static_assert(!(Eigen::SparseMatrix<double>::IsRowMajor), "COLAMDOrdering only supports column major");

                solver.compute(useWeights ? WA : A);

                if (solver.info() != Eigen::Success) {
                    assert(0);
                    std::cout << "computation error" << std::endl;
                    return;
                }

                X = solver.solve(useWeights ? WB : B);
                if (solver.info() != Eigen::Success) {
                    assert(0);
                    std::cout << "solving error" << std::endl;
                    return;
                }
            }
            else{
                Clock clock("solve equations using Eigen::SPQR");

                Eigen::SPQR<Eigen::SparseMatrix<double, Eigen::ColMajor>> solver;
                solver.cholmodCommon()->useGPU = 1;
                solver.compute(useWeights ? WA : A);

                if (solver.info() != Eigen::Success) {
                    assert(0);
                    std::cout << "computation error" << std::endl;
                    return;
                }

                X = solver.solve(useWeights ? WB : B);
                if (solver.info() != Eigen::Success) {
                    assert(0);
                    std::cout << "solving error" << std::endl;
                    return;
                }
            }

            {
                Clock clock("install solved variables");
                for (auto & c : mg.components<RegionData>()){
                    if (!props[c.topo.hd].used)
                        continue;
                    int uhStartPosition = uh2varStartPosition.at(c.topo.hd);
                    for (int i = 0; i < props[c.topo.hd].variables.size(); i++){
                        props[c.topo.hd].variables[i] = X(uhStartPosition + i);
                    }
                }
                for (auto & c : mg.components<LineData>()){
                    if (!props[c.topo.hd].used)
                        continue;
                    int uhStartPosition = uh2varStartPosition.at(c.topo.hd);
                    for (int i = 0; i < props[c.topo.hd].variables.size(); i++){
                        props[c.topo.hd].variables[i] = X(uhStartPosition + i);
                    }
                }
            }
        }

        


        namespace {

            // Generic functor
            template <class InternalFunctorT, class T, int NX = Eigen::Dynamic, int NY = Eigen::Dynamic>
            struct GenericFunctor {
                typedef T Scalar;
                enum {
                    InputsAtCompileTime = NX,
                    ValuesAtCompileTime = NY
                };
                typedef Eigen::Matrix<Scalar, InputsAtCompileTime, 1> InputType;
                typedef Eigen::Matrix<Scalar, ValuesAtCompileTime, 1> ValueType;
                typedef Eigen::Matrix<Scalar, ValuesAtCompileTime, InputsAtCompileTime> JacobianType;

                const int m_inputs, m_values;
                InternalFunctorT m_fun;
                GenericFunctor(const InternalFunctorT & fun) : m_fun(fun), m_inputs(InputsAtCompileTime), m_values(ValuesAtCompileTime) {}
                GenericFunctor(const InternalFunctorT & fun, int inputs, int values) : m_fun(fun), m_inputs(inputs), m_values(values) {}

                int inputs() const { return m_inputs; }
                int values() const { return m_values; }

                // you should define that in the subclass :
                inline int operator() (const InputType& x, ValueType & v) const{
                    m_fun(x, v);
                    return 0;
                }
            };

            template <class T, int NX = Eigen::Dynamic, int NY = Eigen::Dynamic, class InternalFunctorT>
            GenericFunctor<InternalFunctorT, T, NX, NY> MakeGenericFunctor(const InternalFunctorT & fun){
                return GenericFunctor<InternalFunctorT, T, NX, NY>(fun);
            }

            template <class T, int NX = Eigen::Dynamic, int NY = Eigen::Dynamic, class InternalFunctorT>
            GenericFunctor<InternalFunctorT, T, NX, NY> MakeGenericFunctor(const InternalFunctorT & fun, int inputs, int values){
                return GenericFunctor<InternalFunctorT, T, NX, NY>(fun, inputs, values);
            }

            template <class T, int NX = Eigen::Dynamic, int NY = Eigen::Dynamic, class InternalFunctorT>
            Eigen::NumericalDiff<GenericFunctor<InternalFunctorT, T, NX, NY>> MakeGenericNumericDiffFunctor(const InternalFunctorT & fun){
                return Eigen::NumericalDiff<GenericFunctor<InternalFunctorT, T, NX, NY>>(GenericFunctor<InternalFunctorT, T, NX, NY>(fun));
            }

            template <class T, int NX = Eigen::Dynamic, int NY = Eigen::Dynamic, class InternalFunctorT>
            Eigen::NumericalDiff<GenericFunctor<InternalFunctorT, T, NX, NY>> MakeGenericNumericDiffFunctor(const InternalFunctorT & fun, int inputs, int values){
                return Eigen::NumericalDiff<GenericFunctor<InternalFunctorT, T, NX, NY>>(GenericFunctor<InternalFunctorT, T, NX, NY>(fun, inputs, values));
            }



        }

        
        void SolveVariablesUsingNormalDepths(const MixedGraph & mg, MixedGraphPropertyTable & props, bool useWeights){

            THERE_ARE_BUGS_HERE("weights not used yet");

            using namespace Eigen;

            auto uh2varStartPosition = MakeHandledTableForAllComponents<int>(mg);
            auto bh2consStartPosition = MakeHandledTableForAllConstraints<int>(mg);
            auto appliedBinaryAnchors = MakeHandledTableForAllConstraints<std::vector<Vec3>>(mg);
            auto weightsForEachAppliedBinaryAnchor = MakeHandledTableForAllConstraints<double>(mg);

            // initialize vectors
            int varNum, consNum;
            static const bool addAnchor = true;

            varNum = 0;
            for (auto & c : mg.components<LineData>()){
                if (!props[c.topo.hd].used)
                    continue;
                uh2varStartPosition[c.topo.hd] = varNum;
                varNum += props[c.topo.hd].variables.size();
            }
            for (auto & c : mg.components<RegionData>()){
                if (!props[c.topo.hd].used)
                    continue;
                uh2varStartPosition[c.topo.hd] = varNum;
                varNum += props[c.topo.hd].variables.size();
            }

            VectorXd X;
            X.resize(varNum);
            for (auto & c : mg.components<LineData>()){
                if (!props[c.topo.hd].used)
                    continue;
                for (int i = 0; i < props[c.topo.hd].variables.size(); i++){
                    X[uh2varStartPosition.at(c.topo.hd) + i] = props[c.topo.hd].variables.at(i);
                }
            }
            for (auto & c : mg.components<RegionData>()){
                if (!props[c.topo.hd].used)
                    continue;
                for (int i = 0; i < props[c.topo.hd].variables.size(); i++){
                    X[uh2varStartPosition.at(c.topo.hd) + i] = props[c.topo.hd].variables.at(i);
                }
            }

            RegionHandle rhAnchored;

            consNum = 0;
            if (addAnchor){
                consNum++;

                //IMPROVABLE_HERE("find a most connected component to set the anchor");

                // largest plane
                double maxArea = 0.0;
                for (auto & c : mg.components<RegionData>()){
                    if (!props[c.topo.hd].used)
                        continue;
                    if (c.data.area > maxArea){
                        rhAnchored = c.topo.hd;
                        maxArea = c.data.area;
                    }
                }
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

            auto computeDistanceAtAnchors = [consNum, &mg, &props, &uh2varStartPosition, &bh2consStartPosition, 
                &appliedBinaryAnchors, &weightsForEachAppliedBinaryAnchor, rhAnchored](
                const VectorXd & variables, VectorXd & distances){
                
                //std::unordered_set<int> filled;
                if (addAnchor){
                    int varStartPos = uh2varStartPosition.at(rhAnchored);
                    distances(0) = abs(DepthAtDirectionGivenVariables(mg, variables.data() + varStartPos, props, 
                        mg.data(rhAnchored).normalizedCenter, rhAnchored) - 1.0);
                    //filled.insert(0);
                }

                for (auto & c : mg.constraints<RegionBoundaryData>()){
                    if (!props[c.topo.hd].used)
                        continue;
                    int consStartPos = bh2consStartPosition.at(c.topo.hd);
                    auto & anchors = appliedBinaryAnchors.at(c.topo.hd);
                    auto uh1 = c.topo.component<0>();
                    auto uh2 = c.topo.component<1>();
                    int varStartPos1 = uh2varStartPosition.at(uh1);
                    int varStartPos2 = uh2varStartPosition.at(uh2);
                    for (int i = 0; i < anchors.size(); i++){
                        double depth1 = DepthAtDirectionGivenVariables(mg, variables.data() + varStartPos1, props, anchors[i], uh1);
                        double depth2 = DepthAtDirectionGivenVariables(mg, variables.data() + varStartPos2, props, anchors[i], uh2);
                        distances(consStartPos + i) = abs(depth1 - depth2);
                        //filled.insert(consStartPos + i);
                    }
                }

                for (auto & c : mg.constraints<RegionLineConnectionData>()){
                    if (!props[c.topo.hd].used)
                        continue;
                    int consStartPos = bh2consStartPosition.at(c.topo.hd);
                    auto & anchors = appliedBinaryAnchors.at(c.topo.hd);
                    auto uh1 = c.topo.component<0>();
                    auto uh2 = c.topo.component<1>();
                    int varStartPos1 = uh2varStartPosition.at(uh1);
                    int varStartPos2 = uh2varStartPosition.at(uh2);
                    for (int i = 0; i < anchors.size(); i++){
                        double depth1 = DepthAtDirectionGivenVariables(mg, variables.data() + varStartPos1, props, anchors[i], uh1);
                        double depth2 = DepthAtDirectionGivenVariables(mg, variables.data() + varStartPos2, props, anchors[i], uh2);
                        distances(consStartPos + i) = abs(depth1 - depth2);
                        //filled.insert(consStartPos + i);
                    }
                }

                for (auto & c : mg.constraints<LineRelationData>()){
                    if (!props[c.topo.hd].used)
                        continue;
                    int consStartPos = bh2consStartPosition.at(c.topo.hd);
                    auto & anchors = appliedBinaryAnchors.at(c.topo.hd);
                    auto uh1 = c.topo.component<0>();
                    auto uh2 = c.topo.component<1>();
                    int varStartPos1 = uh2varStartPosition.at(uh1);
                    int varStartPos2 = uh2varStartPosition.at(uh2);
                    for (int i = 0; i < anchors.size(); i++){
                        double depth1 = DepthAtDirectionGivenVariables(mg, variables.data() + varStartPos1, props, anchors[i], uh1);
                        double depth2 = DepthAtDirectionGivenVariables(mg, variables.data() + varStartPos2, props, anchors[i], uh2);
                        distances(consStartPos + i) = abs(depth1 - depth2);
                        //filled.insert(consStartPos + i);
                    }
                }

                //for (int i = 0; i < consNum; i++){
                //    assert(Contains(filled, i));
                //}

                distances *= 10000.0;

                std::cout << '.';

            };


            auto functor = MakeGenericNumericDiffFunctor<double>(computeDistanceAtAnchors, varNum, consNum);
            LevenbergMarquardt<decltype(functor)> lm(functor);
            lm.parameters.maxfev = 5000;
            lm.minimize(X);

            // install X
            for (auto & c : mg.components<RegionData>()){
                if (!props[c.topo.hd].used)
                    continue;
                int uhStartPosition = uh2varStartPosition.at(c.topo.hd);
                for (int i = 0; i < props[c.topo.hd].variables.size(); i++){
                    props[c.topo.hd].variables[i] = X(uhStartPosition + i);
                }
            }
            for (auto & c : mg.components<LineData>()){
                if (!props[c.topo.hd].used)
                    continue;
                int uhStartPosition = uh2varStartPosition.at(c.topo.hd);
                for (int i = 0; i < props[c.topo.hd].variables.size(); i++){
                    props[c.topo.hd].variables[i] = X(uhStartPosition + i);
                }
            }

        }


        double ComponentMedianCenterDepth(const MixedGraph & mg, const MixedGraphPropertyTable & props){
            for (auto & c : mg.components<RegionData>()){
                if (!props[c.topo.hd].used)
                    continue;
                for (auto & v : props[c.topo.hd].variables){
                    assert(!IsInfOrNaN(v));
                }
            }
            for (auto & c : mg.components<LineData>()){
                if (!props[c.topo.hd].used)
                    continue;
                for (auto & v : props[c.topo.hd].variables){
                    assert(!IsInfOrNaN(v));
                }
            }

            std::vector<double> centerDepths;
            centerDepths.reserve(mg.internalComponents<RegionData>().size() + mg.internalComponents<LineData>().size());
            for (auto & c : mg.components<RegionData>()){
                if (!props[c.topo.hd].used)
                    continue;
                double d = DepthAt(c.data.normalizedCenter, Instance(mg, props, c.topo.hd));
                if (!IsInfOrNaN(d)){
                    centerDepths.push_back(d);
                }
                assert(!IsInfOrNaN(centerDepths.back()));
            }
            for (auto & c : mg.components<LineData>()){
                if (!props[c.topo.hd].used)
                    continue;
                double d = DepthAt(normalize(c.data.line.center()), Instance(mg, props, c.topo.hd));
                if (!IsInfOrNaN(d)){
                    centerDepths.push_back(d);
                }
                assert(!IsInfOrNaN(centerDepths.back()));
            }
            std::nth_element(centerDepths.begin(), centerDepths.begin() + centerDepths.size() / 2, centerDepths.end());
            double medianCenterDepth = centerDepths[centerDepths.size() / 2];
            return medianCenterDepth;
        }


        double ComputeScore(const MixedGraph & mg, const MixedGraphPropertyTable & props){

            // manhattan fitness
            double sumOfComponentWeightedFitness = 0.0;
            double sumOfComponentWeights = 0.0;

            // lines
            double maxLineSpanAngle = 0.0;
            for (auto & l : mg.components<LineData>()){
                if (!props[l.topo.hd].used)
                    continue;
                double lineSpanAngle = AngleBetweenDirections(l.data.line.first, l.data.line.second);
                if (lineSpanAngle > maxLineSpanAngle)
                    maxLineSpanAngle = lineSpanAngle;
            }

            HandledTable<LineHandle, Line3> lines(mg.internalComponents<LineData>().size());
            HandledTable<RegionHandle, Plane3> planes(mg.internalComponents<RegionData>().size());

            for (auto & l : mg.components<LineData>()){
                if (!props[l.topo.hd].used)
                    continue;
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
                if (!props[r.topo.hd].used)
                    continue;
                if (r.data.area > maxRegionArea){
                    maxRegionArea = r.data.area;
                }
            }

            for (auto & r : mg.components<RegionData>()){
                if (!props[r.topo.hd].used)
                    continue;
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
                if (!props[c.topo.hd].used){
                    sumOfNotUsedConstraintWeights += ElementsNum(c.data.normalizedSampledPoints) * typeWeight;
                    continue;
                }

                auto & plane1 = planes[c.topo.component<0>()];
                auto & plane2 = planes[c.topo.component<1>()];
                for (auto & ss : c.data.normalizedSampledPoints){
                    for (auto & s : ss){
                        double d1 = DepthAt(s, plane1);
                        double d2 = DepthAt(s, plane2);
                        assert(d1 >= 0 && d2 >= 0);
                        
                        double fitness = Gaussian(abs(d1 - d2) / (d1 + d2), 0.1);
                        double weight = typeWeight;

                        sumOfConstraintWeightedFitness += (fitness * weight);
                        sumOfConstraintWeights += weight;
                    }
                }
            }
            //double maxJunctionWeight = 0.0;
            //for (auto & c : mg.constraints<LineRelationData>()){
            //    if (c.data.junctionWeight > maxJunctionWeight){
            //        maxJunctionWeight = c.data.junctionWeight;
            //    }
            //}
            for (auto & c : mg.constraints<LineRelationData>()){
                static const double typeWeight = 8.0;
                if (!props[c.topo.hd].used){
                    sumOfNotUsedConstraintWeights += /* c.data.junctionWeight / maxJunctionWeight **/ typeWeight;
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
                if (!props[c.topo.hd].used){
                    sumOfNotUsedConstraintWeights += ElementsNum(c.data.normalizedAnchors) * typeWeight;
                    continue;
                }

                auto & plane = planes[c.topo.component<0>()];
                auto & line = lines[c.topo.component<1>()];
                for (auto & a : c.data.normalizedAnchors){
                    double d1 = DepthAt(a, plane);
                    double d2 = DepthAt(a, line);
                    assert(d1 >= 0 && d2 >= 0);

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





        void NormalizeVariables(const MixedGraph & mg, MixedGraphPropertyTable & props){
            double medianCenterDepth = ComponentMedianCenterDepth(mg, props);
            
            SetClock();

            std::cout << "median center depth: " << medianCenterDepth << std::endl;

            assert(!IsInfOrNaN(medianCenterDepth));

            // normalize variables
            for (auto & c : mg.components<RegionData>()){
                if (!props[c.topo.hd].used)
                    continue;
                for (int i = 0; i < props[c.topo.hd].variables.size(); i++){
                    props[c.topo.hd].variables[i] *= medianCenterDepth;
                }
            }
            for (auto & c : mg.components<LineData>()){
                if (!props[c.topo.hd].used)
                    continue;
                for (int i = 0; i < props[c.topo.hd].variables.size(); i++){
                    props[c.topo.hd].variables[i] *= medianCenterDepth;
                }
            }
        }









        namespace {
            
            // returns false if confliction occurs
            bool MakeRegionPlaneUsable(RegionHandle rh, bool usable, MixedGraphPropertyTable & props) {
                auto & p = props[rh];
                if (p.used == usable)
                    return true;
                if (p.used && !usable){
                    p.used = false;
                    p.orientationClaz = p.orientationNotClaz = -1;
                    return true;
                }
                p.orientationClaz = p.orientationNotClaz = -1;
                return true;
            }


            // returns false if confliction occurs
            bool MakeRegionPlaneToward(RegionHandle rh, int normalVPId, MixedGraphPropertyTable & props){
                auto & p = props[rh];
                if (!p.used)
                    return true;
                assert(normalVPId != -1);
                if (p.orientationClaz != -1){
                    if (p.orientationClaz != normalVPId)
                        return false;
                    return true;
                }
                if (p.orientationNotClaz == -1) {
                    p.orientationClaz = normalVPId;
                    return true;
                }
                auto & dir = props.vanishingPoints[p.orientationNotClaz];
                if (IsFuzzyPerpendicular(props.vanishingPoints[normalVPId], dir)){
                    p.orientationClaz = normalVPId;
                    p.orientationNotClaz = -1;
                    return true;
                }
                return false;
            }

            // returns false if confliction occurs
            bool MakeRegionPlaneAlsoAlong(RegionHandle rh, int alongVPId, MixedGraphPropertyTable & props){
                auto & p = props[rh];
                if (!p.used)
                    return true;
                assert(alongVPId != -1);
                auto & dir = props.vanishingPoints[alongVPId];
                if (p.orientationClaz != -1){
                    auto & normal = props.vanishingPoints[p.orientationClaz];
                    return IsFuzzyPerpendicular(normal, dir);
                }
                if (p.orientationNotClaz == -1){
                    p.orientationNotClaz = alongVPId;
                    return true;
                }
                if (p.orientationNotClaz == alongVPId)
                    return true;

                auto newNormal = dir.cross(props.vanishingPoints[p.orientationNotClaz]);
                double minAngle = M_PI;
                for (int i = 0; i < props.vanishingPoints.size(); i++){
                    double angle = AngleBetweenUndirectedVectors(props.vanishingPoints[i], newNormal);
                    if (angle < minAngle){
                        p.orientationClaz = i;
                        minAngle = angle;
                    }
                }
                if (p.orientationClaz != -1){
                    p.orientationNotClaz = -1;
                    return true;
                }
                return false;
            }

        }


        static const int regionLineRelatedThreshold = 2;

        void AttachPrincipleDirectionConstraints(const MixedGraph & mg, MixedGraphPropertyTable & props,
            double rangeAngle, bool avoidLineConflictions){

            SetClock();

            // find peaky regions
            std::vector<std::vector<RegionHandle>> peakyRegionHandles(props.vanishingPoints.size());
            for (auto & r : mg.components<RegionData>()){
                auto h = r.topo.hd;
                if (!props[h].used)
                    continue;

                auto & contours = r.data.normalizedContours;
                double radiusAngle = 0.0;
                for (auto & cs : r.data.normalizedContours){
                    for (auto & c : cs){
                        double angle = AngleBetweenDirections(r.data.normalizedCenter, c);
                        if (angle > radiusAngle){
                            radiusAngle = angle;
                        }
                    }
                }

                bool mayCrossAnyVP = false;
                for (auto & vp : props.vanishingPoints){
                    double angle = AngleBetweenUndirectedVectors(vp, r.data.normalizedCenter);
                    if (angle < radiusAngle){
                        mayCrossAnyVP = true;
                        break;
                    }
                }

                if (!mayCrossAnyVP){
                    continue;
                }

                float ppcFocal = 100.0f;
                int ppcSize = 2 * radiusAngle * ppcFocal + 10;
                Vec3 x;
                std::tie(x, std::ignore) = ProposeXYDirectionsFromZDirection(r.data.normalizedCenter);
                PartialPanoramicCamera ppc(ppcSize, ppcSize, ppcFocal, Point3(0, 0, 0), r.data.normalizedCenter, x);
                Imageub mask = Imageub::zeros(ppc.screenSize());

                // project contours to ppc
                std::vector<std::vector<Point2i>> contourProjs(contours.size());
                for (int k = 0; k < contours.size(); k++){
                    auto & contourProj = contourProjs[k];
                    contourProj.reserve(contours[k].size());
                    for (auto & d : contours[k]){
                        contourProj.push_back(vec_cast<int>(ppc.screenProjection(d)));
                    }
                }
                cv::fillPoly(mask, contourProjs, (uint8_t)1);
                
                // intersection test
                for (int i = 0; i < props.vanishingPoints.size(); i++){
                    auto p1 = ToPixelLoc(ppc.screenProjection(props.vanishingPoints[i]));
                    auto p2 = ToPixelLoc(ppc.screenProjection(-props.vanishingPoints[i]));
                    
                    int dilateSize = ppcFocal * rangeAngle;
                    bool intersected = false;
                    for (int x = -dilateSize; x <= dilateSize; x++){
                        if (intersected)
                            break;
                        for (int y = -dilateSize; y <= dilateSize; y++){
                            if (intersected)
                                break;
                            auto pp1 = PixelLoc(p1.x + x, p1.y + y);
                            auto pp2 = PixelLoc(p2.x + x, p2.y + y);
                            if (Contains(mask, pp1) /*Box<int, 2>(Point2i(0, 0), Point2i(mask.cols - 1, mask.rows - 1)).contains(pp1)*/ && mask(pp1)){
                                peakyRegionHandles[i].push_back(r.topo.hd);
                                intersected = true;
                            }
                            else if (Contains(mask, pp2)/*Box<int, 2>(Point2i(0, 0), Point2i(mask.cols - 1, mask.rows - 1)).contains(pp2) */&& mask(pp2)){
                                peakyRegionHandles[i].push_back(r.topo.hd);
                                intersected = true;
                            }
                        }
                    }
                }
            }
       
            for (int i = 0; i < props.vanishingPoints.size(); i++){
                auto & vp = props.vanishingPoints[i];
                auto & rhs = peakyRegionHandles[i];
                for (auto rh : rhs){
                    if (avoidLineConflictions){
                        bool hasLineConflictions = false;
                        for (auto conh : mg.topo(rh).constraints<RegionLineConnectionData>()){
                            if (mg.data(conh).normalizedAnchors.size() < regionLineRelatedThreshold)
                                continue;
                            LineHandle lh = mg.topo(conh).component<1>();
                            auto & lprop = props[lh];
                            if (lprop.orientationClaz == i){
                                hasLineConflictions = true;
                                break;
                            }
                        }
                        if (!hasLineConflictions)
                            MakeRegionPlaneToward(rh, i, props);
                    }
                    else {
                        MakeRegionPlaneToward(rh, i, props);
                        for (auto conh : mg.topo(rh).constraints<RegionLineConnectionData>()){
                            if (mg.data(conh).normalizedAnchors.size() < regionLineRelatedThreshold)
                                continue;
                            LineHandle lh = mg.topo(conh).component<1>();
                            auto & lprop = props[lh];
                            if (lprop.orientationClaz == i){
                                lprop.orientationClaz = -1;
                            }
                        }
                    }
                }
            }

            ForeachMixedGraphConstraintHandle(mg, ConstraintUsableFlagUpdater{ mg, props });
            ForeachMixedGraphComponentHandle(mg, InitializeVariablesForEachHandle{ mg, props });

        }


        void AttachWallConstriants(const MixedGraph & mg, MixedGraphPropertyTable & props,
            double rangeAngle, const Vec3 & verticalSeed){

            SetClock();

            int vertVPId = -1;
            double minAngle = M_PI;
            for (int i = 0; i < props.vanishingPoints.size(); i++){
                double angle = AngleBetweenUndirectedVectors(verticalSeed, props.vanishingPoints[i]);
                if (angle < minAngle){
                    minAngle = angle;
                    vertVPId = i;
                }
            }
            assert(vertVPId != -1);
            auto & vertical = props.vanishingPoints[vertVPId];

            std::vector<RegionHandle> horizontalRegionHandles;
            for (auto & r : mg.components<RegionData>()){
                auto h = r.topo.hd;
                if (!props[h].used)
                    continue;
                auto & contours = r.data.normalizedContours;
                bool intersected = false;
                for (auto & cs : r.data.normalizedContours){
                    if (intersected)
                        break;
                    for (auto & c : cs){
                        double angle = M_PI_2 - AngleBetweenUndirectedVectors(c, vertical);
                        if (angle <= rangeAngle){
                            intersected = true;
                            break;
                        }
                    }
                }

                if (intersected){
                    horizontalRegionHandles.push_back(h);
                }
            }

            for (auto h : horizontalRegionHandles){
                MakeRegionPlaneAlsoAlong(h, vertVPId, props);
            }

            ForeachMixedGraphConstraintHandle(mg, ConstraintUsableFlagUpdater{ mg, props });
            ForeachMixedGraphComponentHandle(mg, InitializeVariablesForEachHandle{ mg, props });
        }



        void AttachFloorAndCeilingConstraints(const MixedGraph & mg, MixedGraphPropertyTable & props,
            double eyeHeightRatioLowerBound, double eyeHeightRatioUpperBound,
            double angleThreshold, const Vec3 & verticalSeed){

            SetClock();

            IMPROVABLE_HERE("we need a better way!");

            int vertVPId = -1;
            double minAngle = M_PI;
            for (int i = 0; i < props.vanishingPoints.size(); i++){
                double angle = AngleBetweenUndirectedVectors(verticalSeed, props.vanishingPoints[i]);
                if (angle < minAngle){
                    minAngle = angle;
                    vertVPId = i;
                }
            }
            assert(vertVPId != -1);
            auto & vertical = props.vanishingPoints[vertVPId];

            int horizVPId = -1;
            double maxAngle = 0;
            for (int i = 0; i < props.vanishingPoints.size(); i++){
                double angle = AngleBetweenUndirectedVectors(vertical, props.vanishingPoints[i]);
                if (angle > maxAngle){
                    maxAngle = angle;
                    horizVPId = i;
                }
            }
            assert(horizVPId != -1);
            auto xDir = normalize(props.vanishingPoints[horizVPId]);
            auto yDir = normalize(vertical.cross(xDir));

            // compute horizontal bounding range of may-be ceilings and floors
            double xmin = std::numeric_limits<double>::max(), ymin = std::numeric_limits<double>::max();
            double xmax = std::numeric_limits<double>::lowest(), ymax = std::numeric_limits<double>::lowest();

            std::vector<Scored<RegionHandle>> maybeCeilinsOrFloors[2];

            auto regionPlanes = MakeHandledTableForComponents<Plane3, RegionData>(mg);
            for (auto & r : mg.components<RegionData>()){
                if (!props[r.topo.hd].used)
                    continue;

                auto plane = Instance(mg, props, r.topo.hd);
                regionPlanes[r.topo.hd] = plane;
                double angleToVert = AngleBetweenUndirectedVectors(plane.normal, vertical);
                if (angleToVert < angleThreshold){
                    int belonging = r.data.normalizedCenter.dot(vertical) < 0 ? 0 : 1;                    
                    double centerZ = IntersectionOfLineAndPlane(Ray3(Point3(0, 0, 0), r.data.normalizedCenter), plane).position.dot(vertical);
                    maybeCeilinsOrFloors[belonging].push_back(ScoreAs(r.topo.hd, centerZ));

                    for (auto & cs : r.data.normalizedContours){
                        for (auto & c : cs){
                            Point3 point = IntersectionOfLineAndPlane(Ray3(Point3(0, 0, 0), c), plane).position;
                            float x = point.dot(xDir), y = point.dot(yDir);
                            if (x < xmin) xmin = x;
                            if (x > xmax) xmax = x;
                            if (y < ymin) ymin = y;
                            if (y > ymax) ymax = y; 
                        }
                    }
                }
            }

            if (maybeCeilinsOrFloors[0].empty() && maybeCeilinsOrFloors[1].empty())
                return;

            // estimate height or ceiling/floor, remove non-ceiling/floor horizontal regions
            double medianCenterDepth = ComponentMedianCenterDepth(mg, props);
            static const double heightAffectRange = medianCenterDepth * 0.01;
            for (int i = 0; i < 2; i++){
                auto & rhsWithHeights = maybeCeilinsOrFloors[i];
                std::sort(rhsWithHeights.begin(), rhsWithHeights.end());
                std::vector<double> votes(rhsWithHeights.size(), 0.0);
                for (int a = 0; a < rhsWithHeights.size(); a++){
                    for (int b = 0; b < rhsWithHeights.size(); b++){
                        votes[a] += Gaussian(rhsWithHeights[a].score - rhsWithHeights[b].score,
                            heightAffectRange);
                    }
                }
                int maxId = std::distance(votes.begin(), std::max_element(votes.begin(), votes.end()));
                
                std::vector<Scored<RegionHandle>> hs; 
                hs.reserve(rhsWithHeights.size());
                for (auto rhh : rhsWithHeights){
                    if (Distance(rhh.score, rhsWithHeights[maxId].score) < heightAffectRange){
                        hs.push_back(rhh);
                    }
                }
                rhsWithHeights = std::move(hs);
            }

            assert(xmax > xmin && ymax > ymin);
            int imSize = 300;
            double scale = imSize / ((xmax + ymax - xmin - ymin) / 2.0);
            int width = static_cast<int>((xmax - xmin) * scale + 1);
            int height = static_cast<int>((ymax - ymin) * scale + 1);
            Imageub ceilingOrFloorHorizontalRegionMasks[] = {
                Imageub::zeros(cv::Size(width, height)),
                Imageub::zeros(cv::Size(width, height))
            };
            for (int i = 0; i < 2; i++){
                for (auto & rhh : maybeCeilinsOrFloors[i]){
                    auto rh = rhh.component;
                    auto & contours = mg.data(rh).normalizedContours;
                    std::vector<std::vector<Point2i>> contourResized;
                    contourResized.reserve(contours.size());
                    for (auto & cs : contours){
                        std::vector<Point2i> csresized;
                        csresized.reserve(cs.size());
                        for (auto & c : cs){
                            Point3 point = IntersectionOfLineAndPlane(Ray3(Point3(0, 0, 0), c), regionPlanes[rh]).position;
                            float x = point.dot(xDir), y = point.dot(yDir);
                            csresized.emplace_back(static_cast<int>((x - xmin) * scale), 
                                static_cast<int>((y - ymin) * scale));
                        }
                        contourResized.push_back(std::move(csresized));
                    }
                    cv::fillPoly(ceilingOrFloorHorizontalRegionMasks[i], contourResized, (uint8_t)1);
                }
            }

            Imagei ceilingOrFloorRegionMasks[] = {
                Imagei(cv::Size(width, height), -1),
                Imagei(cv::Size(width, height), -1)
            };
            for (auto & r : mg.components<RegionData>()){
                const auto & prop = props[r.topo.hd];
                if (!prop.used)
                    continue;
                if (prop.orientationClaz != -1) // already strongly constrained
                    continue;
                if (prop.orientationNotClaz == vertVPId) // constrained to be vertical
                    continue;

                auto & plane = regionPlanes[r.topo.hd];
                auto & contours = r.data.normalizedContours;
                
                std::vector<std::vector<Point2i>> contourResized;
                contourResized.reserve(contours.size());
                
                bool hasPosDot = false, hasNegDot = false;
                for (auto & cs : contours){
                    std::vector<Point2i> csresized;
                    csresized.reserve(cs.size());
                    for (auto & c : cs){
                        Point3 point = IntersectionOfLineAndPlane(Ray3(Point3(0, 0, 0), c), plane).position;
                        float x = point.dot(xDir), y = point.dot(yDir);
                        float z = point.dot(vertical);
                        hasPosDot |= z >= 0;
                        hasNegDot |= z <= 0;
                        csresized.emplace_back(static_cast<int>((x - xmin) * scale),
                            static_cast<int>((y - ymin) * scale));
                    }
                    contourResized.push_back(std::move(csresized));
                }

                if (hasPosDot && hasNegDot){ // this region crosses the horizon
                    continue;
                }
                if (contourResized.empty()){
                    continue;
                }

                auto & mask = ceilingOrFloorRegionMasks[hasNegDot ? 0 : 1];
                cv::fillPoly(mask, contourResized, r.topo.hd.id);
            }

            // match masks
            for (auto & horizontalMask : ceilingOrFloorHorizontalRegionMasks){
                static const int erosionSize = 5;
                cv::erode(horizontalMask, horizontalMask, cv::getStructuringElement(cv::MORPH_ELLIPSE,
                    cv::Size(2 * erosionSize + 1, 2 * erosionSize + 1),
                    cv::Point(erosionSize, erosionSize)));
            }

            //cv::imshow("0", ceilingOrFloorHorizontalRegionMasks[0] * 255);
            //cv::imshow("1", ceilingOrFloorHorizontalRegionMasks[1] * 255);
            //cv::waitKey();

            auto regionAffectedCounts = mg.createComponentTable<RegionData>(0);

            for (int i = 0; i < 2; i++){
                auto & regionMasksHere = ceilingOrFloorRegionMasks[i];
                auto & horizontalMaskThere = ceilingOrFloorHorizontalRegionMasks[1 - i];
                for (auto it = regionMasksHere.begin(); it != regionMasksHere.end(); ++it){
                    if (!horizontalMaskThere(it.pos()))
                        continue;
                    RegionHandle regionHandle(*it);
                    if (regionHandle.invalid())
                        continue;
                    
                    regionAffectedCounts[regionHandle] ++;
                    // add horizontal constraint
                    if (regionAffectedCounts[regionHandle] > Square(imSize * 0.02)){
                        MakeRegionPlaneToward(regionHandle, vertVPId, props);
                    }
                }
            }

            ForeachMixedGraphConstraintHandle(mg, ConstraintUsableFlagUpdater{ mg, props });
            ForeachMixedGraphComponentHandle(mg, InitializeVariablesForEachHandle{ mg, props });

        }



        void AttachGeometricContextConstraints(const MixedGraph & mg, MixedGraphPropertyTable & props,
            const std::vector<GeometricContextEstimator::Feature> & perspectiveGCs,
            const std::vector<PerspectiveCamera> & gcCameras,
            int shrinkRegionOrientationIteration, bool considerGCVerticalConstraint){

            auto orientationVotes = mg.createComponentTable<RegionData, std::unordered_map<OrientationHint, double>>(
                std::unordered_map<OrientationHint, double>((size_t)OrientationHint::Count));

            assert(perspectiveGCs.size() == gcCameras.size());

            // region views
            auto regionMaskViews = mg.createComponentTable<RegionData, View<PartialPanoramicCamera, Imageub>>();
            for (auto & r : mg.components<RegionData>()){
                auto h = r.topo.hd;
                auto & contours = r.data.normalizedContours;
                double radiusAngle = 0.0;
                for (auto & cs : r.data.normalizedContours){
                    for (auto & c : cs){
                        double angle = AngleBetweenDirections(r.data.normalizedCenter, c);
                        if (angle > radiusAngle){
                            radiusAngle = angle;
                        }
                    }
                }
                float ppcFocal = 100.0f;
                int ppcSize = 2 * radiusAngle * ppcFocal;
                Vec3 x;
                std::tie(x, std::ignore) = ProposeXYDirectionsFromZDirection(r.data.normalizedCenter);
                PartialPanoramicCamera ppc(ppcSize, ppcSize, ppcFocal, Point3(0, 0, 0), r.data.normalizedCenter, x);
                Imageub mask = Imageub::zeros(ppc.screenSize());

                // project contours to ppc
                std::vector<std::vector<Point2i>> contourProjs(contours.size());
                for (int k = 0; k < contours.size(); k++){
                    auto & contourProj = contourProjs[k];
                    contourProj.reserve(contours[k].size());
                    for (auto & d : contours[k]){
                        contourProj.push_back(vec_cast<int>(ppc.screenProjection(d)));
                    }
                }
                cv::fillPoly(mask, contourProjs, (uint8_t)1);
                regionMaskViews[h].camera = ppc;
                regionMaskViews[h].image = mask;
            }

            auto regionAreas = mg.createComponentTable<RegionData, double>(0.0);
            for (auto & r : mg.components<RegionData>()){
                auto h = r.topo.hd;
                auto & ppMask = regionMaskViews[h];
                double area = cv::sum(ppMask.image).val[0];
                regionAreas[h] = area;
                for (int i = 0; i < perspectiveGCs.size(); i++){
                    auto sampler = MakeCameraSampler(ppMask.camera, gcCameras[i]);
                    auto occupation = sampler(Imageub::ones(gcCameras[i].screenSize()), cv::BORDER_CONSTANT, 0.0);
                    double areaOccupied = cv::sum(ppMask.image & occupation).val[0];
                    double ratio = areaOccupied / area;
                    assert(ratio <= 1.01);
                    GeometricContextLabel maxLabel = GeometricContextLabel::None;
                    double maxScore = 0.0;
                    for (auto & gcc : perspectiveGCs[i]){
                        auto label = gcc.first;
                        Imaged ppGC = sampler(gcc.second);
                        double score = cv::mean(ppGC, ppMask.image).val[0];
                        if (score > maxScore){
                            maxScore = score;
                            maxLabel = label;
                        }
                    }
                    if (maxLabel != GeometricContextLabel::None){
                        auto & v = orientationVotes[r.topo.hd][ToOrientationHint(maxLabel, considerGCVerticalConstraint)];
                        v = std::max(v, ratio);
                    }
                }
            }

            // optimize
            GCoptimizationGeneralGraph graph(mg.internalComponents<RegionData>().size(), (int)OrientationHint::Count);

            double maxBoundaryLength = 0;
            for (auto & b : mg.constraints<RegionBoundaryData>()){
                maxBoundaryLength = std::max(maxBoundaryLength, b.data.length);
            }
            for (auto & b : mg.constraints<RegionBoundaryData>()){
                graph.setNeighbors(b.topo.component<0>().id, b.topo.component<1>().id,
                    100 * b.data.length / maxBoundaryLength);
            }
            for (auto & r : mg.components<RegionData>()){
                auto & orientationVote = orientationVotes[r.topo.hd];
                // normalize votes
                double votesSum = 0.0;
                for (auto & v : orientationVote){
                    votesSum += v.second;
                }
                assert(votesSum >= 0.0);
                if (votesSum > 0.0){
                    for (auto & v : orientationVote){
                        assert(v.second >= 0.0);
                        v.second /= votesSum;
                    }
                }
                for (int label = 0; label < (int)OrientationHint::Count; label++){
                    double vote = Contains(orientationVote, (OrientationHint)label) ?
                        orientationVote.at((OrientationHint)label) : 0.0;
                    graph.setDataCost(r.topo.hd.id, label,
                        votesSum == 0.0 ? 0 : (10000 * (1.0 - vote)));
                }
            }
            for (int label1 = 0; label1 < (int)OrientationHint::Count; label1++){
                for (int label2 = 0; label2 < (int)OrientationHint::Count; label2++){
                    if (label1 == label2){
                        graph.setSmoothCost(label1, label2, 0);
                    }
                    else{
                        graph.setSmoothCost(label1, label2, 1);
                    }
                }
            }

            graph.expansion();
            graph.swap();

            // get the most vertical vp id
            int vVPId = -1;
            double angleToVert = std::numeric_limits<double>::max();
            for (int i = 0; i < props.vanishingPoints.size(); i++){
                double a = AngleBetweenUndirectedVectors(props.vanishingPoints[i], Vec3(0, 0, 1));
                if (a < angleToVert){
                    vVPId = i;
                    angleToVert = a;
                }
            }
            assert(vVPId != -1);

            // disorient outsided oriented gc labels
            auto regionOrientations = mg.createComponentTable<RegionData, OrientationHint>(OrientationHint::Void);
            for (auto & r : mg.components<RegionData>()){
                regionOrientations[r.topo.hd] = (OrientationHint)(graph.whatLabel(r.topo.hd.id));
            }

            for (int i = 0; i < shrinkRegionOrientationIteration; i++){
                auto shrinked = regionOrientations;
                for (auto & b : mg.constraints<RegionBoundaryData>()){
                    auto l1 = regionOrientations[b.topo.component<0>()];
                    auto l2 = regionOrientations[b.topo.component<1>()];
                    if (l1 == OrientationHint::Horizontal && l2 != OrientationHint::Horizontal){
                        shrinked[b.topo.component<0>()] = OrientationHint::OtherPlanar;
                    }
                    else if (l1 != OrientationHint::Horizontal && l2 == OrientationHint::Horizontal){
                        shrinked[b.topo.component<1>()] = OrientationHint::OtherPlanar;
                    }
                    if (l1 == OrientationHint::Vertical && l2 != OrientationHint::Vertical){
                        shrinked[b.topo.component<0>()] = OrientationHint::OtherPlanar;
                    }
                    else if (l1 != OrientationHint::Vertical && l2 == OrientationHint::Vertical){
                        shrinked[b.topo.component<1>()] = OrientationHint::OtherPlanar;
                    }
                }
                regionOrientations = std::move(shrinked);
            }

            // install gc labels
            for (auto & r : mg.components<RegionData>()){
                OrientationHint oh = regionOrientations[r.topo.hd];
                switch (oh) {
                case OrientationHint::Void:
                    /*props[r.topo.hd].used = false;
                    props[r.topo.hd].orientationClaz = -1;
                    props[r.topo.hd].orientationNotClaz = -1;*/
                    MakeRegionPlaneUsable(r.topo.hd, false, props);
                    break;
                case OrientationHint::Horizontal:
                    /*props[r.topo.hd].used = true;
                    props[r.topo.hd].orientationClaz = vVPId;
                    props[r.topo.hd].orientationNotClaz = -1;*/
                    MakeRegionPlaneToward(r.topo.hd, vVPId, props);
                    break;
                case OrientationHint::Vertical:
                    //props[r.topo.hd].used = true;
                    //props[r.topo.hd].orientationClaz = -1;
                    //props[r.topo.hd].orientationNotClaz = vVPId;
                    MakeRegionPlaneAlsoAlong(r.topo.hd, vVPId, props);
                    break;
                case OrientationHint::OtherPlanar:
                    /*props[r.topo.hd].used = true;
                    props[r.topo.hd].orientationClaz = -1;
                    props[r.topo.hd].orientationNotClaz = -1;*/
                    //MakeRegionPlaneUsable(r.topo.hd, true, props);
                    break;
                case OrientationHint::NonPlanar:
                    /*props[r.topo.hd].used = false;
                    props[r.topo.hd].orientationClaz = -1;
                    props[r.topo.hd].orientationNotClaz = -1;*/
                    MakeRegionPlaneUsable(r.topo.hd, false, props);
                    break;
                default:
                    assert(0);
                    break;
                }
            }

            /*
            // detect big aligned rectangular regions and orient them
            double regionAreaSum = std::accumulate(regionAreas.begin(), regionAreas.end(), 0.0);
            assert(regionAreaSum > 0);

            static const double dotThreshold = 0.001;
            static const double spanAngleThreshold = DegreesToRadians(10);
            static const double wholeDotThreshold = 0.01;
            static const double wholeSpanAngleThreshold = DegreesToRadians(10);

            HandledTable<RegionBoundaryHandle, std::vector<double>> boundaryMaxSpanAnglesForVPs(mg.internalConstraints<RegionBoundaryData>().size());
            for (auto & b : mg.constraints<RegionBoundaryData>()){
            auto & edges = b.data.normalizedEdges;
            std::vector<double> maxSpanAngleForVPs(props.vanishingPoints.size(), 0.0);
            if (b.topo.hd.id == 849){
            std::cout << std::endl;
            }

            for (auto & edge : edges){
            std::vector<std::vector<bool>> edgeVPFlags(props.vanishingPoints.size(),
            std::vector<bool>(edge.size() - 1, false));
            for (int i = 0; i < edge.size() - 1; i++){
            Vec3 n = normalize(edge[i].cross(edge[i + 1]));
            for (int k = 0; k < props.vanishingPoints.size(); k++){
            auto vp = normalize(props.vanishingPoints[k]);
            double dotv = abs(vp.dot(n));
            if (dotv <= dotThreshold){
            edgeVPFlags[k][i] = true;
            }
            }
            }
            // detect aligned longest line spanAngles from edge
            for (int k = 0; k < props.vanishingPoints.size(); k++){
            auto vp = normalize(props.vanishingPoints[k]);
            auto & edgeFlags = edgeVPFlags[k];
            int lastHead = -1, lastTail = -1;
            bool inChain = false;

            double & maxSpanAngle = maxSpanAngleForVPs[k];
            for (int j = 0; j <= edgeFlags.size(); j++){
            if (!inChain && j < edgeFlags.size() && edgeFlags[j]){
            lastHead = j;
            inChain = true;
            }
            else if (inChain && (j == edgeFlags.size() || !edgeFlags[j])){
            lastTail = j;
            inChain = false;
            // examine current chain
            assert(lastHead != -1);
            // compute full span angle
            double spanAngle = 0.0;
            for (int i = lastHead; i < lastTail; i++){
            spanAngle += AngleBetweenDirections(edge[i], edge[i + 1]);
            }
            if (spanAngle < spanAngleThreshold){
            continue;
            }
            // fit line
            const Vec3 & midCorner = edge[(lastHead + lastTail) / 2];
            Vec3 commonNormal = normalize(midCorner.cross(vp));
            std::vector<double> dotsToCommonNormal(lastTail - lastHead + 1, 0.0);
            for (int i = lastHead; i <= lastTail; i++){
            dotsToCommonNormal[i - lastHead] = abs(edge[i].dot(commonNormal));
            }
            if (*std::max_element(dotsToCommonNormal.begin(), dotsToCommonNormal.end()) <= wholeDotThreshold){
            // acceptable!
            if (spanAngle > maxSpanAngle){
            maxSpanAngle = spanAngle;
            }
            }
            }
            }
            }
            }
            if ((b.topo.component<0>().id == 355 || b.topo.component<1>().id == 355) && std::accumulate(maxSpanAngleForVPs.begin(), maxSpanAngleForVPs.end(), 0.0) > 0){
            std::cout << std::endl;
            }
            boundaryMaxSpanAnglesForVPs[b.topo.hd] = std::move(maxSpanAngleForVPs);
            }
            for (auto & r : mg.components<RegionData>()){
            if (!props[r.topo.hd].used)
            continue;
            if (props[r.topo.hd].orientationClaz != -1 || props[r.topo.hd].orientationNotClaz != -1)
            continue;

            if (r.topo.hd.id == 355){
            std::cout << std::endl;
            }

            double area = regionAreas[r.topo.hd];
            if (area / regionAreaSum < 0.005){
            continue;
            }
            std::vector<double> spanAngleSumsForVPs(props.vanishingPoints.size(), 0.0);
            for (auto & h : r.topo.constraints<RegionBoundaryData>()){
            for (int i = 0; i < props.vanishingPoints.size(); i++){
            spanAngleSumsForVPs[i] += boundaryMaxSpanAnglesForVPs[h][i];
            }
            }
            auto maxIter = std::max_element(spanAngleSumsForVPs.begin(), spanAngleSumsForVPs.end());
            int firstMaxVPId = std::distance(spanAngleSumsForVPs.begin(), maxIter);
            double firstMaxSpanAngle = *maxIter;
            *maxIter = -1.0;
            maxIter = std::max_element(spanAngleSumsForVPs.begin(), spanAngleSumsForVPs.end());
            int secondMaxVPId = std::distance(spanAngleSumsForVPs.begin(), maxIter);
            double secondMaxSpanAngle = *maxIter;
            if (secondMaxSpanAngle >= wholeSpanAngleThreshold){
            assert(firstMaxVPId != secondMaxVPId);
            Vec3 normal = props.vanishingPoints[firstMaxVPId].cross(props.vanishingPoints[secondMaxVPId]);
            double minAngle = 0.1;
            int bestMatchedVPId = -1;
            for (int i = 0; i < props.vanishingPoints.size(); i++){
            double angle = AngleBetweenUndirectedVectors(normal, props.vanishingPoints[i]);
            if (angle < minAngle){
            minAngle = angle;
            bestMatchedVPId = i;
            }
            }
            if (bestMatchedVPId != -1){
            std::cout << "vpid:::: " << bestMatchedVPId << std::endl;
            props[r.topo.hd].orientationClaz = bestMatchedVPId;
            props[r.topo.hd].orientationNotClaz = -1;
            }
            }
            else if (firstMaxSpanAngle >= wholeSpanAngleThreshold){
            std::cout << "vpid-along:::: " << firstMaxVPId << std::endl;
            props[r.topo.hd].orientationClaz = -1;
            props[r.topo.hd].orientationNotClaz = firstMaxVPId;
            }
            }*/

            //for (auto & l : mg.internalComponents<LineData>()){
            //    props[l.topo.hd].used = true;
            //    props[l.topo.hd].orientationClaz = l.data.initialClaz;
            //    props[l.topo.hd].orientationNotClaz = -1;
            //}

            ForeachMixedGraphConstraintHandle(mg, ConstraintUsableFlagUpdater{ mg, props });
            ForeachMixedGraphComponentHandle(mg, InitializeVariablesForEachHandle{ mg, props });

        }









        void LooseOrientationConstraintsOnComponents(const MixedGraph & mg, MixedGraphPropertyTable & props,
            double linesLoosableRatio, double regionsLoosableRatio, double distThresRatio){

            SetClock();

            if (linesLoosableRatio <= 0.0 && regionsLoosableRatio <= 0.0)
                return;
            
            double medianCenterDepth = ComponentMedianCenterDepth(mg, props);
            double distThres = medianCenterDepth * distThresRatio;
            double nnSearchRange = distThres;

            bool directApproach = true;
            if (directApproach){

                // collect current instances
                HandledTable<LineHandle, Line3> lines(mg.internalComponents<LineData>().size());
                HandledTable<RegionHandle, Plane3> planes(mg.internalComponents<RegionData>().size());

                for (auto & r : mg.components<RegionData>()){
                    if (!props[r.topo.hd].used)
                        continue;
                    planes[r.topo.hd] = Instance(mg, props, r.topo.hd);
                }
                for (auto & l : mg.components<LineData>()){
                    if (!props[l.topo.hd].used)
                        continue;
                    lines[l.topo.hd] = Instance(mg, props, l.topo.hd);
                }

                if (linesLoosableRatio > 0.0){
                    std::vector<LineHandle> usableLineHandles;
                    usableLineHandles.reserve(mg.internalComponents<LineData>().size());
                    auto distanceToNearestRegions = mg.createComponentTable<LineData>(0.0);
                    for (auto & l : mg.components<LineData>()){
                        if (!props[l.topo.hd].used)
                            continue;
                        usableLineHandles.push_back(l.topo.hd);

                        Point3 points[] = { lines[l.topo.hd].first, lines[l.topo.hd].second };
                        double & ndist = distanceToNearestRegions[l.topo.hd];
                        for (auto & linePoint : points){
                            double distToThis = std::numeric_limits<double>::max();
                            for (auto rlch : l.topo.constraints<RegionLineConnectionData>()){
                                auto & anchors = mg.data(rlch).normalizedAnchors;
                                RegionHandle rh = mg.topo(rlch).component<0>();
                                Plane3 & plane = planes[rh];
                                for (auto & anchor : anchors){
                                    auto planePoint = IntersectionOfLineAndPlane(Ray3(Point3(0, 0, 0), anchor), plane).position;
                                    double distance = Distance(planePoint, linePoint);
                                    if (distance < distToThis){
                                        distToThis = distance;
                                    }
                                }
                            }
                            if (distToThis > ndist){
                                ndist = distToThis;
                            }
                        }
                    }

                    std::sort(usableLineHandles.begin(), usableLineHandles.end(), [&distanceToNearestRegions](LineHandle l1, LineHandle l2){
                        return distanceToNearestRegions[l1] > distanceToNearestRegions[l2];
                    });

                    // loose orientation constraints on these lines
                    for (int i = 0; i < usableLineHandles.size() * linesLoosableRatio; i++){
                        double maxCornerDist = distanceToNearestRegions[usableLineHandles[i]];
                        if (maxCornerDist < distThres)
                            continue;
                        auto & lp = props[usableLineHandles[i]];
                        assert(lp.used);
                        if (lp.orientationClaz >= 0){
                            lp.orientationClaz = -1;
                            lp.variables = { 1.0, 1.0 };
                        }
                        assert(lp.orientationNotClaz == -1);
                        if (lp.orientationClaz == -1){
                            lp.used = false; // disable this line!!!!!!!!!
                            lp.variables.clear();
                        }
                    }
                }

                if (regionsLoosableRatio > 0.0){
                    std::vector<RegionHandle> usableRegionHandles;
                    usableRegionHandles.reserve(mg.internalComponents<RegionData>().size());
                    auto maxDistanceToNearestComponents = mg.createComponentTable<RegionData>(0.0);

                    for (auto & r : mg.components<RegionData>()){
                        if (!props[r.topo.hd].used)
                            continue;
                        usableRegionHandles.push_back(r.topo.hd);
                        Plane3 & thisPlane = planes[r.topo.hd];

                        double & ndist = maxDistanceToNearestComponents[r.topo.hd];

                        // region region
                        for (auto rrch : r.topo.constraints<RegionBoundaryData>()){
                            RegionHandle thatRh = mg.topo(rrch).component<0>();
                            if (thatRh == r.topo.hd){
                                thatRh = mg.topo(rrch).component<1>();
                            }
                            Plane3 & thatPlane = planes[thatRh];
                            for (auto & cs : mg.data(rrch).normalizedSampledPoints){
                                for (auto & c : cs){
                                    Point3 pointHere = IntersectionOfLineAndPlane(Ray3(Point3(0, 0, 0), c), thisPlane).position;
                                    Point3 pointThere = IntersectionOfLineAndPlane(Ray3(Point3(0, 0, 0), c), thatPlane).position;
                                    double dist = Distance(pointHere, pointThere);
                                    if (dist > ndist){
                                        ndist = dist;
                                    }
                                }
                            }
                        }

                        // region line
                        for (auto & rlch : r.topo.constraints<RegionLineConnectionData>()){
                            LineHandle lh = mg.topo(rlch).component<1>();
                            Line3 & thatLine = lines[lh];
                            for (auto & a : mg.data(rlch).normalizedAnchors){
                                Point3 pointHere = IntersectionOfLineAndPlane(Ray3(Point3(0, 0, 0), a), thisPlane).position;
                                Point3 pointThere = DistanceBetweenTwoLines(Ray3(Point3(0, 0, 0), a), thatLine.infiniteLine()).second.first;
                                double dist = Distance(pointHere, pointThere);
                                if (dist > ndist){
                                    ndist = dist;
                                }
                            }
                        }
                    }

                    std::sort(usableRegionHandles.begin(), usableRegionHandles.end(), [&maxDistanceToNearestComponents](RegionHandle l1, RegionHandle l2){
                        return maxDistanceToNearestComponents[l1] > maxDistanceToNearestComponents[l2];
                    });

                    // loose orientation constraints on these lines
                    for (int i = 0; i < usableRegionHandles.size() * regionsLoosableRatio; i++){
                        double maxCornerDist = maxDistanceToNearestComponents[usableRegionHandles[i]];
                        if (maxCornerDist < distThres)
                            continue;
                        auto & rp = props[usableRegionHandles[i]];
                        assert(rp.used);
                        rp.orientationClaz = rp.orientationNotClaz = -1;
                        rp.variables = { 1.0, 1.0, 1.0 };
                        //if (lp.orientationClaz == -1){
                        //    lp.used = false; // disable this line!!!!!!!!!
                        //}
                    }
                }
            }
            else {

                // collect keypoints
                struct ComponentKeyPoint {
                    Point3 point;
                    int handleId;
                    bool isOnRegion;
                    inline ComponentKeyPoint() {}
                    inline ComponentKeyPoint(const Point3 & p, RegionHandle rh) : point(p), isOnRegion(true), handleId(rh.id) {}
                    inline ComponentKeyPoint(const Point3 & p, LineHandle lh) : point(p), isOnRegion(false), handleId(lh.id) {}
                    inline bool matches(RegionHandle rh) const { return isOnRegion && handleId == rh.id; }
                    inline bool matches(LineHandle lh) const { return !isOnRegion && handleId == lh.id; }
                };
                auto getBB = [distThres](const ComponentKeyPoint & p) {return BoundingBox(p.point).expand(distThres); };
                RTreeWrapper<ComponentKeyPoint, decltype(getBB)> keyPointsRTree(getBB);

                // collect current instances
                HandledTable<LineHandle, Line3> lines(mg.internalComponents<LineData>().size());
                HandledTable<RegionHandle, Plane3> planes(mg.internalComponents<RegionData>().size());

                // keypoints on regions
                for (auto & r : mg.components<RegionData>()){
                    if (!props[r.topo.hd].used)
                        continue;
                    auto plane = Instance(mg, props, r.topo.hd);
                    planes[r.topo.hd] = plane;
                    for (auto & cs : r.data.normalizedContours){
                        for (auto & c : cs){
                            double depth = DepthAt(c, plane);
                            keyPointsRTree.insert(ComponentKeyPoint(c * depth, r.topo.hd));
                        }
                    }
                }
                // keypoints on lines
                for (auto & l : mg.components<LineData>()){
                    if (!props[l.topo.hd].used)
                        continue;
                    auto line = Instance(mg, props, l.topo.hd);
                    lines[l.topo.hd] = line;
                    keyPointsRTree.insert(ComponentKeyPoint(line.first, l.topo.hd));
                    keyPointsRTree.insert(ComponentKeyPoint(line.second, l.topo.hd));
                }

                // find most isolated lines
                if (linesLoosableRatio > 0.0){

                    auto lineMaxCornerNNDistances = mg.createComponentTable<LineData, double>();
                    std::vector<LineHandle> usableLineHandles;
                    for (auto & l : mg.components<LineData>()){
                        if (!props[l.topo.hd].used)
                            continue;
                        usableLineHandles.push_back(l.topo.hd);
                        auto & line = lines[l.topo.hd];
                        double nnDistances[] = { std::numeric_limits<double>::max(), std::numeric_limits<double>::max() };
                        keyPointsRTree.search(BoundingBox(line.first).expand(nnSearchRange), [&l, &line, &nnDistances](const ComponentKeyPoint & p){
                            if (p.matches(l.topo.hd))
                                return true;
                            if (!p.isOnRegion)
                                return true;
                            nnDistances[0] = std::min(Distance(line.first, p.point), nnDistances[0]);
                            //nnDistances[1] = std::min(Distance(line.second, p.point), nnDistances[1]);
                            return true;
                        });
                        keyPointsRTree.search(BoundingBox(line.second).expand(nnSearchRange), [&l, &line, &nnDistances](const ComponentKeyPoint & p){
                            if (p.matches(l.topo.hd))
                                return true;
                            //nnDistances[0] = std::min(Distance(line.first, p.point), nnDistances[0]);
                            nnDistances[1] = std::min(Distance(line.second, p.point), nnDistances[1]);
                            return true;
                        });
                        lineMaxCornerNNDistances[l.topo.hd] = std::max(nnDistances[0], nnDistances[1]);
                    }

                    std::sort(usableLineHandles.begin(), usableLineHandles.end(), [&lineMaxCornerNNDistances](LineHandle l1, LineHandle l2){
                        return lineMaxCornerNNDistances[l1] > lineMaxCornerNNDistances[l2];
                    });

                    // loose orientation constraints on these lines
                    for (int i = 0; i < usableLineHandles.size() * linesLoosableRatio; i++){
                        double maxCornerDist = lineMaxCornerNNDistances[usableLineHandles[i]];
                        if (maxCornerDist < distThres)
                            continue;
                        auto & lp = props[usableLineHandles[i]];
                        assert(lp.used);
                        if (lp.orientationClaz >= 0){
                            lp.orientationClaz = -1;
                            lp.variables = { 1.0, 1.0 };
                        }
                        assert(lp.orientationNotClaz == -1);
                        if (lp.orientationClaz == -1){
                            lp.used = false; // disable this line!!!!!!!!!
                            lp.variables.clear();
                        }
                    }
                }

                // find most isolated regions
                if (regionsLoosableRatio > 0.0){

                    auto regionMaxCornerNNDistances = mg.createComponentTable<RegionData, double>();
                    std::vector<RegionHandle> usableRegionHandles;
                    for (auto & r : mg.components<RegionData>()){
                        if (!props[r.topo.hd].used)
                            continue;
                        usableRegionHandles.push_back(r.topo.hd);
                        auto & plane = planes[r.topo.hd];
                        double maxCornerNNDistance = 0.0;
                        for (int i = 0; i < r.data.normalizedContours.size(); i++){
                            auto & contour = r.data.normalizedContours[i];
                            for (int j = 0; j < contour.size(); j++){
                                auto c = contour[j];
                                double depth = DepthAt(c, plane);
                                c *= depth; // the corner position
                                auto id = std::make_pair(i, j); // the corner id
                                double cornerNNDistance = std::numeric_limits<double>::max(); // find nearest neighbor distance to this corner
                                keyPointsRTree.search(BoundingBox(c).expand(nnSearchRange), [&r, &plane, &cornerNNDistance, id, &c](const ComponentKeyPoint & p){
                                    if (p.matches(r.topo.hd))
                                        return true;
                                    double d = Distance(c, p.point);
                                    cornerNNDistance = std::min(cornerNNDistance, d);
                                    return true;
                                });
                                if (cornerNNDistance > maxCornerNNDistance){
                                    maxCornerNNDistance = cornerNNDistance; // get the maximum nn distance of all corners on this region
                                }
                            }
                        }
                        regionMaxCornerNNDistances[r.topo.hd] = maxCornerNNDistance;
                    }

                    std::sort(usableRegionHandles.begin(), usableRegionHandles.end(), [&regionMaxCornerNNDistances](RegionHandle l1, RegionHandle l2){
                        return regionMaxCornerNNDistances[l1] > regionMaxCornerNNDistances[l2];
                    });

                    // loose orientation constraints on these lines
                    for (int i = 0; i < usableRegionHandles.size() * regionsLoosableRatio; i++){
                        double maxCornerDist = regionMaxCornerNNDistances[usableRegionHandles[i]];
                        if (maxCornerDist < distThres)
                            continue;
                        auto & rp = props[usableRegionHandles[i]];
                        assert(rp.used);
                        rp.orientationClaz = rp.orientationNotClaz = -1;
                        rp.variables = { 1.0, 1.0, 1.0 };
                        //if (lp.orientationClaz == -1){
                        //    lp.used = false; // disable this line!!!!!!!!!
                        //}
                    }
                }
            }

            ForeachMixedGraphConstraintHandle(mg, ConstraintUsableFlagUpdater{ mg, props });

        }









        void LooseMaybeOcclusionBoundaryConstraints(const MixedGraph & mg, MixedGraphPropertyTable & props){

            NOT_IMPLEMENTED_YET();

        }









        namespace {

            template <class UhClickHandlerFunT, class UhColorizerFunT = core::ConstantFunctor<gui::Color>>
            void ManuallyOptimizeMixedGraph(const core::Image & panorama,
                const MixedGraph & mg,
                MixedGraphPropertyTable & props,
                UhClickHandlerFunT && uhClicked,
                UhColorizerFunT && uhColorizer = UhColorizerFunT(gui::ColorTag::White),
                bool optimizeInEachIteration = false) {

                bool modified = true;

                struct ComponentID {
                    int handleID;
                    bool isRegion;
                };

                auto sppCallbackFun = [&props, &mg, &modified, &uhClicked](gui::InteractionID iid,
                    const std::pair<ComponentID, gui::Colored<gui::SpatialProjectedPolygon>> & spp) {
                    std::cout << (spp.first.isRegion ? "Region" : "Line") << spp.first.handleID << std::endl;
                    if (spp.first.isRegion){
                        modified = uhClicked(RegionHandle(spp.first.handleID));
                    }
                    else {
                        modified = uhClicked(LineHandle(spp.first.handleID));
                    }
                };

                while (modified){

                    gui::ResourceStore::set("texture", panorama);

                    modified = false;
                    if (optimizeInEachIteration){
                        SolveVariablesUsingInversedDepths(mg, props);
                    }

                    gui::Visualizer viz("mixed graph optimizable");
                    viz.renderOptions.bwColor = 1.0;
                    viz.renderOptions.bwTexColor = 0.0;
                    viz.installingOptions.discretizeOptions.colorTable = gui::ColorTableDescriptor::RGB;
                    std::vector<std::pair<ComponentID, gui::Colored<gui::SpatialProjectedPolygon>>> spps;
                    std::vector<gui::Colored<core::Line3>> lines;

                    for (auto & c : mg.components<RegionData>()){
                        if (!props[c.topo.hd].used)
                            continue;
                        auto uh = c.topo.hd;
                        auto & region = c.data;
                        gui::SpatialProjectedPolygon spp;
                        // filter corners
                        core::ForeachCompatibleWithLastElement(c.data.normalizedContours.front().begin(), c.data.normalizedContours.front().end(),
                            std::back_inserter(spp.corners),
                            [](const core::Vec3 & a, const core::Vec3 & b) -> bool {
                            return core::AngleBetweenDirections(a, b) > M_PI / 1000.0;
                        });
                        if (spp.corners.size() < 3)
                            continue;

                        spp.projectionCenter = core::Point3(0, 0, 0);
                        spp.plane = Instance(mg, props, uh);
                        assert(!HasValue(spp.plane, IsInfOrNaN<double>));
                        spps.emplace_back(ComponentID{ uh.id, true }, std::move(gui::ColorAs(spp, uhColorizer(uh))));
                    }

                    for (auto & c : mg.components<LineData>()){
                        if (!props[c.topo.hd].used)
                            continue;
                        auto uh = c.topo.hd;
                        auto & line = c.data;
                        lines.push_back(gui::ColorAs(Instance(mg, props, uh), uhColorizer(uh)));
                    }

                    viz.begin(spps, sppCallbackFun).shaderSource(gui::OpenGLShaderSourceDescriptor::XPanorama).resource("texture").end();
                    viz.installingOptions.discretizeOptions.color = gui::ColorTag::DarkGray;
                    viz.installingOptions.lineWidth = 5.0;
                    viz.add(lines);

                    std::vector<core::Line3> connectionLines;
                    for (auto & c : mg.constraints<RegionBoundaryData>()){
                        if (!props[c.topo.hd].used)
                            continue;
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
                        if (!props[c.topo.hd].used)
                            continue;
                        auto inst1 = Instance(mg, props, c.topo.component<0>());
                        auto inst2 = Instance(mg, props, c.topo.component<1>());
                        for (auto & s : c.data.normalizedAnchors){
                            double d1 = DepthAt(s, inst1);
                            double d2 = DepthAt(s, inst2);
                            connectionLines.emplace_back(normalize(s) * d1, normalize(s) * d2);
                        }
                    }

                    viz.installingOptions.discretizeOptions.color = gui::ColorTag::Black;
                    viz.installingOptions.lineWidth = 1.0;
                    viz.add(connectionLines);

                    viz.renderOptions.renderMode = gui::RenderModeFlag::Triangles | gui::RenderModeFlag::Lines;
                    viz.renderOptions.backgroundColor = gui::ColorTag::White;
                    viz.renderOptions.bwColor = 0.5;
                    viz.renderOptions.bwTexColor = 0.5;
                    viz.camera(core::PerspectiveCamera(1000, 800, 800, Point3(-1, 1, 1), Point3(0, 0, 0)));
                    viz.show(true, false);

                    gui::ResourceStore::clear();
                }

            }


        }


        struct ComponentHandleColorizer {
            const MixedGraph & mg;
            const MixedGraphPropertyTable & props;
            gui::ColorTable colorTableForVPs;
            gui::ColorTable colorTableForRegionAlongOrientations;
            ComponentHandleColorizer(const MixedGraph & g, const MixedGraphPropertyTable & p) : mg(g), props(p) {
                colorTableForVPs = gui::ColorTable(gui::ColorTableDescriptor::RGB)
                    .appendRandomizedColors(props.vanishingPoints.size()-3);
                colorTableForRegionAlongOrientations = gui::CreateGreyColorTableWithSize(props.vanishingPoints.size());
            }

            inline gui::Color operator()(ComponentHandle<RegionData> rh) const{
                /*double maxPlaneDist = 0.0;
                for (auto & r : mg.components<RegionData>()){
                    if (!props[r.topo.hd].used)
                        continue;
                    double pdist = Instance(mg, props, r.topo.hd).distanceTo(Point3(0, 0, 0));
                    if (pdist > maxPlaneDist)
                        maxPlaneDist = pdist;
                }
                return gui::ColorFromHSV(Instance(mg, props, rh).distanceTo(Point3(0, 0, 0)) / maxPlaneDist, 0.5, 0.8);*/
                auto & prop = props[rh];
                if (!prop.used)
                    return gui::ColorTag::Yellow;
                if (prop.orientationClaz >= 0)
                    return colorTableForVPs[prop.orientationClaz];
                if (prop.orientationNotClaz >= 0)
                    return colorTableForRegionAlongOrientations[prop.orientationNotClaz];
                return gui::ColorTag::Yellow;
            }
            inline gui::Color operator()(ComponentHandle<LineData> rh) const {
                return colorTableForVPs[props[rh].orientationClaz];
            }
        };

        struct ComponentHandleClickHandler {
            const MixedGraph & mg;
            const MixedGraphPropertyTable & props;
            inline bool operator()(ComponentHandle<RegionData> rh) const{
                std::cout << "area: " << mg.data(rh).area
                    << " connected regions: " << mg.topo(rh).constraints<RegionBoundaryData>().size()
                    << " connected lines: " << mg.topo(rh).constraints<RegionLineConnectionData>().size() << std::endl;
                std::cout << "==connected regions info==" << std::endl;
                for (auto & h : mg.topo(rh).constraints<RegionBoundaryData>()){
                    std::cout << "  sample points num: " << ElementsNum(mg.data(h).normalizedSampledPoints) << std::endl;
                }
                std::cout << "==connected lines info==" << std::endl;
                for (auto & h : mg.topo(rh).constraints<RegionLineConnectionData>()){
                    std::cout << "  anchors num: " << ElementsNum(mg.data(h).normalizedAnchors) << std::endl;
                }
                return false;
            }
            inline bool operator()(ComponentHandle<LineData> rh) const {               
                return false;
            }
        };


        void Visualize(const View<PanoramicCamera> & texture,
            const MixedGraph & mg, MixedGraphPropertyTable & props){
            ManuallyOptimizeMixedGraph(texture.image, mg, props, 
                ComponentHandleClickHandler{ mg, props }, 
                ComponentHandleColorizer{ mg, props });
        }
        
        void Visualize(const View<PerspectiveCamera> & texture,
            const MixedGraph & mg, MixedGraphPropertyTable & props) {
            ManuallyOptimizeMixedGraph(texture.sampled(PanoramicCamera(250, Point3(0, 0, 0), Point3(1, 0, 0), Vec3(0, 0, 1))).image, 
                mg, props, ComponentHandleClickHandler{ mg, props },
                ComponentHandleColorizer{ mg, props });
        }



        namespace {

            // graph data
            class Universe {
            public:
                struct Element {
                    int rank;
                    int p; // parent
                    float size;
                };
                explicit Universe(const std::vector<float> & eleSizes)
                    : elements(eleSizes.size()), num(eleSizes.size()) {
                    for (int i = 0; i < eleSizes.size(); i++){
                        elements[i].rank = 0;
                        elements[i].size = eleSizes[i];
                        elements[i].p = i;
                    }
                }
                int find(int x) {
                    int y = x;
                    while (y != elements[y].p)
                        y = elements[y].p;
                    elements[x].p = y;
                    return y;
                }
                void join(int x, int y) {
                    if (elements[x].rank > elements[y].rank) {
                        elements[y].p = x;
                        elements[x].size += elements[y].size;
                    }
                    else {
                        elements[x].p = y;
                        elements[y].size += elements[x].size;
                        if (elements[x].rank == elements[y].rank)
                            elements[y].rank++;
                    }
                    num--;
                }
                float size(int x) const { return elements[x].size; }
                int numSets() const { return num; }
            private:
                int num;
                std::vector<Element> elements;
            };

        }



        HandledTable<RegionHandle, int> ClusterRegions(const MixedGraph & mg, const MixedGraphPropertyTable & props){

            const double medianCenterDepth = ComponentMedianCenterDepth(mg, props);

            auto planes = mg.createComponentTable<RegionData, Plane3>();
            auto regionAreas = mg.createComponentTable<RegionData>(0.0f);
            for (auto & r : mg.components<RegionData>()){
                planes[r.topo.hd] = Instance(mg, props, r.topo.hd);
                regionAreas[r.topo.hd] = r.data.area;
            }

            // graph cut
            std::vector<Scored<RegionBoundaryHandle>> bhs;
            bhs.reserve(mg.internalConstraints<RegionBoundaryData>().size());
            for (auto & b : mg.constraints<RegionBoundaryData>()){
                if (!props[b.topo.hd].used){
                    continue;
                }
                auto rh1 = mg.topo(b.topo.hd).component<0>();
                auto rh2 = mg.topo(b.topo.hd).component<1>();
                assert(props[rh1].used && props[rh2].used);
                
                auto & plane1 = planes[rh1];
                auto & plane2 = planes[rh2];
                double angle = AngleBetweenUndirectedVectors(plane1.normal, plane2.normal);
                double distance = 0.0;
                for (auto & ps : b.data.normalizedSampledPoints){
                    for (auto & p : ps){
                        double d1 = DepthAt(p, plane1);
                        double d2 = DepthAt(p, plane2);
                        distance += abs(d1 - d2);
                    }
                }
                distance /= ElementsNum(b.data.normalizedSampledPoints);
                double weight = angle + distance / medianCenterDepth / 2.0;
                bhs.push_back(ScoreAs(b.topo.hd, weight));
            }

            std::sort(bhs.begin(), bhs.end());


            // segmentation
            //auto regionOnes = mg.createComponentTable<RegionData>(1.0f);
            Universe u(regionAreas.data);
            static const float c = 30.f;
            std::vector<float> threshold(regionAreas.data.size(), c);

            for (auto & bhs : bhs){
                auto bh = bhs.component;
                double weight = bhs.score;
                if (!props[bh].used){
                    continue;
                }

                // components conected by this edge
                auto rh1 = mg.topo(bh).component<0>();
                auto rh2 = mg.topo(bh).component<1>();
                assert(props[rh1].used && props[rh2].used);

                int a = u.find(rh1.id);
                int b = u.find(rh2.id);
                if (a != b){
                    if (weight <= threshold[a] && weight <= threshold[b]){
                        u.join(a, b);
                        a = u.find(a);
                        threshold[a] = weight + c / u.size(a);
                    }
                }
            }

            // merge too small regions?
            double minSize = 0;
            if (minSize > 0){
                for (auto & bhs : bhs) {
                    auto bh = bhs.component;
                    if (!props[bh].used){
                        continue;
                    }
                    auto rh1 = mg.topo(bh).component<0>();
                    auto rh2 = mg.topo(bh).component<1>();
                    assert(props[rh1].used && props[rh2].used);

                    int a = u.find(rh1.id);
                    int b = u.find(rh2.id);
                    if ((a != b) && ((u.size(a) < minSize) || (u.size(b) < minSize)))
                        u.join(a, b);
                }
            }

            // classify regions
            auto planeIds = mg.createComponentTable<RegionData>(-1);
            for (auto & r : mg.components<RegionData>()){
                planeIds[r.topo.hd] = props[r.topo.hd].used ? u.find(r.topo.hd.id) : -1;
            }

            return planeIds;
        }


    }
}


