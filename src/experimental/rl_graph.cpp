
extern "C" {
    #include <gpc.h>
//    #include <mosek.h>
}

//
#include <GCoptimization.h>


#include "../misc/matlab.hpp"
#include "../misc/eigen.hpp"

//
#include "../core/algorithms.hpp"
#include "../core/containers.hpp"
#include "../core/utilities.hpp"
#include "../core/clock.hpp"
#include "../core/homo_graph.hpp"

#include "tools.hpp"
#include "rl_graph.hpp"

#include "../gui/visualizers.hpp"


namespace panoramix {
    namespace experimental {

        template <class FunT>
        inline void ForeachRLGraphComponentHandle(const RLGraph & mg, FunT && fun){
            for (auto & c : mg.components<LineData>()){
                fun(c.topo.hd);
            }
            for (auto & c : mg.components<RegionData>()){
                fun(c.topo.hd);
            }
        }

        template <class FunT>
        inline void ForeachRLGraphConstraintHandle(const RLGraph & mg, FunT && fun){
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

            template <class T>
            inline double ChainLength(const std::vector<T> & v){
                double len = 0;
                for (int i = 1; i < v.size(); i++){
                    len += Distance(v[i - 1], v[i]);
                }
                return len;
            }

            template <class T>
            inline double ChainLength(const std::vector<std::vector<T>> & v){
                double len = 0;
                for (const auto & vv : v){
                    for (int i = 1; i < vv.size(); i++){
                        len += Distance(vv[i - 1], vv[i]);
                    }
                }
                return len;
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

        std::vector<Vec3> EstimateVanishingPointsAndClassifyLines(const PerspectiveCamera & cam,
            std::vector<Classified<Line2>> & lineSegments){
            std::vector<Vec3> lineIntersections;

            int linesNum = 0;
            std::vector<Line2> pureLines(lineSegments.size());
            linesNum += lineSegments.size();
            for (int k = 0; k < pureLines.size(); k++){
                pureLines[k] = lineSegments[k].component;
            }
            auto inters = ComputeLineIntersections(pureLines, nullptr, true, std::numeric_limits<double>::max());
            // insert line intersections
            for (auto & p : inters){
                lineIntersections.push_back(normalize(cam.spatialDirection(p.value())));
            }

            auto vanishingPoints = FindOrthogonalPrinicipleDirections(lineIntersections, 1000, 500, true).unwrap();

            // project lines to space
            std::vector<Classified<Line3>> spatialLineSegments;
            spatialLineSegments.reserve(linesNum);
            for (const auto & line : lineSegments) {
                auto & p1 = line.component.first;
                auto & p2 = line.component.second;
                auto pp1 = cam.spatialDirection(p1);
                auto pp2 = cam.spatialDirection(p2);
                Classified<Line3> cline3;
                cline3.claz = -1;
                cline3.component = Line3{ pp1, pp2 };
                spatialLineSegments.push_back(cline3);
            }

            // classify lines
            ClassifyLines(spatialLineSegments, vanishingPoints, M_PI / 3.0, 0.1, 0.8, M_PI / 18.0);

            int ii = 0;
            for (int j = 0; j < lineSegments.size(); j++){
                lineSegments[j].claz = spatialLineSegments[ii].claz;
                ii++;
            }

            return vanishingPoints;
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
                bool computeJunctionWeights = true){

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


        void AppendLines(RLGraph & mg, const std::vector<Classified<Line2>> & lineSegments,
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
                cam.focal() * incidenceAngleVerticalDirectionThreshold, false, true);

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

            std::vector<RegionHandle> CollectRegionsFromSegmentation(RLGraph & mg, const Imagei & segmentedRegions, const PerspectiveCamera & cam){
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

            std::vector<RegionHandle> CollectRegionsFromSegmentation(RLGraph & mg, const Imagei & segmentedRegions, const PanoramicCamera & cam){
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


            template <class CameraT, bool SampleLineOnImage>
            std::vector<RegionHandle> AppendRegionsTemplate(RLGraph & mg, const Imagei & segmentedRegions, const CameraT & cam,
                double samplingStepAngleOnBoundary, double samplingStepAngleOnLine, 
                int samplerSizeOnBoundary, int samplerSizeOnLine, bool noBoundaryUnderLines){

                auto regionHandles = CollectRegionsFromSegmentation(mg, segmentedRegions, cam);
                int regionNum = regionHandles.size();


                // add region-line connections
                std::map<std::pair<RegionHandle, LineHandle>, RegionLineConnectionData> regionLineConnections;

                for (auto & ld : mg.components<LineData>()){
                    auto & line = ld.data.line;
                    Line2 line2(cam.screenProjection(line.first), cam.screenProjection(line.second));
                    Vec2 vertToLineDir = normalize(PerpendicularDirection(line2.direction())); // vertical to this line
                    double stepOnImage = std::max(cam.focal() * samplingStepAngleOnLine, 0.5);
                    if (SampleLineOnImage){
                        for (double stepLen = 0.0; stepLen <= line2.length(); stepLen += stepOnImage){
                            Point2 sampleP = line2.first + normalize(line2.second - line2.first) * stepLen;
                            PixelLoc originalP = ToPixelLoc(sampleP);
                            std::set<int> connectedRegionIds;
                            for (int x = originalP.x - samplerSizeOnLine; x <= originalP.x + samplerSizeOnLine; x++){
                                for (int y = originalP.y - samplerSizeOnLine; y <= originalP.y + samplerSizeOnLine; y++){
                                    PixelLoc p(x, y);
                                    if (!Contains(segmentedRegions, p))
                                        continue;
                                    int regionId = segmentedRegions(p);
                                    auto rh = regionHandles[regionId];
                                    if (rh.invalid())
                                        continue;
                                    connectedRegionIds.insert(regionId);
                                }
                            }
                            // collect abit farther neighbors for judging whether it can be detachable!
                            std::set<int> abitFarLeftRightRegionIds[2];
                            for (int d = std::min(samplerSizeOnLine, 2); d <= samplerSizeOnLine * 2; d++){
                                for (int dir : {-d, d}){
                                    PixelLoc p = ToPixelLoc(sampleP + vertToLineDir * dir);
                                    if (!Contains(segmentedRegions, p))
                                        continue;
                                    int regionId = segmentedRegions(p);
                                    auto rh = regionHandles[regionId];
                                    if (rh.invalid())
                                        continue;
                                    (dir < 0 ? abitFarLeftRightRegionIds[0] : abitFarLeftRightRegionIds[1]).insert(regionId);
                                }
                            }
                            for (int regionId : connectedRegionIds){
                                RegionLineConnectionData & rd = regionLineConnections[std::make_pair(regionHandles[regionId], ld.topo.hd)];
                                rd.normalizedAnchors.push_back(normalize(cam.spatialDirection(sampleP)));
                                rd.detachable = !(Contains(abitFarLeftRightRegionIds[0], regionId) && Contains(abitFarLeftRightRegionIds[1], regionId));
                            }
                        }
                    }else{
                        double spanAngle = AngleBetweenDirections(line.first, line.second);
                        int stepNum = static_cast<int>(std::ceil(spanAngle / samplingStepAngleOnLine));
                        assert(stepNum >= 1);
                        for (int step = 0; step <= stepNum; step++){
                            double angle = step * samplingStepAngleOnLine;
                            Vec3 sample = RotateDirection(line.first, line.second, angle);
                            if (!cam.isVisibleOnScreen(sample))
                                continue;
                            Point2 sampleP = cam.screenProjection(sample);
                            PixelLoc originalP = ToPixelLoc(sampleP);
                            // collect neighbors
                            std::set<int> connectedRegionIds;
                            for (int x = originalP.x - samplerSizeOnLine; x <= originalP.x + samplerSizeOnLine; x++){
                                for (int y = originalP.y - samplerSizeOnLine; y <= originalP.y + samplerSizeOnLine; y++){
                                    PixelLoc p(x, y);
                                    if (!Contains(segmentedRegions, p))
                                        continue;
                                    int regionId = segmentedRegions(p);
                                    auto rh = regionHandles[regionId];
                                    if (rh.invalid())
                                        continue;
                                    connectedRegionIds.insert(regionId);
                                }
                            }
                            // collect abit farther neighbors for judging whether it can be detachable!
                            std::set<int> abitFarLeftRightRegionIds[2];
                            for (int d = std::min(samplerSizeOnLine, 2); d <= samplerSizeOnLine * 2; d++){
                                for (int dir : {-d, d}){
                                    PixelLoc p = ToPixelLoc(sampleP + vertToLineDir * dir);
                                    if (!Contains(segmentedRegions, p))
                                        continue;
                                    int regionId = segmentedRegions(p);
                                    auto rh = regionHandles[regionId];
                                    if (rh.invalid())
                                        continue;
                                    (dir < 0 ? abitFarLeftRightRegionIds[0] : abitFarLeftRightRegionIds[1]).insert(regionId);
                                }
                            }
                            for (int regionId : connectedRegionIds){
                                RegionLineConnectionData & rd = regionLineConnections[std::make_pair(regionHandles[regionId], ld.topo.hd)];
                                rd.normalizedAnchors.push_back(normalize(sample));
                                rd.detachable = !(Contains(abitFarLeftRightRegionIds[0], regionId) && Contains(abitFarLeftRightRegionIds[1], regionId));
                            }
                        }
                    }
                }                

                for (auto & rlc : regionLineConnections){
                    // compute length;
                    rlc.second.length = AngleBetweenDirections(rlc.second.normalizedAnchors.front(), rlc.second.normalizedAnchors.back());
                    mg.addConstraint(std::move(rlc.second), rlc.first.first, rlc.first.second);
                }


                // add region boundary constraints
                std::map<std::pair<int, int>, std::vector<std::vector<PixelLoc>>> boundaryEdges =
                    FindContoursOfRegionsAndBoundaries(segmentedRegions, samplerSizeOnBoundary, false);

                for (auto & bep : boundaryEdges) {
                    auto & rids = bep.first;
                    auto edges = bep.second;

                    if (regionHandles[rids.first].invalid() || regionHandles[rids.second].invalid())
                        continue;

                    auto rh = regionHandles[rids.first];
                    if (noBoundaryUnderLines){
                        std::vector<std::vector<PixelLoc>> filteredBoundaryPixels;
                        for (auto & e : edges){
                            std::vector<PixelLoc> filteredEdge;
                            for (PixelLoc p : e){
                                bool coveredByLine = false;
                                for (RegionLineConnectionHandle rlcon : mg.topo(rh).constraints<RegionLineConnectionData>()){
                                    LineHandle lh = mg.topo(rlcon).component<1>();
                                    const Line3 & line = mg.data(lh).line;
                                    Line2 line2(cam.screenProjection(line.first), cam.screenProjection(line.second));
                                    double d = DistanceFromPointToLine(vec_cast<double>(p), line2).first;
                                    if (d < samplerSizeOnLine || d < samplerSizeOnBoundary){
                                        coveredByLine = true;
                                        break;
                                    }
                                }
                                if (!coveredByLine){
                                    filteredEdge.push_back(p);
                                }
                            }
                            if (filteredEdge.size() > 2){
                                filteredBoundaryPixels.push_back(std::move(filteredEdge));
                            }
                        }
                        if (filteredBoundaryPixels.size() > 0){
                            edges = std::move(filteredBoundaryPixels);
                        }
                        else{
                            continue;
                        }
                    }

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
                        auto & edge = bd.normalizedEdges[k];
                        if (edge.empty()){
                            continue;
                        }
                        std::vector<Vec3> points = { edge.front() };
                        for (auto & edgeP : edge){
                            double remainedAngle = AngleBetweenDirections(points.back(), edgeP);
                            while (remainedAngle >= samplingStepAngleOnBoundary){
                                points.push_back(normalize(RotateDirection(points.back(), edgeP, samplingStepAngleOnBoundary)));
                                remainedAngle -= samplingStepAngleOnBoundary;
                            }
                        }
                        bd.normalizedSampledPoints[k] = std::move(points);
                    }

                    mg.addConstraint(std::move(bd), regionHandles[rids.first], regionHandles[rids.second]);
                }


                return regionHandles;

            }


        }

        std::vector<RegionHandle> AppendRegions(RLGraph & mg, const Imagei & segmentedRegions, const PerspectiveCamera & cam,
            double samplingStepAngleOnBoundary, double samplingStepAngleOnLine, 
            int samplerSizeOnBoundary, int samplerSizeOnLine, bool noBoundaryUnderLines){
            return AppendRegionsTemplate<PerspectiveCamera, true>(mg, segmentedRegions, cam, samplingStepAngleOnBoundary, samplingStepAngleOnLine, 
                samplerSizeOnBoundary, samplerSizeOnLine, noBoundaryUnderLines);
        }
        std::vector<RegionHandle> AppendRegions(RLGraph & mg, const Imagei & segmentedRegions, const PanoramicCamera & cam,
            double samplingStepAngleOnBoundary, double samplingStepAngleOnLine,
            int samplerSizeOnBoundary, int samplerSizeOnLine, bool noBoundaryUnderLines){
            return AppendRegionsTemplate<PanoramicCamera, false>(mg, segmentedRegions, cam, samplingStepAngleOnBoundary, samplingStepAngleOnLine, 
                samplerSizeOnBoundary, samplerSizeOnLine, noBoundaryUnderLines);
        }




        View<PartialPanoramicCamera, Imageub> PerfectRegionMaskView(const RLGraph & mg, RegionHandle rh, double focal){
            auto h = rh;
            auto & contours = mg.data(rh).normalizedContours;
            double radiusAngle = 0.0;
            for (auto & cs : contours){
                for (auto & c : cs){
                    double angle = AngleBetweenDirections(mg.data(rh).normalizedCenter, c);
                    if (angle > radiusAngle){
                        radiusAngle = angle;
                    }
                }
            }
            int ppcSize = std::ceil(2 * radiusAngle * focal);
            Vec3 x;
            std::tie(x, std::ignore) = ProposeXYDirectionsFromZDirection(mg.data(rh).normalizedCenter);
            PartialPanoramicCamera ppc(ppcSize, ppcSize, focal, Point3(0, 0, 0), mg.data(rh).normalizedCenter, x);
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
            return View<PartialPanoramicCamera, Imageub>{mask, ppc};
        }





        int ConnectedComponents(const RLGraph & mg, RLGraphComponentTable<int> & ccIds) {
            std::vector<RLGraph> ccs;

            struct HComp {
                enum Type { Region, Line };
                Type type;
                int id;
            };
            struct HCons {
                enum Type { RegionBoundary, LineRelation, RegionLine };
                Type type;
                int id;
            };

            std::vector<HComp> comps;
            std::vector<HCons> conss;
            comps.reserve(mg.allComponentsNum());
            conss.reserve(mg.allConstraintsNum());

            auto comp2htable = MakeHandledTableForAllComponents<size_t>(mg);
            auto cons2htable = MakeHandledTableForAllConstraints<size_t>(mg);

            for (auto & r : mg.components<RegionData>()){
                comp2htable[r.topo.hd] = comps.size();
                comps.push_back({ HComp::Region, r.topo.hd.id });
            }
            for (auto & l : mg.components<LineData>()){
                comp2htable[l.topo.hd] = comps.size();
                comps.push_back({ HComp::Line, l.topo.hd.id });
            }

            for (auto & c : mg.constraints<RegionBoundaryData>()){
                cons2htable[c.topo.hd] = conss.size();
                conss.push_back({ HCons::RegionBoundary, c.topo.hd.id });
            }
            for (auto & c : mg.constraints<LineRelationData>()){
                cons2htable[c.topo.hd] = conss.size();
                conss.push_back({ HCons::LineRelation, c.topo.hd.id });
            }
            for (auto & c : mg.constraints<RegionLineConnectionData>()){
                cons2htable[c.topo.hd] = conss.size();
                conss.push_back({ HCons::RegionLine, c.topo.hd.id });
            }

            std::vector<size_t> compIds(comps.size());
            std::iota(compIds.begin(), compIds.end(), 0ull);

            auto getNeighborComps = [&comps, &mg, &comp2htable](size_t compId){
                auto & comp = comps[compId];
                std::vector<size_t> neighbors;
                if (comp.type == HComp::Region){
                    RegionHandle rh(comp.id);
                    for (auto & h : mg.topo(rh).constraints<RegionBoundaryData>()){
                        auto another = mg.topo(h).component<0>();
                        if (another == rh){
                            another = mg.topo(h).component<1>();
                        }
                        neighbors.push_back(comp2htable[another]);
                    }
                    for (auto & h : mg.topo(rh).constraints<RegionLineConnectionData>()){
                        LineHandle another = mg.topo(h).component<1>();
                        neighbors.push_back(comp2htable[another]);
                    }
                }
                else{
                    LineHandle lh(comp.id);
                    for (auto & h : mg.topo(lh).constraints<LineRelationData>()){
                        auto another = mg.topo(h).component<0>();
                        if (another == lh){
                            another = mg.topo(h).component<1>();
                        }
                        neighbors.push_back(comp2htable[another]);
                    }
                    for (auto & h : mg.topo(lh).constraints<RegionLineConnectionData>()){
                        RegionHandle another = mg.topo(h).component<0>();
                        neighbors.push_back(comp2htable[another]);
                    }
                }
                return neighbors;
            };

            return core::ConnectedComponents(compIds.begin(), compIds.end(), std::move(getNeighborComps),
                [&ccIds, &comps](size_t v, int ccid){
                const HComp & hcomp = comps[v];
                if (hcomp.type == HComp::Line){
                    ccIds[LineHandle(hcomp.id)] = ccid;
                }
                else{
                    ccIds[RegionHandle(hcomp.id)] = ccid;
                }
            });
        }

        int ConnectedComponents(const RLGraph & mg, const RLGraphControls & controls, 
            RLGraphComponentTable<int> & ccIds,
            const std::function<bool(const RLGraphConstraintControl &)> & constraintAsConnected) {
            std::vector<RLGraph> ccs;


            struct HComp {
                enum Type { Region, Line };
                Type type;
                int id;
            };
            struct HCons {
                enum Type { RegionBoundary, LineRelation, RegionLine };
                Type type;
                int id;
            };

            std::vector<HComp> comps;
            std::vector<HCons> conss;
            comps.reserve(mg.allComponentsNum());
            conss.reserve(mg.allConstraintsNum());

            auto comp2htable = MakeHandledTableForAllComponents<size_t>(mg);
            auto cons2htable = MakeHandledTableForAllConstraints<size_t>(mg);

            for (auto & r : mg.components<RegionData>()){
                comp2htable[r.topo.hd] = comps.size();
                comps.push_back({ HComp::Region, r.topo.hd.id });
            }
            for (auto & l : mg.components<LineData>()){
                comp2htable[l.topo.hd] = comps.size();
                comps.push_back({ HComp::Line, l.topo.hd.id });
            }

            for (auto & c : mg.constraints<RegionBoundaryData>()){
                cons2htable[c.topo.hd] = conss.size();
                conss.push_back({ HCons::RegionBoundary, c.topo.hd.id });
            }
            for (auto & c : mg.constraints<LineRelationData>()){
                cons2htable[c.topo.hd] = conss.size();
                conss.push_back({ HCons::LineRelation, c.topo.hd.id });
            }
            for (auto & c : mg.constraints<RegionLineConnectionData>()){
                cons2htable[c.topo.hd] = conss.size();
                conss.push_back({ HCons::RegionLine, c.topo.hd.id });
            }

            std::vector<size_t> compIds(comps.size());
            std::iota(compIds.begin(), compIds.end(), 0ull);

            auto getNeighborComps = [&comps, &mg, &comp2htable, &constraintAsConnected, &controls](size_t compId){
                auto & comp = comps[compId];
                std::vector<size_t> neighbors;
                if (comp.type == HComp::Region){
                    RegionHandle rh(comp.id);
                    for (auto & h : mg.topo(rh).constraints<RegionBoundaryData>()){
                        auto & consControl = controls[h];
                        if (!constraintAsConnected){
                            if (!consControl.used || consControl.weight == 0.0)
                                continue;
                        }
                        else{
                            if (!constraintAsConnected(consControl))
                                continue;
                        }
                        auto another = mg.topo(h).component<0>();
                        if (another == rh){
                            another = mg.topo(h).component<1>();
                        }
                        neighbors.push_back(comp2htable[another]);
                    }
                    for (auto & h : mg.topo(rh).constraints<RegionLineConnectionData>()){
                        auto & consControl = controls[h];
                        if (!constraintAsConnected){
                            if (!consControl.used || consControl.weight == 0.0)
                                continue;
                        }
                        else{
                            if (!constraintAsConnected(consControl))
                                continue;
                        }
                        LineHandle another = mg.topo(h).component<1>();
                        neighbors.push_back(comp2htable[another]);
                    }
                }
                else{
                    LineHandle lh(comp.id);
                    for (auto & h : mg.topo(lh).constraints<LineRelationData>()){
                        auto & consControl = controls[h];
                        if (!constraintAsConnected){
                            if (!consControl.used || consControl.weight == 0.0)
                                continue;
                        }
                        else{
                            if (!constraintAsConnected(consControl))
                                continue;
                        }
                        auto another = mg.topo(h).component<0>();
                        if (another == lh){
                            another = mg.topo(h).component<1>();
                        }
                        neighbors.push_back(comp2htable[another]);
                    }
                    for (auto & h : mg.topo(lh).constraints<RegionLineConnectionData>()){
                        auto & consControl = controls[h];
                        if (!constraintAsConnected){
                            if (!consControl.used || consControl.weight == 0.0)
                                continue;
                        }
                        else{
                            if (!constraintAsConnected(consControl))
                                continue;
                        }
                        RegionHandle another = mg.topo(h).component<0>();
                        neighbors.push_back(comp2htable[another]);
                    }
                }
                return neighbors;
            };

            return core::ConnectedComponents(compIds.begin(), compIds.end(), std::move(getNeighborComps),
                [&ccIds, &comps](size_t v, int ccid){
                const HComp & hcomp = comps[v];
                if (hcomp.type == HComp::Line){
                    ccIds[LineHandle(hcomp.id)] = ccid;
                }
                else{
                    ccIds[RegionHandle(hcomp.id)] = ccid;
                }
            });
        }


        template <class T1, class T2>
        std::pair<T2, T1> Inversed(const std::pair<T1, T2> & p){
            return std::make_pair(p.second, p.first);
        }

        std::vector<RLGraph> Decompose(const RLGraph & mg, const RLGraphComponentTable<int> & ccids, int ccnum, 
            RLGraphOldToNew * old2new, RLGraphNewToOld * new2old){
            std::vector<RLGraph> ccs(ccnum);

            RLGraphComponentTable<int> newCompIds = MakeHandledTableForAllComponents<int>(mg);
            for (auto & c : mg.components<RegionData>()){
                int ccid = ccids[c.topo.hd];
                auto h = ccs[ccid].addComponent(c.data);
                newCompIds[c.topo.hd] = h.id;
                auto pair = std::make_pair(c.topo.hd, std::make_pair(ccid, h));
                if (old2new){
                    (*old2new)[pair.first] = pair.second;
                }
                if (new2old){
                    (*new2old)[pair.second] = pair.first;
                }
            }
            for (auto & c : mg.components<LineData>()){
                int ccid = ccids[c.topo.hd];
                auto h = ccs[ccid].addComponent(c.data);
                newCompIds[c.topo.hd] = h.id;
                auto pair = std::make_pair(c.topo.hd, std::make_pair(ccid, h));
                if (old2new){
                    (*old2new)[pair.first] = pair.second;
                }
                if (new2old){
                    (*new2old)[pair.second] = pair.first;
                }
            }

            for (auto & c : mg.constraints<RegionBoundaryData>()){
                int ccid1 = ccids[c.topo.component<0>()];
                int ccid2 = ccids[c.topo.component<1>()];
                if (ccid1 != ccid2){
                    continue;
                }
                int ccid = ccid1;
                auto h = ccs[ccid].addConstraint(c.data, 
                    RegionHandle(newCompIds[c.topo.component<0>()]),
                    RegionHandle(newCompIds[c.topo.component<1>()]));
                auto pair = std::make_pair(c.topo.hd, std::make_pair(ccid, h));
                if (old2new){
                    (*old2new)[pair.first] = pair.second;
                }
                if (new2old){
                    (*new2old)[pair.second] = pair.first;
                }
            }
            for (auto & c : mg.constraints<LineRelationData>()){
                int ccid1 = ccids[c.topo.component<0>()];
                int ccid2 = ccids[c.topo.component<1>()];
                if (ccid1 != ccid2){
                    continue;
                }
                int ccid = ccid1;
                auto h = ccs[ccid].addConstraint(c.data,
                    LineHandle(newCompIds[c.topo.component<0>()]),
                    LineHandle(newCompIds[c.topo.component<1>()]));
                auto pair = std::make_pair(c.topo.hd, std::make_pair(ccid, h));
                if (old2new){
                    (*old2new)[pair.first] = pair.second;
                }
                if (new2old){
                    (*new2old)[pair.second] = pair.first;
                }
            }
            for (auto & c : mg.constraints<RegionLineConnectionData>()){
                int ccid1 = ccids[c.topo.component<0>()];
                int ccid2 = ccids[c.topo.component<1>()];
                if (ccid1 != ccid2){
                    continue;
                }
                int ccid = ccid1;
                auto h = ccs[ccid].addConstraint(c.data,
                    RegionHandle(newCompIds[c.topo.component<0>()]),
                    LineHandle(newCompIds[c.topo.component<1>()]));
                auto pair = std::make_pair(c.topo.hd, std::make_pair(ccid, h));
                if (old2new){
                    (*old2new)[pair.first] = pair.second;
                }
                if (new2old){
                    (*new2old)[pair.second] = pair.first;
                }
            }

            return ccs;
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
                    const RLGraph & mg;
                    const RLGraphControls & controls;
                    RLGraphVars & vars;

                    void operator()(const LineHandle & lh) const {
                        auto & lp = controls[lh];
                        auto & v = vars[lh];
                        if (!lp.used){
                            v.variables = {};
                            return;
                        }
                        if (lp.orientationClaz == -1){
                            // (1/cornerDepth1, 1/cornerDepth2) for LineFree
                            v.variables = { 1.0, 1.0 };
                        }
                        else{
                            // 1/centerDepth for LineOriented,
                            v.variables = { 1.0 };
                        }
                    }

                    void operator()(const RegionHandle & rh) const {
                        auto & rd = mg.data(rh);
                        auto & rp = controls[rh];
                        auto & v = vars[rh];
                        if (!rp.used){
                            v.variables = {};
                            return;
                        }
                        if (rp.orientationClaz == -1 && rp.orientationNotClaz == -1){
                            // (a, b, c) for RegionFree ax+by+c=1, 
                            v.variables = { rd.normalizedCenter[0], rd.normalizedCenter[1], rd.normalizedCenter[2] };
                        }
                        else if (rp.orientationClaz == -1 && rp.orientationNotClaz >= 0){
                            // (a, b), {or (b, c) or (a, c)} for RegionAlongFixedAxis  ax+by+c=1,
                            Vec3 anotherAxis = rd.normalizedCenter.cross(normalize(controls.vanishingPoints[rp.orientationNotClaz]));
                            Vec3 trueNormal = controls.vanishingPoints[rp.orientationNotClaz].cross(anotherAxis);
                            Plane3 plane(rd.normalizedCenter, trueNormal);
                            auto eq = Plane3ToEquation(plane);
                            int c = SwappedComponent(normalize(controls.vanishingPoints[rp.orientationNotClaz]));
                            std::swap(eq[c], eq[2]);
                            v.variables = { eq[0], eq[1] };
                        }
                        else {
                            // 1/centerDepth for RegionWithFixedNormal
                            v.variables = { 1.0 };
                        }
                    }
                };

            

                
            }

        }





        void RLGraphControls::disable(RegionHandle h, const RLGraph & mg){
            componentControls[h].used = false;
            for (auto ch : mg.topo(h).constraints<RegionBoundaryData>()){
                constraintControls[ch].used = false;
            }
            for (auto ch : mg.topo(h).constraints<RegionLineConnectionData>()){
                constraintControls[ch].used = false;
            }
        }

        void RLGraphControls::disable(LineHandle h, const RLGraph & mg){
            componentControls[h].used = false;
            for (auto ch : mg.topo(h).constraints<LineRelationData>()){
                constraintControls[ch].used = false;
            }
            for (auto ch : mg.topo(h).constraints<RegionLineConnectionData>()){
                constraintControls[ch].used = false;
            }
        }

        void RLGraphControls::enable(RegionHandle h, const RLGraph & mg){
            componentControls[h].used = true;
        }

        void RLGraphControls::enable(LineHandle h, const RLGraph & mg){
            componentControls[h].used = true;
        }


        void RLGraphControls::disableAllInvalidConstraints(const RLGraph & mg){
            for (auto & c : mg.constraints<RegionBoundaryData>()){
                if (!componentControls[c.topo.component<0>()].used ||
                    !componentControls[c.topo.component<1>()].used){
                    constraintControls[c.topo.hd].used = false;
                }                
            }
            for (auto & c : mg.constraints<LineRelationData>()){
                if (!componentControls[c.topo.component<0>()].used ||
                    !componentControls[c.topo.component<1>()].used){
                    constraintControls[c.topo.hd].used = false;
                }
            }
            for (auto & c : mg.constraints<RegionLineConnectionData>()){
                if (!componentControls[c.topo.component<0>()].used ||
                    !componentControls[c.topo.component<1>()].used){
                    constraintControls[c.topo.hd].used = false;
                }
            }
        }


        void RLGraphControls::enableAll(){
            for (auto & ct : componentControls.data){
                for (auto & c : ct){
                    c.used = true;
                }
            }
            for (auto & ct : constraintControls.data){
                for (auto & c : ct){
                    c.used = true;
                }
            }
        }


        RLGraphControls::RLGraphControls(const RLGraph & mg, const std::vector<Vec3> & vps) : vanishingPoints(vps) {
            componentControls = MakeHandledTableForAllComponents<RLGraphComponentControl>(mg);
            constraintControls = MakeHandledTableForAllConstraints<RLGraphConstraintControl>(mg);
            for (auto & l : mg.internalComponents<LineData>()){
                componentControls[l.topo.hd].orientationClaz = l.data.initialClaz;
                componentControls[l.topo.hd].orientationNotClaz = -1;
            }
            for (auto & r : mg.internalComponents<RegionData>()){
                componentControls[r.topo.hd].orientationClaz = componentControls[r.topo.hd].orientationNotClaz = -1;
            }
            for (auto & dtable : componentControls.data){
                for (auto & d : dtable){
                    d.used = true;
                }
            }
            for (auto & dtable : constraintControls.data){
                for (auto & d : dtable){
                    d.used = true;
                    d.weight = 1.0;
                }
            }
            ResetWeights(mg, *this);
        }



        std::vector<RLGraphControls> Decompose(const RLGraph & mg,
            const RLGraphControls & controls, const RLGraphComponentTable<int> & ccids, int ccnum){
            std::vector<RLGraphControls> ccs(ccnum);
            RLGraphComponentTable<int> newCompIds = MakeHandledTableForAllComponents<int>(mg);
            for (auto & cc : ccs){
                cc.vanishingPoints = controls.vanishingPoints;
            }
            for (auto & c : mg.components<RegionData>()){
                int ccid = ccids[c.topo.hd];
                auto & regionControls = ccs[ccid].componentControls.dataOfType<RegionHandle>();
                newCompIds[c.topo.hd] = regionControls.size();
                regionControls.push_back(controls[c.topo.hd]);
            }
            for (auto & c : mg.components<LineData>()){
                int ccid = ccids[c.topo.hd];
                auto & lineControls = ccs[ccid].componentControls.dataOfType<LineHandle>();
                newCompIds[c.topo.hd] = lineControls.size();
                lineControls.push_back(controls[c.topo.hd]);
            }
            for (auto & c : mg.constraints<RegionBoundaryData>()){
                int ccid1 = ccids[c.topo.component<0>()];
                int ccid2 = ccids[c.topo.component<1>()];
                if (ccid1 != ccid2){
                    continue;
                }
                int ccid = ccid1;
                auto & ccontrols = ccs[ccid].constraintControls.dataOfType<RegionBoundaryHandle>();
                ccontrols.push_back(controls[c.topo.hd]);
            }
            for (auto & c : mg.constraints<LineRelationData>()){
                int ccid1 = ccids[c.topo.component<0>()];
                int ccid2 = ccids[c.topo.component<1>()];
                if (ccid1 != ccid2){
                    continue;
                }
                int ccid = ccid1;
                auto & ccontrols = ccs[ccid].constraintControls.dataOfType<LineRelationHandle>();
                ccontrols.push_back(controls[c.topo.hd]);
            }
            for (auto & c : mg.constraints<RegionLineConnectionData>()){
                int ccid1 = ccids[c.topo.component<0>()];
                int ccid2 = ccids[c.topo.component<1>()];
                if (ccid1 != ccid2){
                    continue;
                }
                int ccid = ccid1;
                auto & ccontrols = ccs[ccid].constraintControls.dataOfType<RegionLineConnectionHandle>();
                ccontrols.push_back(controls[c.topo.hd]);
            }
            return ccs;
        }






        void ResetWeights(const RLGraph & mg, RLGraphControls & controls){
            auto & consProps = controls.constraintControls;
            double maxChainLen = 0;
            for (auto & b : mg.internalConstraints<RegionBoundaryData>()){
                maxChainLen = std::max(ChainLength(b.data.normalizedSampledPoints), maxChainLen);
            }
            for (auto & c : mg.internalConstraints<RegionLineConnectionData>()){
                maxChainLen = std::max(ChainLength(c.data.normalizedAnchors), maxChainLen);
            }

            for (auto & b : mg.internalConstraints<RegionBoundaryData>()){
                consProps[b.topo.hd].weight = 1.0 * ChainLength(b.data.normalizedSampledPoints) / maxChainLen;
            }
            for (auto & c : mg.internalConstraints<RegionLineConnectionData>()){
                consProps[c.topo.hd].weight = 1.0 * ChainLength(c.data.normalizedAnchors) / maxChainLen;
            }
            for (auto & ll : mg.internalConstraints<LineRelationData>()){
                consProps[ll.topo.hd].weight = 1.0;
            }
            auto & compProps = controls.componentControls;
            for (auto & t : compProps.data){
                for (auto & prop : t){
                    for (auto & wa : prop.weightedAnchors){
                        wa.score = 1.0;
                    }
                }
            }
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



            //enum class OrientationHint {
            //    Void,
            //    Horizontal,
            //    Vertical,
            //    OtherPlanar,
            //    NonPlanar,
            //    Count
            //};

            //inline OrientationHint ToOrientationHint(GeometricContextLabel label, bool leftFrontRightAsVertical){
            //    //THERE_ARE_BUGS_HERE("Only OtherPlanar is returned!");
            //    //return OrientationHint::OtherPlanar;
            //    switch (label){
            //    case GeometricContextLabel::Ceiling: return OrientationHint::Horizontal;
            //    case GeometricContextLabel::Floor: return OrientationHint::Horizontal;
            //    case GeometricContextLabel::Front:
            //    case GeometricContextLabel::Left:
            //    case GeometricContextLabel::Right:
            //        return leftFrontRightAsVertical ? OrientationHint::Vertical : OrientationHint::OtherPlanar;
            //    case GeometricContextLabel::Furniture: return OrientationHint::Void;
            //    case GeometricContextLabel::Ground: return OrientationHint::Horizontal;
            //    case GeometricContextLabel::Sky: return OrientationHint::Void;
            //    case GeometricContextLabel::Vertical: return OrientationHint::OtherPlanar;
            //    case GeometricContextLabel::NotPlanar: return OrientationHint::NonPlanar;
            //    default:
            //        return OrientationHint::Void;
            //    }
            //}

        }


        RLGraphVars MakeVariables(const RLGraph & mg, const RLGraphControls & controls, bool randomized) {
            RLGraphVars vars = MakeHandledTableForAllComponents<RLGraphVar>(mg);
            ForeachRLGraphComponentHandle(mg, InitializeVariablesForEachHandle{ mg, controls, vars });
            if (randomized){
                for (auto & t : vars.data){
                    for (auto & vv : t){
                        for (auto & v : vv.variables){
                            v = ((size_t)std::rand() % 1000 + 1) / 5000.0 + 1.0;
                        }
                    }
                }
            }
            return vars;
        }


        namespace {

            inline double NonZeroize(double d){
                return d == 0.0 ? 1e-6 : d;
            }

        }





        Line3 InstanceGivenVariables(const RLGraph & mg, const double * variables,
            const RLGraphControls & controls, const LineHandle & lh) {
            auto & ld = mg.data(lh);
            auto & c = controls[lh];
            //if (!lp.used)
            //    return Line3();
            if (c.orientationClaz >= 0){
                Ray3 infLine(normalize(ld.line.center()) / NonZeroize(variables[0]), controls.vanishingPoints[c.orientationClaz]);
                return Line3(DistanceBetweenTwoLines(Ray3(Point3(0, 0, 0), normalize(ld.line.first)), infLine).second.second,
                    DistanceBetweenTwoLines(Ray3(Point3(0, 0, 0), normalize(ld.line.second)), infLine).second.second);
            }
            else /*if (line.type == MGUnary::LineFree)*/{
                /*           | sin(theta) | | p | | q |
                len:---------------------------------   => 1/len: [(1/q)sin(phi) - (1/p)sin(phi-theta)] / sin(theta)
                | p sin(phi) - q sin(phi - theta) |
                */
                // variables[0] -> 1/p
                // variables[1] -> 1/q
                return Line3(normalize(ld.line.first) / NonZeroize(variables[0]), normalize(ld.line.second) / NonZeroize(variables[1]));
            }
        }


        Plane3 InstanceGivenVariables(const RLGraph & mg, const double * variables,
            const RLGraphControls & controls, const RegionHandle & rh){
            auto & rd = mg.data(rh);
            auto & c = controls[rh];
            /*if (!rp.used)
            return Plane3();*/
            if (c.orientationClaz >= 0){
                return Plane3(rd.normalizedCenter / NonZeroize(variables[0]), controls.vanishingPoints[c.orientationClaz]);
            }
            else if (c.orientationClaz == -1 && c.orientationNotClaz >= 0){
                double vs[] = { variables[0], variables[1], 0.0 }; // fake vs
                // v1 * o1 + v2 * o2 + v3 * o3 = 0, since region.orientation is orthogonal to normal
                auto orientation = normalize(controls.vanishingPoints[c.orientationNotClaz]);
                int c = SwappedComponent(orientation);
                std::swap(orientation[c], orientation[2]); // now fake orientation
                vs[2] = (-vs[0] * orientation[0] - vs[1] * orientation[1])
                    / orientation[2];
                std::swap(vs[c], vs[2]); // now real vs
                return Plane3FromEquation(vs[0], vs[1], vs[2]);
            }
            else /*if (region.type == MGUnary::RegionWithFixedNormal)*/{
                return Plane3FromEquation(variables[0], variables[1], variables[2]);
            }
        }


        Line3 Instance(const RLGraph & mg, const RLGraphControls & controls, const RLGraphVars & vars, const LineHandle & lh) {
            //auto & ld = mg.data(lh);
            //auto & c = controls[lh];
            //auto & v = vars[lh];
            ////if (!lp.used)
            ////    return Line3();
            //if (c.orientationClaz >= 0){
            //    assert(v.variables.size() == 1);
            //    Ray3 infLine(normalize(ld.line.center()) / NonZeroize(v.variables[0]), controls.vanishingPoints[c.orientationClaz]);
            //    return Line3(DistanceBetweenTwoLines(Ray3(Point3(0, 0, 0), normalize(ld.line.first)), infLine).second.second,
            //        DistanceBetweenTwoLines(Ray3(Point3(0, 0, 0), normalize(ld.line.second)), infLine).second.second);
            //}
            //else /*if (line.type == MGUnary::LineFree)*/{
            //    assert(v.variables.size() == 2);
            //    /*           | sin(theta) | | p | | q |
            //    len:---------------------------------   => 1/len: [(1/q)sin(phi) - (1/p)sin(phi-theta)] / sin(theta)
            //    | p sin(phi) - q sin(phi - theta) |
            //    */
            //    // variables[0] -> 1/p
            //    // variables[1] -> 1/q
            //    return Line3(normalize(ld.line.first) / NonZeroize(v.variables[0]), normalize(ld.line.second) / NonZeroize(v.variables[1]));
            //}
            return InstanceGivenVariables(mg, vars.at(lh).variables.data(), controls, lh);
        }



        Plane3 Instance(const RLGraph & mg, const RLGraphControls & controls, const RLGraphVars & vars, const RegionHandle & rh){
            //auto & rd = mg.data(rh);
            //auto & c = controls[rh];
            //auto & v = vars[rh];
            ///*if (!rp.used)
            //    return Plane3();*/
            //if (c.orientationClaz >= 0){
            //    assert(v.variables.size() == 1);
            //    return Plane3(rd.normalizedCenter / NonZeroize(v.variables[0]), controls.vanishingPoints[c.orientationClaz]);
            //}
            //else if (c.orientationClaz == -1 && c.orientationNotClaz >= 0){
            //    assert(v.variables.size() == 2);
            //    double vs[] = { v.variables[0], v.variables[1], 0.0 }; // fake vs
            //    // v1 * o1 + v2 * o2 + v3 * o3 = 0, since region.orientation is orthogonal to normal
            //    auto orientation = normalize(controls.vanishingPoints[c.orientationNotClaz]);
            //    int c = SwappedComponent(orientation);
            //    std::swap(orientation[c], orientation[2]); // now fake orientation
            //    vs[2] = (-vs[0] * orientation[0] - vs[1] * orientation[1])
            //        / orientation[2];
            //    std::swap(vs[c], vs[2]); // now real vs
            //    return Plane3FromEquation(vs[0], vs[1], vs[2]);
            //}
            //else /*if (region.type == MGUnary::RegionWithFixedNormal)*/{
            //    assert(v.variables.size() == 3);
            //    return Plane3FromEquation(v.variables[0], v.variables[1], v.variables[2]);
            //}
            return InstanceGivenVariables(mg, vars.at(rh).variables.data(), controls, rh);
        }




        std::vector<Polygon3> RegionPolygon(const RLGraph & mg, const RLGraphControls & controls,
            const RLGraphVars & vars, RegionHandle rh) {
            if (!controls[rh].used)
                return std::vector<Polygon3>();
            std::vector<Polygon3> ps;
            Plane3 plane = Instance(mg, controls, vars, rh);
            for (auto & contour : mg.data(rh).normalizedContours){
                Polygon3 polygon;
                polygon.corners.reserve(contour.size());
                for (auto & c : contour){
                    polygon.corners.push_back(PointAt(c, plane));
                }
                polygon.normal = plane.normal;
                ps.push_back(std::move(polygon));
            }
            return ps;
        }

        HandledTable<RegionHandle, std::vector<Polygon3>> RegionPolygons(const RLGraph & mg,
            const RLGraphControls & controls,
            const RLGraphVars & vars) {

            auto polygons = mg.createComponentTable<RegionData, std::vector<Polygon3>>();
            for (auto & r : mg.components<RegionData>()){
                if (!controls[r.topo.hd].used)
                    continue;
                polygons[r.topo.hd] = RegionPolygon(mg, controls, vars, r.topo.hd);
            }
            return polygons;
        }





        std::vector<double> VariableCoefficientsForInverseDepthAtDirection(const RLGraph & mg,
            const RLGraphControls & controls, const RLGraphVars & vars, const Vec3 & direction, const LineHandle & lh){
            auto & ld = mg.data(lh);
            auto & c = controls[lh];
            auto & v = vars[lh];
            if (c.orientationClaz >= 0){
                assert(v.variables.size() == 1);
                Ray3 infLine(normalize(ld.line.center()), controls.vanishingPoints[c.orientationClaz]);
                 // variable is 1.0/centerDepth
                 // corresponding coeff is 1.0/depthRatio
                 // so that 1.0/depth = 1.0/centerDepth * 1.0/depthRatio -> depth = centerDepth * depthRatio
                 double depthRatio = norm(DistanceBetweenTwoLines(Ray3(Point3(0, 0, 0), direction), infLine).second.first);
                 return std::vector<double>{1.0 / depthRatio};
             }
             else /*if(u.type == MGUnary::LineFree)*/{
                 assert(v.variables.size() == 2);
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


        std::vector<double> VariableCoefficientsForInverseDepthAtDirection(const RLGraph & mg,
            const RLGraphControls & controls, const RLGraphVars & vars, const Vec3 & direction, const RegionHandle & rh){
            auto & rd = mg.data(rh);
            auto & c = controls[rh];
            auto & v = vars[rh];
            if (c.orientationClaz >= 0){
                assert(v.variables.size() == 1);
                Plane3 plane(rd.normalizedCenter, controls.vanishingPoints[c.orientationClaz]);
                // variable is 1.0/centerDepth
                // corresponding coeff is 1.0/depthRatio
                // so that 1.0/depth = 1.0/centerDepth * 1.0/depthRatio -> depth = centerDepth * depthRatio
                double depthRatio = norm(IntersectionOfLineAndPlane(Ray3(Point3(0, 0, 0), direction), plane).position);
                return std::vector<double>{1.0 / depthRatio};
            }
            else if (c.orientationClaz == -1 && c.orientationNotClaz >= 0){
                assert(v.variables.size() == 2);
                // depth = 1.0 / (ax + by + cz) where (x, y, z) = direction, (a, b, c) = variables
                // -> 1.0/depth = ax + by + cz
                // -> 1.0/depth = ax + by + (- o1a - o2b)/o3 * z = a(x-o1*z/o3) + b(y-o2*z/o3)
                auto orientation = normalize(controls.vanishingPoints[c.orientationNotClaz]);
                int sc = SwappedComponent(orientation);
                Vec3 forientation = orientation;
                std::swap(forientation[sc], forientation[2]);
                Vec3 fdirection = direction;
                std::swap(fdirection[sc], fdirection[2]);
                return std::vector<double>{
                    fdirection[0] - forientation[0] * fdirection[2] / forientation[2],
                        fdirection[1] - forientation[1] * fdirection[2] / forientation[2]
                };
            }
            else{
                assert(v.variables.size() == 3);
                // depth = 1.0 / (ax + by + cz) where (x, y, z) = direction, (a, b, c) = variables
                // -> 1.0/depth = ax + by + cz
                return std::vector<double>{direction[0], direction[1], direction[2]};
            }
        }





        double DepthAtDirectionGivenVariables(const RLGraph & mg, const double * variables,
            const RLGraphControls & controls, const Vec3 & direction, const LineHandle & lh){
            auto & ld = mg.data(lh);
            auto & lp = controls[lh];
            if (lp.orientationClaz >= 0){
                //assert(lp.variables.size() == 1);
                Ray3 infLine(normalize(ld.line.center()), controls.vanishingPoints[lp.orientationClaz]);
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


        double DepthAtDirectionGivenVariables(const RLGraph & mg, const double * variables,
            const RLGraphControls & controls, const Vec3 & direction, const RegionHandle & rh){
            auto & rd = mg.data(rh);
            auto & rp = controls[rh];
            if (rp.orientationClaz >= 0){
                //assert(rp.variables.size() == 1);
                Plane3 plane(rd.normalizedCenter, controls.vanishingPoints[rp.orientationClaz]);
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
                auto orientation = normalize(controls.vanishingPoints[rp.orientationNotClaz]);
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
                // depth = 1.0 / (ax + by + cz) where (x, y, z) = direction, (a, b, c) = variables
                // -> 1.0/depth = ax + by + cz
                return 1.0 / (direction[0] * variables[0] + direction[1] * variables[1] + direction[2] * variables[2]);
            }
        }




        bool AttachAnchorToCenterOfLargestRegionIfNoAnchorExists(const RLGraph & mg,
            RLGraphControls & controls,
            double depth, double weight){
            std::cout << "calling this feature may cause EYE SKEW!" << std::endl;
            RegionHandle largest;
            double maxArea = 0.0;
            for (auto & r : mg.components<RegionData>()){
                if (!controls[r.topo.hd].used)
                    continue;
                if (controls[r.topo.hd].weightedAnchors.size() > 0)
                    return true;
                if (r.data.area > maxArea){
                    largest = r.topo.hd;
                    maxArea = r.data.area;
                }
            }
            if(largest.valid())
                controls[largest].weightedAnchors.push_back(ScoreAs(mg.data(largest).normalizedCenter * depth, weight));
            return largest.valid();
        }


        int NumberOfAnchors(const RLGraphControls & controls){
            int nanchor = 0;
            for (auto & table : controls.componentControls.data){
                for (auto & compProp : table){
                    if (!compProp.used)
                        continue;
                    nanchor += compProp.weightedAnchors.size();
                }
            }
            return nanchor;
        }

        bool AttachAnchorToCenterOfLargestLineIfNoAnchorExists(const RLGraph & mg,
            RLGraphControls & controls,
            double depth, double weight, bool orientedOnly){
            LineHandle largest;
            double maxArea = 0.0;
            for (auto & r : mg.components<LineData>()){
                if (!controls[r.topo.hd].used)
                    continue;
                if (controls[r.topo.hd].orientationClaz == -1 && orientedOnly)
                    continue;
                if (controls[r.topo.hd].weightedAnchors.size() > 0)
                    return true;
                if (r.data.line.length() > maxArea){
                    largest = r.topo.hd;
                    maxArea = r.data.line.length();
                }
            }
            if(largest.valid())
                controls[largest].weightedAnchors.push_back(ScoreAs(normalize(mg.data(largest).line.center()) * depth, weight));
            return largest.valid();
        }

        bool AttachAnchorToCenterOfLargestRegionIfNoAnchorExists(const RLGraph & mg,
            RLGraphControls & controls,
            double depth, double weight, bool orientedOnly){
            RegionHandle largest;
            double maxArea = 0.0;
            for (auto & r : mg.components<RegionData>()){
                if (!controls[r.topo.hd].used)
                    continue;
                if (controls[r.topo.hd].orientationClaz == -1 && orientedOnly)
                    continue;
                if (controls[r.topo.hd].weightedAnchors.size() > 0)
                    return true;
                if (r.data.area > maxArea){
                    largest = r.topo.hd;
                    maxArea = r.data.area;
                }
            }
            if (largest.valid())
                controls[largest].weightedAnchors.push_back(ScoreAs(normalize(mg.data(largest).normalizedCenter) * depth, weight));
            return largest.valid();
        }


        void ClearAllComponentAnchors(RLGraphControls & controls){
            for (auto & table : controls.componentControls.data){
                for (auto & compProp : table){
                    compProp.weightedAnchors.clear();
                }
            }
        }



        namespace {

            int RegisterVariablePositions(const RLGraph & mg,
                const RLGraphControls & controls, const RLGraphVars & vars,
                RLGraphComponentTable<int> & uh2varStartPosition){
                int varNum = 0;
                for (auto & c : mg.components<LineData>()){
                    if (!controls[c.topo.hd].used)
                        continue;
                    uh2varStartPosition[c.topo.hd] = varNum;
                    varNum += vars[c.topo.hd].variables.size();
                }
                for (auto & c : mg.components<RegionData>()){
                    if (!controls[c.topo.hd].used)
                        continue;
                    uh2varStartPosition[c.topo.hd] = varNum;
                    varNum += vars[c.topo.hd].variables.size();
                }
                return varNum;
            }


            void RegisterVariableValues(const RLGraph & mg,
                const RLGraphControls & controls, const RLGraphVars & vars,
                const RLGraphComponentTable<int> & uh2varStartPosition,
                std::vector<double> & X){
                for (auto & c : mg.components<LineData>()){
                    if (!controls[c.topo.hd].used)
                        continue;
                    for (int i = 0; i < vars[c.topo.hd].variables.size(); i++){
                        X[uh2varStartPosition.at(c.topo.hd) + i] = vars[c.topo.hd].variables.at(i);
                    }
                }
                for (auto & c : mg.components<RegionData>()){
                    if (!controls[c.topo.hd].used)
                        continue;
                    for (int i = 0; i < vars[c.topo.hd].variables.size(); i++){
                        X[uh2varStartPosition.at(c.topo.hd) + i] = vars[c.topo.hd].variables.at(i);
                    }
                }
            }



            struct ExtractNecessaryAnchorsForBinary {
                std::vector<Vec3> operator()(const RLGraph & mg, RegionBoundaryHandle bh) const {
                    size_t n = ElementsNum(mg.data(bh).normalizedSampledPoints);

                    assert(n > 0);
                    const auto & points = mg.data(bh).normalizedSampledPoints;

                    if (n == 1){
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
                    IMPROVABLE_HERE(? );
                    return{ p1, p2, p3 };
                }

                inline std::vector<Vec3> operator()(const RLGraph & mg, LineRelationHandle bh) const {
                    return{ mg.data(bh).normalizedRelationCenter };
                }

                inline std::vector<Vec3> operator()(const RLGraph & mg, RegionLineConnectionHandle bh) const {
                    return{ mg.data(bh).normalizedAnchors.front(), mg.data(bh).normalizedAnchors.back() };
                }
            };

            struct ExtractAllAnchorsForBinary {
                std::vector<Vec3> operator()(const RLGraph & mg, RegionBoundaryHandle bh) const{
                    std::vector<Vec3> anchors;
                    anchors.reserve(ElementsNum(mg.data(bh).normalizedSampledPoints));
                    for (auto & ps : mg.data(bh).normalizedSampledPoints){
                        for (auto & p : ps){
                            anchors.push_back(p);
                        }
                    }
                    return anchors;
                }
                inline std::vector<Vec3> operator()(const RLGraph & mg, LineRelationHandle bh) const {
                    return{ mg.data(bh).normalizedRelationCenter };
                }
                inline std::vector<Vec3> operator()(const RLGraph & mg, RegionLineConnectionHandle bh) const {
                    return mg.data(bh).normalizedAnchors;
                }
            };



            template <class BinaryAnchorExtractorT>
            int RegisterConstraintPositions(const RLGraph & mg,
                const RLGraphControls & controls, 
                const RLGraphComponentTable<int> & uh2varStartPosition,
                RLGraphComponentTable<int> & uh2anchorConsStartPosition,
                RLGraphConstraintTable<int> & bh2consStartPosition,
                RLGraphConstraintTable<std::vector<Vec3>> & appliedBinaryAnchors,
                RLGraphConstraintTable<double> & weightsForEachAppliedBinaryAnchor,
                BinaryAnchorExtractorT && extractBinaryAnchors){

                int consNum = 0;
                // register anchor constraint positions
                for (auto & c : mg.components<RegionData>()){
                    if (!controls[c.topo.hd].used)
                        continue;
                    uh2anchorConsStartPosition[c.topo.hd] = consNum;
                    consNum += controls[c.topo.hd].weightedAnchors.size();
                }
                for (auto & c : mg.components<LineData>()){
                    if (!controls[c.topo.hd].used)
                        continue;
                    uh2anchorConsStartPosition[c.topo.hd] = consNum;
                    consNum += controls[c.topo.hd].weightedAnchors.size();
                }
                // register constraint positions
                for (auto & c : mg.constraints<RegionBoundaryData>()){
                    if (!controls[c.topo.hd].used)
                        continue;
                    if (!(controls[c.topo.component<0>()].used && controls[c.topo.component<1>()].used))
                        continue;
                    bh2consStartPosition[c.topo.hd] = consNum;
                    appliedBinaryAnchors[c.topo.hd] = extractBinaryAnchors(mg, c.topo.hd);
                    weightsForEachAppliedBinaryAnchor[c.topo.hd] = controls[c.topo.hd].weight / sqrt(appliedBinaryAnchors[c.topo.hd].size());
                    consNum += appliedBinaryAnchors[c.topo.hd].size();
                }
                for (auto & c : mg.constraints<LineRelationData>()){
                    if (!controls[c.topo.hd].used)
                        continue;
                    if (!(controls[c.topo.component<0>()].used && controls[c.topo.component<1>()].used))
                        continue;
                    bh2consStartPosition[c.topo.hd] = consNum;
                    appliedBinaryAnchors[c.topo.hd] = extractBinaryAnchors(mg, c.topo.hd);
                    weightsForEachAppliedBinaryAnchor[c.topo.hd] = controls[c.topo.hd].weight / sqrt(appliedBinaryAnchors[c.topo.hd].size());
                    consNum += appliedBinaryAnchors[c.topo.hd].size();
                }
                for (auto & c : mg.constraints<RegionLineConnectionData>()){
                    if (!controls[c.topo.hd].used)
                        continue;
                    if (!(controls[c.topo.component<0>()].used && controls[c.topo.component<1>()].used))
                        continue;
                    bh2consStartPosition[c.topo.hd] = consNum;
                    appliedBinaryAnchors[c.topo.hd] = extractBinaryAnchors(mg, c.topo.hd);
                    weightsForEachAppliedBinaryAnchor[c.topo.hd] = controls[c.topo.hd].weight / sqrt(appliedBinaryAnchors[c.topo.hd].size());
                    consNum += appliedBinaryAnchors[c.topo.hd].size();
                }
                return consNum;
            }




            template <class ComponentDataT, class SparseMatElementT>
            inline void RegisterComponentAnchorEquations(const RLGraph & mg,
                const RLGraphControls & controls,
                const RLGraphVars & vars,
                const RLGraphComponentTable<int> & uh2varStartPosition,
                const RLGraphComponentTable<int> & uh2anchorConsStartPosition,
                std::vector<SparseMatElementT> & Atriplets,
                std::vector<SparseMatElementT> & Wtriplets,
                std::vector<double> & B){
                // add anchors on components
                for (auto & c : mg.components<ComponentDataT>()){
                    if (!controls[c.topo.hd].used)
                        continue;
                    auto uh = c.topo.hd;
                    int uhVarNum = vars[uh].variables.size();
                    int uhVarStartPosition = uh2varStartPosition.at(uh);
                    int eid = uh2anchorConsStartPosition[uh];
                    for (auto & wa : controls[c.topo.hd].weightedAnchors){
                        const Point3 & anchor = wa.component;
                        double weight = wa.score;
                        auto uhVarCoeffsAtAnchorDirection = VariableCoefficientsForInverseDepthAtDirection(mg, controls, vars, anchor, uh);
                        assert(uhVarCoeffsAtAnchorDirection.size() == uhVarNum);
                        for (int i = 0; i < uhVarCoeffsAtAnchorDirection.size(); i++){
                            //A.insert(eid, uhVarStartPosition + i) = uhVarCoeffsAtCenter[i];
                            Atriplets.emplace_back(eid, uhVarStartPosition + i, uhVarCoeffsAtAnchorDirection[i]);
                        }
                        B[eid] = 1.0 / norm(anchor);
                        Wtriplets.emplace_back(eid, eid, weight);
                        eid++;
                    }
                }

            }

            template <class ConstraintDataT, class SparseMatElementT>
            inline void RegisterConstraintEquations(const RLGraph & mg, 
                const RLGraphControls & controls,
                const RLGraphVars & vars,
                const RLGraphComponentTable<int> & uh2varStartPosition,
                const RLGraphConstraintTable<int> & bh2consStartPosition,
                const RLGraphConstraintTable<std::vector<Vec3>> & appliedBinaryAnchors,
                const RLGraphConstraintTable<double> & weightsForEachAppliedBinaryAnchor,
                std::vector<SparseMatElementT> & Atriplets,
                std::vector<SparseMatElementT> & Wtriplets,
                std::vector<double> & B) {

                for (auto & c : mg.constraints<ConstraintDataT>()){
                    if (!controls[c.topo.hd].used)
                        continue;
                    auto bh = c.topo.hd;
                    auto uh1 = mg.topo(bh).component<0>();
                    auto uh2 = mg.topo(bh).component<1>();
                    assert(controls[uh1].used || controls[uh2].used);
                    auto & u1 = mg.data(uh1);
                    auto & u2 = mg.data(uh2);

                    int u1VarStartPosition = uh2varStartPosition.at(uh1);
                    int u1VarNum = vars[uh1].variables.size();

                    int u2VarStartPosition = uh2varStartPosition.at(uh2);
                    int u2VarNum = vars[uh2].variables.size();

                    int eid = bh2consStartPosition[bh];

                    for (auto & a : appliedBinaryAnchors.at(bh)){
                        B[eid] = 0.0;
                        Wtriplets.emplace_back(eid, eid, weightsForEachAppliedBinaryAnchor.at(bh));
                        {
                            auto u1VarCoeffs = VariableCoefficientsForInverseDepthAtDirection(mg, controls, vars, a, uh1);
                            assert(u1VarCoeffs.size() == u1VarNum);
                            for (int i = 0; i < u1VarCoeffs.size(); i++){
                                //A.insert(eid, u1VarStartPosition + i) = u1VarCoeffs[i]; // pos
                                Atriplets.emplace_back(eid, u1VarStartPosition + i, u1VarCoeffs[i]);
                            }
                        }
                        {
                            auto u2VarCoeffs = VariableCoefficientsForInverseDepthAtDirection(mg, controls, vars, a, uh2);
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

        

            template <class VectorT>
            void InstallVariables(const RLGraph & mg, const RLGraphControls & controls,
                const RLGraphComponentTable<int> & uh2varStartPosition, const VectorT & X,
                RLGraphVars & vars){
                for (auto & c : mg.components<RegionData>()){
                    if (!controls[c.topo.hd].used)
                        continue;
                    int uhStartPosition = uh2varStartPosition.at(c.topo.hd);
                    for (int i = 0; i < vars[c.topo.hd].variables.size(); i++){
                        vars[c.topo.hd].variables[i] = X[uhStartPosition + i];
                    }
                }
                for (auto & c : mg.components<LineData>()){
                    if (!controls[c.topo.hd].used)
                        continue;
                    int uhStartPosition = uh2varStartPosition.at(c.topo.hd);
                    for (int i = 0; i < vars[c.topo.hd].variables.size(); i++){
                        vars[c.topo.hd].variables[i] = X[uhStartPosition + i];
                    }
                }
            }

        }

        // inverse depth optimization
        RLGraphVars SolveVariables(const RLGraph & mg, const RLGraphControls & controls, bool useWeights, bool useAllAnchors) {
            RLGraphVars vars = MakeVariables(mg, controls);

            int nanchor = NumberOfAnchors(controls);
            if (nanchor == 0){
                WARNNING("no anchor is given! will ouput very very poor results!!!");
            }
            assert(nanchor > 0);

            SetClock();

            auto uh2varStartPosition = MakeHandledTableForAllComponents<int>(mg);
            auto uh2anchorConsStartPosition = MakeHandledTableForAllComponents<int>(mg);
            auto bh2consStartPosition = MakeHandledTableForAllConstraints<int>(mg);
            auto appliedBinaryAnchors = MakeHandledTableForAllConstraints<std::vector<Vec3>>(mg);
            auto weightsForEachAppliedBinaryAnchor = MakeHandledTableForAllConstraints<double>(mg);

            // register vars
            int varNum = RegisterVariablePositions(mg, controls, vars, uh2varStartPosition);
            std::vector<double> Xdata(varNum, 1.0);
            RegisterVariableValues(mg, controls, vars, uh2varStartPosition, Xdata);

            // register cons
            int consNum = useAllAnchors 
                ? 
                RegisterConstraintPositions(mg, controls, uh2varStartPosition,
                uh2anchorConsStartPosition, bh2consStartPosition, appliedBinaryAnchors,
                weightsForEachAppliedBinaryAnchor, ExtractAllAnchorsForBinary())
                :
                RegisterConstraintPositions(mg, controls, uh2varStartPosition,
                uh2anchorConsStartPosition, bh2consStartPosition, appliedBinaryAnchors,
                weightsForEachAppliedBinaryAnchor, ExtractNecessaryAnchorsForBinary());

            std::vector<double> Bdata(consNum, 0.0);
            std::vector<Eigen::Triplet<double>> Atriplets;
            std::vector<Eigen::Triplet<double>> Wtriplets;
            Atriplets.reserve(consNum * 6);
            Wtriplets.reserve(consNum);

            RegisterComponentAnchorEquations<RegionData>(mg, controls, vars, uh2varStartPosition, uh2anchorConsStartPosition,
                Atriplets, Wtriplets, Bdata);
            RegisterComponentAnchorEquations<LineData>(mg, controls, vars, uh2varStartPosition, uh2anchorConsStartPosition,
                Atriplets, Wtriplets, Bdata);

            RegisterConstraintEquations<RegionBoundaryData>(mg, controls, vars, uh2varStartPosition, bh2consStartPosition,
                appliedBinaryAnchors, weightsForEachAppliedBinaryAnchor, Atriplets, Wtriplets, Bdata);
            RegisterConstraintEquations<LineRelationData>(mg, controls, vars, uh2varStartPosition, bh2consStartPosition,
                appliedBinaryAnchors, weightsForEachAppliedBinaryAnchor, Atriplets, Wtriplets, Bdata);
            RegisterConstraintEquations<RegionLineConnectionData>(mg, controls, vars, uh2varStartPosition, bh2consStartPosition,
                appliedBinaryAnchors, weightsForEachAppliedBinaryAnchor, Atriplets, Wtriplets, Bdata);

            // matrices
            Eigen::SparseMatrix<double> A;
            {
                Clock clock("form matrix A");
                A.resize(consNum, varNum);
                A.setFromTriplets(Atriplets.begin(), Atriplets.end());
            }
            Eigen::Map<const Eigen::VectorXd> B(Bdata.data(), Bdata.size());


            Eigen::SparseMatrix<double> WA;
            Eigen::VectorXd WB;
            if (useWeights){
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
                    SHOULD_NEVER_BE_CALLED();
                }

                X = solver.solve(useWeights ? WB : B);
                if (solver.info() != Eigen::Success) {
                    assert(0);
                    std::cout << "solving error" << std::endl;
                    SHOULD_NEVER_BE_CALLED();
                }
            }
            else{
                Clock clock("solve equations using Eigen::SPQR");

                Eigen::SPQR<Eigen::SparseMatrix<double, Eigen::ColMajor>> solver;
                solver.compute(useWeights ? WA : A);

                if (solver.info() != Eigen::Success) {
                    assert(0);
                    std::cout << "computation error" << std::endl;
                    SHOULD_NEVER_BE_CALLED();
                }

                X = solver.solve(useWeights ? WB : B);
                if (solver.info() != Eigen::Success) {
                    assert(0);
                    std::cout << "solving error" << std::endl;
                    SHOULD_NEVER_BE_CALLED();
                }
            }

            {
                Clock clock("install solved variables");
                InstallVariables(mg, controls, uh2varStartPosition, X, vars);
            }

            return vars;
        }

        


        namespace {

            
            template <class T>
            struct InstanceTableByHandle_ {};
            template <>
            struct InstanceTableByHandle_<RegionHandle> {
                using type = HandledTable<RegionHandle, Plane3>;
            };
            template <>
            struct InstanceTableByHandle_<LineHandle> {
                using type = HandledTable<LineHandle, Line3>;
            };
            template <class T>
            using InstanceTableByHandle = typename InstanceTableByHandle_<T>::type;

            using RLGraphInstanceTable = MetaBind<InstanceTableByHandle, RegionHandle, LineHandle>;


            double DistanceToInstance(const Point3 & p, const Plane3 & plane){
                return plane.distanceTo(p);
            }
            double DistanceToInstance(const Point3 & p, const Line3 & line){
                return DistanceFromPointToLine(p, line).first;
            }


            template <class ComponentDataT>
            void SetComponentAnchorCosts(const Eigen::VectorXd & variables, Eigen::VectorXd & costs,
                const RLGraph & mg, const RLGraphControls & controls,
                const RLGraphComponentTable<int> & compAnchorConsStartPosition,
                const RLGraphInstanceTable & insts,
                bool useWeights){
                for (auto & c : mg.components<ComponentDataT>()){
                    if (!controls[c.topo.hd].used)
                        continue;
                    int anchorStartPosition = compAnchorConsStartPosition[c.topo.hd];
                    const auto & inst = insts[c.topo.hd];
                    int i = 0;
                    for (auto & wa : controls[c.topo.hd].weightedAnchors){
                        const Point3 & anchor = wa.component;
                        double weight = wa.score;
                        double dist = DistanceToInstance(anchor, inst);
                        costs[anchorStartPosition + i] = dist * (useWeights ? weight : 1.0);
                        i++;
                    }
                }
            }

            template <class ConstarintDataT>
            void SetConstraintAnchorCosts(const Eigen::VectorXd & variables, Eigen::VectorXd & costs,
                const RLGraph & mg, const RLGraphControls & controls,
                const RLGraphConstraintTable<int> & consAnchorConsStartPosition,
                const RLGraphConstraintTable<std::vector<Vec3>> & appliedConsAnchors,
                const RLGraphInstanceTable & insts,
                bool useWeights){
                for (auto & c : mg.constraints<ConstarintDataT>()){
                    if (!controls[c.topo.hd].used)
                        continue;
                    // constraint anchors
                    int anchorStartPosition = consAnchorConsStartPosition[c.topo.hd];
                    auto & anchors = appliedConsAnchors.at(c.topo.hd);
                    auto uh1 = c.topo.component<0>();
                    auto uh2 = c.topo.component<1>();
                    auto & inst1 = insts[uh1];
                    auto & inst2 = insts[uh2];
                    for (int i = 0; i < anchors.size(); i++){
                        double depth1 = DepthAt(anchors[i], inst1);
                        double depth2 = DepthAt(anchors[i], inst2);
                        costs(anchorStartPosition + i) = abs(depth1 - depth2) *
                            (useWeights ? (controls[c.topo.hd].weight / sqrt(anchors.size())) : 1.0);
                    }
                }
            }

        }



        void OptimizeVariables(const RLGraph & mg, const RLGraphControls & controls, RLGraphVars & vars, 
            bool useWeights, bool useAllAnchors, 
            const std::function<bool(const RLGraphVars &)> & callback){

            SetClock();

            auto uh2varStartPosition = MakeHandledTableForAllComponents<int>(mg);

            // register vars
            int varNum = RegisterVariablePositions(mg, controls, vars, uh2varStartPosition);
            std::vector<double> Xdata(varNum, 1.0);
            RegisterVariableValues(mg, controls, vars, uh2varStartPosition, Xdata);


            auto compAnchorConsStartPosition = MakeHandledTableForAllComponents<int>(mg);

            auto appliedConsAnchors = MakeHandledTableForAllConstraints<std::vector<Vec3>>(mg);
            auto weightsForEachAppliedConsAnchor = MakeHandledTableForAllConstraints<double>(mg);
            auto consAnchorConsStartPosition = MakeHandledTableForAllConstraints<int>(mg);
            
            //auto regionBoundaryNormalConsStartPosition = MakeHandledTableForConstraints<int, RegionBoundaryData>(mg);


            // register cons
            int consNum = 0;
            // register comp anchors and cons anchors
            consNum += useAllAnchors ?
                RegisterConstraintPositions(mg, controls, uh2varStartPosition, compAnchorConsStartPosition,
                consAnchorConsStartPosition, appliedConsAnchors, weightsForEachAppliedConsAnchor,
                ExtractAllAnchorsForBinary())
                :
                RegisterConstraintPositions(mg, controls, uh2varStartPosition, compAnchorConsStartPosition,
                consAnchorConsStartPosition, appliedConsAnchors, weightsForEachAppliedConsAnchor,
                ExtractNecessaryAnchorsForBinary());
            //// register region boundary normal consistencies
            //for (auto & b : mg.constraints<RegionBoundaryData>()){
            //    regionBoundaryNormalConsStartPosition[b.topo.hd] = consNum;
            //    consNum++;
            //}

            auto costFunctor = misc::MakeGenericNumericDiffFunctor<double>(
                [&](const Eigen::VectorXd & v, Eigen::VectorXd & costs){

                Eigen::VectorXd variables = v;

                // compute current planes and lines
                RLGraphInstanceTable insts;
                insts.container<RegionHandle>() = mg.createComponentTable<RegionData, Plane3>();
                insts.container<LineHandle>() = mg.createComponentTable<LineData, Line3>();
                for (auto & c : mg.components<RegionData>()){
                    insts[c.topo.hd] = InstanceGivenVariables(mg,
                        variables.data() + uh2varStartPosition.at(c.topo.hd), controls, c.topo.hd);
                }
                for (auto & c : mg.components<LineData>()){
                    insts[c.topo.hd] = InstanceGivenVariables(mg,
                        variables.data() + uh2varStartPosition.at(c.topo.hd), controls, c.topo.hd);
                }
                
                SetComponentAnchorCosts<RegionData>(variables, costs, mg, controls, compAnchorConsStartPosition, insts, useWeights);
                SetComponentAnchorCosts<LineData>(variables, costs, mg, controls, compAnchorConsStartPosition, insts, useWeights);

                SetConstraintAnchorCosts<RegionBoundaryData>(variables, costs, mg, controls, 
                    consAnchorConsStartPosition, appliedConsAnchors, insts, useWeights);
                SetConstraintAnchorCosts<LineRelationData>(variables, costs, mg, controls,
                    consAnchorConsStartPosition, appliedConsAnchors, insts, useWeights);
                SetConstraintAnchorCosts<RegionLineConnectionData>(variables, costs, mg, controls,
                    consAnchorConsStartPosition, appliedConsAnchors, insts, useWeights);

                //for (auto & c : mg.constraints<RegionBoundaryData>()){
                //    if (!controls[c.topo.hd].used)
                //        continue;
                //    // constraint anchors
                //    int anchorStartPosition = consAnchorConsStartPosition[c.topo.hd];
                //    auto & anchors = appliedConsAnchors.at(c.topo.hd);
                //    auto uh1 = c.topo.component<0>();
                //    auto uh2 = c.topo.component<1>();
                //    const Plane3 & inst1 = insts[uh1];
                //    const Plane3 & inst2 = insts[uh2];
                //    // region boundary normal consistency
                //    costs(regionBoundaryNormalConsStartPosition[c.topo.hd]) = 
                //        AngleBetweenUndirectedVectors(inst1.normal, inst2.normal) / 10.0 * c.data.length;
                //}
            }, varNum, consNum);


            Eigen::VectorXd X = Eigen::Map<const Eigen::VectorXd>(Xdata.data(), Xdata.size());
            Eigen::LevenbergMarquardt<decltype(costFunctor)> lm(costFunctor);
            lm.parameters.maxfev = 5000;
            //auto status = lm.minimize(X);
            Eigen::LevenbergMarquardtSpace::Status status = lm.minimizeInit(X);
            assert(status != Eigen::LevenbergMarquardtSpace::ImproperInputParameters);
            do {
                status = lm.minimizeOneStep(X);                
                std::cout << "iter: " << lm.iter << "\t\t lm.fnorm = " << lm.fnorm << std::endl;

                // test
                if(0){
                    Eigen::VectorXd costs = Eigen::VectorXd::Zero(consNum);
                    costFunctor(X, costs);
                    std::cout << "norm(costFunctor(X)) = " << costs.norm() << std::endl;
                    costFunctor(X * 2, costs);
                    std::cout << "norm(costFunctor(X * 2)) = " << costs.norm() << std::endl;
                    costFunctor(X * 3, costs);
                    std::cout << "norm(costFunctor(X * 3)) = " << costs.norm() << std::endl;
                    costFunctor(X * 4, costs);
                    std::cout << "norm(costFunctor(X * 4)) = " << costs.norm() << std::endl;
                }

                X.normalize();
                if (callback){
                    InstallVariables(mg, controls, uh2varStartPosition, X, vars);
                    if (!callback(vars))
                        break;
                }
            } while (status == Eigen::LevenbergMarquardtSpace::Running);


            // install X
            X.normalize();
            InstallVariables(mg, controls, uh2varStartPosition, X, vars);

        }



        double MedianCenterDepth(const RLGraph & mg, const RLGraphControls & controls,
            const RLGraphVars & vars){
            for (auto & c : mg.components<RegionData>()){
                if (!controls[c.topo.hd].used)
                    continue;
                for (auto & v : vars[c.topo.hd].variables){
                    assert(!IsInfOrNaN(v));
                }
            }
            for (auto & c : mg.components<LineData>()){
                if (!controls[c.topo.hd].used)
                    continue;
                for (auto & v : vars[c.topo.hd].variables){
                    assert(!IsInfOrNaN(v));
                }
            }

            std::vector<double> centerDepths;
            centerDepths.reserve(mg.internalComponents<RegionData>().size() + mg.internalComponents<LineData>().size());
            for (auto & c : mg.components<RegionData>()){
                if (!controls[c.topo.hd].used)
                    continue;
                double d = DepthAt(c.data.normalizedCenter, Instance(mg, controls, vars, c.topo.hd));
                if (!IsInfOrNaN(d)){
                    centerDepths.push_back(d);
                }
                //assert(!IsInfOrNaN(centerDepths.back()));
            }
            for (auto & c : mg.components<LineData>()){
                if (!controls[c.topo.hd].used)
                    continue;
                double d = DepthAt(normalize(c.data.line.center()), Instance(mg, controls, vars, c.topo.hd));
                if (!IsInfOrNaN(d)){
                    centerDepths.push_back(d);
                }
                //assert(!IsInfOrNaN(centerDepths.back()));
            }
            std::nth_element(centerDepths.begin(), centerDepths.begin() + centerDepths.size() / 2, centerDepths.end());
            double medianCenterDepth = centerDepths[centerDepths.size() / 2];
            return medianCenterDepth;
        }


        double Score(const RLGraph & mg, const RLGraphControls & controls,
            const RLGraphVars & vars){

            // manhattan fitness
            double sumOfComponentWeightedFitness = 0.0;
            double sumOfComponentWeights = 0.0;

            // lines
            double maxLineSpanAngle = 0.0;
            for (auto & l : mg.components<LineData>()){
                if (!controls[l.topo.hd].used)
                    continue;
                double lineSpanAngle = AngleBetweenDirections(l.data.line.first, l.data.line.second);
                if (lineSpanAngle > maxLineSpanAngle)
                    maxLineSpanAngle = lineSpanAngle;
            }

            HandledTable<LineHandle, Line3> lines(mg.internalComponents<LineData>().size());
            HandledTable<RegionHandle, Plane3> planes(mg.internalComponents<RegionData>().size());

            for (auto & l : mg.components<LineData>()){
                if (!controls[l.topo.hd].used)
                    continue;
                lines[l.topo.hd] = Instance(mg, controls, vars, l.topo.hd);

                double fitness = 0.0;
                double weight = 0.0;
                double lineSpanAngle = AngleBetweenDirections(l.data.line.first, l.data.line.second);
                auto & prop = controls[l.topo.hd];
                if (prop.orientationClaz >= 0){
                    fitness = 1.0;
                    weight = lineSpanAngle / maxLineSpanAngle;
                }
                else{
                    static const double angleThreshold = DegreesToRadians(10);
                    double minAngle = angleThreshold;
                    int bestVPId = -1;
                    for (int i = 0; i < controls.vanishingPoints.size(); i++){
                        double angle = AngleBetweenUndirectedVectors(controls.vanishingPoints[i], lines[l.topo.hd].direction());
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
                if (!controls[r.topo.hd].used)
                    continue;
                if (r.data.area > maxRegionArea){
                    maxRegionArea = r.data.area;
                }
            }

            for (auto & r : mg.components<RegionData>()){
                if (!controls[r.topo.hd].used)
                    continue;
                planes[r.topo.hd] = Instance(mg, controls, vars, r.topo.hd);

                double fitness = 0.0;
                double weight = 1.0;
                auto & prop = controls[r.topo.hd];
                if (prop.orientationClaz >= 0){
                    fitness = 1.0;
                    weight = r.data.area / maxRegionArea;
                }
                else{
                    static const double angleThreshold = DegreesToRadians(10);
                    double minAngle = angleThreshold;
                    int bestVPId = -1;
                    for (int i = 0; i < controls.vanishingPoints.size(); i++){
                        double angle = AngleBetweenUndirectedVectors(controls.vanishingPoints[i], planes[r.topo.hd].normal);
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
                if (!controls[c.topo.hd].used){
                    sumOfNotUsedConstraintWeights += ElementsNum(c.data.normalizedSampledPoints) * typeWeight;
                    continue;
                }

                auto & plane1 = planes[c.topo.component<0>()];
                auto & plane2 = planes[c.topo.component<1>()];
                for (auto & ss : c.data.normalizedSampledPoints){
                    for (auto & s : ss){
                        double d1 = DepthAt(s, plane1);
                        double d2 = DepthAt(s, plane2);
                        if (IsInfOrNaN(d1) || IsInfOrNaN(d2)){
                            std::cout << "nan/inf plane!" << std::endl;
                            continue;
                        }
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
                if (!controls[c.topo.hd].used){
                    sumOfNotUsedConstraintWeights += /* c.data.junctionWeight / maxJunctionWeight **/ typeWeight;
                    continue;
                }

                auto & line1 = lines[c.topo.component<0>()];
                auto & line2 = lines[c.topo.component<1>()];
                double d1 = DepthAt(c.data.normalizedRelationCenter, line1);
                double d2 = DepthAt(c.data.normalizedRelationCenter, line2);
                if (IsInfOrNaN(d1) || IsInfOrNaN(d2)){
                    std::cout << "nan/inf line!" << std::endl;
                    continue;
                }

                double fitness = Gaussian(abs(d1 - d2) / (d1 + d2), 0.1);
                double weight = typeWeight;

                sumOfConstraintWeightedFitness += (fitness * weight);
                sumOfConstraintWeights += weight;
            }
            for (auto & c : mg.constraints<RegionLineConnectionData>()){
                static const double typeWeight = 1.0;
                if (!controls[c.topo.hd].used){
                    sumOfNotUsedConstraintWeights += ElementsNum(c.data.normalizedAnchors) * typeWeight;
                    continue;
                }

                auto & plane = planes[c.topo.component<0>()];
                auto & line = lines[c.topo.component<1>()];
                for (auto & a : c.data.normalizedAnchors){
                    double d1 = DepthAt(a, plane);
                    double d2 = DepthAt(a, line);
                    if (IsInfOrNaN(d1) || IsInfOrNaN(d2)){
                        std::cout << "nan/inf regionline!" << std::endl;
                        continue;
                    }
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





        void NormalizeVariables(const RLGraph & mg, const RLGraphControls & controls,
            RLGraphVars & vars){
            double medianCenterDepth = MedianCenterDepth(mg, controls, vars);
            
            SetClock();

            std::cout << "median center depth: " << medianCenterDepth << std::endl;

            assert(!IsInfOrNaN(medianCenterDepth));

            // normalize variables
            for (auto & c : mg.components<RegionData>()){
                if (!controls[c.topo.hd].used)
                    continue;
                for (int i = 0; i < vars[c.topo.hd].variables.size(); i++){
                    vars[c.topo.hd].variables[i] *= medianCenterDepth;
                }
            }
            for (auto & c : mg.components<LineData>()){
                if (!controls[c.topo.hd].used)
                    continue;
                for (int i = 0; i < vars[c.topo.hd].variables.size(); i++){
                    vars[c.topo.hd].variables[i] *= medianCenterDepth;
                }
            }
        }










        namespace {
            
            // returns false if confliction occurs
            bool MakeRegionPlaneUsable(RegionHandle rh, bool usable, RLGraphControls & controls) {
                auto & p = controls[rh];
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
            bool MakeRegionPlaneToward(RegionHandle rh, int normalVPId, RLGraphControls & controls){
                auto & p = controls[rh];
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
                auto & dir = controls.vanishingPoints[p.orientationNotClaz];
                if (IsFuzzyPerpendicular(controls.vanishingPoints[normalVPId], dir)){
                    p.orientationClaz = normalVPId;
                    p.orientationNotClaz = -1;
                    return true;
                }
                return false;
            }

            // returns false if confliction occurs
            bool MakeRegionPlaneAlsoAlong(RegionHandle rh, int alongVPId, RLGraphControls & controls){
                auto & p = controls[rh];
                if (!p.used)
                    return true;
                assert(alongVPId != -1);
                auto & dir = controls.vanishingPoints[alongVPId];
                if (p.orientationClaz != -1){
                    auto & normal = controls.vanishingPoints[p.orientationClaz];
                    return IsFuzzyPerpendicular(normal, dir);
                }
                if (p.orientationNotClaz == -1){
                    p.orientationNotClaz = alongVPId;
                    return true;
                }
                if (p.orientationNotClaz == alongVPId)
                    return true;

                auto newNormal = dir.cross(controls.vanishingPoints[p.orientationNotClaz]);
                double minAngle = M_PI;
                for (int i = 0; i < controls.vanishingPoints.size(); i++){
                    double angle = AngleBetweenUndirectedVectors(controls.vanishingPoints[i], newNormal);
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

        void AttachPrincipleDirectionConstraints(const RLGraph & mg, RLGraphControls & controls,
            double rangeAngle, bool avoidLineConflictions){

            SetClock();

            // find peaky regions
            std::vector<std::vector<RegionHandle>> peakyRegionHandles(controls.vanishingPoints.size());
            for (auto & r : mg.components<RegionData>()){
                auto h = r.topo.hd;
                if (!controls[h].used)
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
                for (auto & vp : controls.vanishingPoints){
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
                for (int i = 0; i < controls.vanishingPoints.size(); i++){
                    auto p1 = ToPixelLoc(ppc.screenProjection(controls.vanishingPoints[i]));
                    auto p2 = ToPixelLoc(ppc.screenProjection(-controls.vanishingPoints[i]));
                    
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
       
            for (int i = 0; i < controls.vanishingPoints.size(); i++){
                auto & vp = controls.vanishingPoints[i];
                auto & rhs = peakyRegionHandles[i];
                for (auto rh : rhs){
                    if (avoidLineConflictions){
                        bool hasLineConflictions = false;
                        for (auto conh : mg.topo(rh).constraints<RegionLineConnectionData>()){
                            if (mg.data(conh).normalizedAnchors.size() < regionLineRelatedThreshold)
                                continue;
                            LineHandle lh = mg.topo(conh).component<1>();
                            auto & lprop = controls[lh];
                            if (lprop.orientationClaz == i){
                                hasLineConflictions = true;
                                break;
                            }
                        }
                        if (!hasLineConflictions)
                            MakeRegionPlaneToward(rh, i, controls);
                    }
                    else {
                        MakeRegionPlaneToward(rh, i, controls);
                        for (auto conh : mg.topo(rh).constraints<RegionLineConnectionData>()){
                            if (mg.data(conh).normalizedAnchors.size() < regionLineRelatedThreshold)
                                continue;
                            LineHandle lh = mg.topo(conh).component<1>();
                            auto & lprop = controls[lh];
                            if (lprop.orientationClaz == i){
                                lprop.orientationClaz = -1;
                            }
                        }
                    }
                }
            }

            controls.disableAllInvalidConstraints(mg);

        }


        void AttachWallConstriants(const RLGraph & mg, RLGraphControls & controls,
            double rangeAngle, const Vec3 & verticalSeed){

            SetClock();

            int vertVPId = -1;
            double minAngle = M_PI;
            for (int i = 0; i < controls.vanishingPoints.size(); i++){
                double angle = AngleBetweenUndirectedVectors(verticalSeed, controls.vanishingPoints[i]);
                if (angle < minAngle){
                    minAngle = angle;
                    vertVPId = i;
                }
            }
            assert(vertVPId != -1);
            auto & vertical = controls.vanishingPoints[vertVPId];

            std::vector<RegionHandle> horizontalRegionHandles;
            for (auto & r : mg.components<RegionData>()){
                auto h = r.topo.hd;
                if (!controls[h].used)
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
                MakeRegionPlaneAlsoAlong(h, vertVPId, controls);
            }

            controls.disableAllInvalidConstraints(mg);
        }



        void AttachFloorAndCeilingConstraints(const RLGraph & mg, 
            RLGraphControls & controls, const RLGraphVars & vars,
            double eyeHeightRatioLowerBound, double eyeHeightRatioUpperBound,
            double angleThreshold, const Vec3 & verticalSeed){

            SetClock();

            IMPROVABLE_HERE("we need a better way!");

            int vertVPId = -1;
            double minAngle = M_PI;
            for (int i = 0; i < controls.vanishingPoints.size(); i++){
                double angle = AngleBetweenUndirectedVectors(verticalSeed, controls.vanishingPoints[i]);
                if (angle < minAngle){
                    minAngle = angle;
                    vertVPId = i;
                }
            }
            assert(vertVPId != -1);
            auto & vertical = controls.vanishingPoints[vertVPId];

            int horizVPId = -1;
            double maxAngle = 0;
            for (int i = 0; i < controls.vanishingPoints.size(); i++){
                double angle = AngleBetweenUndirectedVectors(vertical, controls.vanishingPoints[i]);
                if (angle > maxAngle){
                    maxAngle = angle;
                    horizVPId = i;
                }
            }
            assert(horizVPId != -1);
            auto xDir = normalize(controls.vanishingPoints[horizVPId]);
            auto yDir = normalize(vertical.cross(xDir));

            // compute horizontal bounding range of may-be ceilings and floors
            double xmin = std::numeric_limits<double>::max(), ymin = std::numeric_limits<double>::max();
            double xmax = std::numeric_limits<double>::lowest(), ymax = std::numeric_limits<double>::lowest();

            std::vector<Scored<RegionHandle>> maybeCeilinsOrFloors[2];

            auto regionPlanes = MakeHandledTableForComponents<Plane3, RegionData>(mg);
            for (auto & r : mg.components<RegionData>()){
                if (!controls[r.topo.hd].used)
                    continue;

                auto plane = Instance(mg, controls, vars, r.topo.hd);
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
            double medianCenterDepth = MedianCenterDepth(mg, controls, vars);
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
                const auto & prop = controls[r.topo.hd];
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
                        MakeRegionPlaneToward(regionHandle, vertVPId, controls);
                    }
                }
            }

            controls.disableAllInvalidConstraints(mg);
        }



        //void AttachGeometricContextConstraints(const RLGraph & mg, RLGraphControls & controls,
        //    const std::vector<GeometricContextEstimator::Feature> & perspectiveGCs,
        //    const std::vector<PerspectiveCamera> & gcCameras,
        //    int shrinkRegionOrientationIteration, bool considerGCVerticalConstraint){

        //    auto orientationVotes = mg.createComponentTable<RegionData, std::unordered_map<OrientationHint, double>>(
        //        std::unordered_map<OrientationHint, double>((size_t)OrientationHint::Count));

        //    assert(perspectiveGCs.size() == gcCameras.size());

        //    // region views
        //    auto regionMaskViews = mg.createComponentTable<RegionData, View<PartialPanoramicCamera, Imageub>>();
        //    for (auto & r : mg.components<RegionData>()){
        //        auto h = r.topo.hd;
        //        auto & contours = r.data.normalizedContours;
        //        double radiusAngle = 0.0;
        //        for (auto & cs : r.data.normalizedContours){
        //            for (auto & c : cs){
        //                double angle = AngleBetweenDirections(r.data.normalizedCenter, c);
        //                if (angle > radiusAngle){
        //                    radiusAngle = angle;
        //                }
        //            }
        //        }
        //        float ppcFocal = 100.0f;
        //        int ppcSize = 2 * radiusAngle * ppcFocal;
        //        Vec3 x;
        //        std::tie(x, std::ignore) = ProposeXYDirectionsFromZDirection(r.data.normalizedCenter);
        //        PartialPanoramicCamera ppc(ppcSize, ppcSize, ppcFocal, Point3(0, 0, 0), r.data.normalizedCenter, x);
        //        Imageub mask = Imageub::zeros(ppc.screenSize());

        //        // project contours to ppc
        //        std::vector<std::vector<Point2i>> contourProjs(contours.size());
        //        for (int k = 0; k < contours.size(); k++){
        //            auto & contourProj = contourProjs[k];
        //            contourProj.reserve(contours[k].size());
        //            for (auto & d : contours[k]){
        //                contourProj.push_back(vec_cast<int>(ppc.screenProjection(d)));
        //            }
        //        }
        //        cv::fillPoly(mask, contourProjs, (uint8_t)1);
        //        regionMaskViews[h].camera = ppc;
        //        regionMaskViews[h].image = mask;
        //    }

        //    auto regionAreas = mg.createComponentTable<RegionData, double>(0.0);
        //    for (auto & r : mg.components<RegionData>()){
        //        auto h = r.topo.hd;
        //        auto & ppMask = regionMaskViews[h];
        //        double area = cv::sum(ppMask.image).val[0];
        //        regionAreas[h] = area;
        //        for (int i = 0; i < perspectiveGCs.size(); i++){
        //            auto sampler = MakeCameraSampler(ppMask.camera, gcCameras[i]);
        //            auto occupation = sampler(Imageub::ones(gcCameras[i].screenSize()), cv::BORDER_CONSTANT, 0.0);
        //            double areaOccupied = cv::sum(ppMask.image & occupation).val[0];
        //            double ratio = areaOccupied / area;
        //            assert(ratio <= 1.01);
        //            GeometricContextLabel maxLabel = GeometricContextLabel::None;
        //            double maxScore = 0.0;
        //            for (auto & gcc : perspectiveGCs[i]){
        //                auto label = gcc.first;
        //                Imaged ppGC = sampler(gcc.second);
        //                double score = cv::mean(ppGC, ppMask.image).val[0];
        //                if (score > maxScore){
        //                    maxScore = score;
        //                    maxLabel = label;
        //                }
        //            }
        //            if (maxLabel != GeometricContextLabel::None){
        //                auto & v = orientationVotes[r.topo.hd][ToOrientationHint(maxLabel, considerGCVerticalConstraint)];
        //                v = std::max(v, ratio);
        //            }
        //        }
        //    }

        //    // optimize
        //    GCoptimizationGeneralGraph graph(mg.internalComponents<RegionData>().size(), (int)OrientationHint::Count);

        //    double maxBoundaryLength = 0;
        //    for (auto & b : mg.constraints<RegionBoundaryData>()){
        //        maxBoundaryLength = std::max(maxBoundaryLength, b.data.length);
        //    }
        //    for (auto & b : mg.constraints<RegionBoundaryData>()){
        //        graph.setNeighbors(b.topo.component<0>().id, b.topo.component<1>().id,
        //            100 * b.data.length / maxBoundaryLength);
        //    }
        //    for (auto & r : mg.components<RegionData>()){
        //        auto & orientationVote = orientationVotes[r.topo.hd];
        //        // normalize votes
        //        double votesSum = 0.0;
        //        for (auto & v : orientationVote){
        //            votesSum += v.second;
        //        }
        //        assert(votesSum >= 0.0);
        //        if (votesSum > 0.0){
        //            for (auto & v : orientationVote){
        //                assert(v.second >= 0.0);
        //                v.second /= votesSum;
        //            }
        //        }
        //        for (int label = 0; label < (int)OrientationHint::Count; label++){
        //            double vote = Contains(orientationVote, (OrientationHint)label) ?
        //                orientationVote.at((OrientationHint)label) : 0.0;
        //            graph.setDataCost(r.topo.hd.id, label,
        //                votesSum == 0.0 ? 0 : (10000 * (1.0 - vote)));
        //        }
        //    }
        //    for (int label1 = 0; label1 < (int)OrientationHint::Count; label1++){
        //        for (int label2 = 0; label2 < (int)OrientationHint::Count; label2++){
        //            if (label1 == label2){
        //                graph.setSmoothCost(label1, label2, 0);
        //            }
        //            else{
        //                graph.setSmoothCost(label1, label2, 1);
        //            }
        //        }
        //    }

        //    graph.expansion();
        //    graph.swap();

        //    // get the most vertical vp id
        //    int vVPId = -1;
        //    double angleToVert = std::numeric_limits<double>::max();
        //    for (int i = 0; i < controls.vanishingPoints.size(); i++){
        //        double a = AngleBetweenUndirectedVectors(controls.vanishingPoints[i], Vec3(0, 0, 1));
        //        if (a < angleToVert){
        //            vVPId = i;
        //            angleToVert = a;
        //        }
        //    }
        //    assert(vVPId != -1);

        //    // disorient outsided oriented gc labels
        //    auto regionOrientations = mg.createComponentTable<RegionData, OrientationHint>(OrientationHint::Void);
        //    for (auto & r : mg.components<RegionData>()){
        //        regionOrientations[r.topo.hd] = (OrientationHint)(graph.whatLabel(r.topo.hd.id));
        //    }

        //    for (int i = 0; i < shrinkRegionOrientationIteration; i++){
        //        auto shrinked = regionOrientations;
        //        for (auto & b : mg.constraints<RegionBoundaryData>()){
        //            auto l1 = regionOrientations[b.topo.component<0>()];
        //            auto l2 = regionOrientations[b.topo.component<1>()];
        //            if (l1 == OrientationHint::Horizontal && l2 != OrientationHint::Horizontal){
        //                shrinked[b.topo.component<0>()] = OrientationHint::OtherPlanar;
        //            }
        //            else if (l1 != OrientationHint::Horizontal && l2 == OrientationHint::Horizontal){
        //                shrinked[b.topo.component<1>()] = OrientationHint::OtherPlanar;
        //            }
        //            if (l1 == OrientationHint::Vertical && l2 != OrientationHint::Vertical){
        //                shrinked[b.topo.component<0>()] = OrientationHint::OtherPlanar;
        //            }
        //            else if (l1 != OrientationHint::Vertical && l2 == OrientationHint::Vertical){
        //                shrinked[b.topo.component<1>()] = OrientationHint::OtherPlanar;
        //            }
        //        }
        //        regionOrientations = std::move(shrinked);
        //    }

        //    // install gc labels
        //    for (auto & r : mg.components<RegionData>()){
        //        OrientationHint oh = regionOrientations[r.topo.hd];
        //        switch (oh) {
        //        case OrientationHint::Void:
        //            /*controls[r.topo.hd].used = false;
        //            controls[r.topo.hd].orientationClaz = -1;
        //            controls[r.topo.hd].orientationNotClaz = -1;*/
        //            MakeRegionPlaneUsable(r.topo.hd, false, controls);
        //            break;
        //        case OrientationHint::Horizontal:
        //            /*controls[r.topo.hd].used = true;
        //            controls[r.topo.hd].orientationClaz = vVPId;
        //            controls[r.topo.hd].orientationNotClaz = -1;*/
        //            MakeRegionPlaneToward(r.topo.hd, vVPId, controls);
        //            break;
        //        case OrientationHint::Vertical:
        //            //controls[r.topo.hd].used = true;
        //            //controls[r.topo.hd].orientationClaz = -1;
        //            //controls[r.topo.hd].orientationNotClaz = vVPId;
        //            MakeRegionPlaneAlsoAlong(r.topo.hd, vVPId, controls);
        //            break;
        //        case OrientationHint::OtherPlanar:
        //            /*controls[r.topo.hd].used = true;
        //            controls[r.topo.hd].orientationClaz = -1;
        //            controls[r.topo.hd].orientationNotClaz = -1;*/
        //            //MakeRegionPlaneUsable(r.topo.hd, true, controls);
        //            break;
        //        case OrientationHint::NonPlanar:
        //            /*controls[r.topo.hd].used = false;
        //            controls[r.topo.hd].orientationClaz = -1;
        //            controls[r.topo.hd].orientationNotClaz = -1;*/
        //            MakeRegionPlaneUsable(r.topo.hd, false, controls);
        //            break;
        //        default:
        //            assert(0);
        //            break;
        //        }
        //    }

        //    /*
        //    // detect big aligned rectangular regions and orient them
        //    double regionAreaSum = std::accumulate(regionAreas.begin(), regionAreas.end(), 0.0);
        //    assert(regionAreaSum > 0);

        //    static const double dotThreshold = 0.001;
        //    static const double spanAngleThreshold = DegreesToRadians(10);
        //    static const double wholeDotThreshold = 0.01;
        //    static const double wholeSpanAngleThreshold = DegreesToRadians(10);

        //    HandledTable<RegionBoundaryHandle, std::vector<double>> boundaryMaxSpanAnglesForVPs(mg.internalConstraints<RegionBoundaryData>().size());
        //    for (auto & b : mg.constraints<RegionBoundaryData>()){
        //    auto & edges = b.data.normalizedEdges;
        //    std::vector<double> maxSpanAngleForVPs(props.vanishingPoints.size(), 0.0);
        //    if (b.topo.hd.id == 849){
        //    std::cout << std::endl;
        //    }

        //    for (auto & edge : edges){
        //    std::vector<std::vector<bool>> edgeVPFlags(props.vanishingPoints.size(),
        //    std::vector<bool>(edge.size() - 1, false));
        //    for (int i = 0; i < edge.size() - 1; i++){
        //    Vec3 n = normalize(edge[i].cross(edge[i + 1]));
        //    for (int k = 0; k < props.vanishingPoints.size(); k++){
        //    auto vp = normalize(props.vanishingPoints[k]);
        //    double dotv = abs(vp.dot(n));
        //    if (dotv <= dotThreshold){
        //    edgeVPFlags[k][i] = true;
        //    }
        //    }
        //    }
        //    // detect aligned longest line spanAngles from edge
        //    for (int k = 0; k < props.vanishingPoints.size(); k++){
        //    auto vp = normalize(props.vanishingPoints[k]);
        //    auto & edgeFlags = edgeVPFlags[k];
        //    int lastHead = -1, lastTail = -1;
        //    bool inChain = false;

        //    double & maxSpanAngle = maxSpanAngleForVPs[k];
        //    for (int j = 0; j <= edgeFlags.size(); j++){
        //    if (!inChain && j < edgeFlags.size() && edgeFlags[j]){
        //    lastHead = j;
        //    inChain = true;
        //    }
        //    else if (inChain && (j == edgeFlags.size() || !edgeFlags[j])){
        //    lastTail = j;
        //    inChain = false;
        //    // examine current chain
        //    assert(lastHead != -1);
        //    // compute full span angle
        //    double spanAngle = 0.0;
        //    for (int i = lastHead; i < lastTail; i++){
        //    spanAngle += AngleBetweenDirections(edge[i], edge[i + 1]);
        //    }
        //    if (spanAngle < spanAngleThreshold){
        //    continue;
        //    }
        //    // fit line
        //    const Vec3 & midCorner = edge[(lastHead + lastTail) / 2];
        //    Vec3 commonNormal = normalize(midCorner.cross(vp));
        //    std::vector<double> dotsToCommonNormal(lastTail - lastHead + 1, 0.0);
        //    for (int i = lastHead; i <= lastTail; i++){
        //    dotsToCommonNormal[i - lastHead] = abs(edge[i].dot(commonNormal));
        //    }
        //    if (*std::max_element(dotsToCommonNormal.begin(), dotsToCommonNormal.end()) <= wholeDotThreshold){
        //    // acceptable!
        //    if (spanAngle > maxSpanAngle){
        //    maxSpanAngle = spanAngle;
        //    }
        //    }
        //    }
        //    }
        //    }
        //    }
        //    if ((b.topo.component<0>().id == 355 || b.topo.component<1>().id == 355) && std::accumulate(maxSpanAngleForVPs.begin(), maxSpanAngleForVPs.end(), 0.0) > 0){
        //    std::cout << std::endl;
        //    }
        //    boundaryMaxSpanAnglesForVPs[b.topo.hd] = std::move(maxSpanAngleForVPs);
        //    }
        //    for (auto & r : mg.components<RegionData>()){
        //    if (!props[r.topo.hd].used)
        //    continue;
        //    if (props[r.topo.hd].orientationClaz != -1 || props[r.topo.hd].orientationNotClaz != -1)
        //    continue;

        //    if (r.topo.hd.id == 355){
        //    std::cout << std::endl;
        //    }

        //    double area = regionAreas[r.topo.hd];
        //    if (area / regionAreaSum < 0.005){
        //    continue;
        //    }
        //    std::vector<double> spanAngleSumsForVPs(props.vanishingPoints.size(), 0.0);
        //    for (auto & h : r.topo.constraints<RegionBoundaryData>()){
        //    for (int i = 0; i < props.vanishingPoints.size(); i++){
        //    spanAngleSumsForVPs[i] += boundaryMaxSpanAnglesForVPs[h][i];
        //    }
        //    }
        //    auto maxIter = std::max_element(spanAngleSumsForVPs.begin(), spanAngleSumsForVPs.end());
        //    int firstMaxVPId = std::distance(spanAngleSumsForVPs.begin(), maxIter);
        //    double firstMaxSpanAngle = *maxIter;
        //    *maxIter = -1.0;
        //    maxIter = std::max_element(spanAngleSumsForVPs.begin(), spanAngleSumsForVPs.end());
        //    int secondMaxVPId = std::distance(spanAngleSumsForVPs.begin(), maxIter);
        //    double secondMaxSpanAngle = *maxIter;
        //    if (secondMaxSpanAngle >= wholeSpanAngleThreshold){
        //    assert(firstMaxVPId != secondMaxVPId);
        //    Vec3 normal = props.vanishingPoints[firstMaxVPId].cross(props.vanishingPoints[secondMaxVPId]);
        //    double minAngle = 0.1;
        //    int bestMatchedVPId = -1;
        //    for (int i = 0; i < props.vanishingPoints.size(); i++){
        //    double angle = AngleBetweenUndirectedVectors(normal, props.vanishingPoints[i]);
        //    if (angle < minAngle){
        //    minAngle = angle;
        //    bestMatchedVPId = i;
        //    }
        //    }
        //    if (bestMatchedVPId != -1){
        //    std::cout << "vpid:::: " << bestMatchedVPId << std::endl;
        //    props[r.topo.hd].orientationClaz = bestMatchedVPId;
        //    props[r.topo.hd].orientationNotClaz = -1;
        //    }
        //    }
        //    else if (firstMaxSpanAngle >= wholeSpanAngleThreshold){
        //    std::cout << "vpid-along:::: " << firstMaxVPId << std::endl;
        //    props[r.topo.hd].orientationClaz = -1;
        //    props[r.topo.hd].orientationNotClaz = firstMaxVPId;
        //    }
        //    }*/

        //    //for (auto & l : mg.internalComponents<LineData>()){
        //    //    props[l.topo.hd].used = true;
        //    //    props[l.topo.hd].orientationClaz = l.data.initialClaz;
        //    //    props[l.topo.hd].orientationNotClaz = -1;
        //    //}
        //    
        //    controls.disableAllInvalidConstraints(mg);

        //}







        void LooseOrientationConstraintsOnComponents(const RLGraph & mg,
            RLGraphControls & controls, const RLGraphVars & vars,
            double linesLoosableRatio, double regionsLoosableRatio, double distThresRatio){

            SetClock();

            if (linesLoosableRatio <= 0.0 && regionsLoosableRatio <= 0.0)
                return;
            
            double medianCenterDepth = MedianCenterDepth(mg, controls, vars);
            double distThres = medianCenterDepth * distThresRatio;
            double nnSearchRange = distThres;

            bool directApproach = true;
            if (directApproach){

                // collect current instances
                HandledTable<LineHandle, Line3> lines(mg.internalComponents<LineData>().size());
                HandledTable<RegionHandle, Plane3> planes(mg.internalComponents<RegionData>().size());

                for (auto & r : mg.components<RegionData>()){
                    if (!controls[r.topo.hd].used)
                        continue;
                    planes[r.topo.hd] = Instance(mg, controls, vars, r.topo.hd);
                }
                for (auto & l : mg.components<LineData>()){
                    if (!controls[l.topo.hd].used)
                        continue;
                    lines[l.topo.hd] = Instance(mg, controls, vars, l.topo.hd);
                }

                if (linesLoosableRatio > 0.0){
                    std::vector<LineHandle> usableLineHandles;
                    usableLineHandles.reserve(mg.internalComponents<LineData>().size());
                    auto distanceToNearestRegions = mg.createComponentTable<LineData>(0.0);
                    for (auto & l : mg.components<LineData>()){
                        if (!controls[l.topo.hd].used)
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
                        auto & lp = controls[usableLineHandles[i]];
                        assert(lp.used);
                        if (lp.orientationClaz >= 0){
                            lp.orientationClaz = -1;
                        }
                        assert(lp.orientationNotClaz == -1);
                        if (lp.orientationClaz == -1){
                            lp.used = false; // disable this line!!!!!!!!!
                        }
                    }
                }

                if (regionsLoosableRatio > 0.0){
                    std::vector<RegionHandle> usableRegionHandles;
                    usableRegionHandles.reserve(mg.internalComponents<RegionData>().size());
                    auto maxDistanceToNearestComponents = mg.createComponentTable<RegionData>(0.0);

                    for (auto & r : mg.components<RegionData>()){
                        if (!controls[r.topo.hd].used)
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
                        auto & rp = controls[usableRegionHandles[i]];
                        assert(rp.used);
                        if (rp.orientationClaz == -1 && rp.orientationNotClaz == -1){
                            rp.used = false;
                        }
                        else{
                            rp.orientationClaz = rp.orientationNotClaz = -1;
                        }
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
                    if (!controls[r.topo.hd].used)
                        continue;
                    auto plane = Instance(mg, controls, vars, r.topo.hd);
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
                    if (!controls[l.topo.hd].used)
                        continue;
                    auto line = Instance(mg, controls, vars, l.topo.hd);
                    lines[l.topo.hd] = line;
                    keyPointsRTree.insert(ComponentKeyPoint(line.first, l.topo.hd));
                    keyPointsRTree.insert(ComponentKeyPoint(line.second, l.topo.hd));
                }

                // find most isolated lines
                if (linesLoosableRatio > 0.0){

                    auto lineMaxCornerNNDistances = mg.createComponentTable<LineData, double>();
                    std::vector<LineHandle> usableLineHandles;
                    for (auto & l : mg.components<LineData>()){
                        if (!controls[l.topo.hd].used)
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
                        auto & lp = controls[usableLineHandles[i]];
                        assert(lp.used);
                        if (lp.orientationClaz >= 0){
                            lp.orientationClaz = -1;
                        }
                        assert(lp.orientationNotClaz == -1);
                        if (lp.orientationClaz == -1){
                            lp.used = false; // disable this line!!!!!!!!!
                        }
                    }
                }

                // find most isolated regions
                if (regionsLoosableRatio > 0.0){

                    auto regionMaxCornerNNDistances = mg.createComponentTable<RegionData, double>();
                    std::vector<RegionHandle> usableRegionHandles;
                    for (auto & r : mg.components<RegionData>()){
                        if (!controls[r.topo.hd].used)
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
                        auto & rp = controls[usableRegionHandles[i]];
                        assert(rp.used);
                        rp.orientationClaz = rp.orientationNotClaz = -1;
                        //if (lp.orientationClaz == -1){
                        //    lp.used = false; // disable this line!!!!!!!!!
                        //}
                    }
                }


            }

            controls.disableAllInvalidConstraints(mg);

        }












    }
}


