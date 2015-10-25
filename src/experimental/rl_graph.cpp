
//
#include "../core/algorithms.hpp"
#include "../core/containers.hpp"
#include "../core/utility.hpp"
#include "../core/clock.hpp"
#include "../core/homo_graph.hpp"

#include "rl_graph_vis.hpp"
#include "rl_graph.hpp"
#include "rl_graph_occlusion.hpp"


namespace pano {
    namespace experimental {

      
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
                inline bool operator ()(const Pixel & a, const Pixel & b) const {
                    if (a.x != b.x)
                        return a.x < b.x;
                    return a.y < b.y;
                }
            };


            inline Point2 ToPoint2(const Pixel & p) {
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

            inline Pixel ToPixel(const Point2 & p) {
                return Pixel(static_cast<int>(p[0]), static_cast<int>(p[1]));
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
                Y += v(0, s) * v(1, s) * v(2, s) * DiracDelta(v(0, 1 - s) + v(1, 1 - s) + v(2, 1 - s), 1e-4);
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
                                auto conCenter = HPointFromVector(GetCoeffs(linei.ray())
                                    .cross(GetCoeffs(linej.ray()))).value();

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


                            auto p = ToPixel(lrd.relationCenter);
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
                vps2d[i] = cam.toScreenInHPoint(vps[i]);
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
                ld.line.first = normalize(cam.toSpace(ld2d.line.first));
                ld.line.second = normalize(cam.toSpace(ld2d.line.second));
                lh2dToLh[l2d.topo.hd] = mg.addComponent(std::move(ld));
                newLineHandles.insert(lh2dToLh[l2d.topo.hd]);
            }

            for (auto & r2d : graph2d.elements<1>()){
                auto & rd2d = r2d.data;
                LineRelationData rd;
                rd.type = rd2d.type;
                rd.junctionWeight = rd2d.junctionWeight;
                rd.normalizedRelationCenter = normalize(cam.toSpace(rd2d.relationCenter));
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



        std::vector<LineHandle> AppendLines2(RLGraph & mg, const std::vector<Classified<Line3>> & lines,
            const std::vector<Vec3> & vps,
            double intersectionAngleThreshold,
            double incidenceAngleAlongDirectionThreshold,
            double incidenceAngleVerticalDirectionThreshold) {

            if (lines.empty())
                return std::vector<LineHandle>();

            assert(mg.internalComponents<LineData>().empty() && mg.internalComponents<RegionData>().empty());

            std::vector<LineHandle> lhs;
            lhs.reserve(lines.size());
            mg.internalComponents<LineData>().reserve(lines.size());

            for (auto & line : lines) {
                LineData ld;
                ld.initialClaz = line.claz;
                ld.line = normalize(line.component);
                lhs.push_back(mg.addComponent(std::move(ld)));
            }

            // construct incidence/intersection relations
            auto & linesData = mg.internalComponents<LineData>();
            for (int i = 0; i < linesData.size(); i++) {
                auto & linei = linesData[i].data.line;
                Vec3 ni = normalize(linei.first.cross(linei.second));
                int clazi = linesData[i].data.initialClaz;

                for (int j = i + 1; j < linesData.size(); j++) {
                    auto & linej = linesData[j].data.line;
                    int clazj = linesData[j].data.initialClaz;
                    Vec3 nj = normalize(linej.first.cross(linej.second));

                    auto nearest = DistanceBetweenTwoLines(linei, linej);
                    double d = AngleBetweenDirections(nearest.second.first.position, nearest.second.second.position);

                    if (clazi == clazj && clazi >= 0) { // incidences for classified lines
                        auto conCenter = normalize(nearest.second.first.position + nearest.second.second.position);
                        auto conDir = (nearest.second.first.position - nearest.second.second.position);

                        auto & vp = vps[clazi];

                        if (AngleBetweenDirections(vp, conCenter) < intersectionAngleThreshold)
                            continue;

                        if (d < incidenceAngleAlongDirectionThreshold &&
                            AngleBetweenUndirectedVectors(ni, nj) < incidenceAngleVerticalDirectionThreshold) {
                            LineRelationData lrd;
                            lrd.type = LineRelationData::Type::Incidence;
                            lrd.normalizedRelationCenter = conCenter;
                            lrd.junctionWeight = 5.0;

                            if (HasValue(lrd.normalizedRelationCenter, IsInfOrNaN<double>))
                                continue;

                            mg.addConstraint(std::move(lrd), lhs[i], lhs[j]);
                        }

                    } else if (clazi != clazj && clazi >= 0 && clazj >= 0) { // intersections for classified lines
                        if (d < intersectionAngleThreshold) {
                            auto conCenter = normalize(ni.cross(nj));

                            if (DistanceFromPointToLine(conCenter, linei).first > intersectionAngleThreshold * 4 ||
                                DistanceFromPointToLine(conCenter, linej).first > intersectionAngleThreshold * 4)
                                continue;

                            LineRelationData lrd;
                            lrd.type = LineRelationData::Type::Intersection;
                            lrd.normalizedRelationCenter = conCenter;
                            lrd.junctionWeight = 3.0;

                            if (HasValue(lrd.normalizedRelationCenter, IsInfOrNaN<double>))
                                continue;

                            mg.addConstraint(std::move(lrd), lhs[i], lhs[j]);
                        }
                    }
                }
            }



            // compute junction weights
            static const double angleThreshold = M_PI / 32;
            static const double sigma = 0.1;

            enum LineVotingDirection : int {
                TowardsVanishingPoint = 0,
                TowardsOppositeOfVanishingPoint = 1
            };
            enum class JunctionType : int {
                L, T, Y, W, X
            };

            for (auto & lr : mg.internalConstraints<LineRelationData>()) {
                //std::cout << lr.topo.hd.id << std::endl;
                auto & lrd = lr.data;
                if (lrd.type == LineRelationData::Incidence) {
                    lrd.junctionWeight = IncidenceJunctionWeight(false);
                } else if (lrd.type == LineRelationData::Intersection) {
                    Mat<float, 3, 2> votingData;
                    std::fill(std::begin(votingData), std::end(votingData), 0);

                    for (auto & ld : mg.components<LineData>()) {
                        auto & line = ld.data.line;
                        int claz = ld.data.initialClaz;
                        if (claz == -1 || claz >= 3)
                            continue;

                        auto & vp = vps[claz];
                        Vec3 center = normalize(line.center());

                        Vec3 center2vp = normalize(center.cross(vp));
                        Vec3 center2pos = normalize(center.cross(lrd.normalizedRelationCenter));

                        double angle = AngleBetweenUndirectedVectors(center2vp, center2pos);
                        double angleSmall = angle > M_PI_2 ? (M_PI - angle) : angle;
                        if (IsInfOrNaN(angleSmall))
                            continue;

                        assert(angleSmall >= 0 && angleSmall <= M_PI_2);

                        double angleScore =
                            exp(-(angleSmall / angleThreshold) * (angleSmall / angleThreshold) / sigma / sigma / 2);

                        auto proj = ProjectionOfPointOnLine(lrd.normalizedRelationCenter, line);
                        double projRatio = BoundBetween(proj.ratio, 0.0, 1.0);

                        Vec3 lined = line.first.cross(line.second);
                        double lineSpanAngle = AngleBetweenDirections(line.first, line.second);
                        if (AngleBetweenDirections(center2vp, lined) < M_PI_2) { // first-second-vp
                            votingData(claz, TowardsVanishingPoint) += angleScore * lineSpanAngle * (1 - projRatio);
                            votingData(claz, TowardsOppositeOfVanishingPoint) += angleScore * lineSpanAngle * projRatio;
                        } else { // vp-first-second
                            votingData(claz, TowardsOppositeOfVanishingPoint) += angleScore * lineSpanAngle * (1 - projRatio);
                            votingData(claz, TowardsVanishingPoint) += angleScore * lineSpanAngle * projRatio;
                        }
                    }

                    lrd.junctionWeight = ComputeIntersectionJunctionWeightWithLinesVotes(
                        votingData);
                } else {
                    assert(false && "invalid branch!");
                }
            }
            
            return lhs;
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
                    std::vector<std::vector<Pixel>> contours;

                    cv::findContours(regionMask, contours, CV_RETR_EXTERNAL, CV_CHAIN_APPROX_SIMPLE); // CV_RETR_EXTERNAL: get only the outer contours
                    std::sort(contours.begin(), contours.end(),
                        [](const std::vector<Pixel> & ca, const std::vector<Pixel> & cb){return ca.size() > cb.size(); });

                    if (contours.size() >= 2){
                        //std::cout << "multiple contours for one region in perspective projection!";
                    }

                    auto iter = std::find_if(contours.begin(), contours.end(), [](const std::vector<Pixel> & c){return c.size() <= 2; });
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
                            rd.normalizedContours[k].push_back(normalize(cam.toSpace(p)));
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
                    std::vector<std::vector<Pixel>> contours;
                    cv::findContours(regionMask, contours, CV_RETR_EXTERNAL, CV_CHAIN_APPROX_SIMPLE); // CV_RETR_EXTERNAL: get only the outer contours
                    if (contours.empty()){
                        continue;
                    }

                    Vec3 centerDirection(0, 0, 0);
                    std::vector<Vec3> directions;
                    directions.reserve(ElementsNum(contours));
                    for (auto & cs : contours){
                        for (auto & c : cs){
                            directions.push_back(normalize(cam.toSpace(c)));
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
                        [](const std::vector<Pixel> & ca, const std::vector<Pixel> & cb){return ca.size() > cb.size(); });

                    if (contours.size() >= 2){
                        //std::cout << "multiple contours for one region in projection!";
                    }

                    auto iter = std::find_if(contours.begin(), contours.end(), [](const std::vector<Pixel> & c){return c.size() <= 2; });
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
                            rd.normalizedContours[k].push_back(normalize(sCam.toSpace(p)));
                            center += rd.normalizedContours[k].back();
                        }
                        std::vector<Point2f> contourf(contours[k].size());
                        for (int kk = 0; kk < contours[k].size(); kk++){
                            contourf[kk] = ecast<float>(contours[k][kk]);
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
                    Line2 line2(cam.toScreen(line.first), cam.toScreen(line.second));
                    Vec2 vertToLineDir = normalize(PerpendicularDirection(line2.direction())); // vertical to this line
                    double stepOnImage = std::max(cam.focal() * samplingStepAngleOnLine, 0.5);
                    if (SampleLineOnImage){
                        for (double stepLen = 0.0; stepLen <= line2.length(); stepLen += stepOnImage){
                            Point2 sampleP = line2.first + normalize(line2.second - line2.first) * stepLen;
                            Pixel originalP = ToPixel(sampleP);
                            std::set<int> connectedRegionIds;
                            for (int x = originalP.x - samplerSizeOnLine; x <= originalP.x + samplerSizeOnLine; x++){
                                for (int y = originalP.y - samplerSizeOnLine; y <= originalP.y + samplerSizeOnLine; y++){
                                    Pixel p(x, y);
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
                                    Pixel p = ToPixel(sampleP + vertToLineDir * dir);
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
                                rd.normalizedAnchors.push_back(normalize(cam.toSpace(sampleP)));
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
                            Point2 sampleP = cam.toScreen(sample);
                            Pixel originalP = ToPixel(sampleP);
                            // collect neighbors
                            std::set<int> connectedRegionIds;
                            for (int x = originalP.x - samplerSizeOnLine; x <= originalP.x + samplerSizeOnLine; x++){
                                for (int y = originalP.y - samplerSizeOnLine; y <= originalP.y + samplerSizeOnLine; y++){
                                    Pixel p(x, y);
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
                                    Pixel p = ToPixel(sampleP + vertToLineDir * dir);
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
                std::map<std::pair<int, int>, std::vector<std::vector<Pixel>>> boundaryEdges =
                    FindRegionBoundaries(segmentedRegions, samplerSizeOnBoundary, false);


                for (auto & bep : boundaryEdges) {
                    auto & rids = bep.first;
                    auto edges = bep.second;

                    if (regionHandles[rids.first].invalid() || regionHandles[rids.second].invalid())
                        continue;

                    auto rh = regionHandles[rids.first];
                    if (noBoundaryUnderLines){
                        std::vector<std::vector<Pixel>> filteredBoundaryPixels;
                        for (auto & e : edges){
                            std::vector<Pixel> filteredEdge;
                            for (Pixel p : e){
                                bool coveredByLine = false;
                                for (RegionLineConnectionHandle rlcon : mg.topo(rh).constraints<RegionLineConnectionData>()){
                                    LineHandle lh = mg.topo(rlcon).component<1>();
                                    const Line3 & line = mg.data(lh).line;
                                    Line2 line2(cam.toScreen(line.first), cam.toScreen(line.second));
                                    double d = DistanceFromPointToLine(ecast<double>(p), line2).first;
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
                            bd.normalizedEdges[k].push_back(normalize(cam.toSpace(p)));
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












        std::pair<std::vector<RegionHandle>, std::vector<RegionBoundaryHandle>> AppendRegions(RLGraph & mg, const Imagei & segmentedRegions,
            const std::vector<std::vector<Pixel>> & bndpixels, const std::vector<std::pair<int, int>> & bnd2segs,
            const PanoramicCamera & cam,
            double samplingStepAngleOnBoundary, double samplingStepAngleOnLine,
            int samplerSizeOnBoundary, int samplerSizeOnLine, bool noBoundaryUnderLines) {

            using CameraT = PanoramicCamera;

            auto regionHandles = CollectRegionsFromSegmentation(mg, segmentedRegions, cam);
            int regionNum = regionHandles.size();


            // add region-line connections
            std::map<std::pair<RegionHandle, LineHandle>, RegionLineConnectionData> regionLineConnections;

            for (auto & ld : mg.components<LineData>()) {
                auto & line = ld.data.line;
                Line2 line2(cam.toScreen(line.first), cam.toScreen(line.second));
                Vec2 vertToLineDir = normalize(PerpendicularDirection(line2.direction())); // vertical to this line
                double stepOnImage = std::max(cam.focal() * samplingStepAngleOnLine, 0.5);
                if (std::is_same<CameraT, PerspectiveCamera>::value) {
                    for (double stepLen = 0.0; stepLen <= line2.length(); stepLen += stepOnImage) {
                        Point2 sampleP = line2.first + normalize(line2.second - line2.first) * stepLen;
                        Pixel originalP = ToPixel(sampleP);
                        std::set<int> connectedRegionIds;
                        for (int x = originalP.x - samplerSizeOnLine; x <= originalP.x + samplerSizeOnLine; x++) {
                            for (int y = originalP.y - samplerSizeOnLine; y <= originalP.y + samplerSizeOnLine; y++) {
                                Pixel p(x, y);
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
                        for (int d = std::min(samplerSizeOnLine, 2); d <= samplerSizeOnLine * 2; d++) {
                            for (int dir : {-d, d}) {
                                Pixel p = ToPixel(sampleP + vertToLineDir * dir);
                                if (!Contains(segmentedRegions, p))
                                    continue;
                                int regionId = segmentedRegions(p);
                                auto rh = regionHandles[regionId];
                                if (rh.invalid())
                                    continue;
                                (dir < 0 ? abitFarLeftRightRegionIds[0] : abitFarLeftRightRegionIds[1]).insert(regionId);
                            }
                        }
                        for (int regionId : connectedRegionIds) {
                            RegionLineConnectionData & rd = regionLineConnections[std::make_pair(regionHandles[regionId], ld.topo.hd)];
                            rd.normalizedAnchors.push_back(normalize(cam.toSpace(sampleP)));
                            rd.detachable = !(Contains(abitFarLeftRightRegionIds[0], regionId) && Contains(abitFarLeftRightRegionIds[1], regionId));
                        }
                    }
                } else {
                    double spanAngle = AngleBetweenDirections(line.first, line.second);
                    int stepNum = static_cast<int>(std::ceil(spanAngle / samplingStepAngleOnLine));
                    assert(stepNum >= 1);
                    for (int step = 0; step <= stepNum; step++) {
                        double angle = step * samplingStepAngleOnLine;
                        Vec3 sample = RotateDirection(line.first, line.second, angle);
                        if (!cam.isVisibleOnScreen(sample))
                            continue;
                        Point2 sampleP = cam.toScreen(sample);
                        Pixel originalP = ToPixel(sampleP);
                        // collect neighbors
                        std::set<int> connectedRegionIds;
                        for (int x = originalP.x - samplerSizeOnLine; x <= originalP.x + samplerSizeOnLine; x++) {
                            for (int y = originalP.y - samplerSizeOnLine; y <= originalP.y + samplerSizeOnLine; y++) {
                                Pixel p(x, y);
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
                        for (int d = std::min(samplerSizeOnLine, 2); d <= samplerSizeOnLine * 2; d++) {
                            for (int dir : {-d, d}) {
                                Pixel p = ToPixel(sampleP + vertToLineDir * dir);
                                if (!Contains(segmentedRegions, p))
                                    continue;
                                int regionId = segmentedRegions(p);
                                auto rh = regionHandles[regionId];
                                if (rh.invalid())
                                    continue;
                                (dir < 0 ? abitFarLeftRightRegionIds[0] : abitFarLeftRightRegionIds[1]).insert(regionId);
                            }
                        }
                        for (int regionId : connectedRegionIds) {
                            RegionLineConnectionData & rd = regionLineConnections[std::make_pair(regionHandles[regionId], ld.topo.hd)];
                            rd.normalizedAnchors.push_back(normalize(sample));
                            rd.detachable = !(Contains(abitFarLeftRightRegionIds[0], regionId) && Contains(abitFarLeftRightRegionIds[1], regionId));
                        }
                    }
                }
            }

            for (auto & rlc : regionLineConnections) {
                // compute length;
                rlc.second.length = AngleBetweenDirections(rlc.second.normalizedAnchors.front(), rlc.second.normalizedAnchors.back());
                mg.addConstraint(std::move(rlc.second), rlc.first.first, rlc.first.second);
            }


            // add region boundary constraints
            assert(bndpixels.size() == bnd2segs.size());

            std::vector<RegionBoundaryHandle> boundaryHandles(bndpixels.size());
            std::map<std::pair<int, int>, RegionBoundaryHandle> segs2bh;

            for (int i = 0; i < bndpixels.size(); i++) {
                auto & rids = bnd2segs[i];
                auto edge = bndpixels[i];

                if (regionHandles[rids.first].invalid() || regionHandles[rids.second].invalid())
                    continue;
                if (edge.size() < 2)
                    continue;

                auto rh = regionHandles[rids.first];
                if (noBoundaryUnderLines) {
                    std::vector<Pixel> filteredEdge;
                    for (Pixel p : edge) {
                        bool coveredByLine = false;
                        for (RegionLineConnectionHandle rlcon : mg.topo(rh).constraints<RegionLineConnectionData>()) {
                            LineHandle lh = mg.topo(rlcon).component<1>();
                            const Line3 & line = mg.data(lh).line;
                            Line2 line2(cam.toScreen(line.first), cam.toScreen(line.second));
                            double d = DistanceFromPointToLine(ecast<double>(p), line2).first;
                            if (d < samplerSizeOnLine || d < samplerSizeOnBoundary) {
                                coveredByLine = true;
                                break;
                            }
                        }
                        if (!coveredByLine) {
                            filteredEdge.push_back(p);
                        }
                    }
                    if (filteredEdge.size() <= 2) {
                        continue;
                    }
                    edge = std::move(filteredEdge);
                }

                RegionBoundaryHandle & bh = segs2bh[rids];
                if (bh.invalid()) {
                    RegionBoundaryData bdorigin;
                    bdorigin.length = 0;
                    bh = mg.addConstraint(std::move(bdorigin), regionHandles[rids.first], regionHandles[rids.second]);
                }
                RegionBoundaryData & bd = mg.data(bh);

                std::vector<Vec3> edge3d;
                for (auto & p : edge) {
                    edge3d.push_back(normalize(cam.toSpace(p)));
                }
                bd.normalizedEdges.push_back(std::move(edge3d));
                auto & lastEdge = bd.normalizedEdges.back();

                // get edge point projections
                for (int i = 0; i < int(lastEdge.size()) - 1; i++) {
                    bd.length += AngleBetweenDirections(lastEdge[i], lastEdge[i + 1]);
                }

                // get sample points
                std::vector<Vec3> points = { lastEdge.front() };
                for (auto & edgeP : lastEdge) {
                    double remainedAngle = AngleBetweenDirections(points.back(), edgeP);
                    while (remainedAngle >= samplingStepAngleOnBoundary) {
                        points.push_back(normalize(RotateDirection(points.back(), edgeP, samplingStepAngleOnBoundary)));
                        remainedAngle -= samplingStepAngleOnBoundary;
                    }
                }
                bd.normalizedSampledPoints.push_back(std::move(points));
                assert(bd.normalizedSampledPoints.size() == bd.normalizedEdges.size());

                boundaryHandles[i] = bh;
            }


            return std::make_pair(regionHandles, boundaryHandles);

        }


        std::pair<std::vector<RegionHandle>, std::vector<RegionBoundaryHandle>> AppendRegions2(RLGraph & mg, const Imagei & segmentedRegions,
            const std::vector<std::vector<Pixel>> & bndpixels,
            const std::vector<std::pair<int, int>> & bnd2segs,
            const PanoramicCamera & cam,
            double samplingStepAngleOnBoundary, double samplingStepAngleOnLine) {

            using CameraT = PanoramicCamera;

            auto regionHandles = CollectRegionsFromSegmentation(mg, segmentedRegions, cam);
            int regionNum = regionHandles.size();


            // add region-line connections
            std::map<std::pair<RegionHandle, LineHandle>, RegionLineConnectionData> regionLineConnections;

            for (auto & ld : mg.components<LineData>()) {
                auto & line = ld.data.line;
                Line2 line2(cam.toScreen(line.first), cam.toScreen(line.second));

                double stepOnImage = std::max(cam.focal() * samplingStepAngleOnLine, 0.5);               
                double spanAngle = AngleBetweenDirections(line.first, line.second);
                int stepNum = static_cast<int>(std::ceil(spanAngle / samplingStepAngleOnLine));
                assert(stepNum >= 1);

                for (int step = 0; step <= stepNum; step++) {
                    double angle = step * samplingStepAngleOnLine;
                    Vec3 sample = RotateDirection(line.first, line.second, angle);
                    if (!cam.isVisibleOnScreen(sample))
                        continue;

                    Point2 sampleP = cam.toScreen(sample);
                    Pixel originalP = ToPixel(sampleP);

                    // collect neighbors
                    std::set<int> connectedRegionIds;
                    for (int samplerSizeOnLine = 1; samplerSizeOnLine < 4; samplerSizeOnLine++) {
                        for (int x = originalP.x - samplerSizeOnLine; x <= originalP.x + samplerSizeOnLine; x++) {
                            for (int y = originalP.y - samplerSizeOnLine; y <= originalP.y + samplerSizeOnLine; y++) {
                                if (Distance(Pixel(x, y), originalP) > samplerSizeOnLine) {
                                    continue;
                                }
                                Pixel p(WrapBetween(x, 0, segmentedRegions.cols), WrapBetween(y, 0, segmentedRegions.rows));
                                int regionId = segmentedRegions(p);
                                auto rh = regionHandles[regionId];
                                if (rh.invalid())
                                    continue;
                                connectedRegionIds.insert(regionId);
                            }
                        }
                        if (connectedRegionIds.size() >= 2) {
                            break;
                        }
                    }

                    if (connectedRegionIds.size() < 2) {
                        continue;
                    }

                    for (int regionId : connectedRegionIds) {
                        RegionLineConnectionData & rd = regionLineConnections[std::make_pair(regionHandles[regionId], ld.topo.hd)];
                        rd.normalizedAnchors.push_back(normalize(sample));
                        rd.detachable = false;
                    }
                }
               
            }

            for (auto & rlc : regionLineConnections) {
                // compute length;
                rlc.second.length = AngleBetweenDirections(rlc.second.normalizedAnchors.front(), rlc.second.normalizedAnchors.back());
                mg.addConstraint(std::move(rlc.second), rlc.first.first, rlc.first.second);
            }


            // add region-region
            assert(bndpixels.size() == bnd2segs.size());

            std::vector<RegionBoundaryHandle> boundaryHandles(bndpixels.size());
            std::map<std::pair<int, int>, RegionBoundaryHandle> segs2bh;

            for (int i = 0; i < bndpixels.size(); i++) {
                auto & rids = bnd2segs[i];
                auto edge = bndpixels[i];

                if (regionHandles[rids.first].invalid() || regionHandles[rids.second].invalid())
                    continue;
                if (edge.size() < 2)
                    continue;

                auto rh = regionHandles[rids.first];

                RegionBoundaryHandle & bh = segs2bh[rids];
                if (bh.invalid()) {
                    RegionBoundaryData bdorigin;
                    bdorigin.length = 0;
                    auto rh1 = regionHandles[rids.first];
                    auto rh2 = regionHandles[rids.second];
                    if (rh1.id > rh2.id) {
                        std::swap(rh1, rh2);
                    }
                    bh = mg.addConstraint(std::move(bdorigin), rh1, rh2);
                }

                RegionBoundaryData & bd = mg.data(bh);

                std::vector<Vec3> edge3d;
                for (auto & p : edge) {
                    edge3d.push_back(normalize(cam.toSpace(p)));
                }
                bd.normalizedEdges.push_back(std::move(edge3d));
                auto & lastEdge = bd.normalizedEdges.back();

                // get edge point projections
                for (int i = 0; i < int(lastEdge.size()) - 1; i++) {
                    bd.length += AngleBetweenDirections(lastEdge[i], lastEdge[i + 1]);
                }

                // get sample points
                std::vector<Vec3> points = { lastEdge.front() };
                for (auto & edgeP : lastEdge) {
                    double remainedAngle = AngleBetweenDirections(points.back(), edgeP);
                    while (remainedAngle >= samplingStepAngleOnBoundary) {
                        points.push_back(normalize(RotateDirection(points.back(), edgeP, samplingStepAngleOnBoundary)));
                        remainedAngle -= samplingStepAngleOnBoundary;
                    }
                }
                bd.normalizedSampledPoints.push_back(std::move(points));
                assert(bd.normalizedSampledPoints.size() == bd.normalizedEdges.size());

                boundaryHandles[i] = bh;
            }


            return std::make_pair(regionHandles, boundaryHandles);

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
                    contourProj.push_back(ecast<int>(ppc.toScreen(d)));
                }
            }
            cv::fillPoly(mask, contourProjs, (uint8_t)1);
            return View<PartialPanoramicCamera, Imageub>{mask, ppc};
        }



        HandledTable<RegionBoundaryHandle, int> CollectOcclusionResponseOnBoundaries(const RLGraph & mg, 
            const std::vector<Chain3> & occ, const PanoramicCamera & cam,
            double samplingStepAngleOnOcc, int samplerSizeOnOcc) {

            auto isOccs = HandledTable<RegionBoundaryHandle, int>(mg.internalConstraints<RegionBoundaryData>().size(), 0);

            HandledTable<RegionBoundaryHandle, std::vector<std::vector<int>>> rrSamplePointsTouched(mg.internalConstraints<RegionBoundaryData>().size());
            for (auto & rr : mg.constraints<RegionBoundaryData>()) {
                auto & touched = rrSamplePointsTouched[rr.topo.hd];
                touched.resize(rr.data.normalizedSampledPoints.size());
                for (int i = 0; i < touched.size(); i++) {
                    touched[i].resize(rr.data.normalizedSampledPoints[i].size(), 0);
                }
            }

            RTree<Box3, std::tuple<RegionBoundaryHandle, int, int>> rrSamplingPointsTree;
            for (auto & rr : mg.constraints<RegionBoundaryData>()) {
                for (int i = 0; i < rr.data.normalizedSampledPoints.size(); i++) {
                    auto & ps = rr.data.normalizedSampledPoints[i];
                    for (int j = 0; j < ps.size(); j++) {
                        auto & p = ps[j];
                        rrSamplingPointsTree.insert(BoundingBox(p).expand(2 * samplingStepAngleOnOcc / cam.focal()), std::make_tuple(rr.topo.hd, i, j));
                    }
                }
            }

            for (int i = 0; i < occ.size(); i++) {
                auto & points = occ[i].points;
                std::vector<Vec3> pts; // sampling pts
                if (points.empty()) {
                    continue;
                }
                auto last = points.front();
                for (int j = 1; j < points.size(); j++) {
                    double dist = AngleBetweenDirections(last, points[j]);
                    while (dist >= samplingStepAngleOnOcc) {
                        auto p = RotateDirection(last, points[j], samplingStepAngleOnOcc);
                        pts.push_back(normalize(p));
                        dist -= samplingStepAngleOnOcc;
                        last = p;
                    }
                }
                // for debugging
                if (false) {
                    double originalLen = 0.0;
                    for (int j = 1; j < points.size(); j++) {
                        double dist = AngleBetweenDirections(points[j - 1], points[j]);
                        originalLen += dist;
                    }
                    double len = 0.0;
                    for (int j = 1; j < pts.size(); j++) {
                        double dist = AngleBetweenDirections(pts[j - 1], pts[j]);
                        len += dist;
                    }
                    std::cout << "original len: " << originalLen << " len:" << len << std::endl;
                }
                for (auto & p : pts) {
                    rrSamplingPointsTree.search(BoundingBox(p).expand(2 * samplerSizeOnOcc / cam.focal()),
                        [&p, &rrSamplePointsTouched, &mg, samplerSizeOnOcc, &cam](const std::tuple<RegionBoundaryHandle, int, int> & spId) {
                        const Vec3 & pos = mg.data(Get<0>(spId)).normalizedSampledPoints[Get<1>(spId)][Get<2>(spId)];
                        if (Distance(p, pos) <= samplerSizeOnOcc / cam.focal()) {
                            rrSamplePointsTouched[Get<0>(spId)][Get<1>(spId)][Get<2>(spId)] = 1;
                        }
                        return true;
                    });
                }
            }
            for (auto & rd : mg.constraints<RegionBoundaryData>()) {
                int touchedCount = 0;
                int allCount = 0;
                for (auto & ts : rrSamplePointsTouched[rd.topo.hd]) {
                    for (int t : ts) {
                        allCount++;
                        if (t) {
                            touchedCount++;
                        }
                    }
                }
                if (touchedCount > 0.9 * allCount) {
                    isOccs[rd.topo.hd] = true;
                }
            }

            return isOccs;
        }

        HandledTable<RegionBoundaryHandle, int> CollectOcclusionResponseOnBoundaries(const RLGraph & mg, 
            const std::vector<Scored<Chain3>> & occbnds, const PanoramicCamera & cam,
            double samplingStepAngleOnOcc /*= 1 / 500.0*/, int samplerSizeOnOcc /*= 1*/) {
            std::vector<Chain3> occs;
            for (auto & soc : occbnds) {
                if (soc.score > 0) {
                    occs.push_back(soc.component);
                }
            }
            return CollectOcclusionResponseOnBoundaries(mg, occs, cam, samplingStepAngleOnOcc, samplerSizeOnOcc);
        }


    }
}


