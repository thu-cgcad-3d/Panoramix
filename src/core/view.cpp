
extern "C" {
    #include <gpc.h>
    #include <mosek.h>
}

#include <Eigen/Dense>
#include <Eigen/Sparse>

//#include <ICP.h>

#include "utilities.hpp"
#include "algorithms.hpp"
#include "feature.hpp"
#include "containers.hpp"
#include "view.hpp"

#include "debug.hpp"

namespace panoramix {
    namespace core {

        View<PanoramicCamera> CreatePanoramicView(const Image & panorama) {
            assert(abs(panorama.cols - panorama.rows * 2) < panorama.rows / 10.0f);
            return View<PanoramicCamera>{panorama, PanoramicCamera(panorama.cols / M_PI / 2.0)};
        }


        namespace {

            struct ComparePixelLoc {
                inline bool operator ()(const PixelLoc & a, const PixelLoc & b) const {
                    if (a.x != b.x)
                        return a.x < b.x;
                    return a.y < b.y;
                }
            };

            template <class T>
            inline std::pair<T, T> MakeOrderedPair(const T & a, const T & b) {
                return a < b ? std::make_pair(a, b) : std::make_pair(b, a);
            }

            void FindContoursOfRegionsAndBoundaries(const Imagei & segRegions, int regionNum,
                std::map<std::pair<int, int>, std::vector<std::vector<PixelLoc>>> & boundaryEdges) {

                std::map<std::pair<int, int>, std::set<PixelLoc, ComparePixelLoc>> boundaryPixels;

                int width = segRegions.cols;
                int height = segRegions.rows;
                for (int y = 0; y < height - 1; y++) {
                    for (int x = 0; x < width - 1; x++) {
                        PixelLoc p1(x, y), p2(x + 1, y), p3(x, y + 1), p4(x + 1, y + 1);
                        int rs[] = {
                            segRegions(p1),
                            segRegions(p2),
                            segRegions(p3),
                            segRegions(p4)
                        };

                        if (rs[0] != rs[1]) {
                            boundaryPixels[MakeOrderedPair(rs[0], rs[1])].insert(p1);
                        }
                        if (rs[0] != rs[2]) {
                            boundaryPixels[MakeOrderedPair(rs[0], rs[2])].insert(p1);
                        }
                        if (rs[0] != rs[3]) {
                            boundaryPixels[MakeOrderedPair(rs[0], rs[3])].insert(p1);
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


        RegionsGraph CreateRegionsGraph(const Imagei & segmentedRegions,
            double samplingStepLengthOnBoundary, int dilationSize) {

            RegionsGraph regions;

            double minVal, maxVal;
            std::tie(minVal, maxVal) = MinMaxValOfImage(segmentedRegions);
            assert(minVal == 0.0);

            int regionNum = static_cast<int>(maxVal) + 1;
            regions.internalElements<0>().reserve(regionNum);

            for (int i = 0; i < regionNum; i++){
                RegionData rd;
                rd.area = 0.0;
                rd.center = Point2(0, 0);
                regions.add(std::move(rd));
            }

            // fill basic properties for all region data
            for (auto i = segmentedRegions.begin(); i != segmentedRegions.end(); ++i){
                auto id = *i;
                auto & rd = regions.data(RegionHandle(id));
                rd.area += 1.0;
                PixelLoc p = i.pos();
                rd.boundingBox |= BoundingBox(p);
                rd.center += Vec2(p.x, p.y);
            }
            for (auto & r : regions.elements<0>()){
                r.data.center /= r.data.area;
            }
            
            // calculate contours for each region data
            for (int i = 0; i < regionNum; i++){
                RegionData & vd = regions.data(RegionHandle(i));
                Image regionMask = (segmentedRegions == i);

                // find contour of the region
                cv::Mat regionMaskCopy;
                regionMask.copyTo(regionMaskCopy);
                cv::findContours(regionMaskCopy, vd.contours, CV_RETR_EXTERNAL, CV_CHAIN_APPROX_SIMPLE);
                std::sort(vd.contours.begin(), vd.contours.end(),
                    [](const std::vector<PixelLoc> & ca, const std::vector<PixelLoc> & cb){return ca.size() > cb.size(); });
                assert(!vd.contours.empty() && "no contour? impossible~");
                //for (int k = 1; k < vd.contours.size(); k++)
                //    assert(vd.contours[k].size() <= 2);

                auto dilationElement = cv::getStructuringElement(cv::MORPH_ELLIPSE,
                    cv::Size(2 * dilationSize + 1, 2 * dilationSize + 1),
                    cv::Point(dilationSize, dilationSize));
                regionMask.copyTo(regionMaskCopy);
                cv::dilate(regionMaskCopy, regionMaskCopy, dilationElement);
                cv::findContours(regionMaskCopy, vd.dilatedContours, CV_RETR_EXTERNAL, CV_CHAIN_APPROX_SIMPLE);
                std::sort(vd.dilatedContours.begin(), vd.dilatedContours.end(),
                    [](const std::vector<PixelLoc> & ca, const std::vector<PixelLoc> & cb){return ca.size() > cb.size(); });
                assert(!vd.dilatedContours.empty() && "no contour? impossible~");
                //for (int k = 1; k < vd.dilatedContours.size(); k++)
                //    assert(vd.dilatedContours[k].size() <= 2);
            }



            std::map<std::pair<int, int>, std::vector<std::vector<PixelLoc>>> boundaryEdges;
            FindContoursOfRegionsAndBoundaries(segmentedRegions, regionNum, boundaryEdges);

            for (auto & bep : boundaryEdges) {
                auto & rids = bep.first;
                auto & edges = bep.second;

                RegionBoundaryData bd;
                bd.length = 0;
                bd.edges = edges;
                for (auto & e : edges) {
                    assert(!e.empty() && "edges should never be empty!");
                    for (int i = 0; i < e.size() - 1; i++) {
                        bd.length += Distance(e[i], e[i + 1]);
                    }
                }

                // compute straightness
                double interleavedArea, interleavedLength;
                std::tie(bd.straightness, bd.fittedLine) =
                    ComputeStraightness(edges, &interleavedArea, &interleavedLength);

                // collect sampled points
                double stepLen = samplingStepLengthOnBoundary;
                bd.sampledPoints.resize(edges.size());
                for (int i = 0; i < edges.size(); i++) {
                    auto & e = edges[i];
                    assert(e.size() >= 2 && "invalid point num for an edge");
                    double remLen = 0;
                    auto & s = bd.sampledPoints[i];
                    s.push_back(ToPoint2(e.front()));
                    for (int k = 1; k < e.size(); k++){
                        auto & nextP = ToPoint2(e[k]);
                        remLen += Distance(s.back(), nextP);
                        while (remLen >= stepLen) {
                            auto p = s.back() + (nextP - s.back()) * stepLen / remLen;
                            s.push_back(p);
                            remLen -= stepLen;
                        }
                    }
                }

                regions.add<1>({ RegionHandle(rids.first), RegionHandle(rids.second) }, bd);
            }

            // find connection of end points of boundaries
            struct BoundaryEndPoint {
                RegionBoundaryHandle bh;
                int edgeId;
                bool isFirst;
                PixelLoc position;
            };
            struct BoundaryEndPointBoundingBoxFunctor {
                inline Box2 operator()(const BoundaryEndPoint & p) const {
                    auto bb = BoundingBox(p.position);
                    for (int i = 0; i < bb.Dimension; i++) {
                        bb.maxCorner[i] += 2;
                        bb.minCorner[i] -= 2;
                    }
                    return bb;
                }
            };
            std::map<PixelLoc, std::vector<BoundaryEndPoint>, ComparePixelLoc> mergedBepsTable;
            RTreeWrapper<BoundaryEndPoint, BoundaryEndPointBoundingBoxFunctor> existedBoundaryEndPoints;
            for (auto & boundary : regions.elements<1>()) {
                for (int eid = 0; eid < boundary.data.edges.size(); eid++) {
                    BoundaryEndPoint beps[] = {
                        { boundary.topo.hd, eid, true, boundary.data.edges[eid].front() },
                        { boundary.topo.hd, eid, false, boundary.data.edges[eid].back() }
                    };
                    for (auto & bep : beps) {
                        bool existed = false;
                        existedBoundaryEndPoints.searchNear(bep, [&bep, &mergedBepsTable, &existed](const BoundaryEndPoint & found) {
                            if (Distance(bep.position, found.position) < 4) {
                                mergedBepsTable[found.position].push_back(bep);
                                existed = true;
                                return false;
                            }
                            return true;
                        });
                        if (!existed) {
                            mergedBepsTable[bep.position].push_back(bep);
                            existedBoundaryEndPoints.insert(bep);
                        }
                    }
                }
            }

            // t-junction likelihood for each bep
            for (auto & b : regions.elements<1>()) {
                b.data.tjunctionLikelihood = 0;
            }

            for (auto & beps : mergedBepsTable) {
                std::set<RegionBoundaryHandle> bhs;
                for (auto & bep : beps.second) {
                    bhs.insert(bep.bh);
                }
                if (bhs.size() == 3 && beps.second.size() == 3) {
                    std::vector<std::vector<Point<float, 2>>> sampledPointsOnBoundaries(3);
                    static const size_t tjunctionSampleNum = 5;
                    for (int i = 0; i < 3; i++) {
                        auto & bep = beps.second[i];
                        auto & spts = regions.data(bep.bh).sampledPoints[bep.edgeId];
                        sampledPointsOnBoundaries[i].resize(std::min(tjunctionSampleNum, spts.size()));
                        for (int j = 0; j < sampledPointsOnBoundaries[i].size(); j++) {
                            sampledPointsOnBoundaries[i][j] = ToPoint2f(spts[bep.isFirst ? j : (spts.size() - 1 - j)]);
                        }
                    }
                    double boundaryPairSmoothnesses[3];
                    for (int i = 0; i < 3; i++) {
                        for (int j = i + 1; j < 3; j++) {
                            int k = (0 + 1 + 2) - i - j;
                            std::vector<Point<float, 2>> sptsij(sampledPointsOnBoundaries[i]);
                            sptsij.insert(sptsij.end(), sampledPointsOnBoundaries[j].begin(), sampledPointsOnBoundaries[j].end());
                            cv::Vec4f line;
                            cv::fitLine(sptsij, line, CV_DIST_L2, 0, 0.01, 0.01);
                            auto fittedLine = InfiniteLine2({ line[2], line[3] }, { line[0], line[1] });
                            double interArea = 0;
                            double interLen = 0;
                            for (int ss = 0; ss < sptsij.size() - 1; ss++) {
                                double area, len;
                                std::tie(area, len) = ComputeSpanningArea(
                                    ToPoint2(sptsij[ss]),
                                    ToPoint2(sptsij[ss + 1]),
                                    fittedLine);
                                interArea += area;
                                interLen += len;
                            }
                            boundaryPairSmoothnesses[k] = interArea / interLen;
                        }
                    }
                    double tjunctionLikelihoods[3];
                    for (int i = 0; i < 3; i++) {
                        for (int j = i + 1; j < 3; j++) {
                            int k = (0 + 1 + 2) - i - j;
                            tjunctionLikelihoods[k] = Gaussian(std::max(boundaryPairSmoothnesses[i] / boundaryPairSmoothnesses[k],
                                boundaryPairSmoothnesses[j]), 0.1);
                            //tjunctionLikeliHoodsTable[beps.second[k].bh].push_back(tjunctionLikelihoods[k]);
                            regions.data(beps.second[k].bh).tjunctionLikelihood +=
                                tjunctionLikelihoods[k];
                        }
                    }
                }
            }

            assert(!regions.hasGarbage());
            return regions;

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



        LinesGraph CreateLinesGraph(const std::vector<Classified<Line2>> & lines,
            const std::vector<HPoint2> & vps,
            double intersectionDistanceThreshold,
            double incidenceDistanceAlongDirectionThreshold,
            double incidenceDistanceVerticalDirectionThreshold){

            LinesGraph graph;

            // insert lines
            std::vector<LineHandle> handles;
            handles.reserve(lines.size());
            graph.internalElements<0>().reserve(lines.size());

            for (auto & line : lines){
                if (line.claz == -1)
                    continue;
                LineData ld;
                ld.line = line;
                handles.push_back(graph.add(std::move(ld)));
            }

            // construct incidence/intersection relations
            auto & linesData = graph.internalElements<0>();
            for (int i = 0; i < linesData.size(); i++){
                auto & linei = linesData[i].data.line.component;
                int clazi = linesData[i].data.line.claz;
                assert(clazi != -1);

                for (int j = i + 1; j < linesData.size(); j++){
                    auto & linej = linesData[j].data.line.component;
                    int clazj = linesData[j].data.line.claz;
                    assert(clazj != -1);

                    auto nearest = DistanceBetweenTwoLines(linei, linej);
                    double d = nearest.first;

                    if (clazi == clazj){
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
                            LineRelationData lrd;
                            lrd.type = LineRelationData::Type::Incidence;
                            lrd.relationCenter = conCenter;

                            if (HasValue(lrd.relationCenter, IsInfOrNaN<double>))
                                continue;

                            graph.add<1>({ handles[i], handles[j] }, std::move(lrd));
                        }
                    }
                    else {
                        if (d < intersectionDistanceThreshold){
                            auto conCenter = HPointFromVector(GetCoeffs(linei.infiniteLine())
                                .cross(GetCoeffs(linej.infiniteLine()))).value();

                            assert(conCenter != Point2(0.0, 0.0));

                            if (DistanceFromPointToLine(conCenter, linei).first > intersectionDistanceThreshold * 4 ||
                                DistanceFromPointToLine(conCenter, linej).first > intersectionDistanceThreshold * 4)
                                continue;

                            LineRelationData lrd;
                            lrd.type = LineRelationData::Type::Intersection;
                            lrd.relationCenter = conCenter;

                            if (HasValue(lrd.relationCenter, IsInfOrNaN<double>))
                                continue;

                            graph.add<1>({ handles[i], handles[j] }, std::move(lrd));
                        }
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
                if (lrd.type == lrd.Incidence){
                    lrd.junctionWeight = IncidenceJunctionWeight(false);
                }
                else if (lrd.type == lrd.Intersection){
                    Mat<float, 3, 2> votingData;
                    std::fill(std::begin(votingData.val), std::end(votingData.val), 0);

                    for (auto & ld : graph.elements<0>()){
                        auto & line = ld.data.line.component;
                        int claz = ld.data.line.claz;
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



        namespace {

            std::vector<int> FillInRectangleWithXs(int extendSize){
                std::vector<int> dx;
                dx.reserve(2 * extendSize + 1);
                for (int a = -extendSize; a <= extendSize; a++) {
                    for (int b = -extendSize; b <= extendSize; b++) {
                        dx.push_back(a);
                    }
                }
                return dx;
            }

            std::vector<int> FillInRectangleWithYs(int extendSize){
                std::vector<int> dy;
                dy.reserve(2 * extendSize + 1);
                for (int a = -extendSize; a <= extendSize; a++) {
                    for (int b = -extendSize; b <= extendSize; b++) {
                        dy.push_back(b);
                    }
                }
                return dy;
            }

        }



        std::map<std::pair<RegionHandle, LineHandle>, std::vector<Point2>> 
            RecognizeRegionLineConnections(const Imagei & segmentedRegions,
            const LinesGraph & linesGraph, double samplingStepLengthOnLines){

            // generate sampled points for line-region connections
            std::map<std::pair<RegionHandle, LineHandle>, std::vector<Point2>> regionLineConnections;

            static const int OPT_ExtendSize = 2;
            static const std::vector<int> dx = FillInRectangleWithXs(OPT_ExtendSize);
            static const std::vector<int> dy = FillInRectangleWithYs(OPT_ExtendSize);

            for (auto & ld : linesGraph.elements<0>()) {
                auto & line = ld.data.line.component;
                auto lineDir = normalize(line.direction());
                double sampleStep = samplingStepLengthOnLines;
                int sampledNum = static_cast<int>(std::floor(line.length() / sampleStep));

                for (int i = 0; i < sampledNum; i++) {
                    auto sampledPoint = line.first + lineDir * i * sampleStep;
                    std::set<int32_t> rhids;
                    for (int k = 0; k < dx.size(); k++) {
                        int x = BoundBetween(static_cast<int>(std::round(sampledPoint[0] + dx[k])), 0, segmentedRegions.cols - 1);
                        int y = BoundBetween(static_cast<int>(std::round(sampledPoint[1] + dy[k])), 0, segmentedRegions.rows - 1);
                        PixelLoc p(x, y);
                        rhids.insert(segmentedRegions(p));
                    }
                    for (int32_t rhid : rhids) {
                        regionLineConnections[std::make_pair(RegionHandle(rhid), ld.topo.hd)]
                            .push_back(sampledPoint);
                    }
                }
            }

            return regionLineConnections;
        }






        namespace {

            inline double LatitudeFromLongitudeAndNormalVector(double longitude, const Vec3 & normal) {
                // normal(0)*cos(long)*cos(la) + normal(1)*sin(long)*cos(lat) + normal(2)*sin(la) = 0
                // normal(0)*cos(long) + normal(1)*sin(long) + normal(2)*tan(la) = 0
                return -atan((normal(0)*cos(longitude) + normal(1)*sin(longitude)) / normal(2));
            }

            inline double Longitude1FromLatitudeAndNormalVector(double latitude, const Vec3 & normal) {
                double a = normal(1) * cos(latitude);
                double b = normal(0) * cos(latitude);
                double c = -normal(2) * sin(latitude);
                double sinLong = (a * c + sqrt(Square(a*c) - (Square(a) + Square(b))*(Square(c) - Square(b)))) / (Square(a) + Square(b));
                return asin(sinLong);
            }

            inline double Longitude2FromLatitudeAndNormalVector(double latitude, const Vec3 & normal) {
                double a = normal(1) * cos(latitude);
                double b = normal(0) * cos(latitude);
                double c = -normal(2) * sin(latitude);
                double sinLong = (a * c - sqrt(Square(a*c) - (Square(a) + Square(b))*(Square(c) - Square(b)))) / (Square(a) + Square(b));
                return asin(sinLong);
            }

            inline double UnOrthogonality(const Vec3 & v1, const Vec3 & v2, const Vec3 & v3) {
                return norm(Vec3(v1.dot(v2), v2.dot(v3), v3.dot(v1)));
            }

            std::vector<Vec3> FindVanishingPoints(const std::vector<Vec3>& intersections,
                int longitudeDivideNum = 1000, int latitudeDivideNum = 500) {

                std::vector<Vec3> vps(3);

                Imagef votePanel = Imagef::zeros(longitudeDivideNum, latitudeDivideNum);

                std::cout << "begin voting ..." << std::endl;
                size_t pn = intersections.size();
                for (const Vec3& p : intersections){
                    PixelLoc pixel = PixelLocFromGeoCoord(GeoCoord(p), longitudeDivideNum, latitudeDivideNum);
                    votePanel(pixel.x, pixel.y) += 1.0;
                }
                std::cout << "begin gaussian bluring ..." << std::endl;
                cv::GaussianBlur(votePanel, votePanel, cv::Size((longitudeDivideNum / 50) * 2 + 1, (latitudeDivideNum / 50) * 2 + 1),
                    4, 4, cv::BORDER_REPLICATE);
                std::cout << "done voting" << std::endl;

                double minVal = 0, maxVal = 0;
                int maxIndex[] = { -1, -1 };
                cv::minMaxIdx(votePanel, &minVal, &maxVal, 0, maxIndex);
                cv::Point maxPixel(maxIndex[0], maxIndex[1]);

                vps[0] = GeoCoordFromPixelLoc(maxPixel, longitudeDivideNum, latitudeDivideNum).toVector();
                const Vec3 & vec0 = vps[0];

                // iterate locations orthogonal to vps[0]
                double maxScore = -1;
                for (int x = 0; x < longitudeDivideNum; x++){
                    double longt1 = double(x) / longitudeDivideNum * M_PI * 2 - M_PI;
                    double lat1 = LatitudeFromLongitudeAndNormalVector(longt1, vec0);
                    Vec3 vec1 = GeoCoord(longt1, lat1).toVector();
                    Vec3 vec1rev = -vec1;
                    Vec3 vec2 = vec0.cross(vec1);
                    Vec3 vec2rev = -vec2;
                    Vec3 vecs[] = { vec1, vec1rev, vec2, vec2rev };

                    double score = 0;
                    for (Vec3 & v : vecs){
                        PixelLoc pixel = PixelLocFromGeoCoord(GeoCoord(v), longitudeDivideNum, latitudeDivideNum);
                        score += votePanel(WrapBetween(pixel.x, 0, longitudeDivideNum),
                            WrapBetween(pixel.y, 0, latitudeDivideNum));
                    }
                    if (score > maxScore){
                        maxScore = score;
                        vps[1] = vec1;
                        vps[2] = vec2;
                    }
                }

                if (UnOrthogonality(vps[0], vps[1], vps[2]) < 0.1)
                    return vps;

                // failed, then use y instead of x
                maxScore = -1;
                for (int y = 0; y < latitudeDivideNum; y++){
                    double lat1 = double(y) / latitudeDivideNum * M_PI - M_PI_2;
                    double longt1s[] = { Longitude1FromLatitudeAndNormalVector(lat1, vec0),
                        Longitude2FromLatitudeAndNormalVector(lat1, vec0) };
                    for (double longt1 : longt1s){
                        Vec3 vec1 = GeoCoord(longt1, lat1).toVector();
                        Vec3 vec1rev = -vec1;
                        Vec3 vec2 = vec0.cross(vec1);
                        Vec3 vec2rev = -vec2;
                        Vec3 vecs[] = { vec1, vec1rev, vec2, vec2rev };

                        double score = 0;
                        for (Vec3 & v : vecs){
                            PixelLoc pixel = PixelLocFromGeoCoord(GeoCoord(v), longitudeDivideNum, latitudeDivideNum);
                            score += votePanel(WrapBetween(pixel.x, 0, longitudeDivideNum),
                                WrapBetween(pixel.y, 0, latitudeDivideNum));
                        }
                        if (score > maxScore){
                            maxScore = score;
                            vps[1] = vec1;
                            vps[2] = vec2;
                        }
                    }
                }

                assert(UnOrthogonality(vps[0], vps[1], vps[2]) < 0.1);

                return vps;

            }

            inline Vec3 RotateDirectionTo(const Vec3 & originalDirection, const Vec3 & toDirection, double angle) {
                Vec3 tovec = originalDirection.cross(toDirection).cross(originalDirection);
                Vec3 result3 = originalDirection + tovec * tan(angle);
                return result3 / norm(result3);
            }

        }


        void EstimateVanishingPointsAndBuildLinesGraphs(const std::vector<View<PerspectiveCamera>> & views,
            std::vector<Vec3> & vanishingPoints,
            std::vector<LinesGraph> & linesGraphs,
            double intersectionDistanceThreshold,
            double incidenceDistanceAlongDirectionThreshold,
            double incidenceDistanceVerticalDirectionThreshold){

            // collect lines and intersections
            std::vector<std::vector<Line2>> linesInViews;
            int linesNum = 0;

            linesInViews.reserve(views.size());
            std::vector<Vec3> lineIntersections;
            LineSegmentExtractor lineseg;
            lineseg.params().useLSD = true;
            for (auto & v : views){
                linesInViews.push_back(lineseg(v.image));
                linesNum += linesInViews.back().size();
                auto inters = ComputeLineIntersections(linesInViews.back(), nullptr, true, std::numeric_limits<double>::max());
                // insert line intersections
                for (auto & p : inters){
                    lineIntersections.push_back(normalize(v.camera.spatialDirection(p.value())));
                }
            }

            // find vanishing points using line intersections (projected in space)
            vanishingPoints = FindVanishingPoints(lineIntersections);

            // project lines to space
            std::vector<Classified<Line3>> spatialLineSegments;
            spatialLineSegments.reserve(linesNum);
            for (int i = 0; i < views.size(); i++){
                for (const auto & line : linesInViews[i]) {
                    auto & p1 = line.first;
                    auto & p2 = line.second;
                    auto pp1 = views[i].camera.spatialDirection(p1);
                    auto pp2 = views[i].camera.spatialDirection(p2);
                    Classified<Line3> cline3;
                    cline3.claz = -1;
                    cline3.component = Line3{ pp1, pp2 };
                    spatialLineSegments.push_back(cline3);
                }
            }

            // classify lines
            ClassifyLines3D(spatialLineSegments, vanishingPoints, M_PI / 3.0, 0.1, 0.8);

            // build lines graph
            linesGraphs.clear();
            linesGraphs.reserve(views.size());
            auto spatialLineSegmentsIter = spatialLineSegments.begin();
            for (int i = 0; i < views.size(); i++){
                // set line classes for line2ds
                auto & linesInThisView = linesInViews[i];
                std::vector<Classified<Line2>> clines(linesInThisView.size());
                for (int k = 0; k < linesInThisView.size(); k++){
                    clines[k].component = linesInThisView[k];
                    clines[k].claz = spatialLineSegmentsIter->claz;
                    ++spatialLineSegmentsIter;
                }
                std::vector<HPoint2> vanishingPointsInThisView(vanishingPoints.size());
                for (int k = 0; k < vanishingPoints.size(); k++){
                    vanishingPointsInThisView[k] = views[i].camera.screenProjectionInHPoint(vanishingPoints[k]);
                }
                linesGraphs.push_back(CreateLinesGraph(clines, vanishingPointsInThisView,
                    intersectionDistanceThreshold, 
                    incidenceDistanceAlongDirectionThreshold, incidenceDistanceVerticalDirectionThreshold));
            }
        }



        namespace {

            // polygon conversion
            void ConvertToGPCPolygon(const std::vector<PixelLoc> & pts, gpc_polygon & poly) {
                poly.num_contours = 1;
                poly.contour = new gpc_vertex_list[1];
                poly.contour[0].num_vertices = pts.size();
                poly.contour[0].vertex = new gpc_vertex[pts.size()];
                for (int i = 0; i < pts.size(); i++) {
                    poly.contour[0].vertex[i].x = pts[i].x;
                    poly.contour[0].vertex[i].y = pts[i].y;
                }
                poly.hole = new int[1];
                poly.hole[0] = 0;
            }

            void ConvertToPixelVector(const gpc_polygon & poly, std::vector<PixelLoc> & pts) {
                pts.clear();
                pts.resize(poly.contour[0].num_vertices);
                for (int i = 0; i < pts.size(); i++) {
                    pts[i].x = static_cast<int>(poly.contour[0].vertex[i].x);
                    pts[i].y = static_cast<int>(poly.contour[0].vertex[i].y);
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

            template <class T, int N>
            inline Line<T, N> NormalizeLine(const Line<T, N> & l) {
                return Line<T, N>(normalize(l.first), normalize(l.second));
            }

        }



        template <class CameraT>
        ComponentIndexHashMap<std::pair<RegionIndex, RegionIndex>, double> 
            RecognizeRegionOverlappingsAcrossViewsTemplated(const std::vector<View<CameraT>> & views,
            const std::vector<RegionsGraph> & regionsGraphs){

            assert(views.size() == regionsGraphs.size());
            assert(AllSameInContainer(views, [](const View<CameraT> & v1, const View<CameraT> & v2){
                return v1.camera.eye() == v2.camera.eye(); 
            }));

            ComponentIndexHashMap<std::pair<RegionIndex, RegionIndex>, double> regionOverlappings;

            // compute spatial positions of each region
            ComponentIndexHashMap<RegionIndex, std::vector<Vec3>>
                regionSpatialContours;
            for (int i = 0; i < views.size(); i++) {
                const auto & regions = regionsGraphs[i];
                for (auto & region : regions.elements<0>()) {
                    RegionIndex ri = { i, region.topo.hd };
                    const RegionData & rd = region.data;
                    std::vector<Vec3> spatialContour;
                    if (!rd.dilatedContours.empty()){
                        for (auto & p : rd.dilatedContours.back()) {
                            auto direction = views[i].camera.spatialDirection(p);
                            spatialContour.push_back(direction / norm(direction));
                        }
                    }
                    else{
                        std::cerr << "this region has no dilatedCountour!" << std::endl;
                    }
                    regionSpatialContours[ri] = spatialContour;
                }
            }

            // build spatial rtree for regions
            auto lookupRegionBB = [&regionSpatialContours](const RegionIndex& ri) {
                return BoundingBoxOfContainer(regionSpatialContours[ri]);
            };

            RTreeWrapper<RegionIndex, decltype(lookupRegionBB)> regionsRTree(lookupRegionBB);
            for (auto & region : regionSpatialContours) {
                regionsRTree.insert(region.first);
            }

            // store overlapping ratios between overlapped regions
            for (auto & rip : regionSpatialContours) {
                auto & ri = rip.first;

                auto & riCountours = GetData(ri, regionsGraphs).contours;
                if (riCountours.empty()){
                    std::cerr << "this region has no countour!" << std::endl;
                    continue;
                }

                auto & riContour2d = riCountours.front();
                auto & riCamera = views[ri.viewId].camera;
                double riArea = GetData(ri, regionsGraphs).area;
                //double riArea = cv::contourArea(riContour2d);

                gpc_polygon riPoly;
                ConvertToGPCPolygon(riContour2d, riPoly);

                regionsRTree.search(lookupRegionBB(ri),
                    [&ri, &riContour2d, &riPoly, &riCamera, riArea, &regionOverlappings, &regionSpatialContours](
                    const RegionIndex & relatedRi) {

                    if (ri.viewId == relatedRi.viewId) {
                        return true;
                    }

                    // project relatedRi contour to ri's camera plane
                    auto & relatedRiContour3d = regionSpatialContours[relatedRi];
                    std::vector<core::PixelLoc> relatedRiContour2d(relatedRiContour3d.size());
                    for (int i = 0; i < relatedRiContour3d.size(); i++) {
                        auto p = riCamera.screenProjection(relatedRiContour3d[i]);
                        relatedRiContour2d[i] = PixelLoc(p);
                    }
                    gpc_polygon relatedRiPoly;
                    ConvertToGPCPolygon(relatedRiContour2d, relatedRiPoly);

                    // compute overlapping area ratio
                    gpc_polygon intersectedPoly;
                    gpc_polygon_clip(GPC_INT, &relatedRiPoly, &riPoly, &intersectedPoly);

                    if (intersectedPoly.num_contours > 0 && intersectedPoly.contour[0].num_vertices > 0) {
                        std::vector<core::PixelLoc> intersected;
                        ConvertToPixelVector(intersectedPoly, intersected);
                        double intersectedArea = cv::contourArea(intersected);

                        double overlapRatio = intersectedArea / riArea;
                        //assert(overlapRatio <= 1.5 && "Invalid overlap ratio!");

                        if (overlapRatio > 0.2)
                            regionOverlappings[std::make_pair(relatedRi, ri)] = overlapRatio;
                    }

                    gpc_free_polygon(&relatedRiPoly);
                    gpc_free_polygon(&intersectedPoly);

                    return true;
                });

                gpc_free_polygon(&riPoly);
            }

            return regionOverlappings;

        }



        ComponentIndexHashMap<std::pair<RegionIndex, RegionIndex>, double>
            RecognizeRegionOverlappingsAcrossViews(const std::vector<View<PerspectiveCamera>> & views,
            const std::vector<RegionsGraph> & regionsGraphs){
            return RecognizeRegionOverlappingsAcrossViewsTemplated(views, regionsGraphs);
        }

        ComponentIndexHashMap<std::pair<RegionIndex, RegionIndex>, double>
            RecognizeRegionOverlappingsAcrossViews(const std::vector<View<PanoramicCamera>> & views,
            const std::vector<RegionsGraph> & regionsGraphs){
            return RecognizeRegionOverlappingsAcrossViewsTemplated(views, regionsGraphs);
        }




        ComponentIndexHashMap<std::pair<LineIndex, LineIndex>, Vec3>
            RecognizeLineIncidencesAcrossViews(const std::vector<View<PerspectiveCamera>> & views,
            const std::vector<LinesGraph> & linesGraphs,
            double interViewIncidenceAngleAlongDirectionThreshold,
            double interViewIncidenceAngleVerticalDirectionThreshold){

            ComponentIndexHashMap<std::pair<LineIndex, LineIndex>, Vec3> interViewLineIncidences;

            //// LINES ////
            // compute spatial normal directions for each line
            ComponentIndexHashMap<LineIndex, Classified<Line3>>
                lineSpatialAvatars;
            for (int i = 0; i < views.size(); i++) {
                auto & lines = linesGraphs[i];
                LineIndex li;
                li.viewId = i;
                auto & cam = views[i].camera;
                for (auto & ld : lines.elements<0>()) {
                    li.handle = ld.topo.hd;
                    auto & line = ld.data.line;
                    Classified<Line3> avatar;
                    avatar.claz = line.claz;
                    avatar.component = Line3(
                        cam.spatialDirection(line.component.first),
                        cam.spatialDirection(line.component.second)
                        );
                    lineSpatialAvatars[li] = avatar;
                }
            }

            // build rtree for lines
            auto lookupLineNormal = [&lineSpatialAvatars](const LineIndex & li) -> Box3 {
                auto normal = lineSpatialAvatars[li].component.first.cross(lineSpatialAvatars[li].component.second);
                Box3 b = BoundingBox(normalize(normal));
                static const double s = 0.2;
                b.minCorner = b.minCorner - Vec3(s, s, s);
                b.maxCorner = b.maxCorner + Vec3(s, s, s);
                return b;
            };

            RTreeWrapper<LineIndex, decltype(lookupLineNormal)> linesRTree(lookupLineNormal);
            for (auto & i : lineSpatialAvatars) {
                linesRTree.insert(i.first);
            }

            // recognize incidence constraints between lines of different views
            for (auto & i : lineSpatialAvatars) {
                auto li = i.first;
                auto & lineData = i.second;
                linesRTree.search(lookupLineNormal(li),
                    [interViewIncidenceAngleAlongDirectionThreshold, interViewIncidenceAngleVerticalDirectionThreshold,
                    &li, &lineSpatialAvatars,
                    &views, &linesGraphs, &interViewLineIncidences](const LineIndex & relatedLi) -> bool {
                    if (li.viewId == relatedLi.viewId)
                        return true;
                    if (relatedLi < li) // make sure one relation is stored only once, avoid storing both a-b and b-a
                        return true;

                    auto & line1 = lineSpatialAvatars[li];
                    auto & line2 = lineSpatialAvatars[relatedLi];
                    if (line1.claz != line2.claz) // only incidence relations are recognized here
                        return true;

                    auto normal1 = normalize(line1.component.first.cross(line1.component.second));
                    auto normal2 = normalize(line2.component.first.cross(line2.component.second));

                    if (std::min(std::abs(AngleBetweenDirections(normal1, normal2)),
                        std::abs(AngleBetweenDirections(normal1, -normal2))) <
                        interViewIncidenceAngleVerticalDirectionThreshold) {

                        auto nearest = DistanceBetweenTwoLines(NormalizeLine(line1.component), NormalizeLine(line2.component));
                        if (AngleBetweenDirections(nearest.second.first.position, nearest.second.second.position) >
                            interViewIncidenceAngleAlongDirectionThreshold) // ignore too far-away relations
                            return true;

                        auto relationCenter = (nearest.second.first.position + nearest.second.second.position) / 2.0;
                        relationCenter /= norm(relationCenter);

                        interViewLineIncidences[std::make_pair(li, relatedLi)] = relationCenter;
                    }
                    return true;
                });
            }

            // check whether all interview incidences are valid
            IF_DEBUG_USING_VISUALIZERS{
                double maxDist = 0;
                Line3 farthestLine1, farthestLine2;
                for (auto & lir : interViewLineIncidences) {
                    auto & line1 = lineSpatialAvatars[lir.first.first];
                    auto & line2 = lineSpatialAvatars[lir.first.second];
                    if (line1.claz != line2.claz) {
                        std::cout << "invalid classes!" << std::endl;
                    }
                    auto l1 = NormalizeLine(line1.component);
                    auto l2 = NormalizeLine(line2.component);
                    auto dist = DistanceBetweenTwoLines(l1, l2).first;
                    if (dist > maxDist) {
                        farthestLine1 = l1;
                        farthestLine2 = l2;
                        maxDist = dist;
                    }
                }
                {
                    std::cout << "max dist of interview incidence pair: " << maxDist << std::endl;
                    std::cout << "line1: " << farthestLine1.first << ", " << farthestLine1.second << std::endl;
                    std::cout << "line2: " << farthestLine2.first << ", " << farthestLine2.second << std::endl;
                    auto d = DistanceBetweenTwoLines(farthestLine1, farthestLine2);
                    double angleDist = AngleBetweenDirections(d.second.first.position, d.second.second.position);
                    std::cout << "angle dist: " << angleDist << std::endl;
                }
            }

            return interViewLineIncidences;

        }


        Plane3 PlaneOfMGUnary(const MGUnary & region, const std::vector<Vec3> & vps, const MGUnaryVariable & var) {
            assert(region.type == MGUnary::Region);
            return Plane3(region.normalizedCorners.front() * var.depthOfFirstCorner, vps[var.claz]);
        }

        Line3 LineOfMGUnary(const MGUnary & line, const std::vector<Vec3> & vps, const MGUnaryVariable & var) {
            assert(line.type == MGUnary::Line);
            InfiniteLine3 infLine(line.normalizedCorners.front() * var.depthOfFirstCorner, vps[var.claz]);
            return Line3(line.normalizedCorners.front() * var.depthOfFirstCorner,
                DistanceBetweenTwoLines(InfiniteLine3(Point3(0, 0, 0), line.normalizedCorners.back()), infLine).second.second);
        }

        double DepthRatioOnMGUnary(const Vec3 & direction, const MGUnary & unary, const std::vector<Vec3> & vps, int claz) {
            if (unary.type == MGUnary::Region){
                const auto & region = unary;
                assert(!region.normalizedCorners.empty());
                Plane3 plane(region.normalizedCorners.front(), vps[claz]);
                return norm(IntersectionOfLineAndPlane(InfiniteLine3(Point3(0, 0, 0), direction), plane).position);
            }
            else if (unary.type == MGUnary::Line){
                const auto & line = unary;
                InfiniteLine3 infLine(line.normalizedCorners.front(), vps[claz]);
                return norm(DistanceBetweenTwoLines(InfiniteLine3(Point3(0, 0, 0), direction), infLine).second.first);
            }
        }

        Point3 LocationOnMGUnary(const Vec3 & direction, const MGUnary & unary, const std::vector<Vec3> & vps,
            const MGUnaryVariable & var){
            if (unary.type == MGUnary::Region){
                const auto & region = unary;
                assert(!region.normalizedCorners.empty());
                Plane3 plane(region.normalizedCorners.front() * var.depthOfFirstCorner, vps[var.claz]);
                return IntersectionOfLineAndPlane(InfiniteLine3(Point3(0, 0, 0), direction), plane).position;
            }
            else if (unary.type == MGUnary::Line){
                const auto & line = unary;
                InfiniteLine3 infLine(line.normalizedCorners.front() * var.depthOfFirstCorner, vps[var.claz]);
                return DistanceBetweenTwoLines(InfiniteLine3(Point3(0, 0, 0), direction), infLine).second.first;
            }
        }


        template <class T>
        inline std::vector<T> NumFilter(const std::vector<T> & data, int numLimit){
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


        void UpdateBinaryVars(const MixedGraph & mg, const std::vector<Vec3> & vps, 
            const MGUnaryVarTable & unaryVars, MGBinaryVarTable & binaryVars){
            for (auto & bv : binaryVars){
                auto & uhs = mg.topo(bv.first).lowers;
                auto & samples = mg.data(bv.first).samples;
                for (int i = 0; i < samples.size(); i++){
                    bv.second.sampleDepthsOnRelatedUnaries.front()[i] = unaryVars.at(uhs.front()).depthOfFirstCorner *
                        DepthRatioOnMGUnary(samples[i], mg.data(uhs.front()), vps, unaryVars.at(uhs.front()).claz);
                    bv.second.sampleDepthsOnRelatedUnaries.back()[i] = unaryVars.at(uhs.back()).depthOfFirstCorner *
                        DepthRatioOnMGUnary(samples[i], mg.data(uhs.back()), vps, unaryVars.at(uhs.back()).claz);
                }
            }
        }


        MixedGraph BuildMixedGraph(const std::vector<View<PerspectiveCamera>> & views,
            const std::vector<RegionsGraph> & regionsGraphs, const std::vector<LinesGraph> & linesGraphs,
            const ComponentIndexHashMap<std::pair<RegionIndex, RegionIndex>, double> & regionOverlappingsAcrossViews,
            const ComponentIndexHashMap<std::pair<LineIndex, LineIndex>, Vec3> & lineIncidencesAcrossViews,
            const std::vector<std::map<std::pair<RegionHandle, LineHandle>, std::vector<Point2>>> & regionLineConnections,
            const std::vector<Vec3> & vps,
            MGUnaryVarTable & unaryVars, MGBinaryVarTable & binaryVars,
            double initialDepth){

            MixedGraph mg;
            ComponentIndexHashMap<RegionIndex, MGUnaryHandle> ri2mgh;
            ComponentIndexHashMap<LineIndex, MGUnaryHandle> li2mgh;

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
                    assert(normalizedContour.size() > 2);
                    auto center = cam.spatialDirection(rd.data.center);
                    int initialClaz = std::distance(vps.begin(), std::min_element(vps.begin(), vps.end(),
                        [&center](const Vec3 & vp1, const Vec3 & vp2) -> bool {
                        return AngleBetweenUndirectedVectors(center, vp1) <
                            AngleBetweenUndirectedVectors(center, vp2);
                    }));
                    ri2mgh[ri] = mg.add(MGUnary{
                        MGUnary::Region,
                        std::move(normalizedContour),
                        center
                    });
                    unaryVars[ri2mgh[ri]].claz = initialClaz;
                    unaryVars[ri2mgh[ri]].depthOfFirstCorner = initialDepth;
                }
                // lines
                for (auto & ld : linesGraphs[i].elements<0>()){
                    auto li = LineIndex{ i, ld.topo.hd };
                    li2mgh[li] = mg.add(MGUnary{
                        MGUnary::Line,
                        {
                            normalize(cam.spatialDirection(ld.data.line.component.first)),
                            normalize(cam.spatialDirection(ld.data.line.component.second))
                        },
                        cam.spatialDirection(ld.data.line.component.center())
                    });
                    unaryVars[li2mgh[li]].claz = ld.data.line.claz;
                    unaryVars[li2mgh[li]].depthOfFirstCorner = initialDepth;
                }
                // region-region in each view
                for (auto & bd : regionsGraphs[i].elements<1>()){
                    MGBinary rr;
                    rr.type = MGBinary::RegionRegionConnection;
                    rr.weight = 1.0;
                    rr.samples.reserve(bd.data.sampledPoints.size());
                    for (auto & ps : bd.data.sampledPoints){
                        for (auto & p : ps){
                            rr.samples.push_back(cam.spatialDirection(p));
                        }
                    }
                    auto r1 = RegionIndex{ i, bd.topo.lowers.front() };
                    auto r2 = RegionIndex{ i, bd.topo.lowers.back() };
                    mg.add<1>({ ri2mgh[r1], ri2mgh[r2] }, std::move(rr));
                }
                // region-line
                for (auto & regionLine : regionLineConnections[i]){
                    MGBinary rl;
                    rl.type = MGBinary::RegionLineConnection;
                    rl.weight = 1.0;
                    rl.samples.reserve(regionLine.second.size());
                    for (auto & p : regionLine.second){
                        rl.samples.push_back(cam.spatialDirection(p));
                    }
                    auto ri = RegionIndex{ i, regionLine.first.first };
                    auto li = RegionIndex{ i, regionLine.first.second };
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
                        llinter.samples = { normalize(cam.spatialDirection(rd.data.relationCenter)) };
                        mg.add<1>({ li2mgh[l1], li2mgh[l2] }, std::move(llinter));
                    }
                    else if (rd.data.type == LineRelationData::Incidence){
                        MGBinary llincid;
                        llincid.type = MGBinary::LineLineIncidence;
                        llincid.weight = rd.data.junctionWeight * 10;
                        llincid.samples = { normalize(cam.spatialDirection(rd.data.relationCenter)) };
                        mg.add<1>({ li2mgh[l1], li2mgh[l2] }, std::move(llincid));
                    }
                }
            }
            // add cross view constraints
            // region-region overlappings
            for (auto & regionOverlapping : regionOverlappingsAcrossViews){
                if (regionOverlapping.second < 0.05)
                    continue;
                MGBinary rro;
                rro.type = MGBinary::RegionRegionOverlapping;
                rro.weight = 100;
                auto & r1 = regionOverlapping.first.first;
                auto & r2 = regionOverlapping.first.second;

                // get samples
                Vec3 z = mg.data(ri2mgh[r1]).center + mg.data(ri2mgh[r2]).center;
                Vec3 x, y;
                std::tie(x, y) = ProposeXYDirectionsFromZDirection(z);
                double minx = std::numeric_limits<double>::max(),
                    miny = std::numeric_limits<double>::max();
                double maxx = std::numeric_limits<double>::lowest(),
                    maxy = std::numeric_limits<double>::lowest();
                rro.samples.resize(4, z);
                for (auto & a : mg.data(ri2mgh[r1]).normalizedCorners){
                    double dx = a.dot(x), dy = a.dot(y);
                    if (dx < minx){
                        rro.samples[0] = a;
                        minx = dx;
                    }
                    else if (dx > maxx){
                        rro.samples[1] = a;
                        maxx = dx;
                    }
                    if (dy < miny){
                        rro.samples[2] = a;
                        miny = dy;
                    }
                    else if (dy > maxy){
                        rro.samples[3] = a;
                        maxy = dy;
                    }
                }
                for (auto & a : mg.data(ri2mgh[r2]).normalizedCorners){
                    double dx = a.dot(x), dy = a.dot(y);
                    if (dx < minx){
                        rro.samples[0] = a;
                        minx = dx;
                    }
                    else if (dx > maxx){
                        rro.samples[1] = a;
                        maxx = dx;
                    }
                    if (dy < miny){
                        rro.samples[2] = a;
                        miny = dy;
                    }
                    else if (dy > maxy){
                        rro.samples[3] = a;
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
                llincid.samples = { lineIncidence.second };
                auto & l1 = lineIncidence.first.first;
                auto & l2 = lineIncidence.first.second;
                mg.add<1>({ li2mgh[l1], li2mgh[l2] }, std::move(llincid));
            }

            // compute importance ratios
            std::vector<double> unaryWeightSums(mg.internalElements<0>().size(), 0.0);
            for (auto & b : mg.internalElements<1>()){
                unaryWeightSums[b.topo.lowers.front().id] += b.data.weight * b.data.samples.size();
                unaryWeightSums[b.topo.lowers.back().id] += b.data.weight * b.data.samples.size();
            }
            for (auto & b : mg.internalElements<1>()){
                b.data.importanceRatioInRelatedUnaries.front() = b.data.weight * b.data.samples.size() / 
                    unaryWeightSums[b.topo.lowers.front().id];
                b.data.importanceRatioInRelatedUnaries.back() = b.data.weight * b.data.samples.size() /
                    unaryWeightSums[b.topo.lowers.back().id];
            }

#ifdef _DEBUG
            for (auto & u : mg.elements<0>()){
                double importanceRatioSum = 0.0;
                for (auto & bh : u.topo.uppers){
                    auto & imp = mg.data(bh).importanceRatioInRelatedUnaries;
                    importanceRatioSum += mg.topo(bh).lowers.front() == u.topo.hd ? imp.front() : imp.back();
                }
                assert(FuzzyEquals(importanceRatioSum, 1.0, 0.01));
            }
#endif

            // make overlapped region claz consistent
            std::map<int, int> clazForCCUsingRegionOverlapping;
            std::vector<MGUnaryHandle> uhs;
            for (auto & v : mg.elements<0>()){
                if (v.data.type == MGUnary::Region){
                    uhs.push_back(v.topo.hd);
                }
            }
            core::ConnectedComponents(uhs.begin(), uhs.end(), [&mg](const MGUnaryHandle & uh) {
                std::vector<MGUnaryHandle> neighbors;
                for (auto & bh : mg.topo(uh).uppers){
                    if (mg.data(bh).type == MGBinary::RegionRegionOverlapping){
                        auto anotherUh = mg.topo(bh).lowers.front();
                        if (anotherUh == uh)
                            anotherUh = mg.topo(bh).lowers.back();
                        neighbors.push_back(anotherUh);
                    }
                }
                return neighbors;
            }, [&clazForCCUsingRegionOverlapping, &mg, &unaryVars](const MGUnaryHandle & uh, int ccid){
                if (Contains(clazForCCUsingRegionOverlapping, ccid)){
                    unaryVars[uh].claz = clazForCCUsingRegionOverlapping[ccid];
                }
                else{
                    clazForCCUsingRegionOverlapping[ccid] = unaryVars[uh].claz;
                }
            });

            binaryVars.clear();
            for (auto & b : mg.elements<1>()){
                binaryVars[b.topo.hd].sampleDepthsOnRelatedUnaries.front().resize(b.data.samples.size());
                binaryVars[b.topo.hd].sampleDepthsOnRelatedUnaries.back().resize(b.data.samples.size());
            }
            UpdateBinaryVars(mg, vps, unaryVars, binaryVars);

            return mg;
        }



        int MarkConnectedComponentIds(const MixedGraph & mg, std::unordered_map<MGUnaryHandle, int> & ccids) {
            ccids.reserve(mg.internalElements<0>().size());
            std::vector<MGUnaryHandle> uhs;
            uhs.reserve(mg.internalElements<0>().size());
            for (auto & v : mg.elements<0>()){
                uhs.push_back(v.topo.hd);
            }
            return core::ConnectedComponents(uhs.begin(), uhs.end(), 
                [&mg](const MGUnaryHandle & uh){
                std::vector<MGUnaryHandle> neighborUhs;
                for (auto & bh : mg.topo(uh).uppers){
                    auto anotherUh = mg.topo(bh).lowers.front();
                    if (anotherUh == uh)
                        anotherUh = mg.topo(bh).lowers.back();
                    neighborUhs.push_back(anotherUh);
                }
                return neighborUhs;
            }, [&ccids](const MGUnaryHandle & uh, int ccid){
                ccids[uh] = ccid; 
            });
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
            DepthFirstSearch(uhs.begin(), uhs.end(), [&mg, &patch](const MGUnaryHandle & uh){
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

        std::vector<MGPatch> SplitMixedGraphIntoPatches(const MixedGraph & mg,
            const MGUnaryVarTable & unaryVars, const MGBinaryVarTable & binaryVars){
           
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
                patches[ccid].uhs[uh] = unaryVars.at(uh);
            }
            for (auto & b : mg.elements<1>()){
                auto & uhs = mg.topo(b.topo.hd).lowers;
                assert(ccids[uhs.front()] == ccids[uhs.back()]);
                patches[ccids[uhs.front()]].bhs[b.topo.hd] = binaryVars.at(b.topo.hd);
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
                assert(ccids[uhs.front()] == ccids[uhs.back()]);
                patches[ccids[uhs.front()]].bhs[b.topo.hd] = patch.bhs.at(b.topo.hd);
            }
            return patches;
        }


        double ScoreOfMGPatch(const MixedGraph & mg, const MGPatch & patch) {
            double distanceSum = 0.0;
            int distanceCount = 0;
            for (auto & bhv : patch.bhs){
                for (int i = 0; i < bhv.second.sampleDepthsOnRelatedUnaries.front().size(); i++){
                    distanceSum += 
                        abs(bhv.second.sampleDepthsOnRelatedUnaries.front()[i] - 
                        bhv.second.sampleDepthsOnRelatedUnaries.back()[i]);
                }
                distanceCount += bhv.second.sampleDepthsOnRelatedUnaries.front().size();
            }
            double depthSum = 0.0;
            for (auto & uhv : patch.uhs){
                depthSum += uhv.second.depthOfFirstCorner;
            }
            double violation = (distanceSum / distanceCount) / (depthSum / patch.uhs.size());
            return 1.0 - violation;
        }

        void OptimizeDepthsOfMGPatch(const MixedGraph & mg, MGPatch & patch, 
            const std::vector<Vec3> & vps, double depthLb, double depthUb){

            auto tick = Tick("preparing optimization");

            assert(BinaryHandlesAreValidInPatch(mg, patch));
            int depthNum = patch.uhs.size();
            
            assert(UnariesAreConnectedInPatch(mg, patch));
            int ccConsNum = 1;           
   

            // temp data
            struct UnaryTempData {
                int position;
            };
            std::unordered_map<MGUnaryHandle, UnaryTempData> unaryTemp;
            unaryTemp.reserve(patch.uhs.size());
            int unaryNum = 0;
            for (auto & uhv : patch.uhs){
                unaryTemp[uhv.first].position = unaryNum;
                unaryNum++;
            }

            struct BinaryTempData {
                int firstSamplePosition;
            };
            std::unordered_map<MGBinaryHandle, BinaryTempData> binaryTemp;
            binaryTemp.reserve(patch.bhs.size());

            // count sample connections
            int samplesNum = 0;
            for (auto & bhv : patch.bhs){
                binaryTemp[bhv.first].firstSamplePosition = samplesNum;
                samplesNum += mg.data(bhv.first).samples.size();
            }

            int slackVarNum = samplesNum;
            int varNum = depthNum + slackVarNum; // depths, slackVars
            int consNum = samplesNum * 2;   // depth1 * ratio1 - depth2 * ratio2 < slackVar;  
                                            // depth2 * ratio2 - depth1 * ratio1 < slackVar

            //make as bounds:   depth = defaultValue
            //                  depth \in [depthLb, depthUb]

            MSKenv_t env = nullptr;
            MSKtask_t task = nullptr;

            MSK_makeenv(&env, nullptr);
            MSK_maketask(env, consNum, varNum, &task);

            MSK_appendcons(task, consNum);
            MSK_appendvars(task, varNum);

            // set weights
            int slackVarId = 0;
            for (auto & bhv : patch.bhs){
                auto & bd = mg.data(bhv.first);
                for (int k = 0; k < bd.samples.size(); k++){
                    MSK_putcj(task, slackVarId, bd.weight);
                    slackVarId++;
                }
            }

            // bounds for vars
            {
                int varId = 0;
                MSK_putvarbound(task, varId, MSK_BK_FX,
                    patch.uhs.begin()->second.depthOfFirstCorner,
                    patch.uhs.begin()->second.depthOfFirstCorner);
                varId++;
                for (; varId < depthNum; varId++){ // for depths
                    MSK_putvarbound(task, varId, MSK_BK_RA, depthLb, depthUb);
                }
                for (; varId < depthNum + slackVarNum; varId++){ // for slack vars
                    MSK_putvarbound(task, varId, MSK_BK_LO, 0.0, +MSK_INFINITY);
                }
            }

            // fill constraints related to each vars
            {
                int varId = 0;
                for (auto & uhv : patch.uhs){ // depths
                    auto & uh = uhv.first;

                    auto & relatedBhs = mg.topo(uh).uppers;
                    int relatedAnchorsNum = 0;
                    for (auto & bh : relatedBhs){
                        if (!Contains(patch.bhs, bh))
                            continue;
                        relatedAnchorsNum += mg.data(bh).samples.size();
                    }
                    int relatedConsNum = relatedAnchorsNum * 2;

                    std::vector<MSKint32t> consIds;
                    std::vector<MSKrealt> consValues;

                    consIds.reserve(relatedConsNum);
                    consValues.reserve(relatedConsNum);

                    for (auto & bh : relatedBhs){
                        if (!Contains(patch.bhs, bh))
                            continue;
                        int firstAnchorPosition = binaryTemp[bh].firstSamplePosition;
                        auto & samples = mg.data(bh).samples;
                        for (int k = 0; k < samples.size(); k++){
                            consIds.push_back((firstAnchorPosition + k) * 2); // one for [depth1 * ratio1 - depth2 * ratio2 - slackVar < 0]; 
                            consIds.push_back((firstAnchorPosition + k) * 2 + 1); // another for [- depth1 * ratio1 + depth2 * ratio2 - slackVar < 0];

                            auto & a = samples[k];
                            double ratio = DepthRatioOnMGUnary(a, mg.data(uh), vps, patch.uhs.at(uh).claz);
                            bool isOnLeftSide = uh == mg.topo(bh).lowers.front();
                            if (isOnLeftSide){ // as depth1
                                consValues.push_back(ratio);
                                consValues.push_back(-ratio);
                            }
                            else{
                                consValues.push_back(-ratio);
                                consValues.push_back(ratio);
                            }
                        }
                    }

                    MSK_putacol(task, varId, relatedConsNum, consIds.data(), consValues.data());
                    varId++;
                }

                for (; varId < depthNum + slackVarNum; varId++){ // slack vars
                    MSKint32t consIds[] = {
                        (varId - depthNum) * 2, // one for [depth1 * ratio1 - depth2 * ratio2 - slackVar < 0]; 
                        (varId - depthNum) * 2 + 1  // another for [- depth1 * ratio1 + depth2 * ratio2 - slackVar < 0];
                    };
                    MSKrealt consValues[] = { -1.0, -1.0 };
                    MSK_putacol(task, varId, 2, consIds, consValues);
                }
            }

            // bounds for constraints
            for (int consId = 0; consId < consNum; consId++){
                // all [depth1 * ratio1 - depth2 * ratio2 - slackVar < 0];
                MSK_putconbound(task, consId, MSK_BK_UP, -MSK_INFINITY, 0.0);
            }

            // solve

            Tock(tick);

            tick = Tick("do optimization");

            MSK_putobjsense(task, MSK_OBJECTIVE_SENSE_MINIMIZE);
            MSKrescodee trmcode;
            MSK_optimizetrm(task, &trmcode);
            MSK_solutionsummary(task, MSK_STREAM_LOG);

            MSKsolstae solsta;
            MSK_getsolsta(task, MSK_SOL_BAS, &solsta);

            Tock(tick);

            switch (solsta){
            case MSK_SOL_STA_OPTIMAL:
            case MSK_SOL_STA_NEAR_OPTIMAL:
            {
                                             double *xx = new double[varNum];
                                             MSK_getxx(task, MSK_SOL_BAS, xx);
                                             printf("Optimal primal solution\n");
                                             for (auto & uhv : patch.uhs){ // get resulted depths
                                                 uhv.second.depthOfFirstCorner = xx[unaryTemp[uhv.first].position];
                                             }
                                             delete[] xx;
                                             UpdateBinaryVars(mg, vps, patch.uhs, patch.bhs);
                                             break;
            }
            case MSK_SOL_STA_DUAL_INFEAS_CER:
            case MSK_SOL_STA_PRIM_INFEAS_CER:
            case MSK_SOL_STA_NEAR_DUAL_INFEAS_CER:
            case MSK_SOL_STA_NEAR_PRIM_INFEAS_CER:
                printf("Primal or dual infeasibility certificate found.\n");
                break;
            case MSK_SOL_STA_UNKNOWN:
            {
                                        char symname[MSK_MAX_STR_LEN];
                                        char desc[MSK_MAX_STR_LEN];



                                        MSK_getcodedesc(trmcode,
                                            symname,
                                            desc);

                                        printf("The solution status is unknown.\n");
                                        printf("The optimizer terminitated with code: %s\n", symname);
                                        break;
            }
            default:
                printf("Other solution status.\n");
                break;
            }


            // finalize
            MSK_deletetask(&task);
            MSK_deleteenv(&env);

        }

        struct MosekData {
            MSKenv_t env;
            MSKtask_t task;

            struct UnaryTempData {
                int position;
            };
            std::unordered_map<MGUnaryHandle, UnaryTempData> unaryTemp;
            struct BinaryTempData {
                int firstSamplePosition;
            };
            std::unordered_map<MGBinaryHandle, BinaryTempData> binaryTemp;

            int varNum;
        };


        MGPatchDepthsOptimizer::MGPatchDepthsOptimizer(const MixedGraph & mg, MGPatch & patch,
            const std::vector<Vec3> & vanishingPoints) 
            : _mg(mg), _patch(patch), _vanishingPoints(vanishingPoints){

            MosekData * mosekData = new MosekData;
            _internal = mosekData;
            auto & env = mosekData->env;
            auto & task = mosekData->task;

            auto tick = Tick("preparing optimization");

            assert(BinaryHandlesAreValidInPatch(mg, patch));
            int depthNum = patch.uhs.size();

            assert(UnariesAreConnectedInPatch(mg, patch));
            int ccConsNum = 1;

            // temp data  
            auto & unaryTemp = mosekData->unaryTemp;
            auto & binaryTemp = mosekData->binaryTemp;

            unaryTemp.reserve(patch.uhs.size());
            int unaryNum = 0;
            for (auto & uhv : patch.uhs){
                unaryTemp[uhv.first].position = unaryNum;
                unaryNum++;
            }

           
            binaryTemp.reserve(patch.bhs.size());

            // count sample connections
            int samplesNum = 0;
            for (auto & bhv : patch.bhs){
                binaryTemp[bhv.first].firstSamplePosition = samplesNum;
                samplesNum += mg.data(bhv.first).samples.size();
            }

            int slackVarNum = samplesNum;
            auto & varNum = mosekData->varNum;
            varNum = depthNum + slackVarNum; // depths, slackVars
            int consNum = samplesNum * 2;   // depth1 * ratio1 - depth2 * ratio2 < slackVar;  
            // depth2 * ratio2 - depth1 * ratio1 < slackVar

            //make as bounds:   depth = defaultValue
            //                  depth \in [depthLb, depthUb]

            env = nullptr;
            task = nullptr;

            MSK_makeenv(&env, nullptr);
            MSK_maketask(env, consNum, varNum, &task);

            MSK_appendcons(task, consNum);
            MSK_appendvars(task, varNum);

            // set weights
            int slackVarId = 0;
            for (auto & bhv : patch.bhs){
                auto & bd = mg.data(bhv.first);
                for (int k = 0; k < bd.samples.size(); k++){
                    MSK_putcj(task, slackVarId, bd.weight);
                    slackVarId++;
                }
            }

            // bounds for vars
            {
                int varId = 0;
                MSK_putvarbound(task, varId, MSK_BK_FX,
                    patch.uhs.begin()->second.depthOfFirstCorner,
                    patch.uhs.begin()->second.depthOfFirstCorner);
                varId++;
                for (; varId < depthNum; varId++){ // for depths
                    MSK_putvarbound(task, varId, MSK_BK_RA, 0.1, 10.0);
                }
                for (; varId < depthNum + slackVarNum; varId++){ // for slack vars
                    MSK_putvarbound(task, varId, MSK_BK_LO, 0.0, +MSK_INFINITY);
                }
            }

            // fill constraints related to each vars
            {
                int varId = 0;
                for (auto & uhv : patch.uhs){ // depths
                    auto & uh = uhv.first;

                    auto & relatedBhs = mg.topo(uh).uppers;
                    int relatedAnchorsNum = 0;
                    for (auto & bh : relatedBhs){
                        if (!Contains(patch.bhs, bh))
                            continue;
                        relatedAnchorsNum += mg.data(bh).samples.size();
                    }
                    int relatedConsNum = relatedAnchorsNum * 2;

                    std::vector<MSKint32t> consIds;
                    std::vector<MSKrealt> consValues;

                    consIds.reserve(relatedConsNum);
                    consValues.reserve(relatedConsNum);

                    for (auto & bh : relatedBhs){
                        if (!Contains(patch.bhs, bh))
                            continue;
                        int firstAnchorPosition = binaryTemp[bh].firstSamplePosition;
                        auto & samples = mg.data(bh).samples;
                        for (int k = 0; k < samples.size(); k++){
                            consIds.push_back((firstAnchorPosition + k) * 2); // one for [depth1 * ratio1 - depth2 * ratio2 - slackVar < 0]; 
                            consIds.push_back((firstAnchorPosition + k) * 2 + 1); // another for [- depth1 * ratio1 + depth2 * ratio2 - slackVar < 0];

                            auto & a = samples[k];
                            double ratio = DepthRatioOnMGUnary(a, mg.data(uh), vanishingPoints, patch.uhs.at(uh).claz);
                            bool isOnLeftSide = uh == mg.topo(bh).lowers.front();
                            if (isOnLeftSide){ // as depth1
                                consValues.push_back(ratio);
                                consValues.push_back(-ratio);
                            }
                            else{
                                consValues.push_back(-ratio);
                                consValues.push_back(ratio);
                            }
                        }
                    }

                    MSK_putacol(task, varId, relatedConsNum, consIds.data(), consValues.data());
                    varId++;
                }

                for (; varId < depthNum + slackVarNum; varId++){ // slack vars
                    MSKint32t consIds[] = {
                        (varId - depthNum) * 2, // one for [depth1 * ratio1 - depth2 * ratio2 - slackVar < 0]; 
                        (varId - depthNum) * 2 + 1  // another for [- depth1 * ratio1 + depth2 * ratio2 - slackVar < 0];
                    };
                    MSKrealt consValues[] = { -1.0, -1.0 };
                    MSK_putacol(task, varId, 2, consIds, consValues);
                }
            }

            // bounds for constraints
            for (int consId = 0; consId < consNum; consId++){
                // all [depth1 * ratio1 - depth2 * ratio2 - slackVar < 0];
                MSK_putconbound(task, consId, MSK_BK_UP, -MSK_INFINITY, 0.0);
            }

            Tock(tick);

        }

        MGPatchDepthsOptimizer::~MGPatchDepthsOptimizer(){
            auto mosekData = static_cast<MosekData*>(_internal);
            MSK_deletetask(&mosekData->task);
            MSK_deleteenv(&mosekData->env);
            delete mosekData;
        }

        void MGPatchDepthsOptimizer::setDepthBounds(double depthLb, double depthUb){
            for (int varId = 1; varId < _patch.uhs.size(); varId++){ // for depths
                MSK_putvarbound(static_cast<MosekData*>(_internal)->task, varId, MSK_BK_RA, depthLb, depthUb);
            }
        }

        void MGPatchDepthsOptimizer::setUnaryClass(const MGUnaryHandle & uh, int claz){
            _patch.uhs[uh].claz = claz;

            auto & mg = _mg;
            auto & patch = _patch;
            auto & task = static_cast<MosekData*>(_internal)->task;
            auto & binaryTemp = static_cast<MosekData*>(_internal)->binaryTemp;
            auto & unaryTemp = static_cast<MosekData*>(_internal)->unaryTemp;

            auto & relatedBhs = mg.topo(uh).uppers;
            int relatedAnchorsNum = 0;
            for (auto & bh : relatedBhs){
                if (!Contains(patch.bhs, bh))
                    continue;
                relatedAnchorsNum += mg.data(bh).samples.size();
            }
            int relatedConsNum = relatedAnchorsNum * 2;

            std::vector<MSKint32t> consIds;
            std::vector<MSKrealt> consValues;

            consIds.reserve(relatedConsNum);
            consValues.reserve(relatedConsNum);

            for (auto & bh : relatedBhs){
                if (!Contains(patch.bhs, bh))
                    continue;
                int firstAnchorPosition = binaryTemp[bh].firstSamplePosition;
                auto & samples = mg.data(bh).samples;
                for (int k = 0; k < samples.size(); k++){
                    consIds.push_back((firstAnchorPosition + k) * 2); // one for [depth1 * ratio1 - depth2 * ratio2 - slackVar < 0]; 
                    consIds.push_back((firstAnchorPosition + k) * 2 + 1); // another for [- depth1 * ratio1 + depth2 * ratio2 - slackVar < 0];

                    auto & a = samples[k];
                    double ratio = DepthRatioOnMGUnary(a, mg.data(uh), _vanishingPoints, patch.uhs.at(uh).claz);
                    bool isOnLeftSide = uh == mg.topo(bh).lowers.front();
                    if (isOnLeftSide){ // as depth1
                        consValues.push_back(ratio);
                        consValues.push_back(-ratio);
                    }
                    else{
                        consValues.push_back(-ratio);
                        consValues.push_back(ratio);
                    }
                }
            }

            MSK_putacol(task, unaryTemp[uh].position, relatedConsNum, consIds.data(), consValues.data());
        }

        void MGPatchDepthsOptimizer::optimize() {
            auto tick = Tick("do optimization");

            auto mosekData = static_cast<MosekData*>(_internal);
            auto & task = mosekData->task;
            auto & varNum = mosekData->varNum;
            auto & unaryTemp = mosekData->unaryTemp;

            MSK_putobjsense(task, MSK_OBJECTIVE_SENSE_MINIMIZE);
            MSKrescodee trmcode;
            MSK_optimizetrm(task, &trmcode);
            MSK_solutionsummary(task, MSK_STREAM_LOG);

            MSKsolstae solsta;
            MSK_getsolsta(task, MSK_SOL_BAS, &solsta);

            Tock(tick);

            switch (solsta){
            case MSK_SOL_STA_OPTIMAL:
            case MSK_SOL_STA_NEAR_OPTIMAL:
            {
                                             double *xx = new double[varNum];
                                             MSK_getxx(task, MSK_SOL_BAS, xx);
                                             printf("Optimal primal solution\n");
                                             for (auto & uhv : _patch.uhs){ // get resulted depths
                                                 uhv.second.depthOfFirstCorner = xx[unaryTemp[uhv.first].position];
                                             }
                                             delete[] xx;
                                             UpdateBinaryVars(_mg, _vanishingPoints, _patch.uhs, _patch.bhs);
                                             break;
            }
            case MSK_SOL_STA_DUAL_INFEAS_CER:
            case MSK_SOL_STA_PRIM_INFEAS_CER:
            case MSK_SOL_STA_NEAR_DUAL_INFEAS_CER:
            case MSK_SOL_STA_NEAR_PRIM_INFEAS_CER:
                printf("Primal or dual infeasibility certificate found.\n");
                break;
            case MSK_SOL_STA_UNKNOWN:
            {
                                        char symname[MSK_MAX_STR_LEN];
                                        char desc[MSK_MAX_STR_LEN];
                                        MSK_getcodedesc(trmcode,
                                            symname,
                                            desc);

                                        printf("The solution status is unknown.\n");
                                        printf("The optimizer terminitated with code: %s\n", symname);
                                        break;
            }
            default:
                printf("Other solution status.\n");
                break;
            }
        }

        //struct BinaryAnchors {
        //    double weight;
        //    std::vector<Vec3> anchors;
        //    int positionOfFirstAnchor;
        //    std::vector<double> slackValues;
        //    std::array<std::vector<double>, 2> currentAnchorDepthsOnUnaries;
        //};





        //BinaryAnchors GetBinaryAnchors(const MixedGraph & mg, const MGBinaryHandle & bh){
        //    static const int maxSampleNum = 3;
        //    
        //    if (mg.data(bh).is<MGBinaryRegionRegionBoundary>()){
        //        BinaryAnchors bc;
        //        bc.anchors = NumFilter(mg.data(bh).uncheckedRef<MGBinaryRegionRegionBoundary>().samples, maxSampleNum);
        //        bc.weight = 0.1;
        //        return bc;
        //    }
        //    else if (mg.data(bh).is<MGBinaryRegionRegionOverlapping>()){
        //        BinaryAnchors bc;
        //        bc.weight = 1e2 * mg.data(bh).uncheckedRef<MGBinaryRegionRegionOverlapping>().overlappingRatio;
        //        auto & u1 = mg.data(mg.topo(bh).lowers.front());
        //        auto & u2 = mg.data(mg.topo(bh).lowers.back());
        //        Vec3 z = normalize(u1.uncheckedRef<MGUnaryRegion>().center + u2.uncheckedRef<MGUnaryRegion>().center);
        //        Vec3 x, y;
        //        std::tie(x, y) = ProposeXYDirectionsFromZDirection(z);
        //        double minx = std::numeric_limits<double>::max(), 
        //            miny = std::numeric_limits<double>::max();
        //        double maxx = std::numeric_limits<double>::lowest(), 
        //            maxy = std::numeric_limits<double>::lowest();
        //        //bc.anchors = { u1.uncheckedRef<MGUnaryRegion>().normalizedCorners.front() };
        //        //return bc;
        //        bc.anchors.resize(4, z);
        //        for (auto & a : u1.uncheckedRef<MGUnaryRegion>().normalizedCorners){
        //            double dx = a.dot(x), dy = a.dot(y);
        //            if (dx < minx){
        //                bc.anchors[0] = a;
        //                minx = dx;
        //            }
        //            else if (dx > maxx){
        //                bc.anchors[1] = a;
        //                maxx = dx;
        //            }
        //            if (dy < miny){
        //                bc.anchors[2] = a;
        //                miny = dy;
        //            }
        //            else if (dy > maxy){
        //                bc.anchors[3] = a;
        //                maxy = dy;
        //            }
        //        }
        //        for (auto & a : u2.uncheckedRef<MGUnaryRegion>().normalizedCorners){
        //            double dx = a.dot(x), dy = a.dot(y);
        //            if (dx < minx){
        //                bc.anchors[0] = a;
        //                minx = dx;
        //            }
        //            else if (dx > maxx){
        //                bc.anchors[1] = a;
        //                maxx = dx;
        //            }
        //            if (dy < miny){
        //                bc.anchors[2] = a;
        //                miny = dy;
        //            }
        //            else if (dy > maxy){
        //                bc.anchors[3] = a;
        //                maxy = dy;
        //            }
        //        }
        //        return bc;
        //    }
        //    else if (mg.data(bh).is<MGBinaryRegionLineConnection>()){
        //        BinaryAnchors bc;
        //        bc.anchors = NumFilter(mg.data(bh).uncheckedRef<MGBinaryRegionLineConnection>().samples, maxSampleNum);
        //        bc.weight = 1;
        //        return bc;
        //    }
        //    else if (mg.data(bh).is<MGBinaryLineLineIntersection>()){
        //        BinaryAnchors bc;
        //        auto & llInter = mg.data(bh).uncheckedRef<MGBinaryLineLineIntersection>();
        //        bc.weight = 80;
        //        bc.anchors = { llInter.relationCenter };
        //        return bc;
        //    }
        //    else if (mg.data(bh).is<MGBinaryLineLineIncidence>()){
        //        BinaryAnchors bc;
        //        auto & llIncid = mg.data(bh).uncheckedRef<MGBinaryLineLineIncidence>();
        //        bc.weight = 80;
        //        bc.anchors = { llIncid.relationCenter };
        //        return bc;
        //    }

        //    return BinaryAnchors();
        //}
        //
        //std::vector<BinaryAnchors> GetAllBinaryAnchors(const MixedGraph & mg){
        //    std::vector<BinaryAnchors> binaryAnchors(mg.internalElements<1>().size());
        //    int connectionConsNum = 0;
        //    for (int i = 0; i < mg.internalElements<1>().size(); i++){
        //        binaryAnchors[i] = GetBinaryAnchors(mg, MGBinaryHandle(i));
        //        binaryAnchors[i].positionOfFirstAnchor = connectionConsNum;
        //        connectionConsNum += binaryAnchors[i].anchors.size();
        //        binaryAnchors[i].slackValues.resize(binaryAnchors[i].anchors.size());
        //        binaryAnchors[i].currentAnchorDepthsOnUnaries.front().resize(binaryAnchors[i].anchors.size());
        //        binaryAnchors[i].currentAnchorDepthsOnUnaries.back().resize(binaryAnchors[i].anchors.size());
        //    }
        //    return binaryAnchors;
        //}

        //struct UnaryAnchors{
        //    std::vector<Vec3> anchors;
        //    std::vector<std::vector<double>> anchorDepthRatiosForVPs;
        //    std::vector<double> currentAnchorDepths;
        //};

        //std::vector<UnaryAnchors> GetAllUnaryAnchors(const MixedGraph & mg, const std::vector<Vec3> & vps, 
        //    const std::vector<BinaryAnchors> & binaryAnchors){
        //    std::vector<UnaryAnchors> unaryAnchors(mg.internalElements<0>().size());
        //    for (int i = 0; i < mg.internalElements<0>().size(); i++){
        //        unaryAnchors[i].anchorDepthRatiosForVPs.resize(vps.size());
        //        // collect all anchors
        //        for (auto & bh : mg.topo(MGUnaryHandle(i)).uppers){
        //            unaryAnchors[i].anchors.insert(unaryAnchors[i].anchors.end(),
        //                binaryAnchors[bh.id].anchors.begin(), binaryAnchors[bh.id].anchors.end());
        //        }
        //        // compute anchor depth ratios for each vpid
        //        for (int vpid = 0; vpid < vps.size(); vpid++){
        //            auto & anchorDepthRatios = unaryAnchors[i].anchorDepthRatiosForVPs[vpid];
        //            anchorDepthRatios.reserve(unaryAnchors[i].anchors.size());
        //            for (auto & a : unaryAnchors[i].anchors){
        //                anchorDepthRatios.push_back(DepthRatioOnMGUnaryWithAssignedClass(a, mg.data(MGUnaryHandle(i)), vps, vpid));
        //            }
        //        }
        //        unaryAnchors[i].currentAnchorDepths.resize(unaryAnchors[i].anchors.size());
        //        for (int k = 0; k < unaryAnchors[i].anchors.size(); k++){
        //            unaryAnchors[i].currentAnchorDepths[k] = FirstCornerDepthOfMGUnary(mg.data(MGUnaryHandle(i))) *
        //                DepthRatioOnMGUnary(unaryAnchors[i].anchors[k], mg.data(MGUnaryHandle(i)), vps);
        //        }
        //    }
        //    return unaryAnchors;
        //}


        //static void MSKAPI printstr(void *handle, MSKCONST char str[]) {
        //    printf("%s", str);
        //}

        //class MOSEKTaskForMixedGraph {
        //public:
        //    explicit MOSEKTaskForMixedGraph(MixedGraph & g,
        //        const std::vector<Vec3> & v, 
        //        const std::unordered_map<MGUnaryHandle, int> & c,
        //        double dlb = 0.2, double dub = 0.5)
        //        : mg(g), vps(v), ccids(c), depthLb(dlb), depthUb(dub) {
        //        mg.gc();
        //        initialize();
        //    }

        //    ~MOSEKTaskForMixedGraph() {
        //        finalize();
        //    }

        //public:
        //    double optimizeDepths() {
        //        auto tick = Tick("Solve Equations");
        //        MSKrescodee trmcode;
        //        MSK_optimizetrm(task, &trmcode);
        //        MSK_solutionsummary(task, MSK_STREAM_LOG);

        //        MSKsolstae solsta;
        //        MSK_getsolsta(task, MSK_SOL_BAS, &solsta);

        //        Tock(tick);

        //        switch (solsta){
        //        case MSK_SOL_STA_OPTIMAL:
        //        case MSK_SOL_STA_NEAR_OPTIMAL:
        //        {
        //                                         double *xx = new double[varNum];
        //                                         MSK_getxx(task,
        //                                             MSK_SOL_BAS,    /* Request the basic solution. */
        //                                             xx);
        //                                         printf("Optimal primal solution\n");
        //                                         for (int i = 0; i < mg.internalElements<0>().size(); i++){
        //                                             FirstCornerDepthOfMGUnary(mg.data(MGUnaryHandle(i))) = xx[i];
        //                                         }
        //                                         int slackVarId = mg.internalElements<0>().size();
        //                                         for (auto & ba : binaryAnchors){ // get slack values
        //                                             for (double & slack : ba.slackValues){
        //                                                 slack = xx[slackVarId];
        //                                                 slackVarId++;
        //                                             }
        //                                         }
        //                                         for (int i = 0; i < mg.internalElements<1>().size(); i++){
        //                                             auto & uhs = mg.topo(MGBinaryHandle(i)).lowers;
        //                                             auto & ba = binaryAnchors[i];
        //                                             for (int k = 0; k < ba.anchors.size(); k++){
        //                                                 ba.currentAnchorDepthsOnUnaries.front()[k] =
        //                                                     FirstCornerDepthOfMGUnary(mg.data(uhs.front())) *
        //                                                     DepthRatioOnMGUnary(ba.anchors[k], mg.data(uhs.front()), vps);
        //                                                 ba.currentAnchorDepthsOnUnaries.back()[k] =
        //                                                     FirstCornerDepthOfMGUnary(mg.data(uhs.back())) *
        //                                                     DepthRatioOnMGUnary(ba.anchors[k], mg.data(uhs.back()), vps);
        //                                                 assert(FuzzyEquals(abs(
        //                                                     ba.currentAnchorDepthsOnUnaries.front()[k] - 
        //                                                     ba.currentAnchorDepthsOnUnaries.back()[k]),
        //                                                     ba.slackValues[k], 0.1));
        //                                             }
        //                                         }
        //                                         for (int i = 0; i < mg.internalElements<0>().size(); i++){
        //                                             for (int k = 0; k < unaryAnchors[i].anchors.size(); k++){
        //                                                 unaryAnchors[i].currentAnchorDepths[k] = FirstCornerDepthOfMGUnary(mg.data(MGUnaryHandle(i))) *
        //                                                     unaryAnchors[i].anchorDepthRatiosForVPs[OrientationClassOfMGUnary(mg.data(MGUnaryHandle(i)))][k];
        //                                                 assert(unaryAnchors[i].anchorDepthRatiosForVPs[OrientationClassOfMGUnary(mg.data(MGUnaryHandle(i)))][k]
        //                                                     == DepthRatioOnMGUnary(unaryAnchors[i].anchors[k], mg.data(MGUnaryHandle(i)), vps));
        //                                             }
        //                                         }
        //                                         assert(slackVarId == varNum);
        //                                         delete[] xx;
        //                                         double slackSum = 0.0;
        //                                         for (auto & ba : binaryAnchors){
        //                                             slackSum += std::accumulate(ba.slackValues.begin(), ba.slackValues.end(), 0.0);
        //                                         }
        //                                         return slackSum;
        //                                         //break;
        //        }
        //        case MSK_SOL_STA_DUAL_INFEAS_CER:
        //        case MSK_SOL_STA_PRIM_INFEAS_CER:
        //        case MSK_SOL_STA_NEAR_DUAL_INFEAS_CER:
        //        case MSK_SOL_STA_NEAR_PRIM_INFEAS_CER:
        //            printf("Primal or dual infeasibility certificate found.\n");
        //            break;
        //        case MSK_SOL_STA_UNKNOWN:
        //        {
        //                                    char symname[MSK_MAX_STR_LEN];
        //                                    char desc[MSK_MAX_STR_LEN];

        //                                    /* If the solutions status is unknown, print the termination code
        //                                    indicating why the optimizer terminated prematurely. */

        //                                    MSK_getcodedesc(trmcode,
        //                                        symname,
        //                                        desc);

        //                                    printf("The solution status is unknown.\n");
        //                                    printf("The optimizer terminitated with code: %s\n", symname);
        //                                    break;
        //        }
        //        default:
        //            printf("Other solution status.\n");
        //            break;
        //        }

        //        return std::numeric_limits<double>::max();
        //    }

        //    void adjustRegionOrientations() {
        //        std::vector<double> binarySlackSums(binaryAnchors.size());
        //        for (int i = 0; i < binaryAnchors.size(); i++){
        //            binarySlackSums[i] = std::accumulate(binaryAnchors[i].slackValues.begin(),
        //                binaryAnchors[i].slackValues.end(), 0.0);
        //        }
        //        std::vector<int> unaryRegionIds;
        //        unaryRegionIds.reserve(mg.internalElements<0>().size());
        //        for (int i = 0; i < mg.internalElements<0>().size(); i++){
        //            if (mg.data(MGUnaryHandle(i)).is<MGUnaryLine>())
        //                continue;
        //            unaryRegionIds.push_back(i);
        //        }

        //        std::vector<double> unaryErrorScores(mg.internalElements<0>().size(), 0.0);
        //        //static const double ignoredRatio = 0.25;
        //        for (int i : unaryRegionIds){
        //            auto & relatedBhs = mg.topo(MGUnaryHandle(i)).uppers;
        //            double allSlackSum = 0.0;
        //            double allSlackNum = 0;
        //            for (auto & bh : relatedBhs){
        //                allSlackSum += binarySlackSums[bh.id];
        //                allSlackNum += binaryAnchors[bh.id].anchors.size();
        //            }
        //            //std::vector<double> relatedBinarySlackSumaries;
        //            //relatedBinarySlackSumaries.reserve(relatedBhs.size());
        //            //for (auto & bh : relatedBhs){
        //            //    relatedBinarySlackSumaries.push_back(binarySlackSums[bh.id]/*/binaryAnchors[bh.id].anchors.size()*/);
        //            //}
        //            //std::sort(relatedBinarySlackSumaries.begin(), relatedBinarySlackSumaries.end());
        //            //double errorScore = std::accumulate(relatedBinarySlackSumaries.begin(),
        //            //    relatedBinarySlackSumaries.begin() + (relatedBinarySlackSumaries.size() + 1) * (1.0 - ignoredRatio), 0.0)
        //            //    / ((relatedBinarySlackSumaries.size() + 1) * (1.0 - ignoredRatio));
        //            //unaryErrorScores[i] = errorScore;
        //            unaryErrorScores[i] = allSlackSum / allSlackNum;
        //        }
        //        
        //        std::sort(unaryRegionIds.begin(), unaryRegionIds.end(), [&unaryErrorScores](int uhid1, int uhid2){
        //            return unaryErrorScores[uhid1] > unaryErrorScores[uhid2];
        //        });

        //        // adjust orientations now
        //        int count = 0;

        //        for (int uhid : unaryRegionIds){
        //            if (count++ > unaryRegionIds.size() / 2)
        //                break;

        //            // get neighbor anchor depths
        //            std::vector<double> neighborAnchorDepths;
        //            std::vector<double> neighborAnchorDepthWeights;

        //            neighborAnchorDepths.reserve(unaryAnchorIdsTable[uhid].size() / 2);
        //            neighborAnchorDepthWeights.reserve(unaryAnchorIdsTable[uhid].size() / 2);

        //            auto & relatedBhs = mg.topo(MGUnaryHandle(uhid)).uppers;
        //            for (auto & bh : relatedBhs){
        //                bool isAtFront = uhid == mg.topo(bh).lowers.front().id;
        //                auto & neighborAnchorDepthsOnThisBh = isAtFront ? 
        //                    binaryAnchors[bh.id].currentAnchorDepthsOnUnaries.back() :
        //                    binaryAnchors[bh.id].currentAnchorDepthsOnUnaries.front();
        //                neighborAnchorDepths.insert(neighborAnchorDepths.end(), 
        //                    neighborAnchorDepthsOnThisBh.begin(), neighborAnchorDepthsOnThisBh.end());
        //                neighborAnchorDepthWeights.insert(neighborAnchorDepthWeights.end(),
        //                    neighborAnchorDepthsOnThisBh.size(), binaryAnchors[bh.id].weight);
        //            }

        //            int bestVPId = -1;
        //            double minError = std::numeric_limits<double>::max();

        //            for (int vpid = 0; vpid < vps.size(); vpid++){
        //                const std::vector<double> & thisAnchorDepths = unaryAnchors[uhid].anchorDepthRatiosForVPs[vpid];
        //                assert(thisAnchorDepths.size() == neighborAnchorDepths.size());

        //                // find the min distance between |thisAnchorDepths * x - neighborAnchorDepths|
        //                Eigen::Map<const Eigen::VectorXd> neighborAnchorDepthsVector(neighborAnchorDepths.data(), neighborAnchorDepths.size());
        //                Eigen::Map<const Eigen::VectorXd> neighborAnchorDepthsWeightsVector(neighborAnchorDepthWeights.data(), neighborAnchorDepthWeights.size());
        //                Eigen::Map<const Eigen::VectorXd> thisAnchorDepthsVector(thisAnchorDepths.data(), thisAnchorDepths.size());
        //                double errorWithThisVPId = -Square(neighborAnchorDepthsVector.cwiseProduct(neighborAnchorDepthsWeightsVector)
        //                    .dot(thisAnchorDepthsVector.cwiseProduct(neighborAnchorDepthsWeightsVector).normalized()));

        //                if (errorWithThisVPId < minError){
        //                    bestVPId = vpid;
        //                    minError = errorWithThisVPId;
        //                }
        //            }

        //            setUnaryOrientationClass(MGUnaryHandle(uhid), bestVPId);
        //        }                
        //    }

        //private:
        //    void initialize() {
        //        int depthNum = mg.internalElements<0>().size();
        //        connectionConsNum = 0;

        //        binaryAnchors = GetAllBinaryAnchors(mg);
        //        for (int i = 0; i < binaryAnchors.size(); i++){
        //            connectionConsNum += binaryAnchors[i].anchors.size();
        //        }
        //        unaryAnchors = GetAllUnaryAnchors(mg, vps, binaryAnchors);

        //        // get first uh in each cc
        //        std::map<int, MGUnaryHandle> firstUhsInEachCC;
        //        firstCCUhs.clear();
        //        for (auto & uhCC : ccids){
        //            if (!Contains(firstUhsInEachCC, uhCC.second)){
        //                firstUhsInEachCC[uhCC.second] = uhCC.first;
        //                firstCCUhs.insert(uhCC.first);
        //            }
        //        }

        //        // additional consNum for anchors
        //        int ccConsNum = firstCCUhs.size();

        //        int slackVarNum = connectionConsNum;

        //        varNum = depthNum + slackVarNum; // depths, slackVars
        //        consNum = connectionConsNum * 2; // depth1 * ratio1 - depth2 * ratio2 < slackVar;  
        //        // depth2 * ratio2 - depth1 * ratio1 < slackVar
        //        //ccConsNum +  // depth = defaultValue
        //        //depthNum; // depth \in [depthLb, depthUb]

        //        env = nullptr;
        //        task = nullptr;

        //        MSK_makeenv(&env, nullptr);
        //        MSK_maketask(env, consNum, varNum, &task);
        //        //MSK_linkfunctotaskstream(task, MSK_STREAM_LOG, NULL, printstr);
        //        MSK_appendcons(task, consNum);
        //        MSK_appendvars(task, varNum);

        //        int slackVarId = 0;
        //        for (int i = 0; i < binaryAnchors.size(); i++){
        //            for (auto & a : binaryAnchors[i].anchors){
        //                MSK_putcj(task, slackVarId, binaryAnchors[i].weight);
        //                slackVarId++;
        //            }
        //        }

        //        // bounds for vars
        //        int varId = 0;
        //        for (; varId < depthNum; varId++){
        //            // depth \in  [depthLb, depthUb]
        //            if (!Contains(firstCCUhs, MGUnaryHandle(varId))){
        //                MSK_putvarbound(task, varId, MSK_BK_RA, depthLb, depthUb);
        //            }
        //            else{ // depth = defaultValue
        //                double defaultDepth = FirstCornerDepthOfMGUnary(mg.data(MGUnaryHandle(varId)));
        //                MSK_putvarbound(task, varId, MSK_BK_FX, defaultDepth, defaultDepth);
        //            }
        //        }
        //        for (; varId < depthNum + slackVarNum; varId++){
        //            // slack >= 0
        //            MSK_putvarbound(task, varId, MSK_BK_LO, 0.0, +MSK_INFINITY);
        //        }

        //        // parameters for each var in constraints
        //        unaryDepthsTable = std::vector<std::vector<std::vector<MSKrealt>>>(depthNum);
        //        unaryAnchorIdsTable = std::vector<std::vector<MSKint32t>>(depthNum);

        //        varId = 0;
        //        for (; varId < depthNum; varId++){ // depths
        //            auto & relatedBhs = mg.topo(MGUnaryHandle(varId)).uppers;
        //            int relatedAnchorsNum = 0;
        //            for (auto & bh : relatedBhs){
        //                relatedAnchorsNum += binaryAnchors[bh.id].anchors.size();
        //            }
        //            int relatedConsNum = relatedAnchorsNum * 2;

        //            // collect consIds
        //            std::vector<MSKint32t> & consIds = unaryAnchorIdsTable[varId];
        //            consIds.reserve(relatedConsNum);
        //            for (auto & bh : relatedBhs){
        //                int firstAnchorPosition = binaryAnchors[bh.id].positionOfFirstAnchor;
        //                for (int k = 0; k < binaryAnchors[bh.id].anchors.size(); k++){
        //                    consIds.push_back((firstAnchorPosition + k) * 2); // one for [depth1 * ratio1 - depth2 * ratio2 - slackVar < 0]; 
        //                    consIds.push_back((firstAnchorPosition + k) * 2 + 1); // another for [- depth1 * ratio1 + depth2 * ratio2 - slackVar < 0];
        //                }
        //            }

        //            auto & consValueTable = unaryDepthsTable[varId];

        //            consValueTable.resize(vps.size());

        //            for (int vpid = 0; vpid < vps.size(); vpid++){
        //                auto & consValues = consValueTable[vpid];
        //                consValues.reserve(relatedConsNum);

        //                for (auto & bh : relatedBhs){
        //                    int firstAnchorPosition = binaryAnchors[bh.id].positionOfFirstAnchor;
        //                    for (int k = 0; k < binaryAnchors[bh.id].anchors.size(); k++){
        //                        consIds.push_back((firstAnchorPosition + k) * 2); // one for [depth1 * ratio1 - depth2 * ratio2 - slackVar < 0]; 
        //                        consIds.push_back((firstAnchorPosition + k) * 2 + 1); // another for [- depth1 * ratio1 + depth2 * ratio2 - slackVar < 0];

        //                        auto & a = binaryAnchors[bh.id].anchors[k];
        //                        double ratio = DepthRatioOnMGUnaryWithAssignedClass(a, mg.data(MGUnaryHandle(varId)), vps, vpid);
        //                        bool isOnLeftSide = varId == mg.topo(bh).lowers.front().id;
        //                        if (isOnLeftSide){ // as depth1
        //                            consValues.push_back(ratio);
        //                            consValues.push_back(-ratio);
        //                        }
        //                        else{
        //                            consValues.push_back(-ratio);
        //                            consValues.push_back(ratio);
        //                        }
        //                    }
        //                }
        //            }                 

        //            MSK_putacol(task, varId, relatedConsNum, consIds.data(), 
        //                consValueTable[OrientationClassOfMGUnary(mg.data(MGUnaryHandle(varId)))].data());
        //        }
        //        for (; varId < depthNum + slackVarNum; varId++){ // slack vars
        //            MSKint32t consIds[] = {
        //                (varId - depthNum) * 2, // one for [depth1 * ratio1 - depth2 * ratio2 - slackVar < 0]; 
        //                (varId - depthNum) * 2 + 1  // another for [- depth1 * ratio1 + depth2 * ratio2 - slackVar < 0];
        //            };
        //            MSKrealt consValues[] = { -1.0, -1.0 };
        //            MSK_putacol(task, varId, 2, consIds, consValues);
        //        }

        //        // bounds for constraints
        //        for (int consId = 0; consId < consNum; consId++){
        //            // all [depth1 * ratio1 - depth2 * ratio2 - slackVar < 0];
        //            MSK_putconbound(task, consId, MSK_BK_UP, -MSK_INFINITY, 0.0);
        //        }

        //        MSK_putobjsense(task, MSK_OBJECTIVE_SENSE_MINIMIZE);
        //    }
        //    void finalize(){
        //        // finalize
        //        MSK_deletetask(&task);
        //        MSK_deleteenv(&env);
        //    }
        //    void setUnaryOrientationClass(const MGUnaryHandle & uh, int claz){
        //        if (OrientationClassOfMGUnary(mg.data(uh)) == claz)
        //            return;
        //        OrientationClassOfMGUnary(mg.data(uh)) = claz;
        //        auto & consIds = unaryAnchorIdsTable[uh.id];
        //        auto & consValues = unaryDepthsTable[uh.id][claz];
        //        MSK_putacol(task, uh.id, consValues.size(), consIds.data(), consValues.data());
        //    }

        //private:
        //    MixedGraph & mg;
        //    const std::vector<Vec3> & vps;
        //    const std::unordered_map<MGUnaryHandle, int> & ccids;
        //    double depthLb, depthUb;

        //    std::vector<BinaryAnchors> binaryAnchors;
        //    std::vector<UnaryAnchors> unaryAnchors;

        //    int connectionConsNum;

        //    std::vector<std::vector<std::vector<MSKrealt>>> unaryDepthsTable; // table[unaryId][vpId] -> depths
        //    std::vector<std::vector<MSKint32t>> unaryAnchorIdsTable;

        //    std::set<MGUnaryHandle> firstCCUhs;
        //    
        //    int consNum, varNum;

        //    MSKenv_t env;
        //    MSKtask_t task;
        //};


        //void OptimizeMixedGraph(MixedGraph & mg, const std::vector<Vec3> & vps,
        //    const std::unordered_map<MGUnaryHandle, int> & ccids,
        //    int iterNum,
        //    const std::function<void(const MixedGraph &)> & callbackEachIteration) {

        //    MOSEKTaskForMixedGraph task(mg, vps, ccids);
        //    for (int i = 0; i < iterNum; i++){
        //        std::cout << "iter " << i << std::endl;
        //        double slackSum = task.optimizeDepths();
        //        std::cout << "slack sum is " << slackSum << std::endl;
        //        task.adjustRegionOrientations();
        //        callbackEachIteration(mg);
        //    }

        //}

        /*
        void SolveDepthsInMixedGraph(MixedGraph & mg, const std::vector<Vec3> & vps, 
            const std::unordered_map<MGUnaryHandle, int> & ccids){
            
            using namespace Eigen;
            SparseMatrix<double> A, W;
            VectorXd B;
            VectorXd X;

            mg.gc();
            using namespace Eigen;

            //// PREPARE
            auto tick = Tick("SetUp Equations");
            assert(mg.isDense());

            int varNum = mg.internalElements<0>().size();
            int consNum = 0;

            std::vector<BinaryAnchors> binaryAnchors(mg.internalElements<1>().size());
            for (int i = 0; i < mg.internalElements<1>().size(); i++){
                binaryAnchors[i] = GetBinaryAnchors(mg, MGBinaryHandle(i));
                consNum += binaryAnchors[i].anchors.size();
            }

            // get first uh in each cc
            std::map<int, MGUnaryHandle> firstUhsInEachCC;
            for (auto & uhCC : ccids){
                if (!Contains(firstUhsInEachCC, uhCC.second)){
                    firstUhsInEachCC[uhCC.second] = uhCC.first;
                }
            }
            // additional consNum for anchors
            consNum += firstUhsInEachCC.size();

            A.resize(consNum, varNum);
            W.resize(consNum, consNum);
            B.resize(consNum);

            // write equations
            int eid = 0;
            for (int i = 0; i < binaryAnchors.size(); i++){
                for (auto & a : binaryAnchors[i].anchors){
                    auto uh1 = mg.topo(MGBinaryHandle(i)).lowers.front();
                    auto uh2 = mg.topo(MGBinaryHandle(i)).lowers.back();
                    auto & u1 = mg.data(uh1);
                    auto & u2 = mg.data(uh2);
                    double ratio1 = DepthRatioOnMGUnary(a, u1, vps);
                    double ratio2 = DepthRatioOnMGUnary(a, u2, vps);
                    A.insert(eid, uh1.id) = ratio1;
                    A.insert(eid, uh2.id) = -ratio2;
                    B(eid) = 0.0;
                    W.insert(eid, eid) = binaryAnchors[i].weight;
                    eid++;
                }
            }

            // the anchor equations
            for (auto & uhCCId : firstUhsInEachCC){
                A.insert(eid, 0) = 1.0;
                B(eid) = FirstCornerDepthOfMGUnary(mg.data(uhCCId.second));
                W.insert(eid, eid) = 1.0;
                eid++;
            }
            assert(eid == consNum);
            Tock(tick);


            ///// SOLVE
            X.resize(A.cols());
            for (int i = 0; i < X.size(); i++){
                X(i) = FirstCornerDepthOfMGUnary(mg.data(MGUnaryHandle(i)));
            }
            // solve the equation system
            const bool useWeights = true;
            SparseQR<Eigen::SparseMatrix<double>, COLAMDOrdering<int>> solver;
            static_assert(!(Eigen::SparseMatrix<double>::IsRowMajor), "COLAMDOrdering only supports column major");
            Eigen::SparseMatrix<double> WA = W * A;
            A.makeCompressed();
            WA.makeCompressed();

            tick = Tick("Solve Equations");
            solver.compute(useWeights ? WA : A);
            Tock(tick);

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

            std::cout << "filling back solutions" << std::endl;

            double Xmean = X.mean();
            std::cout << "mean(X) = " << X.mean() << std::endl;
            for (int i = 0; i < X.size(); i++){
                FirstCornerDepthOfMGUnary(mg.data(MGUnaryHandle(i))) = X(i);
            }

        }
        */

        /*
        void SolveDepthsInMixedGraphMOSEKOnce(MixedGraph & mg, const std::vector<Vec3> & vps,
            const std::set<MGUnaryHandle> & firstCCUhs,
            int connectionConsNum,
            std::vector<BinaryAnchors> & binaryAnchors,
            double depthLb, double depthUb){

            assert(mg.isDense());
            int depthNum = mg.internalElements<0>().size();

            // additional consNum for anchors
            int ccConsNum = firstCCUhs.size();

            int slackVarNum = connectionConsNum;

            int varNum = depthNum + slackVarNum; // depths, slackVars
            int consNum = connectionConsNum * 2; // depth1 * ratio1 - depth2 * ratio2 < slackVar;  
            // depth2 * ratio2 - depth1 * ratio1 < slackVar
            //ccConsNum +  // depth = defaultValue
            //depthNum; // depth \in [depthLb, depthUb]

            MSKenv_t env = nullptr;
            MSKtask_t task = nullptr;

            MSK_makeenv(&env, nullptr);
            MSK_maketask(env, consNum, varNum, &task);
            //MSK_linkfunctotaskstream(task, MSK_STREAM_LOG, NULL, printstr);
            MSK_appendcons(task, consNum);
            MSK_appendvars(task, varNum);

            int slackVarId = 0;
            for (int i = 0; i < binaryAnchors.size(); i++){
                for (auto & a : binaryAnchors[i].anchors){
                    MSK_putcj(task, slackVarId, binaryAnchors[i].weight);
                    slackVarId++;
                }
            }

            // bounds for vars
            int varId = 0;
            for (; varId < depthNum; varId++){
                // depth \in  [depthLb, depthUb]
                if (!Contains(firstCCUhs, MGUnaryHandle(varId))){
                    MSK_putvarbound(task, varId, MSK_BK_RA, depthLb, depthUb);
                }
                else{ // depth = defaultValue
                    double defaultDepth = FirstCornerDepthOfMGUnary(mg.data(MGUnaryHandle(varId)));
                    MSK_putvarbound(task, varId, MSK_BK_FX, defaultDepth, defaultDepth);
                }
            }
            for (; varId < depthNum + slackVarNum; varId++){
                // slack >= 0
                MSK_putvarbound(task, varId, MSK_BK_LO, 0.0, +MSK_INFINITY);
            }

            // parameters for each var in constraints
            varId = 0;
            for (; varId < depthNum; varId++){ // depths
                auto & relatedBhs = mg.topo(MGUnaryHandle(varId)).uppers;
                int relatedAnchorsNum = 0;
                for (auto & bh : relatedBhs){
                    relatedAnchorsNum += binaryAnchors[bh.id].anchors.size();
                }
                int relatedConsNum = relatedAnchorsNum * 2;

                std::vector<MSKint32t> consIds;
                std::vector<MSKrealt> consValues;

                consIds.reserve(relatedConsNum);
                consValues.reserve(relatedConsNum);

                for (auto & bh : relatedBhs){
                    int firstAnchorPosition = binaryAnchors[bh.id].positionOfFirstAnchor;
                    for (int k = 0; k < binaryAnchors[bh.id].anchors.size(); k++){
                        consIds.push_back((firstAnchorPosition + k) * 2); // one for [depth1 * ratio1 - depth2 * ratio2 - slackVar < 0]; 
                        consIds.push_back((firstAnchorPosition + k) * 2 + 1); // another for [- depth1 * ratio1 + depth2 * ratio2 - slackVar < 0];

                        auto & a = binaryAnchors[bh.id].anchors[k];
                        double ratio = DepthRatioOnMGUnary(a, mg.data(MGUnaryHandle(varId)), vps);
                        bool isOnLeftSide = varId == mg.topo(bh).lowers.front().id;
                        if (isOnLeftSide){ // as depth1
                            consValues.push_back(ratio);
                            consValues.push_back(-ratio);
                        }
                        else{
                            consValues.push_back(-ratio);
                            consValues.push_back(ratio);
                        }
                    }
                }

                MSK_putacol(task, varId, relatedConsNum, consIds.data(), consValues.data());
            }
            for (; varId < depthNum + slackVarNum; varId++){ // slack vars
                MSKint32t consIds[] = {
                    (varId - depthNum) * 2, // one for [depth1 * ratio1 - depth2 * ratio2 - slackVar < 0]; 
                    (varId - depthNum) * 2 + 1  // another for [- depth1 * ratio1 + depth2 * ratio2 - slackVar < 0];
                };
                MSKrealt consValues[] = { -1.0, -1.0 };

                MSK_putacol(task, varId, 2, consIds, consValues);
            }

            // bounds for constraints
            for (int consId = 0; consId < consNum; consId++){
                // all [depth1 * ratio1 - depth2 * ratio2 - slackVar < 0];
                MSK_putconbound(task, consId, MSK_BK_UP, -MSK_INFINITY, 0.0);
            }

            MSK_putobjsense(task, MSK_OBJECTIVE_SENSE_MINIMIZE);
            MSKrescodee trmcode;
            MSK_optimizetrm(task, &trmcode);
            MSK_solutionsummary(task, MSK_STREAM_LOG);

            MSKsolstae solsta;
            MSK_getsolsta(task, MSK_SOL_BAS, &solsta);

            switch (solsta){
            case MSK_SOL_STA_OPTIMAL:
            case MSK_SOL_STA_NEAR_OPTIMAL:
            {
                                             double *xx = new double[varNum];
                                             MSK_getxx(task,
                                                 MSK_SOL_BAS,    
                                                 xx);
                                             printf("Optimal primal solution\n");
                                             for (int i = 0; i < depthNum; i++){ // get resulted depths
                                                 FirstCornerDepthOfMGUnary(mg.data(MGUnaryHandle(i))) = xx[i];
                                             }
                                             int slackVarId = depthNum;
                                             for (auto & ba : binaryAnchors){ // get slack values
                                                 for (double & slack : ba.slackValues){
                                                     slack = xx[slackVarId];
                                                     slackVarId++;
                                                 }
                                             }
                                             assert(slackVarId == varNum);
                                             delete[] xx;
                                             break;
            }
            case MSK_SOL_STA_DUAL_INFEAS_CER:
            case MSK_SOL_STA_PRIM_INFEAS_CER:
            case MSK_SOL_STA_NEAR_DUAL_INFEAS_CER:
            case MSK_SOL_STA_NEAR_PRIM_INFEAS_CER:
                printf("Primal or dual infeasibility certificate found.\n");
                break;
            case MSK_SOL_STA_UNKNOWN:
            {
                                        char symname[MSK_MAX_STR_LEN];
                                        char desc[MSK_MAX_STR_LEN];

                                     

                                        MSK_getcodedesc(trmcode,
                                            symname,
                                            desc);

                                        printf("The solution status is unknown.\n");
                                        printf("The optimizer terminitated with code: %s\n", symname);
                                        break;
            }
            default:
                printf("Other solution status.\n");
                break;
            }


            // finalize
            MSK_deletetask(&task);
            MSK_deleteenv(&env);

        }
        */
        
        /*void SolveDepthsInMixedGraphMOSEK(MixedGraph & mg, const std::vector<Vec3> & vps,
            const std::unordered_map<MGUnaryHandle, int> & ccids,
            double depthLb, double depthUb) {

            using namespace Eigen;
            mg.gc();

            //// PREPARE
            auto tick = Tick("SetUp Equations");
            assert(mg.isDense());

            int depthNum = mg.internalElements<0>().size();
            int connectionConsNum = 0;

            std::vector<BinaryAnchors> binaryAnchors(mg.internalElements<1>().size());
            std::vector<int> positionOfFirstAnchorInEachBinaryAnchors(mg.internalElements<1>().size());
            for (int i = 0; i < mg.internalElements<1>().size(); i++){
                binaryAnchors[i] = GetBinaryAnchors(mg, MGBinaryHandle(i));
                positionOfFirstAnchorInEachBinaryAnchors[i] = connectionConsNum;
                connectionConsNum += binaryAnchors[i].anchors.size();
            }

            // get first uh in each cc
            std::map<int, MGUnaryHandle> firstUhsInEachCC;
            std::set<MGUnaryHandle> firstCCUhs;
            for (auto & uhCC : ccids){
                if (!Contains(firstUhsInEachCC, uhCC.second)){
                    firstUhsInEachCC[uhCC.second] = uhCC.first;
                    firstCCUhs.insert(uhCC.first);
                }
            }

            // additional consNum for anchors
            int ccConsNum = firstUhsInEachCC.size();            
            
            int slackVarNum = connectionConsNum;

            int varNum = depthNum + slackVarNum; // depths, slackVars
            int consNum = connectionConsNum * 2; // depth1 * ratio1 - depth2 * ratio2 < slackVar;  
                // depth2 * ratio2 - depth1 * ratio1 < slackVar
                //ccConsNum +  // depth = defaultValue
                //depthNum; // depth \in [depthLb, depthUb]

            MSKenv_t env = nullptr;
            MSKtask_t task = nullptr;

            MSK_makeenv(&env, nullptr);
            MSK_maketask(env, consNum, varNum, &task);
            //MSK_linkfunctotaskstream(task, MSK_STREAM_LOG, NULL, printstr);
            MSK_appendcons(task, consNum);
            MSK_appendvars(task, varNum);

            int slackVarId = 0;
            for (int i = 0; i < binaryAnchors.size(); i++){
                for (auto & a : binaryAnchors[i].anchors){
                    MSK_putcj(task, slackVarId, binaryAnchors[i].weight);
                    slackVarId++;
                }
            }

            // bounds for vars
            int varId = 0;
            for (; varId < depthNum; varId++){
                // depth \in  [depthLb, depthUb]
                if (!Contains(firstCCUhs, MGUnaryHandle(varId))){
                    MSK_putvarbound(task, varId, MSK_BK_RA, depthLb, depthUb);
                }
                else{ // depth = defaultValue
                    double defaultDepth = FirstCornerDepthOfMGUnary(mg.data(MGUnaryHandle(varId)));
                    MSK_putvarbound(task, varId, MSK_BK_FX, defaultDepth, defaultDepth);
                }
            }
            for (; varId < depthNum + slackVarNum; varId++){
                // slack >= 0
                MSK_putvarbound(task, varId, MSK_BK_LO, 0.0, +MSK_INFINITY);
            }

            // parameters for each var in constraints
            varId = 0;
            for (; varId < depthNum; varId++){ // depths
                auto & relatedBhs = mg.topo(MGUnaryHandle(varId)).uppers;
                int relatedAnchorsNum = 0;
                for (auto & bh : relatedBhs){
                    relatedAnchorsNum += binaryAnchors[bh.id].anchors.size();
                }
                int relatedConsNum = relatedAnchorsNum * 2;
                
                std::vector<MSKint32t> consIds;
                std::vector<MSKrealt> consValues;
                consIds.reserve(relatedConsNum);
                consValues.reserve(relatedConsNum);

                for (auto & bh : relatedBhs){
                    int firstAnchorPosition = positionOfFirstAnchorInEachBinaryAnchors[bh.id];
                    for (int k = 0; k < binaryAnchors[bh.id].anchors.size(); k++){
                        consIds.push_back((firstAnchorPosition + k) * 2); // one for [depth1 * ratio1 - depth2 * ratio2 - slackVar < 0]; 
                        consIds.push_back((firstAnchorPosition + k) * 2 + 1); // another for [- depth1 * ratio1 + depth2 * ratio2 - slackVar < 0];

                        auto & a = binaryAnchors[bh.id].anchors[k];
                        double ratio = DepthRatioOnMGUnary(a, mg.data(MGUnaryHandle(varId)), vps);
                        bool isOnLeftSide = varId == mg.topo(bh).lowers.front().id;
                        if (isOnLeftSide){ // as depth1
                            consValues.push_back(ratio);
                            consValues.push_back(-ratio);
                        }
                        else{
                            consValues.push_back(-ratio);
                            consValues.push_back(ratio);
                        }
                    }
                }

                MSK_putacol(task, varId, relatedConsNum, consIds.data(), consValues.data());
            }
            for (; varId < depthNum + slackVarNum; varId++){ // slack vars
                MSKint32t consIds[] = { 
                    (varId - depthNum) * 2, // one for [depth1 * ratio1 - depth2 * ratio2 - slackVar < 0]; 
                    (varId - depthNum) * 2 + 1  // another for [- depth1 * ratio1 + depth2 * ratio2 - slackVar < 0];
                };
                MSKrealt consValues[] = { -1.0, -1.0 };

                MSK_putacol(task, varId, 2, consIds, consValues);
            }

            // bounds for constraints
            for (int consId = 0; consId < consNum; consId++){
                // all [depth1 * ratio1 - depth2 * ratio2 - slackVar < 0];
                MSK_putconbound(task, consId, MSK_BK_UP, -MSK_INFINITY, 0.0);
            }

            Tock(tick);

            tick = Tick("Solve Equations");
            MSK_putobjsense(task, MSK_OBJECTIVE_SENSE_MINIMIZE);
            MSKrescodee trmcode;
            MSK_optimizetrm(task, &trmcode);
            MSK_solutionsummary(task, MSK_STREAM_LOG);
            
            MSKsolstae solsta;
            MSK_getsolsta(task, MSK_SOL_BAS, &solsta);

            Tock(tick);

            switch (solsta){
            case MSK_SOL_STA_OPTIMAL:
            case MSK_SOL_STA_NEAR_OPTIMAL:
            {
                                             double *xx = new double[varNum];
                                             MSK_getxx(task,
                                                 MSK_SOL_BAS,    
                                                 xx);
                                             printf("Optimal primal solution\n");
                                             for (int i = 0; i < depthNum; i++){
                                                 FirstCornerDepthOfMGUnary(mg.data(MGUnaryHandle(i))) = xx[i];
                                             }
                                             delete[] xx;
                                             break;
            }
            case MSK_SOL_STA_DUAL_INFEAS_CER:
            case MSK_SOL_STA_PRIM_INFEAS_CER:
            case MSK_SOL_STA_NEAR_DUAL_INFEAS_CER:
            case MSK_SOL_STA_NEAR_PRIM_INFEAS_CER:
                printf("Primal or dual infeasibility certificate found.\n");
                break;
            case MSK_SOL_STA_UNKNOWN:
            {
                                        char symname[MSK_MAX_STR_LEN];
                                        char desc[MSK_MAX_STR_LEN];

                                        MSK_getcodedesc(trmcode,
                                            symname,
                                            desc);

                                        printf("The solution status is unknown.\n");
                                        printf("The optimizer terminitated with code: %s\n", symname);
                                        break;
            }
            default:
                printf("Other solution status.\n");
                break;
            }


            // finalize
            MSK_deletetask(&task);
            MSK_deleteenv(&env);

        }
        */


    }
}