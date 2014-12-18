
extern "C" {
    #include <gpc.h>
    #include <mosek.h>
}

#include <Eigen/Dense>
#include <Eigen/Sparse>
//#include <Eigen/SPQRSupport>

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
                        for (auto & p : rd.dilatedContours.front()) {
                            auto direction = views[i].camera.spatialDirection(p);
                            spatialContour.push_back(direction / norm(direction));
                        }
                    }
                    else{
                        std::cerr << "this region has no dilatedCountour!" << std::endl;
                    }
                    regionSpatialContours[ri] = std::move(spatialContour);
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
            return Plane3(region.normalizedCenter * var.depthOfCenter, vps[var.claz]);
        }

        Line3 LineOfMGUnary(const MGUnary & line, const std::vector<Vec3> & vps, const MGUnaryVariable & var) {
            assert(line.type == MGUnary::Line);
            InfiniteLine3 infLine(line.normalizedCenter * var.depthOfCenter, vps[var.claz]);
            return Line3(DistanceBetweenTwoLines(InfiniteLine3(Point3(0, 0, 0), line.normalizedCorners.front()), infLine).second.second,
                DistanceBetweenTwoLines(InfiniteLine3(Point3(0, 0, 0), line.normalizedCorners.back()), infLine).second.second);
        }

        double DepthRatioOnMGUnary(const Vec3 & direction, const MGUnary & unary, const std::vector<Vec3> & vps, int claz) {
            if (unary.type == MGUnary::Region){
                const auto & region = unary;
                assert(!region.normalizedCorners.empty());
                Plane3 plane(region.normalizedCenter, vps[claz]);
                return norm(IntersectionOfLineAndPlane(InfiniteLine3(Point3(0, 0, 0), direction), plane).position);
            }
            else if (unary.type == MGUnary::Line){
                const auto & line = unary;
                InfiniteLine3 infLine(line.normalizedCenter, vps[claz]);
                return norm(DistanceBetweenTwoLines(InfiniteLine3(Point3(0, 0, 0), direction), infLine).second.first);
            }
        }

        Point3 LocationOnMGUnary(const Vec3 & direction, const MGUnary & unary, const std::vector<Vec3> & vps,
            const MGUnaryVariable & var){
            if (unary.type == MGUnary::Region){
                const auto & region = unary;
                assert(!region.normalizedCorners.empty());
                Plane3 plane(region.normalizedCenter * var.depthOfCenter, vps[var.claz]);
                return IntersectionOfLineAndPlane(InfiniteLine3(Point3(0, 0, 0), direction), plane).position;
            }
            else if (unary.type == MGUnary::Line){
                const auto & line = unary;
                InfiniteLine3 infLine(line.normalizedCenter * var.depthOfCenter, vps[var.claz]);
                return DistanceBetweenTwoLines(InfiniteLine3(Point3(0, 0, 0), direction), infLine).second.first;
            }
        }




        void InitializeUnaryVarDepths(MGUnaryVarTable & unaryVars, double depth){
            for (auto & uhv : unaryVars){
                uhv.second.depthOfCenter = depth;
            }
        }

        void UpdateBinaryVars(const MixedGraph & mg, const std::vector<Vec3> & vps, 
            const MGUnaryVarTable & unaryVars, MGBinaryVarTable & binaryVars){
            for (auto & bv : binaryVars){
                auto & uhs = mg.topo(bv.first).lowers;
                auto & normalizedAnchors = mg.data(bv.first).normalizedAnchors;
                for (int i = 0; i < normalizedAnchors.size(); i++){
                    bv.second.sampleDepthsOnRelatedUnaries.front()[i] = unaryVars.at(uhs.front()).depthOfCenter *
                        DepthRatioOnMGUnary(normalizedAnchors[i], mg.data(uhs.front()), vps, unaryVars.at(uhs.front()).claz);
                    bv.second.sampleDepthsOnRelatedUnaries.back()[i] = unaryVars.at(uhs.back()).depthOfCenter *
                        DepthRatioOnMGUnary(normalizedAnchors[i], mg.data(uhs.back()), vps, unaryVars.at(uhs.back()).claz);
                }
            }
        }

        bool IsBadBinary(const MixedGraph & mg, const MGBinaryHandle & bh, 
            const MGUnaryVarTable & unaryVars, const std::vector<Vec3> & vps){
            if (mg.data(bh).normalizedAnchors.size() == 1)
                return false;
            auto & uhs = mg.topo(bh).lowers;
            if (mg.data(bh).type == MGBinary::RegionLineConnection){
                return !IsFuzzyPerpenducular(vps[unaryVars.at(uhs.front()).claz], 
                    vps[unaryVars.at(uhs.back()).claz], 1e-6);
            }
            if (mg.data(bh).type == MGBinary::RegionRegionOverlapping)
                return unaryVars.at(uhs.front()).claz != unaryVars.at(uhs.back()).claz;
            assert(mg.data(bh).type == MGBinary::RegionRegionConnection);
            if (unaryVars.at(uhs.front()).claz == unaryVars.at(uhs.back()).claz)
                return false;
            // is good only if all anchors are aligned
            Vec3 alignDir = /*mg.data(bh).type == MGBinary::RegionLineConnection ? vps[unaryVars.at(uhs.back()).claz] :*/
                vps[unaryVars.at(uhs.front()).claz].cross(vps[unaryVars.at(uhs.back()).claz]);
            for (int k = 1; k < mg.data(bh).normalizedAnchors.size(); k++){
                Vec3 rotateNorm = mg.data(bh).normalizedAnchors.front().cross(mg.data(bh).normalizedAnchors[k]);
                if (!IsFuzzyPerpenducular(rotateNorm, alignDir, 1e-6))
                    return true;
            }
            return false;
        }


        double ConsistencyOfBinary(const MixedGraph & mg, const MGBinaryHandle & bh,
            const MGUnaryVarTable & unaryVars, const std::vector<Vec3> & vps){
            if (mg.data(bh).normalizedAnchors.size() == 1)
                return 1.0;
            auto & uhs = mg.topo(bh).lowers;
            if (mg.data(bh).type == MGBinary::RegionLineConnection){
                return 1.0 - abs(normalize(vps[unaryVars.at(uhs.front()).claz])
                    .dot(normalize(vps[unaryVars.at(uhs.back()).claz])));
            }
            if (mg.data(bh).type == MGBinary::RegionRegionOverlapping){
                return abs(normalize(vps[unaryVars.at(uhs.front()).claz])
                    .dot(normalize(vps[unaryVars.at(uhs.back()).claz])));
            }
            assert(mg.data(bh).type == MGBinary::RegionRegionConnection);
            if (unaryVars.at(uhs.front()).claz == unaryVars.at(uhs.back()).claz)
                return 1.0;
            
            Vec3 alignDir = normalize(vps[unaryVars.at(uhs.front()).claz].cross(vps[unaryVars.at(uhs.back()).claz]));
            double maxDotProd = 0.0;
            for (int k = 1; k < mg.data(bh).normalizedAnchors.size(); k++){
                Vec3 rotateNorm = mg.data(bh).normalizedAnchors[k-1].cross(mg.data(bh).normalizedAnchors[k]);
                double dotProd = abs(alignDir.dot(rotateNorm));
                if (dotProd > maxDotProd)
                    maxDotProd = dotProd;
            }
            return 1.0 - maxDotProd;
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
            std::unordered_map<MGUnaryHandle, RegionIndex> mgh2ri;
            std::unordered_map<MGUnaryHandle, LineIndex> mgh2li;

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
                    auto center = normalize(cam.spatialDirection(rd.data.center));
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
                    mgh2ri[ri2mgh[ri]] = ri;
                   /* if (ri2mgh[ri].id == 1307){
                        GetVersion();
                    }*/
                    unaryVars[ri2mgh[ri]].claz = initialClaz;
                    unaryVars[ri2mgh[ri]].depthOfCenter = initialDepth;
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
                        normalize(cam.spatialDirection(ld.data.line.component.center()))
                    });
                    mgh2li[li2mgh[li]] = li;
                    unaryVars[li2mgh[li]].claz = ld.data.line.claz;
                    unaryVars[li2mgh[li]].depthOfCenter = initialDepth;
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
                    auto & imp = mg.data(bh).importanceRatioInRelatedUnaries;
                    importanceRatioSum += mg.topo(bh).lowers.front() == u.topo.hd ? imp.front() : imp.back();
                }
                assert(FuzzyEquals(importanceRatioSum, 1.0, 0.01));
            }
#endif

            // make overlapped region claz consistent
            static const bool makeOverlappedRegionClassConsistent = false;
            if (makeOverlappedRegionClassConsistent){
                std::vector<MGUnaryHandle> uhs;
                for (auto & v : mg.elements<0>()){
                    if (v.data.type == MGUnary::Region){
                        uhs.push_back(v.topo.hd);
                    }
                }
                std::map<int, std::vector<MGUnaryHandle>> ccUhs;
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
                }, [&ccUhs](const MGUnaryHandle & uh, int ccid){
                    ccUhs[ccid].push_back(uh);
                });
                for (auto & ccUh : ccUhs){
                    auto & uhs = ccUh.second;
                    Vec3 center(0, 0, 0);
                    for (auto uh : uhs){
                        if (uh.id == 1307){
                            GetVersion();
                        }
                        center += (mg.data(uh).normalizedCenter * GetData(mgh2ri[uh], regionsGraphs).area);
                    }
                    center /= norm(center);
                    int initialClaz = std::distance(vps.begin(), std::min_element(vps.begin(), vps.end(),
                        [&center](const Vec3 & vp1, const Vec3 & vp2) -> bool {
                        return AngleBetweenUndirectedVectors(center, vp1) <
                            AngleBetweenUndirectedVectors(center, vp2);
                    }));
                    for (auto uh : uhs){
                        unaryVars[uh].claz = initialClaz;
                    }
                }
            }

            binaryVars.clear();
            for (auto & b : mg.elements<1>()){
                binaryVars[b.topo.hd].sampleDepthsOnRelatedUnaries.front().resize(b.data.normalizedAnchors.size());
                binaryVars[b.topo.hd].sampleDepthsOnRelatedUnaries.back().resize(b.data.normalizedAnchors.size());
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

        bool IsGoodPatch(const MixedGraph & mg, const MGPatch & patch, const std::vector<Vec3> & vps){
            if (!BinaryHandlesAreValidInPatch(mg, patch) || !UnariesAreConnectedInPatch(mg, patch))
                return false;
            for (auto & bhv : patch.bhs){
                if (!IsGoodBinary(mg, bhv.first, patch.uhs, vps))
                    return false;
            }
            return true;
        }


        double BinaryDistanceOfPatch(const MGBinaryHandle & bh, const MGPatch & patch){
            auto & sampleDepths = patch.bhs.at(bh).sampleDepthsOnRelatedUnaries;
            double distanceSum = 0.0;
            for (int i = 0; i < sampleDepths[0].size(); i++){
                distanceSum += abs(sampleDepths[0][i] - sampleDepths[1][i]);
            }
            return distanceSum / sampleDepths[0].size();
        }

        double AverageBinaryDistanceOfPatch(const MGPatch & patch, int pow){
            double distanceSum = 0.0;
            int distanceCount = 0;
            for (auto & bhv : patch.bhs){
                for (int i = 0; i < bhv.second.sampleDepthsOnRelatedUnaries.front().size(); i++){
                    distanceSum += std::pow(abs(bhv.second.sampleDepthsOnRelatedUnaries.front()[i] -
                        bhv.second.sampleDepthsOnRelatedUnaries.back()[i]), pow);
                }
                distanceCount += bhv.second.sampleDepthsOnRelatedUnaries.front().size();
            }
            return distanceSum / distanceCount;
        }

        double AverageDepthOfPatch(const MGPatch & patch){
            double depthSum = 0.0;
            for (auto & uhv : patch.uhs){
                depthSum += uhv.second.depthOfCenter;
            }
            return depthSum / patch.uhs.size();
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

            for (auto & p : patches){
                assert(UnariesAreConnectedInPatch(mg, p));
                assert(BinaryHandlesAreValidInPatch(mg, p));
            }

            return patches;
        }


        MGPatch MakePatchOnBinary(const MixedGraph & mg, const MGBinaryHandle & bh,
            const MGUnaryVarTable & unaryVars, const MGBinaryVarTable & binaryVars){
            MGPatch patch;
            patch.bhs[bh] = binaryVars.at(bh);
            auto uhs = mg.topo(bh).lowers;
            patch.uhs[uhs[0]] = unaryVars.at(uhs[0]);
            patch.uhs[uhs[1]] = unaryVars.at(uhs[1]);
            assert(BinaryHandlesAreValidInPatch(mg, patch));
            assert(UnariesAreConnectedInPatch(mg, patch));
            return patch;
        }


        MGPatch MakeStarPatchAroundUnary(const MixedGraph & mg, const MGUnaryHandle & uh,
            const MGUnaryVarTable & unaryVars, const MGBinaryVarTable & binaryVars){
            MGPatch patch;
            patch.uhs[uh] = unaryVars.at(uh);
            for (auto & bh : mg.topo(uh).uppers){
                auto anotherUh = mg.topo(bh).lowers.front();
                if (anotherUh == uh)
                    anotherUh = mg.topo(bh).lowers.back();
                patch.bhs[bh] = binaryVars.at(bh);
                patch.uhs[anotherUh] = unaryVars.at(uh);
            }
            assert(BinaryHandlesAreValidInPatch(mg, patch));
            assert(UnariesAreConnectedInPatch(mg, patch));
            return patch;
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


        void ScalePatch(MGPatch & patch, double scale){
            for (auto & uhv : patch.uhs)
                uhv.second.depthOfCenter *= scale;
            for (auto & bhv : patch.bhs){
                for (auto & arr : bhv.second.sampleDepthsOnRelatedUnaries){
                    for (double & d : arr){
                        d *= scale;
                    }
                }
            }
        }


        

      

        void CommitPatchToVariableTable(const MGPatch & patch,
            MGUnaryVarTable & unaryVars, MGBinaryVarTable & binaryVars){
            for (auto & uhv : patch.uhs){
                unaryVars[uhv.first] = uhv.second;
            }
            for (auto & bhv : patch.bhs){
                binaryVars[bhv.first] = bhv.second;
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



        static const int g_NumOfAnchorsUsedForEachBinary = 1;

        struct MGPatchDepthsOptimizerInternalBase {
            virtual void initialize(const MixedGraph & mg, MGPatch & patch, 
                const std::vector<Vec3> & vanishingPoints, bool useWeights) = 0;
            virtual void setDepthBounds(MGPatch & patch, double depthLb, double depthUb) = 0;
            virtual void setDepthsAllGreaterThan(MGPatch & patch, double lob) = 0;
            virtual void setUnaryClass(const MixedGraph & mg, MGPatch & patch,
                const std::vector<Vec3> & vanishingPoints,
                const MGUnaryHandle & uh, int claz) = 0;
            virtual bool optimize(const MixedGraph & mg, MGPatch & patch,
                const std::vector<Vec3> & vanishingPoints) = 0;
            virtual void finalize() = 0;
        };

    

        static struct MosekGlobal {
            MSKenv_t env;
            inline MosekGlobal() { MSK_makeenv(&env, nullptr); }
            inline ~MosekGlobal() { MSK_deleteenv(&env); }
        } mosekGlobal;

        static const bool g_DepthBoundAsConstraint = false;
        static const bool g_BoundSlackPositive = false;

        static void MSKAPI printstr(void *handle,
            MSKCONST char str[]){
            printf("%s", str);
        }


        struct MGPatchDepthsOptimizerInternalMosek : MGPatchDepthsOptimizerInternalBase {
            MSKtask_t task;

            std::unordered_map<MGUnaryHandle, int> uhPositions;
            std::unordered_map<MGBinaryHandle, std::pair<std::vector<Vec3>, int>> bhAnchors;

            int varNum;
            int consNum;

            virtual void initialize(const MixedGraph & mg, MGPatch & patch,
                const std::vector<Vec3> & vanishingPoints, bool useWeights) override;
            virtual void setDepthBounds(MGPatch & patch, double depthLb, double depthUb) override;
            virtual void setDepthsAllGreaterThan(MGPatch & patch, double lob) override;
            virtual void setUnaryClass(const MixedGraph & mg, MGPatch & patch,
                const std::vector<Vec3> & vanishingPoints,
                const MGUnaryHandle & uh, int claz) override;
            virtual bool optimize(const MixedGraph & mg, MGPatch & patch,
                const std::vector<Vec3> & vanishingPoints) override;
            virtual void finalize() override;
        };

        void MGPatchDepthsOptimizerInternalMosek::initialize(const MixedGraph & mg, MGPatch & patch,
            const std::vector<Vec3> & vanishingPoints, bool useWeights){
            auto internalData = this;

            auto & task = internalData->task;

            //auto tick = Tick("preparing optimization");

            assert(BinaryHandlesAreValidInPatch(mg, patch));
            int depthNum = patch.uhs.size();

            assert(UnariesAreConnectedInPatch(mg, patch));

            // temp data  
            auto & uhPositions = internalData->uhPositions;
            auto & bhAnchors = internalData->bhAnchors;

            uhPositions.reserve(patch.uhs.size());
            int unaryNum = 0;
            for (auto & uhv : patch.uhs){
                uhPositions[uhv.first] = unaryNum;
                unaryNum++;
            }

            bhAnchors.reserve(patch.bhs.size());

            // count sample connections
            int samplesNum = 0;
            for (auto & bhv : patch.bhs){
                bhAnchors[bhv.first].first = SelectSubset(mg.data(bhv.first).normalizedAnchors, g_NumOfAnchorsUsedForEachBinary);
                bhAnchors[bhv.first].second = samplesNum;
                samplesNum += bhAnchors[bhv.first].first.size();
            }

            int slackVarNum = samplesNum;
            auto & varNum = internalData->varNum;
            varNum = depthNum + slackVarNum; // depths, slackVars
            auto & consNum = internalData->consNum;
            consNum = samplesNum * 2;   // depth1 * ratio1 - depth2 * ratio2 < slackVar;  
            // depth2 * ratio2 - depth1 * ratio1 < slackVar

            //make as bounds:   depth = defaultValue
            //                  depth \in [depthLb, depthUb]
            if (g_DepthBoundAsConstraint){
                consNum += depthNum;
            }

            //env = nullptr;
            task = nullptr;
            auto & env = mosekGlobal.env;

            MSK_maketask(env, consNum, varNum, &task);
            //MSK_linkfunctotaskstream(task, MSK_STREAM_LOG, NULL, printstr);

            MSK_appendcons(task, consNum);
            MSK_appendvars(task, varNum);

            // set weights
            int slackVarId = 0;
            for (auto & bhv : patch.bhs){
                auto & bd = mg.data(bhv.first);
                for (int k = 0; k < bhAnchors[bhv.first].first.size(); k++){
                    MSK_putcj(task, slackVarId, useWeights ? bd.weight : 1.0);
                    slackVarId++;
                }
            }

            // bounds for vars
            {
                int varId = 0;
                if (!g_DepthBoundAsConstraint){
                    for (; varId < depthNum; varId++){ // for depths
                        MSK_putvarbound(task, varId, MSK_BK_LO, 1.0, +MSK_INFINITY);
                    }
                }
                if (g_BoundSlackPositive){
                    for (; varId < depthNum + slackVarNum; varId++){ // for slack vars 
                        MSK_putvarbound(task, varId, MSK_BK_LO, 0.0, +MSK_INFINITY);
                    }
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
                        relatedAnchorsNum += bhAnchors[bh].first.size();
                    }

                    int relatedConsNum = relatedAnchorsNum * 2;
                    if (g_DepthBoundAsConstraint){
                        relatedConsNum += 1;
                    }

                    std::vector<MSKint32t> consIds;
                    std::vector<MSKrealt> consValues;

                    consIds.reserve(relatedConsNum);
                    consValues.reserve(relatedConsNum);

                    for (auto & bh : relatedBhs){
                        if (!Contains(patch.bhs, bh))
                            continue;

                        int firstAnchorPosition = bhAnchors[bh].second;
                        auto & samples = bhAnchors[bh].first;
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

                    if (g_DepthBoundAsConstraint){ // add depth bound
                        consIds.push_back(samplesNum * 2 + varId);
                        consValues.push_back(1.0); // (1.0) * depth > 1.0;
                    }

                    assert(varId == uhPositions[uh]);
                    MSK_putacol(task, uhPositions[uh], relatedConsNum, consIds.data(), consValues.data());
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
            int consId = 0;
            for (; consId < samplesNum * 2; consId++){
                // all [depth1 * ratio1 - depth2 * ratio2 - slackVar < 0];
                MSK_putconbound(task, consId, MSK_BK_UP, -MSK_INFINITY, 0.0);
            }
            if (g_DepthBoundAsConstraint){
                for (; consId < consNum; consId++){
                    // depth > 1.0
                    MSK_putconbound(task, consId, MSK_BK_LO, 1.0, +MSK_INFINITY);
                }
            }
            else{
                assert(consId == consNum);
            }


        }

        void MGPatchDepthsOptimizerInternalMosek::setDepthBounds(MGPatch & patch, double depthLb, double depthUb) {
            for (int varId = 0; varId < patch.uhs.size(); varId++){ // for depths
                MSK_putvarbound(task,
                    varId, MSK_BK_RA, depthLb, depthUb);
            }
        }

        void MGPatchDepthsOptimizerInternalMosek::setDepthsAllGreaterThan(MGPatch & patch, double lob) {
            for (int varId = 0; varId < patch.uhs.size(); varId++){ // for depths
                MSK_putvarbound(task,
                    varId, MSK_BK_LO, lob, +MSK_INFINITY);
            }
        }

        void MGPatchDepthsOptimizerInternalMosek::setUnaryClass(const MixedGraph & mg, MGPatch & patch,
            const std::vector<Vec3> & vanishingPoints,
            const MGUnaryHandle & uh, int claz) {

            auto & relatedBhs = mg.topo(uh).uppers;
            int relatedAnchorsNum = 0;
            for (auto & bh : relatedBhs){
                if (!Contains(patch.bhs, bh))
                    continue;
                relatedAnchorsNum += bhAnchors[bh].first.size();
            }
            int relatedConsNum = relatedAnchorsNum * 2;

            std::vector<MSKint32t> consIds;
            std::vector<MSKrealt> consValues;

            consIds.reserve(relatedConsNum);
            consValues.reserve(relatedConsNum);

            for (auto & bh : relatedBhs){
                if (!Contains(patch.bhs, bh))
                    continue;
                int firstAnchorPosition = bhAnchors.at(bh).second;
                auto & samples = bhAnchors.at(bh).first;
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

            MSK_putacol(task, uhPositions.at(uh), relatedConsNum, consIds.data(), consValues.data());
        }

        bool MGPatchDepthsOptimizerInternalMosek::optimize(const MixedGraph & mg, MGPatch & patch,
            const std::vector<Vec3> & vanishingPoints){

            MSK_putobjsense(task, MSK_OBJECTIVE_SENSE_MINIMIZE);
            MSKrescodee trmcode;
            MSK_optimizetrm(task, &trmcode);
            MSK_solutionsummary(task, MSK_STREAM_LOG);

            MSKsolstae solsta;
            MSK_getsolsta(task, MSK_SOL_BAS, &solsta);

            //Tock(tick);

            switch (solsta){
            case MSK_SOL_STA_OPTIMAL:
            case MSK_SOL_STA_NEAR_OPTIMAL:
            {
                                             double *xx = new double[varNum];
                                             MSK_getxx(task, MSK_SOL_BAS, xx);
                                             /*std::cout << "results: ";
                                             for (int i = 0; i < varNum; i++)
                                             std::cout << xx[i] << ' ';
                                             std::cout << std::endl;*/
                                             //printf("Optimal primal solution\n");
                                             for (auto & uhv : patch.uhs){ // get resulted depths
                                                 uhv.second.depthOfCenter = xx[uhPositions.at(uhv.first)];
                                             }
                                             delete[] xx;
                                             UpdateBinaryVars(mg, vanishingPoints, patch.uhs, patch.bhs);
                                             return true;
            }
            case MSK_SOL_STA_DUAL_INFEAS_CER:
            case MSK_SOL_STA_PRIM_INFEAS_CER:
            case MSK_SOL_STA_NEAR_DUAL_INFEAS_CER:
            case MSK_SOL_STA_NEAR_PRIM_INFEAS_CER:
                printf("Primal or dual infeasibility certificate found.\n");
                return false;
            case MSK_SOL_STA_UNKNOWN:
            {
                                        char symname[MSK_MAX_STR_LEN];
                                        char desc[MSK_MAX_STR_LEN];
                                        MSK_getcodedesc(trmcode,
                                            symname,
                                            desc);

                                        printf("The solution status is unknown.\n");
                                        printf("The optimizer terminitated with code: %s\n", symname);
                                        return false;
            }
            default:
                printf("Other solution status.\n");
                return false;
            }
        }

        void MGPatchDepthsOptimizerInternalMosek::finalize(){
            MSK_deletetask(&task);
        }




        struct MGPatchDepthsOptimizerInternalEigen : MGPatchDepthsOptimizerInternalBase {
            Eigen::SparseMatrix<double> A, W;
            Eigen::VectorXd B;
            bool useWeights;
            std::unordered_map<MGUnaryHandle, int> uhPositions;
            std::unordered_map<MGBinaryHandle, std::vector<Vec3>> bhAnchors;

            virtual void initialize(const MixedGraph & mg, MGPatch & patch,
                const std::vector<Vec3> & vanishingPoints, bool useWeights) override;
            virtual void setDepthBounds(MGPatch & patch, double depthLb, double depthUb) override;
            virtual void setDepthsAllGreaterThan(MGPatch & patch, double lob) override;
            virtual void setUnaryClass(const MixedGraph & mg, MGPatch & patch,
                const std::vector<Vec3> & vanishingPoints,
                const MGUnaryHandle & uh, int claz) override;
            virtual bool optimize(const MixedGraph & mg, MGPatch & patch,
                const std::vector<Vec3> & vanishingPoints) override;
            virtual void finalize() override;
        };

        void MGPatchDepthsOptimizerInternalEigen::initialize(const MixedGraph & mg, MGPatch & patch,
            const std::vector<Vec3> & vanishingPoints, bool useWeights) {

            this->useWeights = useWeights;

            assert(BinaryHandlesAreValidInPatch(mg, patch));
            assert(UnariesAreConnectedInPatch(mg, patch));

            int uhid = 0;
            for (auto & uhv : patch.uhs){
                uhPositions[uhv.first] = uhid++;
            }

            int varNum = patch.uhs.size();
            int consNum = 0;
            
            consNum++;
            for (auto & bhv : patch.bhs){
                bhAnchors[bhv.first] = SelectSubset(mg.data(bhv.first).normalizedAnchors, g_NumOfAnchorsUsedForEachBinary);
                consNum += bhAnchors.at(bhv.first).size();
            }
            
            A.resize(consNum, varNum);
            W.resize(consNum, consNum);
            B.resize(consNum);
            
            // write equations
            int eid = 0;
            A.insert(eid, 0) = 1.0;
            B(eid) = 1.0;
            W.insert(eid, eid) = 1.0;
            eid++;

            for (auto & bhv : patch.bhs){
                auto & bh = bhv.first;
                for (auto & a : bhAnchors.at(bh)){
                    auto uh1 = mg.topo(bh).lowers.front();
                    auto uh2 = mg.topo(bh).lowers.back();
                    auto & u1 = mg.data(uh1);
                    auto & u2 = mg.data(uh2);
                    double ratio1 = DepthRatioOnMGUnary(a, u1, vanishingPoints, patch.uhs.at(uh1).claz);
                    double ratio2 = DepthRatioOnMGUnary(a, u2, vanishingPoints, patch.uhs.at(uh2).claz);
                    A.insert(eid, uhPositions.at(uh1)) = ratio1;
                    A.insert(eid, uhPositions.at(uh2)) = -ratio2;
                    B(eid) = 0.0;
                    W.insert(eid, eid) = mg.data(bh).weight;
                    eid++;
                }
            }

            assert(eid == consNum);
        }

        void MGPatchDepthsOptimizerInternalEigen::setDepthBounds(MGPatch & patch, double depthLb, double depthUb) {
            NOT_IMPLEMENTED_YET();
        }

        void MGPatchDepthsOptimizerInternalEigen::setDepthsAllGreaterThan(MGPatch & patch, double lob) {
            NOT_IMPLEMENTED_YET();
        }

        void MGPatchDepthsOptimizerInternalEigen::setUnaryClass(const MixedGraph & mg, MGPatch & patch,
            const std::vector<Vec3> & vanishingPoints,
            const MGUnaryHandle & uh, int claz) {
            int eid = 1;
            for (auto & bhv : patch.bhs){
                auto & bh = bhv.first;
                auto uh1 = mg.topo(bh).lowers.front();
                auto uh2 = mg.topo(bh).lowers.back();
                if (uh1 == uh || uh2 == uh){
                    for (auto & a : bhAnchors.at(bh)){
                        auto & u1 = mg.data(uh1);
                        auto & u2 = mg.data(uh2);
                        double ratio1 = DepthRatioOnMGUnary(a, u1, vanishingPoints, claz);
                        double ratio2 = DepthRatioOnMGUnary(a, u2, vanishingPoints, claz);
                        A.insert(eid, uhPositions.at(uh1)) = ratio1;
                        A.insert(eid, uhPositions.at(uh2)) = -ratio2;
                        B(eid) = 0.0;
                        W.insert(eid, eid) = mg.data(bh).weight;
                    }
                }
                eid += bhAnchors.at(bh).size();
            }
        }

        bool MGPatchDepthsOptimizerInternalEigen::optimize(const MixedGraph & mg, MGPatch & patch,
            const std::vector<Vec3> & vanishingPoints) {
            
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
                uhv.second.depthOfCenter = X(uhPositions.at(uhv.first));
            }

            UpdateBinaryVars(mg, vanishingPoints, patch.uhs, patch.bhs);

            return true;

        }

        void MGPatchDepthsOptimizerInternalEigen::finalize() {}


        



        MGPatchDepthsOptimizer::MGPatchDepthsOptimizer(const MixedGraph & mg, MGPatch & patch,
            const std::vector<Vec3> & vanishingPoints, bool useWeights, AlgorithmType at)
            : _mg(mg), _patch(patch), _vanishingPoints(vanishingPoints), _at(at){

            if (_at == AlgorithmType::MosekLinearProgramming){
                _internal = new MGPatchDepthsOptimizerInternalMosek;
            }
            else if (_at == AlgorithmType::EigenSparseQR){
                _internal = new MGPatchDepthsOptimizerInternalEigen;
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

        void MGPatchDepthsOptimizer::setDepthBounds(double depthLb, double depthUb){
            static_cast<MGPatchDepthsOptimizerInternalBase*>(_internal)->setDepthBounds(_patch, depthLb, depthUb);
        }

        void MGPatchDepthsOptimizer::setDepthsAllGreaterThan(double lob){
            static_cast<MGPatchDepthsOptimizerInternalBase*>(_internal)->setDepthsAllGreaterThan(_patch, lob);
        }

        void MGPatchDepthsOptimizer::setUnaryClass(const MGUnaryHandle & uh, int claz){
            static_cast<MGPatchDepthsOptimizerInternalBase*>(_internal)->setUnaryClass(_mg, _patch, _vanishingPoints, uh, claz);
        }

        bool MGPatchDepthsOptimizer::optimize() {
            return static_cast<MGPatchDepthsOptimizerInternalBase*>(_internal)->optimize(_mg, _patch, _vanishingPoints);
        }

  

        void GrowPatch(const MixedGraph & mg, MGPatch & root, 
            const MGUnaryVarTable & unaryVars, const MGBinaryVarTable & binaryVars){

            MaxHeap<MGUnaryHandle> Q;
            for (auto & uhv : root.uhs){
                Q.push(uhv.first, 0.0);
            }
            for (auto & bhv : root.bhs){
                auto & uhs = mg.topo(bhv.first).lowers;
                Q.setScore(uhs[0], Q.at(uhs[0]) + mg.data(bhv.first).importanceRatioInRelatedUnaries[0]);
                Q.setScore(uhs[1], Q.at(uhs[1]) + mg.data(bhv.first).importanceRatioInRelatedUnaries[1]);
            }

            while (!Q.empty()){
                // todo
            }

        }

    }
}