#include "regions_net.hpp"

#include "../vis/visualize2d.hpp"

#include "../core/debug.hpp"

namespace panoramix {
    namespace rec {

        RegionsNet::Params::Params() : samplingStepLengthOnBoundary(15.0) {
            SegmentationExtractor::Params segmenterParams;
            segmenterParams.c = 80;
            segmenterParams.minSize = 300;
            segmenterParams.sigma = 1.0;
            segmenter = SegmentationExtractor(segmenterParams);
        }

        RegionsNet::RegionsNet(const Image & image, const Params & params) 
            : _image(image), _params(params)/*, _regionsRTree(RegionDataBoundingBoxFunctor(_regions))*/{}

        namespace {

            void ComputeRegionProperties(const Image & regionMask, const Image & image,
                Vec2 & center, double & area, Box2 & boundingBox) {
                assert(regionMask.depth() == CV_8U && regionMask.channels() == 1);
                center = { 0, 0 };
                area = 0;
                boundingBox = Box2();
                for (auto i = regionMask.begin<uint8_t>(); i != regionMask.end<uint8_t>(); ++i){
                    if (*i) { // masked
                        area += 1.0;
                        PixelLoc p = i.pos();
                        boundingBox |= BoundingBox(p);
                        center += Vec2(p.x, p.y);
                    }
                }
                center /= area;
            }

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

            void FindContoursOfRegionsAndBoundaries(const SegmentationExtractor::Feature & segRegions, int regionNum,
                std::map<std::pair<int, int>, std::vector<std::vector<PixelLoc>>> & boundaryEdges) {

                std::map<std::pair<int, int>, std::set<PixelLoc, ComparePixelLoc>> boundaryPixels;

                int width = segRegions.cols;
                int height = segRegions.rows;
                for (int y = 0; y < height - 1; y++) {
                    for (int x = 0; x < width - 1; x++) {
                        PixelLoc p1(x, y), p2(x + 1, y), p3(x, y + 1), p4(x + 1, y + 1);
                        int rs[] = { 
                            segRegions.at<int32_t>(p1),
                            segRegions.at<int32_t>(p2),
                            segRegions.at<int32_t>(p3),
                            segRegions.at<int32_t>(p4) 
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
                    
                    static const int xdirs[] = { 1, 0, -1, 0, -1, 1, 1, -1, 0, 0,  2, -2};
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
                            } else {
                                bool closed = Distance(edges.back().front(), edges.back().back()) <= 1.5;
                                cv::approxPolyDP(edges.back(), edges.back(), 2, closed);
                                if (edges.back().size() <= 1)
                                    edges.pop_back();
                            }

                            if (pixels.empty()) { // no more pixels
                                break;
                            } else { // more pixels
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

        void RegionsNet::buildNetAndComputeGeometricFeatures() {
            _segmentedRegions = _params.segmenter(_image);
            int regionNum = static_cast<int>(MinMaxValOfImage(_segmentedRegions).second) + 1;
            _regions.internalElements<0>().reserve(regionNum);
            for (int i = 0; i < regionNum; i++){
                RegionData vd;
                vd.regionMask = (_segmentedRegions == i);
                ComputeRegionProperties(vd.regionMask, _image, vd.center, vd.area, vd.boundingBox);

                // find contour of the region
                cv::Mat regionMaskCopy;
                vd.regionMask.copyTo(regionMaskCopy);
                cv::findContours(regionMaskCopy, vd.contours, CV_RETR_LIST, CV_CHAIN_APPROX_SIMPLE);
                assert(!vd.contours.empty() && "no contour? impossible~");

                int dilationSize = 3;
                auto dilationElement = cv::getStructuringElement(cv::MORPH_ELLIPSE,
                    cv::Size(2 * dilationSize + 1, 2 * dilationSize + 1),
                    cv::Point(dilationSize, dilationSize));
                vd.regionMask.copyTo(regionMaskCopy);
                cv::dilate(regionMaskCopy, regionMaskCopy, dilationElement);
                cv::findContours(regionMaskCopy, vd.dilatedContours, CV_RETR_LIST, CV_CHAIN_APPROX_SIMPLE);
                assert(!vd.dilatedContours.empty() && "no contour? impossible~");
               
                auto vh = _regions.add(vd);
                //_regionsRTree.insert(vh);  
            }

            std::map<std::pair<int, int>, std::vector<std::vector<PixelLoc>>> boundaryEdges;
            FindContoursOfRegionsAndBoundaries(_segmentedRegions, regionNum, boundaryEdges);

            for (auto & bep : boundaryEdges) {
                auto & rids = bep.first;
                auto & edges = bep.second;

                BoundaryData bd;
                bd.edges = edges;
                bd.length = 0;
                for (auto & e : edges) {
                    assert(!e.empty() && "edges should never be empty!");
                    for (int i = 0; i < e.size()-1; i++) {
                        bd.length += Distance(e[i], e[i + 1]);
                    }
                } 

                // compute straightness
                std::vector<Point<float, 2>> points;
                for (auto & e : edges) {
                    for (auto & p : e) {
                        points.push_back(Point<float, 2>(p.x, p.y));
                    }
                }

                cv::Vec4f line;
                cv::fitLine(points, line, CV_DIST_L2, 0, 0.01, 0.01);
                bd.fittedLine = InfiniteLine2({ line[2], line[3] }, { line[0], line[1] });
                double interArea = 0;
                double interLen = 0;
                for (auto & e : edges) {
                    for (int i = 0; i < e.size() - 1; i++) {
                        double area, len;
                        std::tie(area, len) = ComputeSpanningArea(
                            ToPoint2(e[i]),
                            ToPoint2(e[i + 1]),
                            bd.fittedLine);
                        interArea += area;
                        interLen += len;
                    }
                }

                bd.straightness = Gaussian(interArea / interLen, 1.0);
                if (edges.size() == 1 && edges.front().size() == 2) {
                    assert(FuzzyEquals(bd.straightness, 1.0, 0.01) && "simple line should has the best straightness..");
                }

                // collect sampled points
                double stepLen = _params.samplingStepLengthOnBoundary;
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

                _regions.add<1>({ RegionHandle(rids.first), RegionHandle(rids.second) }, bd);
            }

            // find connection of end points of boundaries
            struct BoundaryEndPoint {
                BoundaryHandle bh;
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
            for (auto & boundary : _regions.elements<1>()) {
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
            for (auto & b : _regions.elements<1>()) {
                b.data.tjunctionLikelihood = 0;
            }

            for (auto & beps : mergedBepsTable) {
                std::set<BoundaryHandle> bhs;
                for (auto & bep : beps.second) {
                    bhs.insert(bep.bh);
                }
                if (bhs.size() == 3 && beps.second.size() == 3) {
                    std::vector<std::vector<Point<float, 2>>> sampledPointsOnBoundaries(3);
                    static const size_t tjunctionSampleNum = 5;
                    for (int i = 0; i < 3; i ++) {
                        auto & bep = beps.second[i];
                        auto & spts = _regions.data(bep.bh).sampledPoints[bep.edgeId];
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
                            _regions.data(beps.second[k].bh).tjunctionLikelihood += 
                                tjunctionLikelihoods[k];
                        }
                    }
                }
            }

           

            // visualize region contours
            IF_DEBUG_USING_VISUALIZERS {
                Image regionVis(_image.rows, _image.cols, CV_8UC3, vis::Color(100, 100, 100));
                //Image regionVis = _params.segmenter(_image, true);
                std::vector<std::vector<std::vector<PixelLoc>>> boundaries;
                std::vector<std::vector<std::vector<Point2>>> sampledPoints;
                std::vector<double> straightnesses, tjunctionlikelihoods;
                boundaries.reserve(_regions.internalElements<1>().size());
                sampledPoints.reserve(_regions.internalElements<1>().size());
                straightnesses.reserve(_regions.internalElements<1>().size());
                tjunctionlikelihoods.reserve(_regions.internalElements<1>().size());
                for (auto & b : _regions.elements<1>()) {
                    boundaries.push_back(b.data.edges);
                    sampledPoints.push_back(b.data.sampledPoints);
                    straightnesses.push_back(b.data.straightness);
                    tjunctionlikelihoods.push_back(b.data.tjunctionLikelihood);
                }
                auto & ctable = vis::PredefinedColorTable(vis::ColorTableDescriptor::AllColors);
                //for (auto & r : _regions.elements<0>()) {
                //    auto color = vis::Color(0, 0, 0, 1);
                //    for (auto & polyline : r.data.dilatedContours) {
                //        for (int j = 0; j < polyline.size() - 1; j++) {
                //            cv::line(regionVis, polyline[j], polyline[j + 1], color, 1);
                //        }
                //    }
                //}
                for (int i = 0; i < boundaries.size(); i++) {
                    auto pcolor = ctable[i % ctable.size()];
                    auto color = vis::Color(255, 255, 255, 1) * straightnesses[i];
                    auto tjcolor = vis::Color(255, 255, 255, 1) * tjunctionlikelihoods[i];
                    for (auto & polyline : boundaries[i]) {
                        for (int j = 0; j < polyline.size() - 1; j++) {
                            cv::line(regionVis, polyline[j], polyline[j + 1], tjcolor);
                        }
                    }
                    for (auto & s : sampledPoints[i]) {
                        for (auto & p : s) {
                            cv::circle(regionVis, ToPixelLoc(p), 1, pcolor);
                        }
                    }
                }
                //for (auto & bep : mergedBepsTable) {
                //    auto color = vis::ColorFromTag(vis::ColorTag::White);
                //    for (auto & p : bep.second) {
                //        //regionVis.at<Vec<uint8_t, 3>>(p.position) = Vec<uint8_t, 3>(255, 255, 255);
                //        cv::circle(regionVis, p.position, 1, color);
                //    }
                //}
                vis::Visualizer2D(regionVis) << vis::manip2d::Show();
            }
        }

        void RegionsNet::computeImageFeatures() {
            //NOT_IMPLEMENTED_YET();
        }

    }
}