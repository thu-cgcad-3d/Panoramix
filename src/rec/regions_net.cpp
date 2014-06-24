#include "regions_net.hpp"

#include "../vis/visualize2d.hpp"

namespace panoramix {
    namespace rec {

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
                std::map<std::pair<int, int>, std::vector<std::vector<PixelLoc>>> & boundaryEdges,
                std::map<std::tuple<int, int, int>, std::set<PixelLoc, ComparePixelLoc>> & triPixels) {

                std::map<std::pair<int, int>, std::set<PixelLoc, ComparePixelLoc>> boundaryPixels;

                int width = segRegions.cols;
                int height = segRegions.rows;
                for (int y = 0; y < height - 1; y++) {
                    for (int x = 0; x < width - 1; x++) {
                        PixelLoc p1(x, y), p2(x + 1, y), p3(x, y + 1), p4(x + 1, y + 1);
                        int rs[] = { segRegions.at<int32_t>(p1),
                            segRegions.at<int32_t>(p2),
                            segRegions.at<int32_t>(p3),
                            segRegions.at<int32_t>(p4) };

                        if (rs[0] != rs[1]) {
                            boundaryPixels[MakeOrderedPair(rs[0], rs[1])].insert(p1);
                        }
                        if (rs[0] != rs[2]) {
                            boundaryPixels[MakeOrderedPair(rs[0], rs[2])].insert(p1);
                        }
                        if (rs[0] != rs[3]) {
                            boundaryPixels[MakeOrderedPair(rs[0], rs[3])].insert(p1);
                        }

                        std::sort(std::begin(rs), std::end(rs));
                        auto last = std::unique(std::begin(rs), std::end(rs));
                        if (std::distance(std::begin(rs), last) >= 3) { // a tri-pixel
                            auto & triPixel = triPixels[std::make_tuple(rs[0], rs[1], rs[2])];
                            bool hasClosePixels = false;
                            for (auto & p : triPixel) {
                                if (Distance(p, p1) < 2) {
                                    hasClosePixels = true;
                                    break;
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
                    
                    static const int xdirs[] = { 1, 0, -1, 0, -1, 1, 1, -1 };
                    static const int ydirs[] = { 0, 1, 0, -1, 1, -1, 1, -1 };

                    std::vector<std::vector<PixelLoc>> edges;
                    edges.push_back({ p });
                    pixels.erase(p);

                    while (true) {
                        auto & curEdge = edges.back();
                        auto & curTail = curEdge.back();

                        bool foundMore = false;
                        for (int i = 0; i < 8; i++) {
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

            inline PixelLoc ToPixelLoc(const Point2 & p) {
                return PixelLoc(static_cast<int>(p[0]), static_cast<int>(p[1]));
            }

            double ComputeSpanningArea(const Point2 & a, const Point2 & b, const InfiniteLine2 & line) {
                auto ad = SignedDistanceFromPointToLine(a, line);
                auto bd = SignedDistanceFromPointToLine(b, line);
                auto ap = DistanceFromPointToLine(a, line).second;
                auto bp = DistanceFromPointToLine(b, line).second;
                if (ad * bd >= 0) {
                    return norm(ap - bp) * std::abs(ad + bd) / 2.0;
                }
                ad = abs(ad);
                bd = abs(bd);
                auto len = norm(ap - bp);
                return (ad * ad + bd * bd) * len / (ad + bd) / 2.0;
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
                std::vector<std::vector<PixelLoc>> contours;
                cv::findContours(regionMaskCopy, contours, CV_RETR_EXTERNAL, CV_CHAIN_APPROX_SIMPLE);
                vd.contour = contours.front();
                assert(!vd.contour.empty() && "no contour? impossible~");
               
                auto vh = _regions.add(vd);
                //_regionsRTree.insert(vh);  
            }

            std::map<std::pair<int, int>, std::vector<std::vector<PixelLoc>>> boundaryEdges;
            std::map<std::tuple<int, int, int>, std::set<PixelLoc, ComparePixelLoc>> triPixels;
            FindContoursOfRegionsAndBoundaries(_segmentedRegions, regionNum, boundaryEdges, triPixels);

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
                for (auto & e : edges) {
                    for (int i = 0; i < e.size() - 1; i++) {
                        interArea += ComputeSpanningArea(
                            Point2(e[i].x, e[i].y),
                            Point2(e[i + 1].x, e[i + 1].y),
                            bd.fittedLine);
                    }
                }

                bd.straightness = Gaussian(interArea / bd.length, 1.0);
                if (edges.size() == 1 && edges.front().size() == 2) {
                    assert(FuzzyEquals(bd.straightness, 1.0, 0.01) && "simple line should has the best straightness..");
                }

                // collect sampled points
                static const double stepLen = 8;
                bd.sampledPoints.resize(bd.edges.size());
                for (int i = 0; i < bd.edges.size(); i++) {
                    auto & e = bd.edges[i];
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

            for (auto & boundary : _regions.elements<1>()) {
                // compute tjunction likelihood
                auto & bd = boundary.data;
                auto rh1 = boundary.topo.lowers[0];
                auto rh2 = boundary.topo.lowers[1];
                
                for (auto & e : bd.edges) {
                    if (e.size() <= 1)
                        continue;
                    auto head = e.front();
                    auto tail = e.back();

                }
            }


            // visualize region contours
            IF_DEBUG_USING_VISUALIZERS {
                Image regionVis(_image.rows, _image.cols, CV_8UC3, vis::Color(100, 100, 100));
                std::vector<std::vector<std::vector<PixelLoc>>> boundaries;
                std::vector<std::vector<std::vector<Point2>>> sampledPoints;
                std::vector<double> straightnesses;
                boundaries.reserve(_regions.internalElements<1>().size());
                sampledPoints.reserve(_regions.internalElements<1>().size());
                straightnesses.reserve(_regions.internalElements<1>().size());
                for (auto & b : _regions.elements<1>()) {
                    boundaries.push_back(b.data.edges);
                    sampledPoints.push_back(b.data.sampledPoints);
                    straightnesses.push_back(b.data.straightness);
                }
                auto & ctable = vis::PredefinedColorTable(vis::ColorTableDescriptor::AllColors);
                for (auto & r : _regions.elements<0>()) {
                    auto color = vis::Color(0, 0, 0, 1);
                    auto & polyline = r.data.contour;
                    for (int j = 0; j < polyline.size() - 1; j++) {
                        cv::line(regionVis, polyline[j], polyline[j + 1], color, 1);
                    }
                }
                for (int i = 0; i < boundaries.size(); i++) {
                    auto color = vis::Color(255, 255, 255, 1)*straightnesses[i];
                    for (auto & polyline : boundaries[i]) {
                        for (int j = 0; j < polyline.size() - 1; j++) {
                            cv::line(regionVis, polyline[j], polyline[j + 1], color);
                        }
                    }
                    auto pcolor = ctable[i % ctable.size()];
                    for (auto & s : sampledPoints[i]) {
                        for (auto & p : s) {
                            cv::circle(regionVis, ToPixelLoc(p), 1, pcolor);
                        }
                    }
                }
                vis::Visualizer2D(regionVis) << vis::manip2d::Show();
            }
        }

        void RegionsNet::computeImageFeatures() {
            //NOT_IMPLEMENTED_YET();
        }

    }
}