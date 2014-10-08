#include <iostream>
#include <unordered_map>

#include <Eigen/Dense>
#include <Eigen/Sparse>

extern "C" {
    #include <lsd.h>
}

#include <SLIC.h>

#include "misc.hpp"
#include "feature.hpp"
#include "utilities.hpp"

#include "debug.hpp"

namespace panoramix {
    namespace core {

        void NonMaximaSuppression(const Image & src, Image & dst, int sz, std::vector<PixelLoc> * pixels,
            const ImageWithType<bool> & mask) {

            const int M = src.rows;
            const int N = src.cols;
            const bool masked = !mask.empty();

            cv::Mat block = 255 * cv::Mat_<uint8_t>::ones(Size(2 * sz + 1, 2 * sz + 1));
            dst = cv::Mat::zeros(src.size(), src.type());

            // iterate over image blocks
            for (int m = 0; m < M; m += sz + 1) {
                for (int n = 0; n < N; n += sz + 1) {
                    cv::Point ijmax;
                    double vcmax, vnmax;

                    // get the maximal candidate within the block
                    cv::Range ic(m, std::min(m + sz + 1, M));
                    cv::Range jc(n, std::min(n + sz + 1, N));
                    cv::minMaxLoc(src(ic, jc), NULL, &vcmax, NULL, &ijmax, masked ? mask(ic, jc) : cv::noArray());
                    cv::Point cc = ijmax + cv::Point(jc.start, ic.start);

                    // search the neighbours centered around the candidate for the true maxima
                    cv::Range in(std::max(cc.y - sz, 0), std::min(cc.y + sz + 1, M));
                    cv::Range jn(std::max(cc.x - sz, 0), std::min(cc.x + sz + 1, N));

                    // mask out the block whose maxima we already know
                    cv::Mat_<uint8_t> blockmask;
                    block(cv::Range(0, in.size()), cv::Range(0, jn.size())).copyTo(blockmask);
                    cv::Range iis(ic.start - in.start, std::min(ic.start - in.start + sz + 1, in.size()));
                    cv::Range jis(jc.start - jn.start, std::min(jc.start - jn.start + sz + 1, jn.size()));
                    blockmask(iis, jis) = cv::Mat_<uint8_t>::zeros(Size(jis.size(), iis.size()));
                    minMaxLoc(src(in, jn), NULL, &vnmax, NULL, &ijmax, masked ? mask(in, jn).mul(blockmask) : blockmask);
                    cv::Point cn = ijmax + cv::Point(jn.start, in.start);

                    // if the block centre is also the neighbour centre, then it's a local maxima
                    if (vcmax > vnmax) {
                        std::memcpy(dst.ptr(cc.y, cc.x), src.ptr(cc.y, cc.x), src.elemSize());
                        if (pixels){
                            pixels->push_back(cc);
                        }
                    }
                }
            }
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


        namespace {

            void ExtractLines(const cv::Mat& im, std::vector<Line2> & lines,
                int minlen, int xborderw, int yborderw, int numDir) {

                std::cout << "image processing..." << std::endl;

                cv::Mat gim;
                cv::cvtColor(im, gim, CV_BGR2GRAY);
                int h = gim.rows;
                int w = gim.cols;

                cv::Mat dx, dy;
                cv::Mat ggim;
                cv::GaussianBlur(gim, ggim, cv::Size(7, 7), 1.5);
                cv::Sobel(ggim, dx, CV_64F, 1, 0);
                cv::Sobel(ggim, dy, CV_64F, 0, 1);

                cv::Mat imCanny;
                cv::Canny(gim, imCanny, 5, 20);

                std::cout << "gradient binning..." << std::endl;

                cv::Mat imBinIds(im.size(), CV_32SC1); // int32_t
                std::vector<std::set<int>> pixelIdsSet(numDir);

                for (int x = 0; x < imCanny.cols; x++){
                    for (int y = 0; y < imCanny.rows; y++){
                        if (imCanny.at<uchar>(y, x) > 0){
                            double a = atan(dy.at<double>(y, x) / dx.at<double>(y, x));
                            if (isnan(a)){ // NaN
                                continue;
                            }
                            // compute bin id
                            int binId = int((a / M_PI + 0.5) * numDir);
                            if (binId == -1) binId = 0;
                            if (binId == numDir) binId = numDir - 1;

                            int pixelId = y + imCanny.rows * x;
                            
                            imBinIds.at<int32_t>(y, x) = binId;

                            pixelIdsSet[(binId + numDir - 1) % numDir].insert(pixelId);
                            pixelIdsSet[binId].insert(pixelId);
                            pixelIdsSet[(binId + 1) % numDir].insert(pixelId);
                        }
                        else{
                            imBinIds.at<int32_t>(y, x) = -1;
                        }
                    }
                }

                std::vector<int> xs, ys, ids;
                xs.reserve(512);
                ys.reserve(512);
                ids.reserve(512);

                std::cout << "collecting pixels.." << std::endl;

                for (int binId = 0; binId < numDir; binId++){
                    // search in bins
                    auto pixelIdNotSearchedYet = pixelIdsSet[binId];

                    while (true){
                        if (pixelIdNotSearchedYet.empty())
                            break;
                        int rootId = *pixelIdNotSearchedYet.begin();

                        // BFS
                        xs.clear(); ys.clear(); ids.clear();
                        int x = rootId / imCanny.rows;
                        int y = rootId - imCanny.rows * x;
                        xs.push_back(x);
                        ys.push_back(y);
                        ids.push_back(rootId);

                        // used for search only
                        pixelIdNotSearchedYet.erase(rootId);
                        int head = 0;

                        static const int xdirs[] = { 1, 1, 0, -1, -1, -1, 0, 1 };
                        static const int ydirs[] = { 0, 1, 1, 1, 0, -1, -1, -1 };
                        while (true){
                            if (head == xs.size())
                                break;
                            x = xs[head];
                            y = ys[head];
                            for (int k = 0; k < 8; k++){
                                int nx = x + xdirs[k];
                                int ny = y + ydirs[k];
                                int npixelId = ny + imCanny.rows * nx;
                                if (pixelIdNotSearchedYet.find(npixelId) != pixelIdNotSearchedYet.end()){
                                    xs.push_back(nx);
                                    ys.push_back(ny);
                                    ids.push_back(npixelId);
                                    pixelIdNotSearchedYet.erase(npixelId);
                                }
                            }
                            head++;
                        }

                        int edgeSize = (int)xs.size();
                        if (edgeSize < minlen)
                            continue;

                        cv::Mat xsmat(xs), ysmat(ys);
                        double meanx = cv::mean(xsmat).val[0], meany = cv::mean(ysmat).val[0];
                        cv::Mat zmx = xsmat - meanx, zmy = ysmat - meany;

                        cv::Mat v, lambda;
                        Mat<double, 2, 2> D;
                        D(0, 0) = cv::sum(zmx.mul(zmx)).val[0];
                        D(0, 1) = D(1, 0) = cv::sum(zmx.mul(zmy)).val[0];
                        D(1, 1) = cv::sum(zmy.mul(zmy)).val[0];
                        cv::eigen(D, true, lambda, v);

                        double theta = atan2(v.at<double>(0, 1), v.at<double>(0, 0));
                        double confidence = std::numeric_limits<double>::max();
                        if (lambda.at<double>(1) > 0){
                            confidence = lambda.at<double>(0) / lambda.at<double>(1);
                        }

                        // build line
                        if (confidence >= 200){ ////
                            for (int pid : ids){
                                pixelIdsSet[binId].erase(pid);
                                pixelIdsSet[(binId - 1 + numDir) % numDir].erase(pid);
                                pixelIdsSet[(binId + 1) % numDir].erase(pid);
                            }

                            auto xends = std::minmax_element(xs.begin(), xs.end());
                            auto yends = std::minmax_element(ys.begin(), ys.end());
                            double minx = *xends.first, maxx = *xends.second;
                            double miny = *yends.first, maxy = *yends.second;

                            if (maxx <= xborderw || minx >= w - xborderw || maxy <= yborderw || miny >= h - yborderw)
                                continue;

                            double len = sqrt((maxx - minx)*(maxx - minx) + (maxy - miny)*(maxy - miny));
                            double x1 = meanx - cos(theta) * len / 2;
                            double x2 = meanx + cos(theta) * len / 2;
                            double y1 = meany - sin(theta) * len / 2;
                            double y2 = meany + sin(theta) * len / 2;

                            lines.push_back({ Vec2(x1, y1), Vec2(x2, y2) });
                        }
                    }
                }

                std::cout << "done" << std::endl;

            }


            void ExtractLinesUsingLSD(const cv::Mat & im, std::vector<Line2> & lines, double minlen, int xbwidth, int ybwidth,
                std::vector<double> * lineWidths = nullptr,
                std::vector<double> * anglePrecisions = nullptr,
                std::vector<double> * negLog10NFAs = nullptr){

                cv::Mat gim;
                cv::cvtColor(im, gim, CV_BGR2GRAY);
                gim.convertTo(gim, CV_64FC1);
                //cv::imshow("gim", gim);
                //cv::waitKey();

                int h = gim.rows;
                int w = gim.cols;

                double * imgData = new double[h * w];
                for (int x = 0; x < w; x++){
                    for (int y = 0; y < h; y++){
                        imgData[x + y * w] = gim.at<double>(cv::Point(x, y));
                    }
                }
                
                int nOut;

                // x1,y1,x2,y2,width,p,-log10(NFA)
                double * linesData = lsd(&nOut, imgData, w, h);

                lines.clear();
                lines.reserve(nOut);
                if (lineWidths){
                    lineWidths->clear();
                    lineWidths->reserve(nOut);
                }
                if (anglePrecisions){
                    anglePrecisions->clear();
                    anglePrecisions->reserve(nOut);
                }
                if (negLog10NFAs){
                    negLog10NFAs->clear();
                    negLog10NFAs->reserve(nOut);
                }
                for (int i = 0; i < nOut; i++){
                    Line2 line = { Point2(linesData[7 * i + 0], linesData[7 * i + 1]), Point2(linesData[7 * i + 2], linesData[7 * i + 3]) };
                    if (line.length() < minlen)
                        continue;
                    if (line.first[0] <= xbwidth && line.second[0] <= xbwidth || 
                        line.first[0] >= w - xbwidth && line.second[0] >= w - xbwidth || 
                        line.first[1] <= ybwidth && line.second[1] <= ybwidth || 
                        line.first[1] >= h - ybwidth && line.second[1] >= h - ybwidth){
                        continue;
                    }
                    lines.emplace_back(line);
                    if (lineWidths){
                        lineWidths->push_back(linesData[7 * i + 4]);
                    }
                    if (anglePrecisions){
                        anglePrecisions->push_back(linesData[7 * i + 5]);
                    }
                    if (negLog10NFAs){
                        negLog10NFAs->push_back(linesData[7 * i + 6]);
                    }
                }

                delete[] linesData;
                delete[] imgData;
            }
        }


        LineSegmentExtractor::Feature LineSegmentExtractor::operator() (const Image & im) const {
            LineSegmentExtractor::Feature lines;
            if (_params.useLSD) {
                ExtractLinesUsingLSD(im, lines, _params.minLength, _params.xBorderWidth, _params.yBorderWidth);
            }
            else{
                ExtractLines(im, lines, _params.minLength, _params.xBorderWidth, _params.yBorderWidth, _params.numDirs);
            }
            return lines;
        }


        std::vector<HPoint2> ComputeLineIntersections(const std::vector<Line2> & lines,
            std::vector<std::pair<int, int>> * lineids,
            bool suppresscross, 
            double minDistanceOfLinePairs) {

            std::vector<HPoint2> hinterps;
            
            size_t lnum = lines.size();
            for (int i = 0; i < lnum; i++){
                Vec3 eqi = cv::Vec3d(lines[i].first[0], lines[i].first[1], 1)
                    .cross(cv::Vec3d(lines[i].second[0], lines[i].second[1], 1));
                for (int j = i + 1; j < lnum; j++){
                    if (minDistanceOfLinePairs < std::numeric_limits<double>::max()){
                        if (DistanceBetweenTwoLines(lines[i], lines[j]).first < minDistanceOfLinePairs)
                            continue;
                    }

                    Vec3 eqj = cv::Vec3d(lines[j].first[0], lines[j].first[1], 1)
                        .cross(cv::Vec3d(lines[j].second[0], lines[j].second[1], 1));
                    Vec3 interp = eqi.cross(eqj);
                    if (interp[0] == 0 && interp[1] == 0 && interp[2] == 0){ // lines overlapped
                        interp[0] = -eqi[1];
                        interp[1] = eqi[0];
                    }
                    interp /= norm(interp);

                    if (suppresscross){
                        auto& a1 = lines[i].first;
                        auto& a2 = lines[i].second;
                        auto& b1 = lines[j].first;
                        auto& b2 = lines[j].second;
                        double q = a1[0] * b1[1] - a1[1] * b1[0] - a1[0] * b2[1] + a1[1] * b2[0] -
                            a2[0] * b1[1] + a2[1] * b1[0] + a2[0] * b2[1] - a2[1] * b2[0];
                        double t = (a1[0] * b1[1] - a1[1] * b1[0] - a1[0] * b2[1] +
                            a1[1] * b2[0] + b1[0] * b2[1] - b1[1] * b2[0]) / q;
                        if (t > 0 && t < 1 && t == t)
                            continue;
                    }
                    hinterps.push_back(HPointFromVector(interp));
                    if (lineids)
                        lineids->emplace_back(i, j);
                }
            }

            return hinterps;

        }


        namespace {

            inline Point2 ToPoint2(const PixelLoc & p) {
                return Point2(p.x, p.y);
            }

            template <class T>
            inline Point2 ToPoint2(const Point<T, 2> & p) {
                return Point2(static_cast<double>(p[0]), static_cast<double>(p[1]));
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


        std::pair<double, InfiniteLine2> ComputeStraightness(const std::vector<std::vector<PixelLoc>> & edges,
            double * interleavedArea, double * interleavedLen) {
            
            std::vector<Point<float, 2>> points;
            for (auto & e : edges) {
                for (auto & p : e) {
                    points.push_back(Point<float, 2>(p.x, p.y));
                }
            }

            cv::Vec4f line;
            cv::fitLine(points, line, CV_DIST_L2, 0, 0.01, 0.01);
            auto fittedLine = InfiniteLine2({ line[2], line[3] }, { line[0], line[1] });
            double interArea = 0;
            double interLen = 0;
            for (auto & e : edges) {
                for (int i = 0; i < e.size() - 1; i++) {
                    double area, len;
                    std::tie(area, len) = ComputeSpanningArea(
                        ToPoint2(e[i]),
                        ToPoint2(e[i + 1]),
                        fittedLine);
                    interArea += area;
                    interLen += len;
                }
            }

            if (interleavedArea)
                *interleavedArea = interArea;
            if (interleavedLen)
                *interleavedLen = interLen;
            double straightness = Gaussian(interArea / interLen, 1.0);
            if (edges.size() == 1 && edges.front().size() == 2) {
                assert(FuzzyEquals(straightness, 1.0, 0.01) && "simple line should has the best straightness..");
            }

            return std::make_pair(straightness, fittedLine);

        }




        namespace {

            struct Edge {
                float w;
                int a, b;
            };

            class Universe {
            public:
                struct Element {
                    int rank;
                    int p;
                    int size;
                };
                inline Universe(int eleNum) : elements(eleNum), num(eleNum) {
                    for (int i = 0; i < eleNum; i++){
                        elements[i].rank = 0;
                        elements[i].size = 1;
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
                    } else {
                        elements[x].p = y;
                        elements[y].size += elements[x].size;
                        if (elements[x].rank == elements[y].rank)
                            elements[y].rank++;
                    }
                    num--;
                }
                inline int size(int x) const { return elements[x].size; }
                inline int numSets() const { return num; }
            private:
                int num;
                std::vector<Element> elements;
            };

            inline float Threshold(int size, float c) {
                return c / size;
            }

            Universe SegmentGraph(int numVertices, std::vector<Edge> & edges, float c) {
                std::sort(edges.begin(), edges.end(), [](const Edge & e1, const Edge & e2){
                    return e1.w < e2.w;
                });

                Universe u(numVertices);
                std::vector<float> threshold(numVertices);

                for (int i = 0; i < numVertices; i++)
                    threshold[i] = Threshold(1, c);

                for (int i = 0; i < edges.size(); i++) {
                    const Edge & edge = edges[i];

                    // components conected by this edge
                    int a = u.find(edge.a);
                    int b = u.find(edge.b);
                    if (a != b) {
                        if ((edge.w <= threshold[a]) &&
                            (edge.w <= threshold[b])) {
                            u.join(a, b);
                            a = u.find(a);
                            threshold[a] = edge.w + Threshold(u.size(a), c);
                        }
                    }
                }

                return u;
            }

            inline float PixelDiff(const Image & im, const cv::Point & p1, const cv::Point & p2) {
                assert(im.depth() == CV_8U && im.channels() == 3);
                Vec3 c1 = im.at<cv::Vec<uint8_t, 3>>(p1);
                Vec3 c2 = im.at<cv::Vec<uint8_t, 3>>(p2);
                return static_cast<float>(norm(c1 - c2));
            }

            // first return is CV_32SC1, the second is CV_8UC3 (for display)
            std::pair<Image, Image> SegmentImage(const Image & im, float sigma, float c, int minSize, 
                int & numCCs, bool returnColoredResult = false) {

                assert(im.depth() == CV_8U && im.channels() == 3);
                //std::cout << "depth: " << ImageDepth2Str(im.depth()) << std::endl;
                //std::cout << "channels: " << im.channels() << std::endl;

                int width = im.cols;
                int height = im.rows;
                Image smoothed;
                cv::GaussianBlur(im, smoothed, cv::Size(5, 5), sigma);

                // build pixel graph
                std::vector<Edge> edges;
                edges.reserve(width * height * 4);
                for (int y = 0; y < height; y++) {
                    for (int x = 0; x < width; x++) {
                        if (x < width - 1) {
                            Edge edge;
                            edge.a = y * width + x;
                            edge.b = y * width + (x + 1);
                            edge.w = PixelDiff(smoothed, { x, y }, { x + 1, y });
                            edges.push_back(edge);
                        }

                        if (y < height - 1) {
                            Edge edge;
                            edge.a = y * width + x;
                            edge.b = (y + 1) * width + x;
                            edge.w = PixelDiff(smoothed, { x, y }, { x, y + 1 });
                            edges.push_back(edge);
                        }

                        if ((x < width - 1) && (y < height - 1)) {
                            Edge edge;
                            edge.a = y * width + x;
                            edge.b = (y + 1) * width + (x + 1);
                            edge.w = PixelDiff(smoothed, { x, y }, { x + 1, y + 1 });
                            edges.push_back(edge);
                        }

                        if ((x < width - 1) && (y > 0)) {
                            Edge edge;
                            edge.a = y * width + x;
                            edge.b = (y - 1) * width + (x + 1);
                            edge.w = PixelDiff(smoothed, { x, y }, { x + 1, y - 1 });
                            edges.push_back(edge);
                        }
                    }
                }

                int num = (int)edges.size();
                Universe u = SegmentGraph(width * height, edges, c);

                for (int i = 0; i < num; i++) {
                    int a = u.find(edges[i].a);
                    int b = u.find(edges[i].b);
                    if ((a != b) && ((u.size(a) < minSize) || (u.size(b) < minSize)))
                        u.join(a, b);
                }

                numCCs = u.numSets();
                std::unordered_map<int, int> compIntSet;
                Image output(im.size(), CV_32SC1);
                for (int y = 0; y < height; y++) {
                    for (int x = 0; x < width; x++) {
                        int comp = u.find(y * width + x);
                        if (compIntSet.find(comp) == compIntSet.end()){
                            compIntSet.insert(std::make_pair(comp, (int)compIntSet.size()));
                        }
                        output.at<int32_t>(cv::Point(x, y)) = compIntSet[comp];
                    }
                }
                assert(compIntSet.size() == numCCs);

                if (!returnColoredResult){
                    return std::make_pair(output, Image());
                }

                Image coloredOutput(im.size(), CV_8UC3);
                std::vector<cv::Vec<uint8_t, 3>> colors(numCCs);
                std::generate(colors.begin(), colors.end(), [](){
                    return cv::Vec<uint8_t, 3>(uint8_t(std::rand() % 256), 
                        uint8_t(std::rand() % 256), 
                        uint8_t(std::rand() % 256));
                });
                for (int y = 0; y < height; y++) {
                    for (int x = 0; x < width; x++) {
                        coloredOutput.at<cv::Vec<uint8_t, 3>>(cv::Point(x, y)) = 
                            colors[output.at<int32_t>(cv::Point(x, y))];
                    }
                }

                //cv::imshow("result", coloredOutput);
                //cv::waitKey();

                return std::make_pair(output, coloredOutput);
            }


            std::pair<ImageWithType<int32_t>, int> SegmentImageUsingSLIC(const Image & im, int spsize, int spnum){
                assert(im.depth() == CV_8U && im.channels() == 3);

                SLIC slic;
                unsigned int * ubuff = new unsigned int[im.cols * im.rows];
                // fill buffer
                for (int x = 0; x < im.cols; x++){
                    for (int y = 0; y < im.rows; y++){
                        auto & pixel = ubuff[y * im.cols + x];
                        auto & color = im.at<Vec3b>(PixelLoc(x, y));
                        pixel = (color[2] << 16) + (color[1] << 8) + color[0];
                    }
                }

                int * klabels = nullptr;
                int nlabels = 0;
                if (spsize > 0){
                    slic.DoSuperpixelSegmentation_ForGivenSuperpixelSize(ubuff, im.cols, im.rows, klabels, nlabels, spsize, 50);
                }
                else{
                    slic.DoSuperpixelSegmentation_ForGivenNumberOfSuperpixels(ubuff, im.cols, im.rows, klabels, nlabels, spnum, 50);
                }

                // fill labels back
                Image labels = Image::zeros(im.size(), CV_32SC1);
                for (int x = 0; x < im.cols; x++){
                    for (int y = 0; y < im.rows; y++){
                        labels.at<int32_t>(PixelLoc(x, y)) = klabels[y * im.cols + x];
                    }
                }

                delete[] ubuff;
                delete[] klabels;

                return std::make_pair(labels, nlabels);
            }

        }


        std::pair<ImageWithType<int32_t>, int> SegmentationExtractor::operator() (const Image & im) const {
            if (_params.useSLIC){
                return SegmentImageUsingSLIC(im, _params.superpixelSizeSuggestion, _params.superpixelNumberSuggestion);
            }
            else{
                int numCCs;
                Image segim = SegmentImage(im, _params.sigma, _params.c, _params.minSize, numCCs, false).first;
                return std::make_pair(segim, numCCs);
            }
        }






        namespace {

            struct LineVPScoreFunctor {
                inline LineVPScoreFunctor(double angleThres = M_PI / 3.0, double s = 0.1) 
                : angleThreshold(angleThres), sigma(s) {}
                inline double operator()(double angle, bool liesOnLine) const {
                    if (angle >= angleThreshold)
                        return 0;
                    if (liesOnLine)
                        return 0;
                    double vote = (1 - (1 / angleThreshold) * angle);
                    vote = exp(-Square(1 - vote) / sigma / sigma / 2);
                    return vote;
                }
                double angleThreshold, sigma;
            };

            template <class HPoint2ContainerT, class LineVPScoreFunctorT = LineVPScoreFunctor>
            ImageWithType<double> LinesVotesToPoints(const HPoint2ContainerT & points, const std::vector<Line2> & lines,
                LineVPScoreFunctorT && scoreFun = LineVPScoreFunctorT()) {
                size_t nlines = lines.size();
                size_t npoints = points.size();
                ImageWithType<double> votes = ImageWithType<double>::zeros(nlines, npoints);
                for (int i = 0; i < nlines; i++) {
                    auto & line = lines[i];
                    for (int j = 0; j < npoints; j++) {
                        const HPoint2 & point = points[j];
                        Vec2 mid2vp = (point - HPoint2(line.center())).value();
                        // project point on line
                        double proj = mid2vp.dot(normalize(line.direction()));
                        bool liesOnLine = abs(proj) <= line.length() / 2.0;
                        double angle = AngleBetweenUndirectedVectors(mid2vp, line.direction());
                        double score = scoreFun(angle, liesOnLine);
                        if (std::isinf(score) || std::isnan(score)) {
                            score = 0.0;
                        }
                        votes(i, j) = score;
                    }
                }
                return votes;
            }

            std::vector<int> ClassifyLines(const ImageWithType<double> & votes, double scoreThreshold = 0.5) {
                std::vector<int> lineClasses(votes.rows, -1);
                int nlines = votes.rows;
                int npoints = votes.cols;
                for (size_t i = 0; i < nlines; i++) {
                    // classify lines
                    lineClasses[i] = -1;
                    double curscore = scoreThreshold;
                    for (int j = 0; j < npoints; j++) {
                        if (votes(i, j) >= curscore) {
                            lineClasses[i] = j;
                            curscore = votes(i, j);
                        }
                    }
                }
                return lineClasses;
            }


            template <class T>
            inline Vec<T, 2> PerpendicularDirection(const Vec<T, 2> & d) {
                return Vec<T, 2>(-d(1), d(0));
            }

            template <class T>
            inline Vec<T, 3> PerpendicularRootOfLineEquation(const Vec<T, 3> & lineeq) {
                auto & a = lineeq[0];
                auto & b = lineeq[1];
                auto & c = lineeq[2];
                return Vec<T, 3>(-a * c, -b * c, a * a + b * b);
            }

            std::pair<Point2, double> ComputeProjectionCenterAndFocalLength(const Point2 & vp1, const Point2 & vp2, const Point2 & vp3) {
                /* lambda = dot(vp1 - vp3, vp2 - vp3, 2) . / ...
                ((vp1(:, 1) - vp2(:, 1)).*(vp1(:, 2) - vp3(:, 2)) - ...
                (vp1(:, 1) - vp3(:, 1)).*(vp1(:, 2) - vp2(:, 2)));
                principalPoint(infcounts == 0, :) = vp3 + (vp1(:, [2 1]) - vp2(:, [2 1])).*...
                [-lambda lambda];
                focalLength(infcounts == 0) = ...
                sqrt(-dot(vp1 - principalPoint(infcounts == 0, :), ...
                vp2 - principalPoint(infcounts == 0, :), 2));*/
                auto lambda = (vp1 - vp3).dot(vp2 - vp3) /
                    ((vp1(0) - vp2(0))*(vp1(1) - vp3(1)) -
                    (vp1(0) - vp3(0))*(vp1(1) - vp2(1)));
                Point2 pp = vp3 + PerpendicularDirection(vp1 - vp2) * lambda;
                double focalLength = sqrt(-(vp1 - pp).dot(vp2 - pp));
                return std::make_pair(pp, focalLength);
            }


            std::vector<std::pair<Point2, double>> ComputeProjectionCenterAndFocalLength(
                const std::vector<Point2> & vp1s, const std::vector<Point2> & vp2s, const Point2 & vp3) {
                assert(vp1s.size() == vp2s.size());
                std::vector<std::pair<Point2, double>> ppAndFocals(vp1s.size());
                
                using namespace Eigen;
                Array<double, Dynamic, 2> vp1m, vp2m;
                vp1m.resize(vp1s.size(), 2);
                vp2m.resize(vp2s.size(), 2);

                for (int i = 0; i < vp1s.size(); i++){
                    vp1m(i, 0) = vp1s[i][0];
                    vp1m(i, 1) = vp1s[i][1];
                    vp2m(i, 0) = vp2s[i][0];
                    vp2m(i, 1) = vp2s[i][1];
                }

                Array<double, 1, 2> vp3m;
                vp3m(0, 0) = vp3[0];
                vp3m(0, 1) = vp3[1];

                auto lambdaUppers = (vp1m.rowwise() - vp3m).cwiseProduct(vp2m.rowwise() - vp3m).rowwise().sum();
                auto lambdaLowers = (vp1m.col(0) - vp2m.col(0)).cwiseProduct(vp1m.col(1) - vp3m(1))
                    - (vp1m.col(0) - vp3m(0)).cwiseProduct(vp1m.col(1) - vp2m.col(1));
                auto lambdas = lambdaUppers / lambdaLowers;

                Matrix<double, 2, 2> perpendicular;
                perpendicular.setZero();
                perpendicular(0, 1) = 1;
                perpendicular(1, 0) = -1;

                Array<double, Dynamic, 2> pps = (((vp1m - vp2m).matrix() * perpendicular).array().colwise() * lambdas).rowwise() + vp3m;
                auto vp1_pp = vp1m - pps;
                auto vp2_pp = vp2m - pps;
                Array<double, Dynamic, 1> focalLengths = sqrt(-(vp1_pp.col(0) * vp2_pp.col(0) + vp1_pp.col(1) + vp2_pp.col(1)));

                for (int i = 0; i < vp1s.size(); i++){
                    ppAndFocals[i].first = { pps(i, 0), pps(i, 1) };
                    ppAndFocals[i].second = focalLengths(i);
                }
                return ppAndFocals;
            }

            
            // make infinite points finite
            // merge close points
            void RefineIntersections(std::vector<HPoint2> & intersections, std::vector<std::pair<int, int>> & intersectionMakerLineIds, 
                double distanceThreshold = 2.0) {
                for (auto & hp : intersections) {
                    if (hp.denominator == 0.0)
                        hp.denominator = 1e-5;
                }

                std::vector<HPoint2> mergedIntersections;
                mergedIntersections.reserve(intersections.size());
                std::vector<std::pair<int, int>> mergedIntersectionMakerLineIds;
                mergedIntersectionMakerLineIds.reserve(intersections.size());

                RTreeWrapper<HPoint2> rtreeRecorder;
                for (int i = 0; i < intersections.size(); i++) {
                    if (rtreeRecorder.contains(intersections[i], [&distanceThreshold](const HPoint2 & a, const HPoint2 & b) {
                        return Distance(a, b) < distanceThreshold;
                    })) {
                        continue;
                    }
                    rtreeRecorder.insert(intersections[i]);
                    mergedIntersections.push_back(intersections[i]);
                    mergedIntersectionMakerLineIds.push_back(intersectionMakerLineIds[i]);
                }

                intersections = std::move(mergedIntersections);
                intersectionMakerLineIds = std::move(mergedIntersectionMakerLineIds);
            }

            std::vector<Vec3> RefineIntersectionsAndProjectToSpace(const std::vector<HPoint2> & intersections, 
                double fakeFocal, double angleThres) {
                std::vector<Vec3> dirs;
                dirs.reserve(intersections.size());

                RTreeWrapper<Vec3> rtreeRecorder;
                for (int i = 0; i < intersections.size(); i++){
                    Vec3 inter = normalize(VectorFromHPoint(intersections[i], fakeFocal));
                    if (rtreeRecorder.contains(inter, [&angleThres](const Vec3 & a, const Vec3 & b){
                        return AngleBetweenDirections(a, b) < angleThres;
                    })){
                        continue;
                    }
                    rtreeRecorder.insert(inter);
                    dirs.push_back(inter);
                }

                return dirs;
            }

            ImageWithType<double> GetLineLengthRatios(const std::vector<Line2> & lines){
                // get max line length
                double maxLineLen = 0;
                for (auto & line : lines) {
                    if (line.length() > maxLineLen)
                        maxLineLen = line.length();
                }
                ImageWithType<double> lineLengthRatios(lines.size(), 1);
                for (int i = 0; i < lines.size(); i++){
                    lineLengthRatios(i) = lines[i].length() / maxLineLen;
                }
                return lineLengthRatios;
            }


            inline HPoint2 ProjectOnToImagePlane(const Vec3 & d, const Point2 & pp, double focal) {
                return HPointFromVector(d, focal) + HPoint2(pp);
            }


            int AppendTheBestNPPAndFocalData(
                const std::vector<Point2> & vp2cands, const std::vector<int> & vp2candIdInRemainedIntersections,
                const std::vector<Point2> & vp3cands, const std::vector<int> & vp3candIdInRemainedIntersections,
                const Point2 & vp1p, int vp1id,
                const ImageWithType<double> & votesPanel, const ImageWithType<double> & votesRemainedPanel,
                const std::vector<Line2> & lines,
                double maxPrinciplePointOffset, double minFocalLength, double maxFocalLength,
                std::vector<std::pair<Point2, double>> & ppAndFocals, std::vector<float> & scores,
                int N){

                auto ppAndFocalsThisTime = ComputeProjectionCenterAndFocalLength(vp2cands, vp3cands, vp1p);

                // compute scores
                std::cout << ".";
                std::vector<float> scoresThisTime(ppAndFocalsThisTime.size(), 0);

                for (int i = 0; i < ppAndFocalsThisTime.size(); i++){
                    int vp2candId = vp2candIdInRemainedIntersections[i];
                    int vp3candId = vp3candIdInRemainedIntersections[i];

                    auto & principlePoint = ppAndFocalsThisTime[i].first;
                    auto & focalLength = ppAndFocalsThisTime[i].second;
                    if (std::isnan(focalLength) || std::isinf(focalLength)){
                        continue;
                    }
                    if (!(norm(principlePoint) < maxPrinciplePointOffset &&
                        IsBetween(focalLength, minFocalLength, maxFocalLength))){
                        continue;
                    }
                    // get votes from each line
                    float score = 0.0;
                    for (int lineId = 0; lineId < lines.size(); lineId++) {
                        double voteForVP1 = votesPanel(lineId, vp1id);
                        double voteForVP2 = votesRemainedPanel(lineId, vp2candId);
                        double voteForVP3 = votesRemainedPanel(lineId, vp3candId);
                        score += std::max({ voteForVP1, voteForVP2, voteForVP3 });
                    }
                    scoresThisTime[i] = score;
                }

                // sort
                std::vector<int> sortedPPAndFocalIds(scoresThisTime.size());
                for (int i = 0; i < scoresThisTime.size(); i++)
                    sortedPPAndFocalIds[i] = i;

                std::sort(sortedPPAndFocalIds.begin(), sortedPPAndFocalIds.end(), [&scoresThisTime](int a, int b){
                    return scoresThisTime[a] > scoresThisTime[b];
                });

                int keptSize = std::min<size_t>(N, ppAndFocalsThisTime.size());
                scores.reserve(scores.size() + scoresThisTime.size());
                ppAndFocals.reserve(ppAndFocals.size() + ppAndFocalsThisTime.size());

                for (int i = 0; i < keptSize; i++){
                    scores.push_back(scoresThisTime[sortedPPAndFocalIds[i]]);
                    ppAndFocals.push_back(ppAndFocalsThisTime[sortedPPAndFocalIds[i]]);
                }
                return scoresThisTime.size();
            }

        }


        std::tuple<std::array<HPoint2, 3>, double, std::vector<int>> VanishingPointsDetector::estimateWithProjectionCenterAtOrigin (
            const std::vector<Line2> & lines) const {

            // get max line length
            auto lineLengthRatios = GetLineLengthRatios(lines);

            std::vector<int> lineClasses(lines.size(), -1);

            // get all intersections
            std::vector<std::pair<int, int>> intersectionMakerLineIds;
            auto intersections = ComputeLineIntersections(lines, &intersectionMakerLineIds, true);

            //std::cout << "intersection num: " << intersections.size() << std::endl;        
            RefineIntersections(intersections, intersectionMakerLineIds);
            //std::cout << "intersection num: " << intersections.size() << std::endl;

            // nlines x npoints (without consideration of line length ratios)
            ImageWithType<double> votesPanel = LinesVotesToPoints(intersections, lines);

            // vote all lines for all intersections (with consideration of line length ratios)
            std::vector<double> votesForIntersections(intersections.size(), 0.0);            
            for (int i = 0; i < intersections.size(); i++) {
                votesForIntersections[i] = votesPanel.col(i).dot(lineLengthRatios);
            }

            // get vp1
            auto intersectionIdWithMaxVotes = std::distance(votesForIntersections.begin(),
                std::max_element(votesForIntersections.begin(), votesForIntersections.end()));
            const HPoint2 & vp1 = intersections[intersectionIdWithMaxVotes];
            std::cout << "vp1: " << vp1.value() << std::endl;
            std::cout << "score: " << votesForIntersections[intersectionIdWithMaxVotes] << std::endl;


            // classify lines for vp1 to be 0, and collect remained lines
            std::vector<Line2> remainedLines;
            remainedLines.reserve(lines.size() / 2);
            for (int i = 0; i < lines.size(); i++) {
                if (votesPanel(i, intersectionIdWithMaxVotes) > 0.8) {
                    lineClasses[i] = 0;
                } else {
                    remainedLines.push_back(lines[i]);
                }
            }
            //std::cout << "remained lines num: " << remainedLines.size() << std::endl;




            // get remained intersections
            std::vector<std::pair<int, int>> remainedIntersectionMakerLineIds;
            auto remainedIntersections = ComputeLineIntersections(remainedLines, &remainedIntersectionMakerLineIds, true);

            //std::cout << "remained intersection num: " << remainedLines.size() << std::endl;
            RefineIntersections(remainedIntersections, remainedIntersectionMakerLineIds);
            //std::cout << "remained intersection num: " << remainedLines.size() << std::endl;

            // vote remained lines for remained intersections
            std::vector<double> votesForRemainedIntersections(remainedIntersections.size(), 0.0);
            // nlines x npoints
            ImageWithType<double> votesRemainedPanel = LinesVotesToPoints(remainedIntersections, lines); // all lines participated!!!!!
            for (int i = 0; i < remainedIntersections.size(); i++) {
                votesForRemainedIntersections[i] = cv::sum(votesRemainedPanel.col(i)).val[0];
            }


            //// traverse other vp pairs
            double curMaxScore = 0.0;
            HPoint2 vp2, vp3;
            double curFocal;
            Point2 curPrinciplePoint;

            Point2 vp1ccenter = vp1.value() / 2.0;
            double vp1cdist = norm(vp1).value();

            int maxNum = 5000;
            int count = 0;

            for (int i = 0; i < remainedIntersections.size(); i++) {
                if (count >= maxNum)
                    break;

                auto & vp2cand = remainedIntersections[i];

                if (Distance(vp2cand, vp1) < _params.minFocalLength)
                    continue;
                if (Distance(vp2cand.value(), vp1ccenter) < vp1cdist / 2.0 - _params.minFocalLength)
                    continue;

                Point2 vp12center = (vp1 + vp2cand).value() / 2.0;
                double vp12dist = norm(vp1 - vp2cand).value();

                for (int j = i + 1; j < remainedIntersections.size(); j++) {
                    if (count >= maxNum)
                        break;

                    auto & vp3cand = remainedIntersections[j];
                    if (Distance(vp3cand, vp1) < _params.minFocalLength ||
                        Distance(vp2cand, vp3cand) < _params.minFocalLength)
                        continue;
                    if (Distance(vp3cand.value(), vp12center) < vp12dist / 2.0)
                        continue;

                    double focalLength;
                    Point2 principlePoint;
                    std::tie(principlePoint, focalLength) =
                        ComputeProjectionCenterAndFocalLength(vp1.value(), vp2cand.value(), vp3cand.value());

                    if (std::isnan(focalLength) || std::isinf(focalLength))
                        continue;

                    if (norm(principlePoint) < _params.maxPrinciplePointOffset &&
                        IsBetween(focalLength, _params.minFocalLength, _params.maxFocalLength)) {
                        count++;

                        // get votes from each line
                        double score = 0.0;
                        for (int lineId = 0; lineId < lines.size(); lineId++) {
                            double voteForVP1 = votesPanel(lineId, intersectionIdWithMaxVotes);
                            double voteForVP2 = votesRemainedPanel(lineId, i);
                            double voteForVP3 = votesRemainedPanel(lineId, j);
                            score += std::max({ voteForVP1, voteForVP2, voteForVP3 });
                        }

                        if (score > curMaxScore) {
                            curMaxScore = score;
                            vp2 = vp2cand;
                            vp3 = vp3cand;
                            curFocal = focalLength;
                            curPrinciplePoint = principlePoint;
                        }
                    }
                }
            }

            std::array<HPoint2, 3> vps = { { vp1, vp2, vp3 } };

            lineClasses = ClassifyLines(LinesVotesToPoints(vps, lines), 0.5);
            return std::make_tuple(vps, curFocal, lineClasses);
        }

        std::tuple<std::array<HPoint2, 3>, double, std::vector<int>> VanishingPointsDetector::operator() (
            const std::vector<Line2> & lines, const Point2 & projCenter) const {
            std::vector<Line2> offsetedLines = lines;
            for (auto & line : offsetedLines) {
                line.first -= projCenter;
                line.second -= projCenter;
            }
            auto results = estimateWithProjectionCenterAtOrigin(offsetedLines);
            for (HPoint2 & vp : std::get<0>(results)) {
                vp = vp + HPoint2(projCenter);
            }
            return results;
        }

       



        LocalManhattanVanishingPointsDetector::Result LocalManhattanVanishingPointsDetector::estimateWithProjectionCenterAtOriginIII(
            const std::vector<Line2> & lines) const {

            //Point2 offseted(_params.image.cols / 2.0, _params.image.rows / 2.0);

            auto lineLengthRatios = GetLineLengthRatios(lines);

            std::cout << "line num: " << lines.size() << std::endl;

            Box2 bbox = BoundingBoxOfContainer(lines);
            //auto diag = bbox.maxCorner - bbox.minCorner;
            double scale = norm(bbox.minCorner, bbox.maxCorner);

            Result result;
            result.lineClasses.clear();
            result.lineClasses.resize(lines.size(), -1);


            //// find vertical vp (vp1)

            // get all intersections
            std::vector<std::pair<int, int>> intersectionMakerLineIds;
            auto intersections = ComputeLineIntersections(lines, &intersectionMakerLineIds, true);

            std::cout << "intersection num: " << intersections.size() << std::endl;
            RefineIntersections(intersections, intersectionMakerLineIds);
            std::cout << "intersection num: " << intersections.size() << std::endl;

            // vote all lines for all intersections
            std::vector<std::pair<int, double>> intersectionIdsWithVotes(intersections.size());

            // nlines x npoints
            ImageWithType<double> votesPanel = LinesVotesToPoints(intersections, lines); //// TODO: bottleneck!!!
            for (int i = 0; i < intersections.size(); i++) {
                intersectionIdsWithVotes[i].first = i;
                intersectionIdsWithVotes[i].second = votesPanel.col(i).dot(lineLengthRatios);
            }

            // sort all intersectoins in votings descending order
            std::sort(intersectionIdsWithVotes.begin(), intersectionIdsWithVotes.end(),
                [](const std::pair<int, double> & a, const std::pair<int, double> & b) -> bool {
                return a.second > b.second;
            });

            // get vertical vp with max votes
            std::pair<int, double> verticalIntersectionIdWithMaxVotes = { -1, 0.0 };
            for (auto & idWithVotes : intersectionIdsWithVotes){
                auto & direction = intersections[idWithVotes.first];
                double angle = AngleBetweenUndirectedVectors(Vec2(0, 1), direction.numerator);
                if (angle < _params.verticalVPAngleRange &&
                    norm(direction.value()) >  scale / 2.0 * _params.verticalVPMinDistanceRatioToCenter){
                    verticalIntersectionIdWithMaxVotes = idWithVotes;
                    break;
                }
            }

            {
                auto & direction = intersections[verticalIntersectionIdWithMaxVotes.first];
                double angle = AngleBetweenUndirectedVectors(Vec2(0, 1), direction.numerator);
                std::cout << "angle: " << angle << std::endl;
            }

            assert(verticalIntersectionIdWithMaxVotes.first != -1 &&
                "failed to find vertical vp! "
                "try "
                "1) increasing verticalVPAngleRange or "
                "2) decreasing verticalVPMinDistanceRatioToCenter");

            HPoint2 vp1 = intersections[verticalIntersectionIdWithMaxVotes.first];
            std::cout << "vp1: " << vp1.value() << std::endl;

            // refine vp1
            double thresholdVotesForVP1 = 0.7;

            for (int i = 0; i < lines.size(); i++) {
                if (votesPanel(i, verticalIntersectionIdWithMaxVotes.first) > thresholdVotesForVP1) {
                    result.lineClasses[i] = 0;
                }
                else{
                    result.lineClasses[i] = -1;
                }
            }

            //static const int refineRepTime = 5;
            //for (int k = 0; k < refineRepTime; k++){
            //    // recalculate vp1 using classified vp1 lines
            //    Point2 vp1IntersectionsSum(0, 0);
            //    double vp1IntersectionsWeightSum = 0;
            //    for (int i = 0; i < intersectionMakerLineIds.size(); i++){
            //        auto & p = intersectionMakerLineIds[i];
            //        if (result.lineClasses[p.first] == 0 && result.lineClasses[p.second] == 0){
            //            auto v = intersections[i].value();
            //            double w = 1.0 - Gaussian(lines[p.first].length() + lines[p.second].length(), 20.0);
            //            assert(!HasValue(v, IsInfOrNaN<double>) && !HasValue(w, IsInfOrNaN<double>));
            //            vp1IntersectionsSum += v * w;
            //            vp1IntersectionsWeightSum += w;
            //        }
            //    }
            //    HPoint2 nvp1 = HPoint2(vp1IntersectionsSum / vp1IntersectionsWeightSum);
            //    vp1 = nvp1;

            //    std::cout << "refined vp1: " << vp1.value() << std::endl;

            //    // reclassify lines
            //    ImageWithType<double> refinedVotesForVP1 = LinesVotesToPoints(std::array<HPoint2, 1>{{ vp1 }}, lines);
            //    for (int i = 0; i < lines.size(); i++) {
            //        if (refinedVotesForVP1(i) > thresholdVotesForVP1) {
            //            result.lineClasses[i] = 0;
            //        }
            //        else{
            //            result.lineClasses[i] = -1;
            //        }
            //    }
            //}

            // collect remained lines
            std::vector<Line2> remainedLines;
            remainedLines.reserve(lines.size() / 2);
            for (int i = 0; i < lines.size(); i++) {
                if (result.lineClasses[i] == -1) {
                    remainedLines.push_back(lines[i]);
                }
            }
            std::cout << "remained lines num: " << remainedLines.size() << std::endl;


            // get remained intersections
            std::vector<std::pair<int, int>> remainedIntersectionMakerLineIds;
            auto remainedIntersections = ComputeLineIntersections(remainedLines, &remainedIntersectionMakerLineIds, true);

            //std::cout << "remained intersection num: " << remainedLines.size() << std::endl;
            RefineIntersections(remainedIntersections, remainedIntersectionMakerLineIds);
            //std::cout << "remained intersection num: " << remainedLines.size() << std::endl;

            // vote remained lines for remained intersections
            std::vector<double> votesForRemainedIntersections(remainedIntersections.size(), 0.0);
            // nlines x npoints
            ImageWithType<double> votesRemainedPanel = LinesVotesToPoints(remainedIntersections, lines); // all lines participated!!!!!
            for (int i = 0; i < remainedIntersections.size(); i++) {
                votesForRemainedIntersections[i] = votesRemainedPanel.col(i).dot(lineLengthRatios);
            }


            // test pps and focals
            Point2 vp1p = vp1.value();
            Point2 vp1ccenter = vp1.value() / 2.0;
            double vp1cdist = norm(vp1).value();

            int maxNumPerTime = 5e6;
            // keep the best N
            static const int N = 1e2;

            std::vector<Point2> vp2cands, vp3cands;
            vp2cands.reserve(maxNumPerTime);
            vp3cands.reserve(maxNumPerTime);
            std::vector<int> vp2candIdInRemainedIntersections, vp3candIdInRemainedIntersections;
            vp2candIdInRemainedIntersections.reserve(maxNumPerTime);
            vp3candIdInRemainedIntersections.reserve(maxNumPerTime);

            std::vector<float> scores;
            std::vector<std::pair<Point2, double>> ppAndFocals;
            scores.reserve(maxNumPerTime * 4);
            ppAndFocals.reserve(maxNumPerTime * 4);

            std::cout << "collecting all vp2-vp3 combinations ..." << std::endl;

            for (int i = 0; i < remainedIntersections.size(); i++) {
                /*if (vp2cands.size() >= maxNum)
                    break;*/

                auto & vp2cand = remainedIntersections[i];
                Point2 vp2candp = vp2cand.value();

                if (Distance(vp2cand, vp1) < _params.minFocalLength)
                    continue;
                if (Distance(vp2candp, vp1ccenter) < vp1cdist / 2.0 - _params.minFocalLength)
                    continue;

                Point2 vp12center = (vp1p + vp2candp) / 2.0;
                double vp12dist = norm(vp1p - vp2candp);

                for (int j = i + 1; j < remainedIntersections.size(); j++) {
                    /*if (vp2cands.size() >= maxNum)
                        break;*/

                    auto & vp3cand = remainedIntersections[j];
                    Point2 vp3candp = vp3cand.value();

                    if (Distance(vp3cand, vp1) < _params.minFocalLength ||
                        Distance(vp2cand, vp3cand) < _params.minFocalLength)
                        continue;
                    if (Distance(vp3candp, vp12center) < vp12dist / 2.0)
                        continue;

                    vp2cands.push_back(vp2candp);
                    vp3cands.push_back(vp3candp);
                    vp2candIdInRemainedIntersections.push_back(i);
                    vp3candIdInRemainedIntersections.push_back(j);

                    if (vp2cands.size() >= maxNumPerTime){
                        std::cout << ".";
                        AppendTheBestNPPAndFocalData(vp2cands, vp2candIdInRemainedIntersections, vp3cands, vp3candIdInRemainedIntersections,
                            vp1p, verticalIntersectionIdWithMaxVotes.first, votesPanel, votesRemainedPanel,
                            lines, _params.maxPrinciplePointOffset, _params.minFocalLength, _params.maxFocalLength, 
                            ppAndFocals, scores, N);

                        vp2cands.clear();
                        vp3cands.clear();
                        vp2candIdInRemainedIntersections.clear();
                        vp3candIdInRemainedIntersections.clear();
                    }
                }
            }

            if(vp2cands.size() > 0) { // finalize
                std::cout << "." << std::endl;
                AppendTheBestNPPAndFocalData(vp2cands, vp2candIdInRemainedIntersections, vp3cands, vp3candIdInRemainedIntersections,
                    vp1p, verticalIntersectionIdWithMaxVotes.first, votesPanel, votesRemainedPanel,
                    lines, _params.maxPrinciplePointOffset, _params.minFocalLength, _params.maxFocalLength,
                    ppAndFocals, scores, N);

                vp2cands.clear();
                vp3cands.clear();
                vp2candIdInRemainedIntersections.clear();
                vp3candIdInRemainedIntersections.clear();
            }

            std::cout << "done collecting all valid vp2-vp3 combination scores and configs, total num: " << ppAndFocals.size() << std::endl;

            // sort and ...
            // only keep the best N
            {
                // sort
                std::vector<int> sortedPPAndFocalIds(scores.size());
                for (int i = 0; i < scores.size(); i++)
                    sortedPPAndFocalIds[i] = i;

                std::sort(sortedPPAndFocalIds.begin(), sortedPPAndFocalIds.end(), [&scores](int a, int b){
                    return scores[a] > scores[b];
                });
                
                int keptSize = std::min<size_t>(N, ppAndFocals.size());
                std::vector<float> keptScores;
                std::vector<std::pair<Point2, double>> keptPPAndFocals;
                keptScores.reserve(keptSize);
                keptPPAndFocals.reserve(keptSize);
                for (int i = 0; i < keptSize; i++){
                    keptScores.push_back(scores[sortedPPAndFocalIds[i]]);
                    keptPPAndFocals.push_back(ppAndFocals[sortedPPAndFocalIds[i]]);
                }
                scores = std::move(keptScores);
                ppAndFocals = std::move(keptPPAndFocals);
            }


            IF_DEBUG_USING_VISUALIZERS{
                for (int i = 0; i < std::min<size_t>(30, scores.size()); i++){
                    auto & principlePoint = ppAndFocals[i].first;
                    auto & focalLength = ppAndFocals[i].second;
                    std::cout << "focal: " << focalLength << "  pp: " << principlePoint << "   score: " << scores[i] << std::endl;
                }
            }


            // collect remained close lines
            static const double closeLinePairDistanceThreshold = 40;
            std::map<std::pair<int, int>, double> remainedCloseLineIdPairs;
            for (int i = 0; i < remainedLines.size(); i++){
                for (int j = i + 1; j < remainedLines.size(); j++){
                    double distance = DistanceBetweenTwoLines(remainedLines[i], remainedLines[j]).first;
                    if (distance < closeLinePairDistanceThreshold){
                        remainedCloseLineIdPairs[std::make_pair(i, j)] = distance;
                    }
                }
            }

            std::cout << "close line pair num: " << remainedCloseLineIdPairs.size() << std::endl;
            std::cout << "testing local manhattan consistency ..." << std::endl;
            std::map<int, double> lmanScores;

            for (int id = 0; id < scores.size(); id++){
                auto & principlePoint = ppAndFocals[id].first;
                auto & focalLength = ppAndFocals[id].second;

                if (std::isnan(focalLength) || std::isinf(focalLength))
                    continue;

                if (scores[id] <= 0.0)
                    continue;

                Vec3 vp1v = Concat(vp1p - principlePoint, focalLength);
                double lmanScore = 0;

                for (auto & p : remainedCloseLineIdPairs) {
                    auto & line1 = remainedLines[p.first.first];
                    auto & line2 = remainedLines[p.first.second];
                    double distance = p.second;

                    Vec3 line1eq = Concat(line1.first - principlePoint, focalLength)
                        .cross(Concat(line1.second - principlePoint, focalLength));
                    Vec3 line2eq = Concat(line2.first - principlePoint, focalLength)
                        .cross(Concat(line2.second - principlePoint, focalLength));

                    Vec3 inter1 = normalize(line1eq.cross(vp1v));
                    Vec3 inter2 = normalize(line2eq.cross(vp1v));

                    lmanScore += Gaussian(inter1.dot(inter2), 0.1) * Gaussian(distance, 20.0);
                }

                static const double lmFactor = 0.5;
                lmanScores[id] = lmanScore * lmFactor + scores[id] * (1 - lmFactor);
            }

            {
                // sort considering lmanScores
                std::vector<int> sortedPPAndFocalIds(scores.size());
                for (int i = 0; i < scores.size(); i++)
                    sortedPPAndFocalIds[i] = i;

                std::sort(sortedPPAndFocalIds.begin(), sortedPPAndFocalIds.end(), [&lmanScores](int a, int b){
                    return lmanScores[a] > lmanScores[b];
                });

                int keptSize = scores.size();
                std::vector<float> keptScores;
                std::vector<std::pair<Point2, double>> keptPPAndFocals;
                keptScores.reserve(keptSize);
                keptPPAndFocals.reserve(keptSize);
                for (int i = 0; i < keptSize; i++){
                    keptScores.push_back(scores[sortedPPAndFocalIds[i]]);
                    keptPPAndFocals.push_back(ppAndFocals[sortedPPAndFocalIds[i]]);
                }
                scores = std::move(keptScores);
                ppAndFocals = std::move(keptPPAndFocals);
            }


            IF_DEBUG_USING_VISUALIZERS{
                for (int i = 0; i < std::min<size_t>(30, scores.size()); i++){
                    auto & principlePoint = ppAndFocals[i].first;
                    auto & focalLength = ppAndFocals[i].second;
                    std::cout << "focal: " << focalLength << "  pp: " << principlePoint
                        << "   score: " << scores[i] << "  lmscore: " << lmanScores[i] << std::endl;
                }
            }
            
            result.principlePoint = ppAndFocals.front().first;
            result.focalLength = ppAndFocals.front().second;

            // get horizontal vanishing points
            result.vanishingPoints.clear();
            result.vanishingPoints.push_back(vp1);
            result.verticalVanishingPointId = 0;
            result.horizontalVanishingPointIds.clear();

            Vec3 vp1v = Concat(vp1p - result.principlePoint, result.focalLength);
            std::vector<std::pair<Vec3, Vec3>> hvpvs;
            for (auto & p : remainedCloseLineIdPairs) {
                auto & line1 = remainedLines[p.first.first];
                auto & line2 = remainedLines[p.first.second];
                double distance = p.second;

                Vec3 line1eq = Concat(line1.first - result.principlePoint, result.focalLength)
                    .cross(Concat(line1.second - result.principlePoint, result.focalLength));
                Vec3 line2eq = Concat(line2.first - result.principlePoint, result.focalLength)
                    .cross(Concat(line2.second - result.principlePoint, result.focalLength));

                assert(Distance(ProjectOnToImagePlane(Concat(line1.first - result.principlePoint, result.focalLength), 
                    result.principlePoint, result.focalLength).value(), line1.first) < 0.5);
                assert(Distance(ProjectOnToImagePlane(Concat(line2.first - result.principlePoint, result.focalLength), 
                    result.principlePoint, result.focalLength).value(), line2.first) < 0.5);

                Vec3 inter1 = normalize(line1eq.cross(vp1v));
                Vec3 inter2 = normalize(line2eq.cross(vp1v));

                if (inter1.dot(inter2) < 0.03 && distance < 30){
                    hvpvs.emplace_back(inter1, inter2);
                }
            }
            std::cout << "initial horizontal direction num: " << (hvpvs.size() * 2) << std::endl;

            std::vector<HPoint2> horizonVPs;
            horizonVPs.push_back(vp1);
            std::vector<std::pair<int, int>> orthoPairs;
            std::vector<int> hvpvIdOfVP; // hvpvIdOfVP[{id in horizonVPs}] = {id in hvpvs}
            hvpvIdOfVP.push_back(-1);
            std::vector<bool> vpIsAtFirstInHVPV;
            vpIsAtFirstInHVPV.push_back(false);
            for (int i = 0; i < hvpvs.size(); i++){
                const Vec3 & v1 = hvpvs[i].first;
                const Vec3 & v2 = hvpvs[i].second;
                HPoint2 hinter1 = ProjectOnToImagePlane(v1, result.principlePoint, result.focalLength);
                HPoint2 hinter2 = ProjectOnToImagePlane(v2, result.principlePoint, result.focalLength);
                horizonVPs.push_back(hinter1);
                hvpvIdOfVP.push_back(i);
                vpIsAtFirstInHVPV.push_back(true);
                horizonVPs.push_back(hinter2);
                hvpvIdOfVP.push_back(i);
                vpIsAtFirstInHVPV.push_back(false);
                orthoPairs.emplace_back(static_cast<int>(horizonVPs.size() - 2), static_cast<int>(horizonVPs.size() - 1));
            }

            auto lineClasses = ClassifyLines(LinesVotesToPoints(horizonVPs, lines), 0.5);

            // merge close horizontal vanishing point directions
            std::vector<std::pair<Vec3, Vec3>> mergedHVPVs;
            std::vector<int> oldHVPVIdToMergedHVPVId(hvpvs.size(), -1);
            std::vector<bool> oldHVPVIdToMergedHVPVIdSwapped(hvpvs.size(), -1);
            for (int i = 0; i < hvpvs.size(); i++){
                auto & hvpvPair = hvpvs[i];
                int nearestMergedHVPVId = -1;
                bool swapped = false; // whether old hvpvPair has to be swapped to match the nearestMergedHVPV
                double minDist = 0.03;
                for (int j = 0; j < mergedHVPVs.size(); j++){
                    auto & hvpvPairRecorded = mergedHVPVs[j];
                    double dist1 = std::min({
                        AngleBetweenUndirectedVectors(hvpvPair.first, hvpvPairRecorded.first),
                        AngleBetweenUndirectedVectors(hvpvPair.second, hvpvPairRecorded.second)
                    });
                    double dist2 = std::min({
                        AngleBetweenUndirectedVectors(hvpvPair.first, hvpvPairRecorded.second),
                        AngleBetweenUndirectedVectors(hvpvPair.second, hvpvPairRecorded.first)
                    });
                    double dist = std::min(dist1, dist2);
                    if (dist <= minDist){
                        minDist = dist;
                        nearestMergedHVPVId = j;
                        swapped = dist1 > dist2;
                    }                    
                }
                if (nearestMergedHVPVId == -1){
                    mergedHVPVs.push_back(hvpvPair);
                    nearestMergedHVPVId = mergedHVPVs.size() - 1;
                }
                oldHVPVIdToMergedHVPVId[i] = nearestMergedHVPVId;
                oldHVPVIdToMergedHVPVIdSwapped[i] = swapped;
            }

            for (int i = 0; i < mergedHVPVs.size(); i++){
                const Vec3 & v1 = mergedHVPVs[i].first;
                const Vec3 & v2 = mergedHVPVs[i].second;
                HPoint2 hinter1 = ProjectOnToImagePlane(v1, result.principlePoint, result.focalLength);
                HPoint2 hinter2 = ProjectOnToImagePlane(v2, result.principlePoint, result.focalLength);
                result.vanishingPoints.push_back(hinter1);
                result.vanishingPoints.push_back(hinter2);
                result.horizontalVanishingPointIds.emplace_back(
                    static_cast<int>(result.vanishingPoints.size() - 2), 
                    static_cast<int>(result.vanishingPoints.size() - 1));
            }

            // reassign line classes
            for (int i = 0; i < lineClasses.size(); i++){
                auto & line = lines[i];
                int c = lineClasses[i];
                if (c == 0)
                    continue;
                if (c == -1)
                    continue;

                int hvpvId = hvpvIdOfVP[c];
                if (hvpvId != -1){
                    bool isAtFirst = vpIsAtFirstInHVPV[lineClasses[i]];
                    if (oldHVPVIdToMergedHVPVIdSwapped[hvpvId]){
                        isAtFirst = !isAtFirst;
                    }
                    auto & orthoPair = result.horizontalVanishingPointIds[oldHVPVIdToMergedHVPVId[hvpvId]];
                    lineClasses[i] = isAtFirst ? orthoPair.first : orthoPair.second;
                }
            }
            result.lineClasses = lineClasses;
            
            std::cout << "merged horizontal direction num: " << (result.vanishingPoints.size() - 1) << std::endl;


            return result;
        }




        LocalManhattanVanishingPointsDetector::Result LocalManhattanVanishingPointsDetector::operator()(
            const std::vector<Line2> & lines, const Point2 & projCenter) const {
            std::vector<Line2> offsetedLines = lines;
            for (auto & line : offsetedLines) {
                line.first -= projCenter;
                line.second -= projCenter;
            }
            Result result = estimateWithProjectionCenterAtOriginIII(lines);
            
            for (HPoint2 & vp : result.vanishingPoints) {
                vp = vp + HPoint2(projCenter);
            }
            
            result.horizon.anchor += projCenter;
            for (auto & hline : result.hlineCands){
                hline.anchor += projCenter;
            }

            result.principlePoint += projCenter;
            
            return result;
        }





        std::pair<double, double> ComputeFocalsFromHomography(const Mat3 & H, std::pair<bool, bool> * ok) {
            /// from OpenCV code

            const double* h = reinterpret_cast<const double*>(H.val);

            double d1, d2; // Denominators
            double v1, v2; // Focal squares value candidates
            double f0, f1;

            if (ok)
                ok->second = true;
            d1 = h[6] * h[7];
            d2 = (h[7] - h[6]) * (h[7] + h[6]);
            v1 = -(h[0] * h[1] + h[3] * h[4]) / d1;
            v2 = (h[0] * h[0] + h[3] * h[3] - h[1] * h[1] - h[4] * h[4]) / d2;
            if (v1 < v2) std::swap(v1, v2);
            if (v1 > 0 && v2 > 0) f1 = sqrt(std::abs(d1) > std::abs(d2) ? v1 : v2);
            else if (v1 > 0) f1 = sqrt(v1);
            else if(ok) 
                ok->second = false;

            if (ok)
                ok->first = true;
            d1 = h[0] * h[3] + h[1] * h[4];
            d2 = h[0] * h[0] + h[1] * h[1] - h[3] * h[3] - h[4] * h[4];
            v1 = -h[2] * h[5] / d1;
            v2 = (h[5] * h[5] - h[2] * h[2]) / d2;
            if (v1 < v2) std::swap(v1, v2);
            if (v1 > 0 && v2 > 0) f0 = sqrt(std::abs(d1) > std::abs(d2) ? v1 : v2);
            else if (v1 > 0) f0 = sqrt(v1);
            else if(ok)
                ok->first = false;

            return std::make_pair(f0, f1);
        }




    }
}

    
