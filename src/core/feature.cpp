#include <iostream>
#include <unordered_map>
#include <thread>

extern "C" {
    #include <lsd.h>
}
#include <quickshift_common.h>

#include <SLIC.h>

#include <VPSample.h>
#include <VPCluster.h>

#include "MSAC.h"

#include "cameras.hpp"
#include "feature.hpp"
#include "utilities.hpp"
#include "containers.hpp"
#include "clock.hpp"

#include "../misc/matlab_engine.hpp"
#include "../misc/eigen.hpp"




namespace panoramix {
    namespace core {

        using misc::MatlabEngine;

#pragma region NonMaximaSuppression

        void NonMaximaSuppression(const Image & src, Image & dst, int sz, std::vector<PixelLoc> * pixels,
            const Imageb & mask) {

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

#pragma endregion NonMaximaSuppression



#pragma region LineSegmentExtractor

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
            if (_params.algorithm == LSD) {
                ExtractLinesUsingLSD(im, lines, _params.minLength, _params.xBorderWidth, _params.yBorderWidth);
            }
            else if(_params.algorithm == GradientGrouping){
                ExtractLines(im, lines, _params.minLength, _params.xBorderWidth, _params.yBorderWidth, _params.numDirs);
            }
            return lines;
        }

        LineSegmentExtractor::Feature LineSegmentExtractor::operator() (const Image & im, int pyramidHeight, int minSize) const{
            Feature lines;
            Image image = im.clone();
            for (int i = 0; i < pyramidHeight; i++){
                if (image.cols < minSize || image.rows < minSize)
                    break;
                Feature ls = (*this)(image);
                for (auto & l : ls){
                    lines.push_back(l * (double(im.cols) / image.cols));
                }
                cv::pyrDown(image, image);
            }
            return lines;
        }

#pragma endregion LineSegmentExtractor






        std::vector<HPoint2> ComputeLineIntersections(const std::vector<Line2> & lines,
            std::vector<std::pair<int, int>> * lineids,
            bool suppresscross, 
            double minDistanceOfLinePairs) {

            SetClock();

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




        void ClassifyLines(std::vector<Classified<Line2>> &lines, const std::vector<HPoint2> & vps,
            double angleThreshold, double sigma, double scoreThreshold, double avoidVPDistanceThreshold){

            for (auto & line : lines){
                // classify lines
                line.claz = -1;

                // classify
                std::vector<double> linescores(vps.size());

                // get score based on angle
                for (int j = 0; j < vps.size(); j++){
                    auto & point = vps[j];
                    double angle = std::min(
                        AngleBetweenDirections(line.component.direction(), (point - HPoint2(line.component.center())).numerator),
                        AngleBetweenDirections(-line.component.direction(), (point - HPoint2(line.component.center())).numerator));
                    double score = exp(-(angle / angleThreshold) * (angle / angleThreshold) / sigma / sigma / 2);
                    if (avoidVPDistanceThreshold >= 0.0 && 
                        DistanceFromPointToLine(point.value(), line.component).first < avoidVPDistanceThreshold){
                        linescores[j] = -1.0;
                    }
                    else{
                        linescores[j] = (angle > angleThreshold) ? 0 : score;
                    }
                }

                double curscore = scoreThreshold;
                for (int j = 0; j < vps.size(); j++){
                    if (linescores[j] > curscore){
                        line.claz = j;
                        curscore = linescores[j];
                    }
                }
            }
        }


        void ClassifyLines(std::vector<Classified<Line3>> & lines, const std::vector<Vec3> & vps,
            double angleThreshold, double sigma, double scoreThreshold, double avoidVPAngleThreshold) {

            size_t nlines = lines.size();
            size_t npoints = vps.size();

            for (size_t i = 0; i < nlines; i++){
                const Vec3 & a = lines[i].component.first;
                const Vec3 & b = lines[i].component.second;
                Vec3 normab = a.cross(b);
                normab /= norm(normab);

                std::vector<double> linescores(npoints);

                // get score based on angle
                for (int j = 0; j < npoints; j++){
                    const Vec3 & point = vps[j];
                    double angle = abs(asin(normab.dot(point)));
                    double score = exp(-(angle / angleThreshold) * (angle / angleThreshold) / sigma / sigma / 2);
                    if (avoidVPAngleThreshold >= 0.0 && std::min(
                            DistanceFromPointToLine(normalize(point), normalize(lines[i].component)).first,
                            DistanceFromPointToLine(normalize(-point), normalize(lines[i].component)).first
                        ) < 2.0 * sin(avoidVPAngleThreshold / 2.0)){ // avoid that a line belongs to its nearby vp
                        linescores[j] = -1.0;
                    }
                    else{
                        linescores[j] = (angle > angleThreshold) ? 0 : score;
                    }
                }

                // classify lines
                lines[i].claz = -1;
                double curscore = scoreThreshold;
                for (int j = 0; j < npoints; j++){

                    if (linescores[j] > curscore){
                        lines[i].claz = j;
                        curscore = linescores[j];
                    }
                }
            }
        }




        namespace {

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


        std::pair<double, Ray2> ComputeStraightness(const std::vector<std::vector<PixelLoc>> & edges,
            double * interleavedArea, double * interleavedLen) {
            
            std::vector<Point<float, 2>> points;
            for (auto & e : edges) {
                for (auto & p : e) {
                    points.push_back(Point<float, 2>(p.x, p.y));
                }
            }

            cv::Vec4f line;
            cv::fitLine(points, line, CV_DIST_L2, 0, 0.01, 0.01);
            auto fittedLine = Ray2({ line[2], line[3] }, { line[0], line[1] });
            double interArea = 0;
            double interLen = 0;
            for (auto & e : edges) {
                for (int i = 0; i < e.size() - 1; i++) {
                    double area, len;
                    std::tie(area, len) = ComputeSpanningArea(
                        vec_cast<double>(e[i]),
                        vec_cast<double>(e[i + 1]),
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
        }

        Failable<std::vector<Vec3>> FindOrthogonalPrinicipleDirections(const std::vector<Vec3> & intersections,
            int longitudeDivideNum, int latitudeDivideNum, bool allowMoreThan2HorizontalVPs, const Vec3 & verticalSeed) {

            std::vector<Vec3> vps(3);

            // collect votes of intersection directions
            Imagef votePanel = Imagef::zeros(longitudeDivideNum, latitudeDivideNum);
            size_t pn = intersections.size();
            for (const Vec3& p : intersections){
                PixelLoc pixel = PixelLocFromGeoCoord(GeoCoord(p), longitudeDivideNum, latitudeDivideNum);
                votePanel(pixel.x, pixel.y) += 1.0;
            }
            cv::GaussianBlur(votePanel, votePanel, cv::Size((longitudeDivideNum / 50) * 2 + 1, (latitudeDivideNum / 50) * 2 + 1),
                4, 4, cv::BORDER_REPLICATE);

            // set the direction with the max votes as the first vanishing point
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

                double score = 0;
                for (const Vec3 & v : { vec1, vec1rev, vec2, vec2rev }){
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

            if (UnOrthogonality(vps[0], vps[1], vps[2]) >= 0.1){
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
            }

            if (UnOrthogonality(vps[0], vps[1], vps[2]) >= 0.1){
                return nullptr; // failed
            }

            // make vps[0] the vertical vp
            {
                int vertVPId = -1;
                double minAngle = std::numeric_limits<double>::max();
                for (int i = 0; i < vps.size(); i++){
                    double a = AngleBetweenUndirectedVectors(vps[i], verticalSeed);
                    if (a < minAngle){
                        vertVPId = i;
                        minAngle = a;
                    }
                }
                std::swap(vps[0], vps[vertVPId]);
            }

            if (allowMoreThan2HorizontalVPs){            
                // find more horizontal vps
                double nextMaxScore = maxScore * 0.5; // threshold
                double minAngleToCurHorizontalVPs = DegreesToRadians(30);
                for (int x = 0; x < longitudeDivideNum; x++){
                    double longt1 = double(x) / longitudeDivideNum * M_PI * 2 - M_PI;
                    double lat1 = LatitudeFromLongitudeAndNormalVector(longt1, vps[0]);
                    Vec3 vec1 = GeoCoord(longt1, lat1).toVector();
                    Vec3 vec1rev = -vec1;
                    Vec3 vec2 = vps[0].cross(vec1);
                    Vec3 vec2rev = -vec2;

                    double score = 0;
                    for (const Vec3 & v : { vec1, vec1rev, vec2, vec2rev }){
                        PixelLoc pixel = PixelLocFromGeoCoord(GeoCoord(v), longitudeDivideNum, latitudeDivideNum);
                        score += votePanel(WrapBetween(pixel.x, 0, longitudeDivideNum),
                            WrapBetween(pixel.y, 0, latitudeDivideNum));
                    }

                    bool tooCloseToExistingVP = false;
                    for (int i = 0; i < vps.size(); i++){
                        if (AngleBetweenUndirectedVectors(vps[i], vec1) < minAngleToCurHorizontalVPs){
                            tooCloseToExistingVP = true;
                            break;
                        }
                    }
                    if (tooCloseToExistingVP)
                        continue;

                    if (score > nextMaxScore){
                        nextMaxScore = score;
                        vps.push_back(vec1);
                        vps.push_back(vec2);
                    }
                }
            }

            return std::move(vps);

        }

        int NearestDirectionId(const std::vector<Vec3> & directions,
            const Vec3 & verticalSeed){
            int vid = -1;
            double minAngle = M_PI;
            for (int i = 0; i < directions.size(); i++){
                double angle = AngleBetweenUndirectedVectors(directions[i], verticalSeed);
                if (angle < minAngle){
                    vid = i;
                    minAngle = angle;
                }
            }
            return vid;
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
                lineIntersections.push_back(normalize(cam.toSpace(p.value())));
            }

            auto vanishingPoints = FindOrthogonalPrinicipleDirections(lineIntersections, 1000, 500, true).unwrap();

            // project lines to space
            std::vector<Classified<Line3>> spatialLineSegments;
            spatialLineSegments.reserve(linesNum);
            for (const auto & line : lineSegments) {
                auto & p1 = line.component.first;
                auto & p2 = line.component.second;
                auto pp1 = cam.toSpace(p1);
                auto pp2 = cam.toSpace(p2);
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
                    lineIntersections.push_back(normalize(cams[i].toSpace(p.value())));
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
                    auto pp1 = cams[i].toSpace(p1);
                    auto pp2 = cams[i].toSpace(p2);
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





#pragma region VanishingPointsDetector

        namespace {

            template <class T>
            inline Vec<T, 3> PerpendicularRootOfLineEquation(const Vec<T, 3> & lineeq) {
                auto & a = lineeq[0];
                auto & b = lineeq[1];
                auto & c = lineeq[2];
                return Vec<T, 3>(-a * c, -b * c, a * a + b * b);
            }

        }

        std::pair<Point2, double> ComputePrinciplePointAndFocalLength(const Point2 & vp1, const Point2 & vp2, const Point2 & vp3) {
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
            double focalLength = sqrt(abs(-(vp1 - pp).dot(vp2 - pp)));
            return std::make_pair(pp, focalLength);
        }

        namespace {

            //struct LineVPScoreFunctor {
            //    inline LineVPScoreFunctor(double angleThres = M_PI / 3.0, double s = 0.1) 
            //        : angleThreshold(angleThres), sigma(s) {}
            //    inline double operator()(double angle, bool liesOnLine) const {
            //        if (angle >= angleThreshold)
            //            return 0;
            //        if (liesOnLine)
            //            return 0;
            //        double vote = (1 - (1 / angleThreshold) * angle);
            //        vote = exp(-Square(1 - vote) / sigma / sigma / 2);
            //        return vote;
            //    }
            //    double angleThreshold, sigma;
            //};

            Imaged LinesVotesToPoints(const std::vector<HPoint2> & points, const std::vector<Line2> & lines) {

                SetClock();

                static const double angleThreshold = M_PI / 3.0;
                static const double sigma = 0.1;

                size_t nlines = lines.size();
                size_t npoints = points.size();
                Imaged votes = Imaged::zeros(nlines, npoints);
                for (auto it = votes.begin(); it != votes.end(); ++it) {
                    auto & line = lines[it.pos().y];
                    const HPoint2 & point = points[it.pos().x];
                    Vec2 mid2vp = (point - HPoint2(line.center())).value();
                    // project point on line
                    double proj = mid2vp.dot(normalize(line.direction()));
                    bool liesOnLine = abs(proj) <= line.length() / 2.0;
                    double angle = AngleBetweenUndirectedVectors(mid2vp, line.direction());

                    double score = 0.0;
                    if (angle >= angleThreshold)
                        continue;
                    if (liesOnLine)
                        continue;
                    score = exp(-Square(angle / angleThreshold) / sigma / sigma / 2.0);
                    if (std::isinf(score) || std::isnan(score)) {
                        score = 0.0;
                    }
                    *it = score;
                }
                return votes;
            }

            std::vector<int> ClassifyLines(const Imaged & votes, double scoreThreshold = 0.5) {

                SetClock();

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


            std::vector<std::pair<Point2, double>> ComputeProjectionCenterAndFocalLength(
                const std::vector<Point2> & vp1s, const std::vector<Point2> & vp2s, const Point2 & vp3) {
                assert(vp1s.size() == vp2s.size());
                std::vector<std::pair<Point2, double>> ppAndFocals(vp1s.size());                

                THERE_ARE_BOTTLENECKS_HERE("use Eigen::Map!");

                using namespace Eigen;
                Array<double, Eigen::Dynamic, 2> vp1m, vp2m;
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

                Array<double, Eigen::Dynamic, 2> pps = (((vp1m - vp2m).matrix() * perpendicular).array().colwise() * lambdas).rowwise() + vp3m;
                auto vp1_pp = vp1m - pps;
                auto vp2_pp = vp2m - pps;
                Array<double, Eigen::Dynamic, 1> focalLengths = sqrt(-(vp1_pp.col(0) * vp2_pp.col(0) + vp1_pp.col(1) + vp2_pp.col(1)));

                for (int i = 0; i < vp1s.size(); i++){
                    ppAndFocals[i].first = { pps(i, 0), pps(i, 1) };
                    ppAndFocals[i].second = focalLengths(i);
                }
                return ppAndFocals;
            }

            
            // make infinite points finite
            // merge close points
            void RefineIntersections(std::vector<HPoint2> & intersections, 
                std::vector<std::pair<int, int>> * intersectionMakerLineIds = nullptr, 
                double distanceThreshold = 1.0) {

                SetClock();

                static const bool useVotes = false;

                if (useVotes){

                    THERE_ARE_BUGS_HERE("returns nothing!!!!");

                    static const double fakeFocal = 500;
                    PanoramicCamera cam(fakeFocal);
                    core::Imaged votes = core::Imaged::zeros(cam.screenSize());

                    for (auto & inter : intersections){
                        Vec3 dir(inter.numerator[0], inter.numerator[1], inter.denominator * fakeFocal);
                        for (auto & p : { cam.toScreen(dir), cam.toScreen(-dir) }){
                            int x = p[0], y = p[1];
                            if (x < 0) x = 0;
                            if (x >= votes.cols) x = votes.cols - 1;
                            if (y < 0) y = 0;
                            if (y >= votes.rows) y = votes.rows - 1;
                            votes(y, x) += 1.0;
                        }
                    }

                    cv::GaussianBlur(votes, votes, cv::Size(5, 5), 0.1);
                    std::vector<PixelLoc> points;
                    NonMaximaSuppression(votes, votes, 50, &points);
                    intersections.clear();
                    for (auto & p : points){
                        Vec3 dir = cam.toSpace(p);
                        intersections.emplace_back(Point2(dir[0], dir[1]), dir[2] / fakeFocal);
                    }

                    assert(!intersectionMakerLineIds);
                }
                else{
                    for (auto & hp : intersections) {
                        if (hp.denominator == 0.0)
                            hp.denominator = 1e-5;
                    }

                    std::vector<HPoint2> mergedIntersections;
                    mergedIntersections.reserve(intersections.size());
                    std::vector<std::pair<int, int>> mergedIntersectionMakerLineIds;
                    mergedIntersectionMakerLineIds.reserve(intersections.size());

                    RTreeWrapper<HPoint2, DefaultInfluenceBoxFunctor<double>> rtreeRecorder(DefaultInfluenceBoxFunctor<double>(distanceThreshold * 2.0));
                    for (int i = 0; i < intersections.size(); i++) {
                        if (rtreeRecorder.contains(intersections[i], [&distanceThreshold](const HPoint2 & a, const HPoint2 & b) {
                            return Distance(a.value(), b.value()) < distanceThreshold;
                        })) {
                            continue;
                        }
                        rtreeRecorder.insert(intersections[i]);
                        mergedIntersections.push_back(intersections[i]);
                        if (intersectionMakerLineIds)
                            mergedIntersectionMakerLineIds.push_back((*intersectionMakerLineIds)[i]);
                    }

                    intersections = std::move(mergedIntersections);
                    if (intersectionMakerLineIds)
                        *intersectionMakerLineIds = std::move(mergedIntersectionMakerLineIds);
                }
            }

            std::vector<Vec3> RefineIntersectionsAndProjectToSpace(const std::vector<HPoint2> & intersections, 
                double fakeFocal, double angleThres) {
                std::vector<Vec3> dirs;
                dirs.reserve(intersections.size());

                RTreeWrapper<Vec3, DefaultInfluenceBoxFunctor<double>> rtreeRecorder(DefaultInfluenceBoxFunctor<double>(angleThres * 2.0));
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

            Imaged GetLineLengthRatios(const std::vector<Line2> & lines){
                // get max line length
                double maxLineLen = 0;
                for (auto & line : lines) {
                    if (line.length() > maxLineLen)
                        maxLineLen = line.length();
                }
                Imaged lineLengthRatios(lines.size(), 1);
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
                const Imaged & votesPanel, const Imaged & votesRemainedPanel,
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

            std::tuple<std::vector<HPoint2>, double, std::vector<int>> EstimateVanishingPointsWithProjectionCenterAtOrigin(
                const std::vector<Line2> & lines, double minFocalLength, double maxFocalLength, double maxPrinciplePointOffset,
                bool & succeed) {

                SetClock();

                // get max line length
                auto lineLengthRatios = GetLineLengthRatios(lines);

                std::vector<int> lineClasses(lines.size(), -1);

                // get all intersections
                std::vector<std::pair<int, int>> intersectionMakerLineIds;
                auto intersections = ComputeLineIntersections(lines, &intersectionMakerLineIds, true, std::numeric_limits<double>::max());

                std::cout << "intersection num: " << intersections.size() << std::endl;        
                RefineIntersections(intersections, &intersectionMakerLineIds);
                std::cout << "intersection num: " << intersections.size() << std::endl;

                // nlines x npoints (without consideration of line length ratios)
                Imaged votesPanel = LinesVotesToPoints(intersections, lines);

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
                    }
                    else {
                        remainedLines.push_back(lines[i]);
                    }
                }

                if (remainedLines.empty()){
                    std::cout << "no remained lines to locate other two vps, failed" << std::endl;
                    succeed = false;
                    return std::tuple<std::vector<HPoint2>, double, std::vector<int>>();
                }

                //std::cout << "remained lines num: " << remainedLines.size() << std::endl;

                // get remained intersections
                std::vector<std::pair<int, int>> remainedIntersectionMakerLineIds;
                auto remainedIntersections = ComputeLineIntersections(remainedLines, &remainedIntersectionMakerLineIds, true, std::numeric_limits<double>::max());

                //std::cout << "remained intersection num: " << remainedLines.size() << std::endl;
                RefineIntersections(remainedIntersections, &remainedIntersectionMakerLineIds);
                //std::cout << "remained intersection num: " << remainedLines.size() << std::endl;

                // vote remained lines for remained intersections
                std::vector<double> votesForRemainedIntersections(remainedIntersections.size(), 0.0);
                // remainedlines x remainedpoints
                Imaged votesRemainedPanel = LinesVotesToPoints(remainedIntersections, lines); // all lines participated!!!!!
                for (int i = 0; i < remainedIntersections.size(); i++) {
                    votesForRemainedIntersections[i] = //cv::sum(votesRemainedPanel.col(i)).val[0];
                        votesRemainedPanel.col(i).dot(lineLengthRatios);
                }
                // sort remained intersections by votes
                std::vector<int> orderedRemainedIntersectionIDs(remainedIntersections.size());
                std::iota(orderedRemainedIntersectionIDs.begin(), orderedRemainedIntersectionIDs.end(), 0);
                std::sort(orderedRemainedIntersectionIDs.begin(), orderedRemainedIntersectionIDs.end(), 
                    [&votesForRemainedIntersections](int id1, int id2){
                    return votesForRemainedIntersections[id1] > votesForRemainedIntersections[id2];
                });

                //// traverse other vp pairs
                double curMaxScore = 0.0;
                HPoint2 vp2, vp3;
                double curFocal;
                Point2 curPrinciplePoint;

                Point2 vp1ccenter = vp1.value() / 2.0;
                double vp1cdist = norm(vp1).value();

                int maxNum = 500000;
                int count = 0;

                for (int k = 0; k < orderedRemainedIntersectionIDs.size() * 2; k++) {

                    //Clock clock(std::to_string(k) + "-th guessing of vp triplets");

                    if (count >= maxNum)
                        break;

                    int lb = std::max(0, k - (int)orderedRemainedIntersectionIDs.size());
                    int ub = std::min((int)orderedRemainedIntersectionIDs.size(), k);

                    for (int i = lb; i < ub; i++){
                        int j = k - 1 - i;
                        if (j <= i)
                            continue;

                        assert(IsBetween(i, 0, orderedRemainedIntersectionIDs.size()));
                        assert(IsBetween(j, 0, orderedRemainedIntersectionIDs.size()));

                        if (count >= maxNum)
                            break;

                        auto & vp2cand = remainedIntersections[orderedRemainedIntersectionIDs[i]];
                        auto & vp3cand = remainedIntersections[orderedRemainedIntersectionIDs[j]];

                        if (Distance(vp2cand.value(), vp1.value()) < minFocalLength)
                            continue;
                        if (Distance(vp2cand.value(), vp1ccenter) < vp1cdist / 2.0 - minFocalLength)
                            continue;

                        HPoint2 vp12center = (vp1 + vp2cand) / 2.0;
                        double vp12dist = norm(vp1 - vp2cand).value();

                        if (Distance(vp3cand.value(), vp1.value()) < minFocalLength ||
                            Distance(vp2cand.value(), vp3cand.value()) < minFocalLength)
                            continue;
                        if (Distance(vp3cand.value(), vp12center.value()) < vp12dist / 2.0)
                            continue;


                        double focalLength;
                        Point2 principlePoint;
                        std::tie(principlePoint, focalLength) =
                            ComputePrinciplePointAndFocalLength(vp1.value(), vp2cand.value(), vp3cand.value());

                        if (std::isnan(focalLength) || std::isinf(focalLength))
                            continue;

                        if (norm(principlePoint) < maxPrinciplePointOffset &&
                            IsBetween(focalLength, minFocalLength, maxFocalLength)) {
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

                if (count == 0){
                    std::cout << "no valid vps, failed" << std::endl;
                    succeed = false;
                    return std::tuple<std::vector<HPoint2>, double, std::vector<int>>();
                }

                std::vector<HPoint2> vps = { { vp1, vp2, vp3 } };

                std::cout << "vp2: " << vp2.value() << "   " << "vp3: " << vp3.value() << std::endl;

                lineClasses = ClassifyLines(LinesVotesToPoints(vps, lines), 0.8);
                succeed = true;
                return std::make_tuple(std::move(vps), curFocal, std::move(lineClasses));
            }

        }


        namespace {




            // estimate vanishing points using classified lines
            HPoint2 EstimateVanishingPointsFromLines(const std::vector<Line2> & lines, double * score = nullptr){
                assert(!lines.empty());
                if (lines.size() == 1){
                    if (score){
                        *score = -0.1;
                    }
                    return HPoint2(lines.front().direction(), 0.0);
                }
                if (lines.size() == 2){
                    if (score){
                        *score = 0.0;
                    }
                    auto eq1 = Concat(lines[0].first, 1.0).cross(Concat(lines[0].second, 1.0));
                    auto eq2 = Concat(lines[1].first, 1.0).cross(Concat(lines[1].second, 1.0));
                    return HPointFromVector(eq1.cross(eq2));
                }

                using namespace Eigen;

                //// 0    1   2   3
                //// [e1x e1y e2x e2y]
                Map<const Array<double, 4, Eigen::Dynamic>> lineMat((const double*)lines.data(), 4, lines.size());
                auto e1x = lineMat.row(0).eval();
                auto e1y = lineMat.row(1).eval();
                auto e2x = lineMat.row(2).eval();
                auto e2y = lineMat.row(3).eval();

                int lineNum = lines.size();
                auto costFunction = [lineNum, &e1x, &e1y, &e2x, &e2y](const VectorXd & x, VectorXd & v){
                    // x: input, v: output
                    assert(x.size() == 2);
                    double vx(cos(x(0))*sin(x(1))), vy(sin(x[0])*sin(x[1])), vz(cos(x[1]));
                    //                                                             2 
                    //(e1y vx - e2y vx - e1x vy + e2x vy + e1x e2y vz - e2x e1y vz) 
                    //-------------------------------------------------------------- 
                    //                            2                           2 
                    //    (e1x vz - 2 vx + e2x vz)  + (e1y vz - 2 vy + e2y vz)
                    auto top = abs(e1y * vx - e2y * vx - e1x * vy + e2x * vy + e1x * e2y * vz - e2x * e1y * vz);
                    auto bottom = sqrt((e1x * vz - vx * 2.0 + e2x * vz).square() + (e1y * vz - 2.0 * vy + e2y * vz).square());
                    v = (top / bottom).matrix().transpose();
                };

                auto functor = misc::MakeGenericNumericDiffFunctor<double>(costFunction, 2, lineNum);
                LevenbergMarquardt<decltype(functor)> lm(std::move(functor));

                // calc initial vp
                Vec3 initialVP = Concat(lines[0].first, 1.0).cross(Concat(lines[0].second, 1.0))
                    .cross(Concat(lines[1].first, 1.0).cross(Concat(lines[1].second, 1.0)));
                initialVP /= norm(initialVP);
                VectorXd x(2);
                x << atan2(initialVP[1], initialVP[0]), acos(initialVP[2]);
                lm.minimize(x);

                if (score){
                    VectorXd scores;
                    costFunction(x, scores);
                    *score = scores.squaredNorm();
                }

                double vx(cos(x(0))*sin(x(1))), vy(sin(x[0])*sin(x[1])), vz(cos(x[1]));
                return HPoint2(Point2(vx, vy), vz);
            }

        }


        Failable<std::tuple<std::vector<HPoint2>, double, std::vector<int>>> VanishingPointsDetector::operator() (
            const std::vector<Line2> & lines, const SizeI & imSize) const {
            
            Point2 projCenter(imSize.width / 2.0, imSize.height / 2.0);
            double imScale = sqrt(imSize.area());
            double minFocalLength = imScale * _params.minFocalLengthRatio;
            double maxFocalLength = imScale * _params.maxFocalLengthRatio;
            double maxPrinciplePointOffset = imScale * _params.maxPrinciplePointOffsetRatio;

            if (_params.algorithm == Naive){
                std::vector<Line2> offsetedLines = lines;
                for (auto & line : offsetedLines) {
                    line.first -= projCenter;
                    line.second -= projCenter;
                }
                bool ok = false;
                auto results = EstimateVanishingPointsWithProjectionCenterAtOrigin(offsetedLines, 
                    minFocalLength, maxFocalLength, maxPrinciplePointOffset, ok);
                if (!ok){
                    return nullptr;
                }
                for (HPoint2 & vp : std::get<0>(results)) {
                    vp = vp + HPoint2(projCenter);
                }
                return std::move(results);
            }
            else if (_params.algorithm == MATLAB_PanoContext){
                assert(MatlabEngine::IsBuilt());
                MatlabEngine matlab;
                // install lines
                Imaged linesData(lines.size(), 4);
                for (int i = 0; i < lines.size(); i++){
                    linesData(i, 0) = lines[i].first[0];
                    linesData(i, 1) = lines[i].first[1];
                    linesData(i, 2) = lines[i].second[0];
                    linesData(i, 3) = lines[i].second[1];
                }
                MatlabEngine::PutVariable("linesData", linesData);
                MatlabEngine::PutVariable("projCenter", projCenter);
                // convert to struct array
                matlab << "[vp, f, lineclasses] = panoramix_wrapper_vpdetection(linesData, projCenter');";
                Imaged vpData;
                Imaged focal;
                Imaged lineClassesData;
                MatlabEngine::GetVariable("vp", vpData, false);
                MatlabEngine::GetVariable("f", focal, false);
                MatlabEngine::GetVariable("lineclasses", lineClassesData, false);
                if (!(vpData.cols == 2 && vpData.rows == 3)){
                    return nullptr;
                }
                assert(lineClassesData.cols * lineClassesData.rows == lines.size());
                std::vector<HPoint2> vps = {
                    HPoint2({ vpData(0, 0), vpData(0, 1) }, 1.0),
                    HPoint2({ vpData(1, 0), vpData(1, 1) }, 1.0),
                    HPoint2({ vpData(2, 0), vpData(2, 1) }, 1.0)
                };
                std::vector<int> lineClasses(lines.size(), -1);
                for (int i = 0; i < lines.size(); i++){
                    lineClasses[i] = (int)std::round(lineClassesData(i) - 1);
                }
                return std::make_tuple(std::move(vps), focal(0), std::move(lineClasses));
            }
            else if (_params.algorithm == MATLAB_Tardif){
                assert(MatlabEngine::IsBuilt());
                // install lines
                Imaged linesData(lines.size(), 5);
                for (int i = 0; i < lines.size(); i++){
                    linesData(i, 0) = lines[i].first[0];
                    linesData(i, 1) = lines[i].first[1];
                    linesData(i, 2) = lines[i].second[0];
                    linesData(i, 3) = lines[i].second[1];
                    linesData(i, 4) = 1.0;
                }
                //MatlabEngine::PutVariable("linesData", linesData);
                NOT_IMPLEMENTED_YET();
            }
            else /*if (_params.algorithm == TardifSimplified)*/{

                std::vector<std::vector<Line2>> lineClusters;
                std::vector<double> clusterInitialScores;
                std::vector<int> intLabels;
                int classNum = 0;

                {
                    std::vector< std::vector<float> *> pts(lines.size());
                    for (int i = 0; i < lines.size(); i++){
                        pts[i] = new std::vector<float>{
                            (float)lines[i].first[0], (float)lines[i].first[1], (float)lines[i].second[0], (float)lines[i].second[1]
                        };
                    }

                    std::vector<std::vector<float> *> *mModels =
                        VPSample::run(&pts, 5000, 2, 0, 3);
                    std::vector<unsigned int> labels;
                    std::vector<unsigned int> labelCount;
                    classNum = VPCluster::run(labels, labelCount, &pts, mModels, 2, 2);

                    assert(classNum >= 3);

                    for (unsigned int i = 0; i < mModels->size(); ++i)
                        delete (*mModels)[i];
                    delete mModels;
                    for (auto & p : pts){
                        delete p;
                    }

                    // estimate vps
                    assert(lines.size() == labels.size());
                    lineClusters.resize(classNum);
                    clusterInitialScores.resize(classNum, 0.0);
                    intLabels.resize(labels.size());
                    for (int i = 0; i < labels.size(); i++){
                        lineClusters[labels[i]].push_back(lines[i]);
                        clusterInitialScores[labels[i]] += lines[i].length();
                        intLabels[i] = labels[i];
                    }
                }

                // initial vps
                std::vector<HPoint2> vps(classNum);
                for (int i = 0; i < lineClusters.size(); i++){
                    if (lineClusters[i].empty())
                        continue;
                    vps[i] = EstimateVanishingPointsFromLines(lineClusters[i]);
                }



                // find class with maximum initial score as the first class
                int firstClass = std::distance(clusterInitialScores.begin(), 
                    std::max_element(clusterInitialScores.begin(), clusterInitialScores.end()));


                int finalClasses[3] = { firstClass, -1, -1 };
                double finalFocal = 1.0;
                Point2 finalPP;
                double minViolation = std::numeric_limits<double>::max();

                double logMinFocal = std::log10(minFocalLength);
                double logMaxFocal = std::log10(maxFocalLength);
                double logMiddleFocal = (logMinFocal + logMaxFocal) / 2.0;

                // select the other two classes by traversing
                using namespace Eigen;
                Map<const Array<double, 4, Eigen::Dynamic>> lineMat1((const double*)lineClusters[firstClass].data(), 4, lineClusters[firstClass].size());
                for (int i = 0; i < classNum; i++){
                    if (i == firstClass)
                        continue;
                    Map<const Array<double, 4, Eigen::Dynamic>> lineMat2((const double*)lineClusters[i].data(), 4, lineClusters[i].size());
                    for (int j = i + 1; j < classNum; j++){
                        if (j == firstClass)
                            continue;
                        Map<const Array<double, 4, Eigen::Dynamic>> lineMat3((const double*)lineClusters[j].data(), 4, lineClusters[j].size());

                        int ids[] = { firstClass, i, j };
                        
                        // calc initial vp
                        VectorXd x(6);
                        for (int k = 0; k < 3; k++){
                            auto & lineCluster = lineClusters[ids[k]];
                            Vec3 initialVP = normalize(VectorFromHPoint(vps[ids[k]]));
                            x(k * 2) = atan2(initialVP[1], initialVP[0]);
                            x(k * 2 + 1) = acos(initialVP[2]);
                        }
                        
                        int startPoses[] = { 0, lineMat1.cols(), lineMat1.cols() + lineMat2.cols() };
                        auto costFunction = [&lineMat1, &lineMat2, &lineMat3, &startPoses, logMiddleFocal, &projCenter, &vps, &ids](const VectorXd & x, VectorXd & v){
                            // x: input, v: output
                            assert(x.size() == 6);
                            v.resize(lineMat1.cols() + lineMat2.cols() + lineMat3.cols());

                            Map<const Array<double, 4, Eigen::Dynamic>> * lineMats[] = { &lineMat1, &lineMat2, &lineMat3 };
                            Vec3 currentVPs[3];
                            double maxOutAngle = 0.0;
                            for (int k = 0; k < 3; k++){
                                currentVPs[k] = Vec3(cos(x(k * 2))*sin(x(k * 2 + 1)), sin(x[k * 2])*sin(x[k * 2 + 1]), cos(x[k * 2 + 1]));
                                maxOutAngle = std::max(maxOutAngle, AngleBetweenUndirectedVectors(VectorFromHPoint(vps[ids[k]]), currentVPs[k]));
                            }
                            for (int k = 0; k < 3; k++){
                                auto & lineMat = *lineMats[k];
                                double vx = currentVPs[k][0], vy = currentVPs[k][1], vz = currentVPs[k][2];
                                auto e1x = lineMat.row(0);
                                auto e1y = lineMat.row(1);
                                auto e2x = lineMat.row(2);
                                auto e2y = lineMat.row(3);

                                //                                                             2 
                                //(e1y vx - e2y vx - e1x vy + e2x vy + e1x e2y vz - e2x e1y vz) 
                                //-------------------------------------------------------------- 
                                //                            2                           2 
                                //    (e1x vz - 2 vx + e2x vz)  + (e1y vz - 2 vy + e2y vz)
                                auto top = abs(e1y * vx - e2y * vx - e1x * vy + e2x * vy + e1x * e2y * vz - e2x * e1y * vz);
                                auto bottom = sqrt((e1x * vz - vx * 2.0 + e2x * vz).square() + (e1y * vz - 2.0 * vy + e2y * vz).square());
                                v.middleRows(startPoses[k], lineMats[k]->cols()) = (top / bottom).matrix().transpose();
                            }

                            Point2 pp;
                            double focal;
                            std::tie(pp, focal) = ComputePrinciplePointAndFocalLength(
                                Vec2(currentVPs[0][0], currentVPs[0][1]) / currentVPs[0][2],
                                Vec2(currentVPs[1][0], currentVPs[1][1]) / currentVPs[1][2],
                                Vec2(currentVPs[2][0], currentVPs[2][1]) / currentVPs[2][2]);
                           
                            if (std::isinf(focal) || std::isnan(focal)){
                                //v.bottomRows(1)(0) = 1;
                                v *= 1e6;
                            }
                            else{
                                static const double angleThreshold = M_PI / 20.0;
                                //v.bottomRows(1)(0) = 1;
                                //v.bottomRows(1)(0) = 10 * (Distance(pp, projCenter)) + 20.0 * abs(log10(focal) - logMiddleFocal);
                                /*if (maxOutAngle > angleThreshold){
                                    v *= (maxOutAngle / angleThreshold);
                                }*/
                            }
                        };

                        auto functor = misc::MakeGenericNumericDiffFunctor<double>(costFunction, 6, lineMat1.cols() + lineMat2.cols() + lineMat3.cols());
                        LevenbergMarquardt<decltype(functor)> lm(functor);

                    
                        VectorXd oldx = x;
                        lm.minimize(x);
                        /*if (oldx == x){
                            std::cout << "NO OPTIMIZATION AT ALL !!!!!!!!!!!!!!!!!!!!!!!" << std::endl;
                        }*/

                        // get the optimized vps
                        Vec3 currentVPs[3];
                        for (int k = 0; k < 3; k++){
                            currentVPs[k] = Vec3(cos(x(k * 2))*sin(x(k * 2 + 1)), sin(x[k * 2])*sin(x[k * 2 + 1]), cos(x[k * 2 + 1]));
                        }

                        Point2 pp;
                        double focal;
                        std::tie(pp, focal) = ComputePrinciplePointAndFocalLength(
                            Vec2(currentVPs[0][0], currentVPs[0][1]) / currentVPs[0][2],
                            Vec2(currentVPs[1][0], currentVPs[1][1]) / currentVPs[1][2],
                            Vec2(currentVPs[2][0], currentVPs[2][1]) / currentVPs[2][2]);
                        
                        if (std::isnan(focal) || std::isinf(focal))
                            continue;

                        if (!IsBetween(focal, minFocalLength, maxFocalLength))
                            continue;

                        double ppViolation = Distance(pp, projCenter);

                        auto minmaxScore = std::minmax({ clusterInitialScores[i], clusterInitialScores[j] });
                        double violation = ( ppViolation > maxPrinciplePointOffset ? (80.0 * (ppViolation - maxPrinciplePointOffset) / maxPrinciplePointOffset) : 0.0 )
                            - 100 * (minmaxScore.first + minmaxScore.second) / clusterInitialScores[firstClass];

                        if (violation < minViolation){
                            minViolation = violation;
                            for (int k = 0; k < 3; k++){
                                vps[ids[k]] = HPointFromVector(currentVPs[k]);
                            }
                            finalClasses[1] = i;
                            finalClasses[2] = j;
                            finalFocal = focal;
                            finalPP = pp;
                        }
                    }
                }

                if (finalClasses[1] == -1 || finalClasses[2] == -1){
                    return nullptr;
                }

                for (int k = 0; k < 3; k++){
                    std::swap(vps[k], vps[finalClasses[k]]);
                    std::cout << "vp[" << k << "] = " << vps[k].value() << std::endl;
                }
                std::cout << "focal = " << finalFocal << std::endl;

                for (int & l : intLabels){
                    for (int k = 0; k < 3; k++){
                        if (l == k){
                            l = finalClasses[k];
                            break;
                        }
                        else if (l == finalClasses[k]){
                            l = k;
                            break;
                        }
                    }
                }
                return std::make_tuple(std::move(vps), finalFocal, std::move(intLabels));
            }

        }

       
        Failable<std::tuple<std::vector<HPoint2>, double>> VanishingPointsDetector::operator() (
            std::vector<Classified<Line2>> & lines, const SizeI & imSize) const {
            std::vector<HPoint2> vps;
            double focalLength = 0.0;
            std::vector<int> lineClassies;
            std::vector<Line2> plines(lines.size());
            for (int i = 0; i < lines.size(); i++)
                plines[i] = lines[i].component;
            auto result = (*this)(plines, imSize);
            if (result.null()){
                return nullptr;
            }
            std::tie(vps, focalLength, lineClassies) = result.unwrap();
            assert(lineClassies.size() == lines.size());
            for (int i = 0; i < lines.size(); i++)
                lines[i].claz = lineClassies[i];
            return std::make_tuple(std::move(vps), focalLength);
        }



#pragma endregion VanishingPointsDetector








        std::pair<Failable<double>, Failable<double>> ComputeFocalsFromHomography(const Mat3 & H) {
            /// from OpenCV code

            const double* h = reinterpret_cast<const double*>(H.val);

            double d1, d2; // Denominators
            double v1, v2; // Focal squares value candidates
            double f0 = 0.0, f1 = 0.0;

            std::pair<Failable<double>, Failable<double>> results;

            d1 = h[6] * h[7];
            d2 = (h[7] - h[6]) * (h[7] + h[6]);
            v1 = -(h[0] * h[1] + h[3] * h[4]) / d1;
            v2 = (h[0] * h[0] + h[3] * h[3] - h[1] * h[1] - h[4] * h[4]) / d2;
            if (v1 < v2) std::swap(v1, v2);
            if (v1 > 0 && v2 > 0) 
                results.second = sqrt(std::abs(d1) > std::abs(d2) ? v1 : v2);
            else if (v1 > 0)
                results.second = sqrt(v1);

            d1 = h[0] * h[3] + h[1] * h[4];
            d2 = h[0] * h[0] + h[1] * h[1] - h[3] * h[3] - h[4] * h[4];
            v1 = -h[2] * h[5] / d1;
            v2 = (h[5] * h[5] - h[2] * h[2]) / d2;
            if (v1 < v2) std::swap(v1, v2);
            if (v1 > 0 && v2 > 0) 
                results.first = sqrt(std::abs(d1) > std::abs(d2) ? v1 : v2);
            else if (v1 > 0) 
                results.first = sqrt(v1);

            return results;
        }







#pragma region SegmentationExtractor


        namespace {

            // the edge of graph
            struct Edge {
                float w;
                int a, b;
            };

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

            inline float DefaultThreshold(float size, float c) {
                return size == 0.0 ? 1e8 : c / size;
            }

            // merge similar nodes in graph
            template <class ThresholdFunT = decltype(DefaultThreshold)>
            Universe SegmentGraph(const std::vector<float> & verticesSizes, std::vector<Edge> & edges, float c, 
                ThresholdFunT && thresholdFun = DefaultThreshold) {
                
                std::sort(edges.begin(), edges.end(), [](const Edge & e1, const Edge & e2){
                    return e1.w < e2.w;
                });

                int numVertices = verticesSizes.size();
                Universe u(verticesSizes);
                std::vector<float> threshold(numVertices);

                for (int i = 0; i < numVertices; i++)
                    threshold[i] = thresholdFun(1, c);

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
                            threshold[a] = edge.w + thresholdFun(u.size(a), c);
                        }
                    }
                }

                return u;
            }


            inline float ColorDistance(const Vec3 & a, const Vec3 & b, bool useYUV){
                static const Mat3 RGB2YUV(
                    0.299, 0.587, 0.114,
                    -0.14713, -0.28886, 0.436,
                    0.615, -0.51499, -0.10001
                );
                static const Mat3 BGR2VUY = RGB2YUV.t();
                return useYUV ? norm(BGR2VUY * (a - b)) * 3.0 : norm(a - b);
            }

            // for edge weight computation
            // measure pixel distance
            inline float PixelDiff(const Image & im, const cv::Point & p1, const cv::Point & p2, bool useYUV) {
                assert(im.depth() == CV_8U && im.channels() == 3);
                Vec3 c1 = im.at<cv::Vec<uint8_t, 3>>(p1);
                Vec3 c2 = im.at<cv::Vec<uint8_t, 3>>(p2);
                return ColorDistance(c1, c2, useYUV);
            }

            // measure pixel distance with 2d line cuttings
            inline float PixelDiff(const Image & im, const Imagei & linesOccupation,
                const PixelLoc & p1, const PixelLoc & p2,
                const std::vector<Line2> & lines, bool useYUV){

                assert(im.depth() == CV_8U && im.channels() == 3);
                for (int lineId : { linesOccupation(p1), linesOccupation(p2) }){
                    if (lineId >= 0){
                        auto & line = lines[lineId];
                        double p1OnLeftFlag = (p1 - ToPixelLoc(line.first)).cross(ToPixelLoc(line.direction()));
                        double p2OnLeftFlag = (p2 - ToPixelLoc(line.first)).cross(ToPixelLoc(line.direction()));
                        if (p1OnLeftFlag * p2OnLeftFlag < 0){
                            return 1e5;
                        }
                    }
                }
                Vec3 c1 = im.at<cv::Vec<uint8_t, 3>>(p1);
                Vec3 c2 = im.at<cv::Vec<uint8_t, 3>>(p2);
                return ColorDistance(c1, c2, useYUV);
            }

            inline float PixelDiff(const Image & im, const Imagei & linesOccupation,
                const PixelLoc & p1, const PixelLoc & p2,
                const std::vector<Line3> & lines, const PanoramicCamera & cam, bool useYUV){

                assert(im.depth() == CV_8U && im.channels() == 3);
                auto direction1 = cam.toSpace(p1);
                auto direction2 = cam.toSpace(p2);
                for (int lineId : { linesOccupation(p1), linesOccupation(p2) }){
                    if (lineId >= 0){
                        auto & line = lines[lineId];
                        auto lineNormal = line.first.cross(line.second);
                        double p1OnLeftFlag = lineNormal.dot(direction1);
                        double p2OnLeftFlag = lineNormal.dot(direction2);
                        if (p1OnLeftFlag * p2OnLeftFlag < 0){
                            return 1e5;
                        }
                    }
                }
                Vec3 c1 = im.at<cv::Vec<uint8_t, 3>>(p1);
                Vec3 c2 = im.at<cv::Vec<uint8_t, 3>>(p2);
                return ColorDistance(c1, c2, useYUV);
            }



            template <class PixelDiffFuncT>
            std::vector<Edge> ComposeGraphEdges(int width, int height, bool isPanorama,
                const Image & smoothed, const PixelDiffFuncT & pixelDiff){

                std::vector<Edge> edges;
                edges.reserve(width * height * 4);
                for (int y = 0; y < height; y++) {
                    for (int x = 0; x < width; x++) {
                        if (x < width - 1) {
                            Edge edge;
                            edge.a = y * width + x;
                            edge.b = y * width + ((x + 1) % width);
                            edge.w = pixelDiff(smoothed, { x, y }, { (x + 1) % width, y });
                            edges.push_back(edge);
                        }

                        if (y < height - 1) {
                            Edge edge;
                            edge.a = y * width + x;
                            edge.b = ((y + 1) % height) * width + x;
                            edge.w = pixelDiff(smoothed, { x, y }, { x, (y + 1) % height });
                            edges.push_back(edge);
                        }

                        if ((x < width - 1) && (y < height - 1)) {
                            Edge edge;
                            edge.a = y * width + x;
                            edge.b = ((y + 1) % height) * width + ((x + 1) % width);
                            edge.w = pixelDiff(smoothed, { x, y }, { (x + 1) % width, (y + 1) % height });
                            edges.push_back(edge);
                        }

                        if ((x < width - 1) && (y > 0)) {
                            Edge edge;
                            edge.a = y * width + x;
                            edge.b = ((y + height - 1) % height) * width + ((x + 1) % width);
                            edge.w = pixelDiff(smoothed, { x, y }, { (x + 1) % width, (y + height - 1) % height });
                            edges.push_back(edge);
                        }
                    }
                }
                if (isPanorama){ // collect panorama borders
                    for (int y = 0; y < height; y++){
                        Edge edge;
                        edge.a = y * width + 0;
                        edge.b = y * width + width - 1;
                        edge.w = pixelDiff(smoothed, cv::Point{ 0, y }, cv::Point{ width - 1, y });
                        edges.push_back(edge);
                        if (y < height - 1){
                            edge.b = (y + 1) * width + width - 1;
                            edge.w = pixelDiff(smoothed, cv::Point{ 0, y }, cv::Point{ width - 1, y + 1});
                            edges.push_back(edge);
                        }
                        if (y > 0){
                            edge.b = (y - 1) * width + width - 1;
                            edge.w = pixelDiff(smoothed, cv::Point{ 0, y }, cv::Point{ width - 1, y - 1 });
                            edges.push_back(edge);
                        }
                    }
                    for (int x1 = 0; x1 < width; x1++){
                        for (int x2 = x1; x2 < width; x2++){
                            Edge edge;
                            edge.a = 0 * width + x1;
                            edge.b = 0 * width + x2;
                            edge.w = pixelDiff(smoothed, cv::Point{ x1, 0 }, cv::Point{ x2, 0 });
                            edges.push_back(edge);
                            edge.a = (height - 1) * width + x2;
                            edge.b = (height - 1) * width + x1;
                            edge.w = pixelDiff(smoothed, cv::Point{ x2, height - 1 }, cv::Point{ x1, height - 1 });
                            edges.push_back(edge);
                        }
                    }
                }
                return edges;
            }


            std::vector<float> ComposeGraphVerticesSizes(int width, int height, bool isPanorama){
                if (!isPanorama){
                    return std::vector<float>(width * height, 1.0f);
                }
                std::vector<float> vsizes(width * height);
                float radius = width / 2.0 / M_PI;
                for (int y = 0; y < height; y++){
                    float longitude = (y - height/2.0f) / radius;
                    float scale = cos(longitude);
                    if (scale < 0){
                        scale = 0.0;
                    }
                    std::fill_n(vsizes.data() + y * width, width, scale);
                    /*for (int x = 0; x < width; x++){
                        vsizes[y * width + x] = scale;
                    }*/
                }
                return vsizes;
            }


            std::pair<Imagei, Image> PerformSegmentation(const std::vector<float> & verticesSizes, 
                std::vector<Edge> & edges, int width, int height, 
                float sigma, float c, int minSize, int & numCCs, bool returnColoredResult = false){

                int num = (int)edges.size();
                Universe u = SegmentGraph(verticesSizes, edges, c);

                for (int i = 0; i < num; i++) {
                    int a = u.find(edges[i].a);
                    int b = u.find(edges[i].b);
                    if ((a != b) && ((u.size(a) < minSize) || (u.size(b) < minSize)))
                        u.join(a, b);
                }

                numCCs = u.numSets();
                std::unordered_map<int, int> compIntSet;
                Imagei output(cv::Size2i(width, height));
                for (int y = 0; y < height; y++) {
                    for (int x = 0; x < width; x++) {
                        int comp = u.find(y * width + x);
                        if (compIntSet.find(comp) == compIntSet.end()){
                            compIntSet.insert(std::make_pair(comp, (int)compIntSet.size()));
                        }
                        output(cv::Point(x, y)) = compIntSet[comp];
                    }
                }
                assert(compIntSet.size() == numCCs);

                if (!returnColoredResult){
                    return std::make_pair(output, Image());
                }

                Image coloredOutput(cv::Size2i(width, height), CV_8UC3);
                std::vector<cv::Vec<uint8_t, 3>> colors(numCCs);
                std::generate(colors.begin(), colors.end(), [](){
                    return cv::Vec<uint8_t, 3>(uint8_t(std::rand() % 256),
                        uint8_t(std::rand() % 256),
                        uint8_t(std::rand() % 256));
                });
                for (int y = 0; y < height; y++) {
                    for (int x = 0; x < width; x++) {
                        coloredOutput.at<cv::Vec<uint8_t, 3>>(cv::Point(x, y)) =
                            colors[output(cv::Point(x, y))];
                    }
                }
                return std::make_pair(output, coloredOutput);
            }




            // first return is CV_32SC1, the second is CV_8UC3 (for display)
            std::pair<Imagei, Image> SegmentImage(const Image & im, float sigma, float c, int minSize, bool isPanorama, 
                int & numCCs, bool returnColoredResult = false, bool useYUV = true) {

                assert(im.depth() == CV_8U && im.channels() == 3);

                int width = im.cols;
                int height = im.rows;
                Image smoothed;
                cv::GaussianBlur(im, smoothed, cv::Size(5, 5), sigma);

                // build pixel graph
                std::vector<float> vSizes = ComposeGraphVerticesSizes(width, height, isPanorama);
                //float(*pixelDiff)(const Image & im, const cv::Point & p1, const cv::Point & p2, bool useYUV) = PixelDiff;
                std::vector<Edge> edges = ComposeGraphEdges(width, height, isPanorama, smoothed, 
                    [useYUV](const Image & im, const cv::Point & p1, const cv::Point & p2){
                    return PixelDiff(im, p1, p2, useYUV);
                });

                return PerformSegmentation(vSizes, edges, width, height, sigma, c, minSize, numCCs, returnColoredResult);
            }

            // first return is CV_32SC1, the second is CV_8UC3 (for display)
            std::pair<Imagei, Image> SegmentImage(const Image & im, float sigma, float c, int minSize,
                const std::vector<Line2> & lines,
                int & numCCs, bool returnColoredResult = false, bool useYUV = true){

                assert(im.depth() == CV_8U && im.channels() == 3);

                int width = im.cols;
                int height = im.rows;
                Image smoothed;
                cv::GaussianBlur(im, smoothed, cv::Size(5, 5), sigma);

                Imagei linesOccupation(im.size(), -1);
                for (int i = 0; i < lines.size(); i++){
                    auto & l = lines[i];
                    cv::line(linesOccupation, ToPixelLoc(l.first), ToPixelLoc(l.second), i, 2);
                }

                // build pixel graph
                std::vector<float> vSizes = ComposeGraphVerticesSizes(width, height, false);
                std::vector<Edge> edges = ComposeGraphEdges(width, height, false, smoothed, 
                    [&linesOccupation, &lines, useYUV](const Image & im, const cv::Point & p1, const cv::Point & p2){
                    return PixelDiff(im, linesOccupation, p1, p2, lines, useYUV);
                });

                return PerformSegmentation(vSizes, edges, width, height, sigma, c, minSize, numCCs, returnColoredResult);
            }

            // first return is CV_32SC1, the second is CV_8UC3 (for display)
            std::pair<Imagei, Image> SegmentImage(const Image & im, float sigma, float c, int minSize,
                const std::vector<Line3> & lines, const PanoramicCamera & cam,
                int & numCCs, bool returnColoredResult = false, bool useYUV = true){

                assert(im.depth() == CV_8U && im.channels() == 3);

                int width = im.cols;
                int height = im.rows;
                Image smoothed;
                cv::GaussianBlur(im, smoothed, cv::Size(5, 5), sigma);

                Imagei linesOccupation(im.size(), -1);
                for (int i = 0; i < lines.size(); i++){
                    auto & l = lines[i];
                    double spanAngle = AngleBetweenDirections(l.first, l.second);
                    std::vector<std::vector<PixelLoc>> pline(1);
                    for (double a = 0.0; a <= spanAngle; a += 0.01){
                        auto direction = RotateDirection(l.first, l.second, a);
                        pline.front().push_back(ToPixelLoc(cam.toScreen(direction)));
                    }                    
                    cv::polylines(linesOccupation, pline, false, i, 2);
                }

                // build pixel graph
                std::vector<float> vSizes = ComposeGraphVerticesSizes(width, height, true);
                std::vector<Edge> edges = ComposeGraphEdges(width, height, true, smoothed,
                    [&linesOccupation, &lines, &cam, useYUV](const Image & im, const cv::Point & p1, const cv::Point & p2){
                    return PixelDiff(im, linesOccupation, p1, p2, lines, cam, useYUV);
                });

                return PerformSegmentation(vSizes, edges, width, height, sigma, c, minSize, numCCs, returnColoredResult);
            }





            std::pair<Imagei, int> SegmentImageUsingSLIC(const Image & im, int spsize, int spnum){
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
                Imagei labels = Imagei::zeros(im.size());
                for (int x = 0; x < im.cols; x++){
                    for (int y = 0; y < im.rows; y++){
                        labels(PixelLoc(x, y)) = klabels[y * im.cols + x];
                    }
                }

                delete[] ubuff;
                delete[] klabels;

                return std::make_pair(labels, nlabels);
            }


            inline void ToQuickShiftImage(const Image & IMG, image_t & im){
                im.N1 = IMG.rows;
                im.N2 = IMG.cols;
                im.K = IMG.channels();
                assert(im.K == 3);
                im.I = (float *)calloc(im.N1*im.N2*im.K, sizeof(float));
                for (int k = 0; k < im.K; k++)
                for (int col = 0; col < im.N2; col++)
                for (int row = 0; row < im.N1; row++)
                {
                    auto & pt = IMG.at<cv::Vec<uint8_t, 3>>(/*im.N1 - 1 - row*/row, col);
                    im.I[row + col*im.N1 + k*im.N1*im.N2] = 32. * pt[k] / 255.; // Scale 0-32
                }
            }

            int * map_to_flatmap(float * map, unsigned int size) {
                /********** Flatmap **********/
                int *flatmap = (int *)malloc(size*sizeof(int));
                for (unsigned int p = 0; p < size; p++)
                {
                    flatmap[p] = map[p];
                }

                bool changed = true;
                while (changed)
                {
                    changed = false;
                    for (unsigned int p = 0; p < size; p++)
                    {
                        changed = changed || (flatmap[p] != flatmap[flatmap[p]]);
                        flatmap[p] = flatmap[flatmap[p]];
                    }
                }

                /* Consistency check */
                for (unsigned int p = 0; p < size; p++)
                    assert(flatmap[p] == flatmap[flatmap[p]]);

                return flatmap;
            }

            image_t imseg(image_t im, int * flatmap) {
                /********** Mean Color **********/
                float * meancolor = (float *)calloc(im.N1*im.N2*im.K, sizeof(float));
                float * counts = (float *)calloc(im.N1*im.N2, sizeof(float));

                for (int p = 0; p < im.N1*im.N2; p++)
                {
                    counts[flatmap[p]]++;
                    for (int k = 0; k < im.K; k++)
                        meancolor[flatmap[p] + k*im.N1*im.N2] += im.I[p + k*im.N1*im.N2];
                }

                int roots = 0;
                for (int p = 0; p < im.N1*im.N2; p++)
                {
                    if (flatmap[p] == p)
                        roots++;
                }
                printf("Roots: %d\n", roots);

                int nonzero = 0;
                for (int p = 0; p < im.N1*im.N2; p++)
                {
                    if (counts[p] > 0)
                    {
                        nonzero++;
                        for (int k = 0; k < im.K; k++)
                            meancolor[p + k*im.N1*im.N2] /= counts[p];
                    }
                }
                if (roots != nonzero)
                    printf("Nonzero: %d\n", nonzero);
                assert(roots == nonzero);


                /********** Create output image **********/
                image_t imout = im;
                imout.I = (float *)calloc(im.N1*im.N2*im.K, sizeof(float));
                for (int p = 0; p < im.N1*im.N2; p++)
                for (int k = 0; k < im.K; k++)
                    imout.I[p + k*im.N1*im.N2] = meancolor[flatmap[p] + k*im.N1*im.N2];

                free(meancolor);
                free(counts);

                return imout;
            }


            std::pair<Imagei, int> SegmentImageUsingQuickShiftCPU(const Image & originalIm, float sigma, float tau) {
                image_t im;
                ToQuickShiftImage(originalIm, im);
                float *map, *E, *gaps;
                int * flatmap;
                image_t imout;

                map = (float *)calloc(im.N1*im.N2, sizeof(float));
                gaps = (float *)calloc(im.N1*im.N2, sizeof(float));
                E = (float *)calloc(im.N1*im.N2, sizeof(float));

                quickshift(im, sigma, tau, map, gaps, E);

                /* Consistency check */
                for (int p = 0; p < im.N1*im.N2; p++)
                if (map[p] == p) assert(gaps[p] == INF);

                flatmap = map_to_flatmap(map, im.N1*im.N2);
                imout = imseg(im, flatmap);
                Imagei segmented(im.N1, im.N2);
                for (int col = 0; col < im.N2; col++)
                for (int row = 0; row < im.N1; row++)
                {
                    segmented(row, col) = flatmap[row + col*im.N1];
                }
                int segnum = *std::max_element(flatmap, flatmap + im.N1 * im.N2) + 1;

                free(im.I);
                free(imout.I);
                free(map);
                free(gaps);
                free(E);
                free(flatmap);
                return std::make_pair(segmented, segnum);
            }

            std::pair<Imagei, int> SegmentImageUsingQuickShiftGPU(const Image & originalIm, float sigma, float tau) {
                image_t im;
                ToQuickShiftImage(originalIm, im);
                float *map, *E, *gaps;
                int * flatmap;
                image_t imout;

                map = (float *)calloc(im.N1*im.N2, sizeof(float));
                gaps = (float *)calloc(im.N1*im.N2, sizeof(float));
                E = (float *)calloc(im.N1*im.N2, sizeof(float));

                quickshift_gpu(im, sigma, tau, map, gaps, E);

                /* Consistency check */
                for (int p = 0; p < im.N1*im.N2; p++)
                if (map[p] == p) assert(gaps[p] == INF);

                flatmap = map_to_flatmap(map, im.N1*im.N2);
                imout = imseg(im, flatmap);
                Imagei segmented(im.N1, im.N2);
                for (int col = 0; col < im.N2; col++)
                for (int row = 0; row < im.N1; row++)
                {
                    segmented(row, col) = flatmap[row + col*im.N1];
                }
                int segnum = *std::max_element(flatmap, flatmap + im.N1 * im.N2) + 1;

                free(im.I);
                free(imout.I);
                free(map);
                free(gaps);
                free(E);
                free(flatmap);
                return std::make_pair(segmented, segnum);
            }

        }


        std::pair<Imagei, int> SegmentationExtractor::operator() (const Image & im, bool isPanorama) const {
            if (_params.algorithm == SLIC){
                assert(!isPanorama);
                return SegmentImageUsingSLIC(im, _params.superpixelSizeSuggestion, _params.superpixelNumberSuggestion);
            }
            else if (_params.algorithm == GraphCut){                
                int numCCs;
                Imagei segim = SegmentImage(im, _params.sigma, _params.c, _params.minSize, isPanorama, numCCs, false, _params.useYUVColorSpace).first;
                return std::make_pair(segim, numCCs);
            }
            else if (_params.algorithm == QuickShiftCPU){
                assert(!isPanorama);
                return SegmentImageUsingQuickShiftCPU(im, 6, 10);
            }
            else if (_params.algorithm == QuickShiftGPU){
                assert(!isPanorama);
                return SegmentImageUsingQuickShiftGPU(im, 6, 10);
            }
            else{
                SHOULD_NEVER_BE_CALLED();
            }
        }

        std::pair<Imagei, int>  SegmentationExtractor::operator() (const Image & im, const std::vector<Line2> & lines, double extensionLength) const {
            assert(_params.algorithm == GraphCut);
            int numCCs;
            if (extensionLength == 0.0){
                Imagei segim = SegmentImage(im, _params.sigma, _params.c, _params.minSize, lines, numCCs, false, _params.useYUVColorSpace).first;
                return std::make_pair(segim, numCCs);
            }
            else{
                auto extLines = lines;
                for (auto & line : extLines){
                    auto d = normalize(line.direction());
                    line.first -= (d * extensionLength);
                    line.second += (d * extensionLength);
                }
                Imagei segim = SegmentImage(im, _params.sigma, _params.c, _params.minSize, extLines, numCCs, false, _params.useYUVColorSpace).first;
                return std::make_pair(segim, numCCs);
            }
        }

        std::pair<Imagei, int> SegmentationExtractor::operator() (const Image & im, const std::vector<Line3> & lines, const PanoramicCamera & cam, double extensionAngle) const {
            assert(_params.algorithm == GraphCut);
            int numCCs;
            if (extensionAngle == 0.0){
                Imagei segim = SegmentImage(im, _params.sigma, _params.c, _params.minSize, lines, cam, numCCs, false, _params.useYUVColorSpace).first;
                return std::make_pair(segim, numCCs);
            }
            else{
                auto extLines = lines;
                for (auto & line : extLines){
                    auto p1 = RotateDirection(line.first, line.second, -extensionAngle);
                    auto p2 = RotateDirection(line.second, line.first, -extensionAngle);
                    line.first = p1;
                    line.second = p2;
                }
                Imagei segim = SegmentImage(im, _params.sigma, _params.c, _params.minSize, extLines, cam, numCCs, false, _params.useYUVColorSpace).first;
                return std::make_pair(segim, numCCs);
            }
        }

#pragma endregion SegmentationExtractor


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

        }

        bool IsDenseSegmentation(const Imagei & segRegions){
            int minv, maxv;
            std::tie(minv, maxv) = MinMaxValOfImage(segRegions);
            if (minv != 0) return false;
            std::vector<bool> stored(maxv + 1, false);
            for (int i : segRegions){
                stored[i] = true;
            }
            for (auto b : stored){
                if (!b)
                    return false;
            }
            return true;
        }

        std::map<std::pair<int, int>, std::vector<std::vector<PixelLoc>>> FindContoursOfRegionsAndBoundaries(
            const Imagei & segRegions, int connectionExtendSize, bool simplifyStraightEdgePixels){

            std::map<std::pair<int, int>, std::vector<std::vector<PixelLoc>>> boundaryEdges;
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
                            if (simplifyStraightEdgePixels){
                                bool closed = Distance(edges.back().front(), edges.back().back()) <= 1.5;
                                cv::approxPolyDP(edges.back(), edges.back(), 2, closed);
                            }
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
       
            return boundaryEdges;
        }


        std::map<std::set<int>, std::vector<PixelLoc>> ExtractBoundaryJunctions(const Imagei & regions){
            std::map<std::set<int>, std::vector<PixelLoc>> junctions;
            for (auto it = regions.begin(); it != regions.end(); ++it){
                auto p = it.pos();
                if (p.x == regions.cols - 1)
                    continue;
                if (p.y == regions.rows - 1)
                    continue;
                std::set<int> regionIds = {
                    regions(p), regions(PixelLoc(p.x + 1, p.y)),
                    regions(PixelLoc(p.x, p.y + 1)), regions(PixelLoc(p.x + 1, p.y + 1))
                };
                if (regionIds.size() >= 3){
                    junctions[regionIds].push_back(p);
                }
            }
            return junctions;
        }













        std::tuple<std::vector<Point2>, std::vector<Point2>, std::vector<Point2>, std::vector<int>, int>  sample_line(const std::vector<Line2> &lines, const std::vector<int> &lineClasses){


            int sample_rate = 10; //sample every 5 pixel on line
            std::vector<Point2> ls[3];
            int count = 0;
            std::vector<int> lsclass;

            for (int i = 0; i < lines.size(); i++)
            {
                int n_sample = ceil(norm(lines[i].first - lines[i].second) / sample_rate);
                int lclass = lineClasses[i];

                /*ls(i).sample = [...
                linspace(lines(i).point1(1), lines(i).point2(1), n_sample)' ...
                linspace(lines(i).point1(2), lines(i).point2(2), n_sample)' ];*/
                Point2 eps;
                eps[0] = (lines[i].first - lines[i].second).val[0] / n_sample;
                eps[1] = (lines[i].first - lines[i].second).val[1] / n_sample;
                Point2 curr = lines[i].second;

                if (lclass == -1)
                {
                    continue;
                }


                for (int j = 0; j < n_sample; j++)
                {
                    ls[lclass].push_back(curr);
                    lsclass.push_back(lclass);
                    count++;
                    curr = curr + eps;
                }

                /*ls(i).lineclass = repmat(lines(i).lineclass, n_sample, 1);*/
            }


            return std::make_tuple(ls[0], ls[1], ls[2], lsclass, count);
        }


        Point2 move_line_towards_vp(Point2 curp1, Point2 curp2, Point2 vp, int amount, bool &atvp, Point2 &newp2){
            Point2 vec1 = curp1 - vp;
            Point2 vec2 = curp2 - vp;
            double len1 = sqrt(vec1[0] * vec1[0] + vec1[1] * vec1[1]);
            double len2 = sqrt(vec2[0] * vec2[0] + vec2[1] * vec2[1]);
            Point2 norm1;
            norm1[0] = vec1[0] / len1*amount;
            norm1[1] = vec1[1] / len1*amount;
            Point2 newp1;
            Point2 norm2;
            double ratio21 = len2 / len1;
            norm2[0] = vec2[0] / len2  * amount* ratio21;
            norm2[1] = vec2[1] / len2 * amount * ratio21;



            if (len1 < abs(amount)){
                newp1 = curp1;
                newp2 = curp2;
                atvp = 1;
            }
            else
            {

                newp1 = curp1 + norm1;
                newp2 = curp2 + norm2;
                atvp = 0;
            }

            return newp1;
        };


        std::vector<cv::Point2f> getpoly(const Line2 & line,
            const Point2 & vanishingPoints, int to_away, const cv::Size & imageSize, std::vector<Point2> sample){
            cv::Point2f p1;
            p1.x = line.first[0];
            p1.y = line.first[1];
            cv::Point2f p2;
            p2.x = line.second[0];
            p2.y = line.second[1];



            cv::Point2f curp1 = p1;
            cv::Point2f curp2 = p2;
            int moveAmount = 64;



            std::vector<cv::Point2f> result;


            //	Imagei oo = Imagei::ones(imageSize) * -1;;

            while (moveAmount >= 1)
            {
                bool atvp = 0;
                Point2 hcurp1;
                Point2 hcurp2;

                hcurp1[0] = curp1.x;
                hcurp1[1] = curp1.y;
                hcurp2[0] = curp2.x;
                hcurp2[1] = curp2.y;


                Point2 hnewp2;
                Point2 hnewp1 = move_line_towards_vp(hcurp1, hcurp2, vanishingPoints, to_away * moveAmount, atvp, hnewp2);
                cv::Point2f newp1;
                cv::Point2f newp2;
                newp1.x = hnewp1[0];
                newp1.y = hnewp1[1];
                newp2.x = hnewp2[0];
                newp2.y = hnewp2[1];




                bool failed = 0;
                if (atvp == 1)
                {

                    failed = 1;
                }
                else if ((hnewp1[0] > imageSize.width || hnewp1[0]<1 || hnewp1[1]>imageSize.height || hnewp1[1]<1) || (hnewp2[0]>imageSize.width || hnewp2[0]<1 || hnewp2[1]>imageSize.height || hnewp2[1] < 1))
                {

                    failed = 1;
                }
                else
                {
                    std::vector<cv::Point2f> poly;
                    poly.push_back(p1);
                    poly.push_back(p2);
                    cv::Point2f pt;
                    for (int i = 0; i < sample.size(); i++)
                    {
                        pt.x = (sample[i])[0];
                        pt.y = (sample[i])[1];
                        poly.push_back(newp2);
                        poly.push_back(newp1);
                        poly.push_back(p1);
                        double isstop = cv::pointPolygonTest(poly, pt, false);

                        if (isstop>0){
                            failed = 1;
                            break;
                        }
                    }
                }

                if (failed == 1)
                {

                    moveAmount = moveAmount / 2;
                }
                else{
                    curp1 = newp1;
                    curp2 = newp2;
                }
            }
            //cv::Point** pts;
            //int* npts;
            //int ncontours;
            //cv::Mat IMG = cv::Mat::zeros(imageSize, CV_32SC1);
            //
            //cv::Point offset = cv::Point();
            ///*cv::fillPoly(IMG, pts, npts, ncontours, 225, 8, 0, offset);
            //cv::fillPoly()*/

            result.push_back(p1);
            result.push_back(p2);

            result.push_back(curp2);
            result.push_back(curp1);

            return result;

        }


        void myfillPoly(cv::Mat& img, std::vector<std::vector<  Point<int32_t, 2>  > > & pts, int label)
        {
            for (auto & ps : pts)
            {
                cv::fillConvexPoly(img, ps, label);
            }
        }


        Imagei ExtractOrientationMaps(
            const cv::Size & imageSize,
            const std::vector<Line2> & lines,
            const std::array<HPoint2, 3> & vanishingPoints,
            const std::vector<int> & lineClasses) {

            // fill with -1 as default
            Imagei omap = Imagei::ones(imageSize) * -1;
            std::vector<std::vector<  Point<int32_t, 2>  > > lineextimg[3][3][2];

            std::vector<Point2> ls[3];
            std::vector<int> lsclass;
            int lsSize;

            std::tie(ls[0], ls[1], ls[2], lsclass, lsSize) = sample_line(lines, lineClasses);

            size_t lnum = lines.size();
            //size_t lnum = 40;
            for (int i = 0; i < lnum; i++) {
                int lclass = lineClasses[i];
                if (lclass == -1)
                    continue;
                for (int j = 0; j < 3; j++)
                {
                    if (j == lclass)
                    {
                        ;
                    }
                    else{
                        Point2 Vp = vanishingPoints[j].value();
                        //Vp[0];
                        int orient;
                        for (int k = 0; k < 3; k++)
                        {
                            if (k != lclass && k != j)
                                orient = k;
                        }

                        std::vector<cv::Point2f> poly1 = getpoly(lines[i], Vp, -1, imageSize, ls[orient]);
                        std::vector<Point<int32_t, 2>> ppoly1;
                        std::vector<Point<int32_t, 2>> ppoly2;
                        for (int ii = 0; ii < poly1.size(); ii++)
                        {
                            Point<int32_t, 2> temppoint;
                            temppoint[0] = poly1.at(ii).x;
                            temppoint[1] = poly1.at(ii).y;
                            ppoly1.push_back(temppoint);
                        }
                        lineextimg[lclass][j][0].push_back(ppoly1);
                        std::vector<cv::Point2f> poly2 = getpoly(lines[i], Vp, 1, imageSize, ls[orient]);
                        for (int ii = 0; ii < poly1.size(); ii++)
                        {
                            Point<int32_t, 2> temppoint;
                            temppoint[0] = poly2.at(ii).x;
                            temppoint[1] = poly2.at(ii).y;
                            ppoly2.push_back(temppoint);
                        }
                        lineextimg[lclass][j][1].push_back(ppoly2);
                        //std::cout << "i =  " << i << "," << j << std::endl;
                        /*if (i == 4 && j == 0)
                        int ass = 0;*/


                    }
                }


            }


            Imagei omap1 = Imagei::ones(imageSize) * -1;
            Imagei omap2 = Imagei::ones(imageSize) * -1;
            Imagei omap3 = Imagei::ones(imageSize) * -1;
            Imagei omap12 = Imagei::ones(imageSize) * -1;
            Imagei omap21 = Imagei::ones(imageSize) * -1;
            Imagei omap10 = Imagei::ones(imageSize) * -1;
            Imagei omap01 = Imagei::ones(imageSize) * -1;
            Imagei omap02 = Imagei::ones(imageSize) * -1;
            Imagei omap20 = Imagei::ones(imageSize) * -1;
            myfillPoly(omap12, lineextimg[1][2][0], 0);
            myfillPoly(omap12, lineextimg[1][2][1], 0);
            myfillPoly(omap21, lineextimg[2][1][0], 0);
            myfillPoly(omap21, lineextimg[2][1][1], 0);
            myfillPoly(omap10, lineextimg[1][0][0], 2);
            myfillPoly(omap10, lineextimg[1][0][1], 2);
            myfillPoly(omap01, lineextimg[0][1][0], 2);
            myfillPoly(omap01, lineextimg[0][1][1], 2);
            myfillPoly(omap20, lineextimg[2][0][0], 1);
            myfillPoly(omap20, lineextimg[2][0][1], 1);
            myfillPoly(omap02, lineextimg[0][2][0], 1);
            myfillPoly(omap02, lineextimg[0][2][1], 1);

            for (int j = 0; j < imageSize.height; j++){
                for (int i = 0; i < imageSize.width; i++)
                {
                    if (omap12[j][i] == 0 && omap21[j][i] == 0)
                    {
                        omap1[j][i] = 0;
                    }
                }
            }
            for (int j = 0; j < imageSize.height; j++){
                for (int i = 0; i < imageSize.width; i++)
                {
                    if (omap02[j][i] == 1 && omap20[j][i] == 1)
                    {
                        omap3[j][i] = 1;
                    }
                }
            }
            for (int j = 0; j < imageSize.height; j++){
                for (int i = 0; i < imageSize.width; i++)
                {
                    if (omap10[j][i] == 2 && omap01[j][i] == 2)
                    {
                        omap2[j][i] = 2;
                    }
                }
            }





            Image3d omapImage1 = Image3d::zeros(omap.size());
            Image3d omapImage2 = Image3d::zeros(omap.size());
            Image3d omapImage3 = Image3d::zeros(omap.size());

            for (int x = 0; x < omap.cols; x++) {
                for (int y = 0; y < omap.rows; y++) {
                    Vec3 color(0, 0, 0);
                    int omapValue = omap1(y, x);
                    if (omapValue >= 0 && omapValue < 3)
                        color[omapValue] = 255;
                    omapImage1(y, x) = color;
                }
            }
            for (int x = 0; x < omap.cols; x++) {
                for (int y = 0; y < omap.rows; y++) {
                    Vec3 color(0, 0, 0);
                    int omapValue = omap2(y, x);
                    if (omapValue >= 0 && omapValue < 3)
                        color[omapValue] = 255;
                    omapImage2(y, x) = color;
                }
            }
            for (int x = 0; x < omap.cols; x++) {
                for (int y = 0; y < omap.rows; y++) {
                    Vec3 color(0, 0, 0);
                    int omapValue = omap3(y, x);
                    if (omapValue >= 0 && omapValue < 3)
                        color[omapValue] = 255;
                    omapImage3(y, x) = color;
                }
            }
            //cv::imshow("Orientation Map1", omapImage1);
            //cv::imshow("Orientation Map2", omapImage2);
            //cv::imshow("Orientation Map3", omapImage3);
            //cv::waitKey();
            for (int j = 0; j < imageSize.height; j++){
                for (int i = 0; i < imageSize.width; i++)
                {
                    if (omap1[j][i] == 0 && omap2[j][i] == -1 && omap3[j][i] == -1)
                    {
                        omap[j][i] = 0;
                    }
                    if (omap1[j][i] == -1 && omap2[j][i] == 2 && omap3[j][i] == -1)
                    {
                        omap[j][i] = 2;
                    }
                    if (omap1[j][i] == -1 && omap2[j][i] == -1 && omap3[j][i] == 1)
                    {
                        omap[j][i] = 1;
                    }

                }

            }

            //Imagei ao23 = (lineextimg[1][2][0] | lineextimg[1][2][1]);
            //Imagei ao32 = (lineextimg[2][1][0] | lineextimg[2][1][1]);
            //Imagei ao13 = (lineextimg[0][2][0] | lineextimg[0][2][1]);
            //Imagei ao31 = (lineextimg[2][0][0] | lineextimg[2][0][1]);
            //Imagei ao12 = (lineextimg[0][1][0] | lineextimg[0][1][1]);
            //Imagei ao21 = (lineextimg[1][0][0] | lineextimg[1][0][1]);
            //Imagei aa23 = (lineextimg[1][2][0] & lineextimg[1][2][1]);
            //Imagei aa32 = (lineextimg[2][1][0] & lineextimg[2][1][1]);
            //Imagei aa13 = (lineextimg[0][2][0] & lineextimg[0][2][1]);
            //Imagei aa31 = (lineextimg[2][0][0] & lineextimg[2][0][1]);
            //Imagei aa12 = (lineextimg[0][1][0] & lineextimg[0][1][1]);
            //Imagei aa21 = (lineextimg[1][0][0] & lineextimg[1][0][1]);
            //
            //
            //
            //Imagei Aa1  = ao23 & ao32;
            //Imagei Aa2 = ao13 & ao31;
            //Imagei Aa3  = ao12 & ao21;

            //Imagei b1  = Aa1 &~Aa2 &~Aa3;
            //Imagei b2 = ~Aa1 &Aa2 &~Aa3;
            //Imagei b3 = ~Aa1 &~Aa2 &Aa3;

            ////omap = cat(3, b{ 1 }, b{ 2 }, b{ 3 });


            //Imagei testmap = Imagei::ones(imageSize) * -1;
            //std::vector<std::vector<Point<int32_t, 2>>> polys = { 
            //	{
            //	{ 0, 0 },
            //	{ 20, 100 },
            //	{ 50, 100 },
            //	{ 150, 20 }
            //},
            //{
            //	{ 0, 15 },
            //	{ 20, 100 },
            //	{ 50, 100 },
            //	{ 150, 20 }
            //} };
            //// fill with 1
            //myfillPoly(testmap, polys, 1);
            //Image3d testmapImage = Image3d::zeros(omap.size());
            //for (int x = 0; x < omap.cols; x++) {
            //	for (int y = 0; y < omap.rows; y++) {
            //		Vec3 color(0, 0, 0);
            //		int omapValue = testmap(y, x);
            //		if (omapValue >= 0 && omapValue < 3)
            //			color[omapValue] = 255;
            //		testmapImage(y, x) = color;
            //	}
            //}
            //cv::imshow("Orientation Map test", testmapImage);




            return omap;
        }





        Imagei ComputeOrientationMaps(const std::vector<Classified<Line2>> & lines,
            const std::vector<HPoint2> & vps, const SizeI & imSize){

            std::array<HPoint2, 3> vanishingPoints = { { vps[0], vps[1], vps[2] } };
            std::vector<Line2> ls; ls.reserve(lines.size());
            std::vector<int> lcs; lcs.reserve(lines.size());
            for (int i = 0; i < lines.size(); i++){
                if (lines[i].claz == -1 || lines[i].claz >= 3)
                    continue;
                ls.push_back(lines[i].component);
                lcs.push_back(lines[i].claz);
            }
            return ExtractOrientationMaps(imSize, ls, vanishingPoints, lcs);

        }







        ImageOfType<Vec<double, 7>> ComputeGeometricContext(const Image & im, SceneClass sceneClass, bool useHedauForIndoor){
            MatlabEngine matlab;
            matlab.PutVariable("im", im);
            if (sceneClass == SceneClass::Indoor && useHedauForIndoor){
                matlab << "[~, ~, slabelConfMap] = calcGC(im);";
            }
            else{
                matlab << (std::string("slabelConfMap = panoramix_wrapper_gc(im, ") + (sceneClass == SceneClass::Outdoor ? "true" : "false") + ");");
            }
            ImageOfType<Vec<double, 7>> gc;
            matlab.GetVariable("slabelConfMap", gc);
            assert(gc.channels() == 7);
            assert(gc.size() == im.size());
            return gc;
        }

        ImageOfType<Vec<double, 5>> GeometricContextEstimator::operator()(const ImageOfType<Vec<double, 7>> & rawgc,
            SceneClass sceneClass) const{
            ImageOfType<Vec<double, 5>> result(rawgc.size(), Vec<double, 5>());
            for (auto it = result.begin(); it != result.end(); ++it){
                auto & p = rawgc(it.pos());
                auto & resultv = *it;
                if (sceneClass == SceneClass::Indoor){
                    // 0: front, 1: left, 2: right, 3: floor, 4: ceiling, 5: clutter, 6: unknown
                    /*resultv[II_FrontVerticalPlanarFace] += p[0];
                    resultv[II_SideVerticalPlanarFace] += p[1];
                    resultv[II_SideVerticalPlanarFace] += p[2];
                    resultv[II_HorizontalPlanarFace] += p[3];
                    resultv[II_HorizontalPlanarFace] += p[4];
                    resultv[II_Clutter] += p[5];
                    resultv[II_Other] += p[6];*/
                    resultv[II_HorizontalPlanarFace] += p[3];
                    resultv[II_HorizontalPlanarFace] += p[4];
                    resultv[II_VerticalPlanarFace] += p[0];
                    resultv[II_VerticalPlanarFace] += p[1];
                    resultv[II_VerticalPlanarFace] += p[2];
                    resultv[II_Clutter] += p[5];
                    resultv[II_Other] += p[6];
                }
                else{
                    // 0: ground, 1,2,3: vertical, 4:clutter, 5:poros, 6: sky
                    resultv[OI_Ground] += p[0];
                    resultv[OI_VerticalPlanarFace] += p[1];
                    resultv[OI_VerticalPlanarFace] += p[2];
                    resultv[OI_VerticalPlanarFace] += p[3];
                    resultv[OI_Clutter] += p[4];
                    resultv[OI_Porous] += p[5];
                    resultv[OI_Sky] += p[6];
                }
            }
            return result;
        }

        ImageOfType<Vec<double, 5>> GeometricContextEstimator::operator() (const Image & im, SceneClass sceneClass) const{
            if (im.channels() == 7){
                return (*this)((const ImageOfType<Vec<double, 7>>&)im, sceneClass);
            }
            ImageOfType<Vec<double, 7>> gc = ComputeGeometricContext(im, sceneClass);
            return (*this)(gc, sceneClass);
        }

        std::pair<ImageOfType<Vec<double, 5>>, Imagei> GeometricContextEstimator::operator() (const Image & image,
            const PanoramicCamera & camera, const std::vector<Vec3> & allvps, SceneClass sceneClass) const {

            NOT_IMPLEMENTED_YET(); 

            ImageOfType<Vec<double, 5>> result = ImageOfType<Vec<double, 7>>::zeros(image.size());
            Imagei votes = Imagei::zeros(image.size());

            assert(allvps.size() >= 3);
            std::vector<Vec3> vps = { allvps[0], allvps[1], allvps[2] };

            int vertVPid = NearestDirectionId(vps, camera.up());
            std::vector<PerspectiveCamera> hcams;
            for (int i = -1; i <= 1; i++){
                std::vector<Vec3> vvps = vps;
                for (auto & d : vvps){
                    d += normalize(camera.up()) * i * 0.5;
                }
                auto cs = CreateHorizontalPerspectiveCameras(camera, vvps, 600, 400, 200);
                hcams.insert(hcams.end(), cs.begin(), cs.end());
            }
            // 0: front, 1: left, 2: right, 3: floor, 4: ceiling, 5: clutter, 6: unknown
            for (int i = 0; i < hcams.size(); i++){
                auto gc = ComputeGeometricContext(PanoramicView{ image, camera }.sampled(hcams[i]).image, sceneClass);
                int frontVPid = NearestDirectionId(vps, hcams[i].forward());
                int sideVPid = NearestDirectionId(vps, hcams[i].leftward());

                for (auto it = result.begin(); it != result.end(); ++it){
                    auto gcPos = ToPixelLoc(hcams[i].toScreen(camera.toSpace(it.pos())));
                    if (!Contains(gc, gcPos))
                        continue;
                    auto & resultv = *it;
                    auto & p = gc(gcPos);

                    if (sceneClass == SceneClass::Indoor){
                       /* resultv[frontVPid] += p[0];
                        resultv[sideVPid] += p[1];
                        resultv[sideVPid] += p[2];
                        resultv[vertVPid] += p[3];
                        resultv[vertVPid] += p[4];
                        resultv[3] += p[5];
                        resultv[4] += p[6];*/
                        resultv[II_VerticalPlanarFace] += p[0];
                        resultv[II_VerticalPlanarFace] += p[1];
                        resultv[II_VerticalPlanarFace] += p[2];
                        resultv[II_HorizontalPlanarFace] += p[3];
                        resultv[II_HorizontalPlanarFace] += p[4];
                        resultv[II_Clutter] += p[5];
                        resultv[II_Other] += p[6];
                    }
                    else{
                        // 0: ground, 1,2,3: vertical, 4:clutter, 5:poros, 6: sky
                        resultv[OI_Ground] += p[0];
                        resultv[OI_VerticalPlanarFace] += p[1];
                        resultv[OI_VerticalPlanarFace] += p[2];
                        resultv[OI_VerticalPlanarFace] += p[3];
                        resultv[OI_Clutter] += p[4];
                        resultv[OI_Porous] += p[5];
                        resultv[OI_Sky] += p[6];
                    }
                    votes(it.pos())++;
                }
            }

            return std::make_pair(result, votes);
        }


    }
}

    
