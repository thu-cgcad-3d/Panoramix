#include <iostream>
#include <unordered_map>

#include <Eigen/Dense>

#include "misc.hpp"
#include "feature.hpp"
#include "utilities.hpp"

namespace panoramix {
    namespace core {

        PerspectiveCamera::PerspectiveCamera(int w, int h, double focal, const Vec3 & eye,
            const Vec3 & center, const Vec3 & up, double near, double far) 
        : _screenW(w), _screenH(h), _focal(focal), _eye(eye), _center(center), _up(up), _near(near), _far(far) {
            updateMatrices();
        }

        void PerspectiveCamera::updateMatrices() {
            _viewMatrix = MakeMat4LookAt(_eye, _center, _up);

            double verticalViewAngle = atan(_screenH / 2.0 / _focal) * 2;
            double aspect = double(_screenW) / double(_screenH);
            _projectionMatrix = MakeMat4Perspective(verticalViewAngle, aspect, _near, _far);

            _viewProjectionMatrix = _projectionMatrix * _viewMatrix;
            _viewProjectionMatrixInv = _viewProjectionMatrix.inv();
        }

        Vec2 PerspectiveCamera::screenProjection(const Vec3 & p3) const {
            Vec4 p4(p3(0), p3(1), p3(2), 1);
            Vec4 position = _viewProjectionMatrix * p4;
            double xratio = position(0) / position(3) / 2;
            double yratio = position(1) / position(3) / 2;
            double x = (xratio + 0.5) * _screenW;
            double y = _screenH - (yratio + 0.5) * _screenH;
            return Vec2(x, y);
        }

        bool PerspectiveCamera::isVisibleOnScreen(const Vec3 & p3d) const {
            Vec4 p4(p3d(0), p3d(1), p3d(2), 1);
            Vec4 position = _viewProjectionMatrix * p4;
            return position(3) > 0 && position(2) > 0;
        }

        HPoint2 PerspectiveCamera::screenProjectionInHPoint(const Vec3 & p3) const {
            Vec4 p4(p3(0), p3(1), p3(2), 1);
            Vec4 position = _viewProjectionMatrix * p4;
            double xratio = position(0) / 2;
            double yratio = position(1) / 2;
            double zratio = position(3);

            double x = (xratio + 0.5 * zratio) * _screenW;
            double y = _screenH * zratio - (yratio + 0.5 * zratio) * _screenH;
            return HPoint2({x, y}, zratio);
        }

        Vec3 PerspectiveCamera::spatialDirection(const Vec2 & p2d) const {
            double xratio = (p2d(0) / _screenW - 0.5) * 2;
            double yratio = ((_screenH - p2d(1)) / _screenH - 0.5) * 2;
            Vec4 position(xratio, yratio, 1, 1);
            Vec4 realPosition = _viewProjectionMatrixInv * position;
            return Vec3(realPosition(0) / realPosition(3), 
                realPosition(1) / realPosition(3), 
                realPosition(2) / realPosition(3));
        }

        void PerspectiveCamera::resizeScreen(const Size & sz, bool updateMat) {
            if (_screenH == sz.height && _screenW == sz.width)
                return;
            _screenH = sz.height;
            _screenW = sz.width;
            if (updateMat)
                updateMatrices();
        }

        void PerspectiveCamera::setFocal(double f, bool updateMat) {
            if (f == _focal)
                return;
            _focal = f;
            if (updateMat)
                updateMatrices();
        }

        void PerspectiveCamera::setEye(const Vec3 & e, bool updateMat) {
            if (_eye == e)
                return;
            _eye = e;
            if (updateMat)
                updateMatrices();
        }

        void PerspectiveCamera::setCenter(const Vec3 & c, bool updateMat) {
            if (_center == c)
                return;
            _center = c;
            if (updateMat)
                updateMatrices();
        }

        void PerspectiveCamera::setUp(const Vec3 & up, bool updateMat) {
            if (_up == up)
                return;
            _up = up;
            if (updateMat)
                updateMatrices();
        }

        void PerspectiveCamera::setNearAndFarPlanes(double near, double far, bool updateMat) {
            if (_near == near && _far == far)
                return;
            _near = near;
            _far = far;
            if (updateMat)
                updateMatrices();
        }



        PanoramicCamera::PanoramicCamera(double focal, const Vec3 & eye,
            const Vec3 & center, const Vec3 & up)
            : _focal(focal), _eye(eye), _center(center), _up(up) {
            _xaxis = (_center - _eye); _xaxis /= core::norm(_xaxis);
            _yaxis = _up.cross(_xaxis); _yaxis /= core::norm(_yaxis);
            _zaxis = _xaxis.cross(_yaxis);
        }

        Vec2 PanoramicCamera::screenProjection(const Vec3 & p3) const {
            double xx = p3.dot(_xaxis);
            double yy = p3.dot(_yaxis);
            double zz = p3.dot(_zaxis);
            GeoCoord pg = core::Vec3(xx, yy, zz);
            auto sz = screenSize();
            double x = (pg.longitude + M_PI) / 2.0 / M_PI * sz.width;
            double y = (pg.latitude + M_PI_2) / M_PI * sz.height;
            return Vec2(x, y);
        }

        Vec3 PanoramicCamera::spatialDirection(const Vec2 & p2d) const {
            auto sz = screenSize();
            double longi = p2d(0) / double(sz.width) * 2 * M_PI - M_PI;
            double lati = p2d(1) / double(sz.height) * M_PI - M_PI_2;
            Vec3 dd = (GeoCoord(longi, lati).toVector());
            return dd(0) * _xaxis + dd(1) * _yaxis + dd(2) * _zaxis;
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
                        if (confidence >= 400){
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

           
            //void ExtractLinesExperimental(const cv::Mat& im, std::vector<Line2> & lines,
            //    int minlen, int xborderw, int yborderw, int numDir) {

            //    std::cout << "image processing..." << std::endl;

            //    cv::Mat gim;
            //    cv::cvtColor(im, gim, CV_BGR2GRAY);
            //    int h = gim.rows;
            //    int w = gim.cols;

            //    cv::Mat dx, dy;
            //    cv::Mat ggim;
            //    cv::GaussianBlur(gim, ggim, cv::Size(7, 7), 1.5);
            //    cv::Sobel(ggim, dx, CV_64F, 1, 0);
            //    cv::Sobel(ggim, dy, CV_64F, 0, 1);

            //    cv::Mat imCanny;
            //    cv::Canny(gim, imCanny, 5, 20);

            //    std::cout << "gradient binning..." << std::endl;

            //    cv::Mat imBinIds(im.size(), CV_32SC1); // int32_t
            //    std::set<int> edgePixels; // store all edge pixels in a set
            //    for (int x = 0; x < imCanny.cols; x++){
            //        for (int y = 0; y < imCanny.rows; y++){
            //            if (imCanny.at<uchar>(y, x) > 0){
            //                double a = atan(dy.at<double>(y, x) / dx.at<double>(y, x));
            //                if (isnan(a)){ // NaN
            //                    continue;
            //                }
            //                // compute bin id
            //                int binId = int((a / M_PI + 0.5) * numDir);
            //                imBinIds.at<int32_t>(y, x) = WrapBetween(binId, 0, numDir);
            //                edgePixels.insert(y * w + x);
            //            }
            //            else{
            //                imBinIds.at<int32_t>(y, x) = -1;
            //            }
            //        }
            //    }

            //    std::cout << "collecting pixels.." << std::endl;

            //    cv::Mat imComponentIds(im.size(), CV_32SC1); // component ids int32_t
            //    cv::Mat imParentsInTraversal(im.size(), CV_32SC2); // 
            //    std::map<int, bool> componentUsed;
            //    std::fill(imComponentIds.begin<int32_t>(), imComponentIds.end<int32_t>(), -1); // init to -1

            //    std::vector<int> xs, ys;
            //    xs.reserve(2048);
            //    ys.reserve(2048);

            //    int componentIdGenerator = 0;
            //    while (true){
            //        bool allPixelsAreUsed = true;
            //        PixelLoc unusedPixel;
            //        for (int index : edgePixels){
            //            int y = index / w;
            //            int x = index - y * w;
            //            unusedPixel = {x, y};
            //            if (imBinIds.at<int32_t>(unusedPixel) == -1) // this is not an edge pixel
            //                continue;
            //            int32_t componentId = imComponentIds.at<int32_t>(unusedPixel);
            //            if (componentUsed.find(componentId)!=componentUsed.end() && componentUsed[componentId]) 
            //                // this component is already used to compose a line
            //                continue;
            //            allPixelsAreUsed = false;
            //            break;
            //        }
            //        if (allPixelsAreUsed)
            //            break;

            //        const PixelLoc & fistPixel = unusedPixel;
            //        int activeComponentId = componentIdGenerator++;
            //        int activeBinId = imBinIds.at<int32_t>(fistPixel);
            //        int activeSimiliarBinIds[] = { WrapBetween(activeBinId + 1, 0, numDir), WrapBetween(activeBinId - 1, 0, numDir) };
            //        int pixelsNum = 0;
            //        xs.clear(); 
            //        ys.clear();

            //        static const int xdirs[] = { 1, 1, 0, -1, -1, -1, 0, 1 };
            //        static const int ydirs[] = { 0, 1, 1, 1, 0, -1, -1, -1 };
            //        PixelLoc detector = fistPixel;
            //        PixelLoc detectorPrev = fistPixel;
            //        while (true){ // depth first search
            //            imComponentIds.at<int32_t>(detector) = activeComponentId;
            //            xs.push_back(detector.x);
            //            ys.push_back(detector.y);
            //            pixelsNum++;
            //            
            //            bool allConnectedPixelsAreFound = false;
            //            while (true){
            //                bool noWayOut = true;
            //                for (int k = 0; k < 8; k++){
            //                    PixelLoc next = detector;
            //                    next.x += xdirs[k];
            //                    next.y += ydirs[k];
            //                    if (!IsBetween(next.x, 0, w - 1) || !IsBetween(next.y, 0, h - 1)) // outsided
            //                        continue;

            //                    int32_t binId = imBinIds.at<int32_t>(next);
            //                    if (binId == -1)
            //                        continue;
            //                    if (binId != activeBinId && binId != activeSimiliarBinIds[0] && binId != activeSimiliarBinIds[1])
            //                        // bin ids are not even similiar
            //                        continue;
            //                    int32_t componentId = imComponentIds.at<int32_t>(next);
            //                    if (componentId == activeComponentId) // already visited
            //                        continue;
            //                    if (componentUsed.find(componentId) != componentUsed.end() && componentUsed[componentId])
            //                        // this component is already used to compose a line
            //                        continue;

            //                    detectorPrev = detector;
            //                    detector = next;
            //                    imParentsInTraversal.at<Vec<int32_t, 2>>(detector) = { detectorPrev.x, detectorPrev.y };
            //                    noWayOut = false;
            //                    break;
            //                }
            //                if (!noWayOut){
            //                    break; // a way is found
            //                }else if (detector == fistPixel){ // no way at root pixel
            //                    allConnectedPixelsAreFound = true;
            //                    break;
            //                }else{ // not root pixel yet, traverse back
            //                    auto parent = imParentsInTraversal.at<Vec<int32_t, 2>>(detector);
            //                    detector = { parent[0], parent[1] };
            //                }
            //            }
            //            if (allConnectedPixelsAreFound)
            //                break;
            //        }

            //        std::cout << "pixels num: " << pixelsNum << std::endl;
            //        if (pixelsNum < minlen){
            //            // declare that this component is abandoned
            //            componentUsed[activeComponentId] = true;
            //            continue;
            //        }

            //        // check confidence
            //        cv::Mat xsmat(xs), ysmat(ys);
            //        double meanx = cv::mean(xsmat).val[0], meany = cv::mean(ysmat).val[0];
            //        cv::Mat zmx = xsmat - meanx, zmy = ysmat - meany;

            //        cv::Mat v, lambda;
            //        Mat<double, 2, 2> D;
            //        D(0, 0) = cv::sum(zmx.mul(zmx)).val[0];
            //        D(0, 1) = D(1, 0) = cv::sum(zmx.mul(zmy)).val[0];
            //        D(1, 1) = cv::sum(zmy.mul(zmy)).val[0];
            //        cv::eigen(D, true, lambda, v);

            //        double theta = atan2(v.at<double>(0, 1), v.at<double>(0, 0));
            //        double confidence = std::numeric_limits<double>::max();
            //        if (lambda.at<double>(1) > 0){
            //            confidence = lambda.at<double>(0) / lambda.at<double>(1);
            //        }

            //        if (confidence < 400){
            //            continue;
            //        }

            //        // declare that this component is used
            //        componentUsed[activeComponentId] = true;

            //        // install line
            //        auto xends = std::minmax_element(xs.begin(), xs.end());
            //        auto yends = std::minmax_element(ys.begin(), ys.end());
            //        double minx = *xends.first, maxx = *xends.second;
            //        double miny = *yends.first, maxy = *yends.second;

            //        // ignore boundary line
            //        if (maxx <= xborderw || minx >= w - xborderw || maxy <= yborderw || miny >= h - yborderw)
            //            continue;

            //        double len = sqrt((maxx - minx)*(maxx - minx) + (maxy - miny)*(maxy - miny));
            //        double x1 = meanx - cos(theta) * len / 2;
            //        double x2 = meanx + cos(theta) * len / 2;
            //        double y1 = meany - sin(theta) * len / 2;
            //        double y2 = meany + sin(theta) * len / 2;

            //        lines.push_back({ Vec2(x1, y1), Vec2(x2, y2) });
            //    }



                /*
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
                        if (confidence >= 400){
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
                */
                /*std::cout << "done" << std::endl;

            }*/

        }


        LineSegmentExtractor::Feature LineSegmentExtractor::operator() (const Image & im) const {
            LineSegmentExtractor::Feature lines;
            lines.reserve(1000);
            if (_params.useExperimentalAlgorithm) {
                /*ExtractLinesExperimental(im, lines, _params.minLength,
                    _params.xBorderWidth, _params.yBorderWidth, _params.numDirs);*/
            } else
                ExtractLines(im, lines, _params.minLength, 
                    _params.xBorderWidth, _params.yBorderWidth, _params.numDirs);
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

            std::string ImageDepth2Str(int depth)
            {
                switch (depth)
                {
                case CV_8U: return "CV_8U";
                case CV_8S: return "CV_8S";
                case CV_16U: return "CV_16U";
                case CV_16S: return "CV_16S";
                case CV_32S: return "CV_32S";
                case CV_32F: return "CV_32F";
                case CV_64F: return "CV_64F";
                default:
                    return "unknown depth type";
                }
            }

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
        }


        SegmentationExtractor::Feature SegmentationExtractor::operator() (const Image & im, 
            bool forVisualization) const {
            int numCCs;
            return forVisualization ? SegmentImage(im, _params.sigma, _params.c, _params.minSize, numCCs, true).second 
                : SegmentImage(im, _params.sigma, _params.c, _params.minSize, numCCs, false).first;
        }






        namespace {

            struct LineVPScoreFunctor {
                inline LineVPScoreFunctor(double angleThres = M_PI_4, double s = 0.1) : angleThreshold(angleThres), sigma(s) {}
                inline double operator()(double angle, double lengthRatio) const {
                    return exp(-(angle / angleThreshold) * (angle / angleThreshold) / sigma / sigma / 2) +
                        lengthRatio;
                }
                double angleThreshold, sigma;
            };

            template <class HPoint2ContainerT, class LineVPScoreFunctorT = LineVPScoreFunctor>
            ImageWithType<double> LinesVotesToPoints(const HPoint2ContainerT & points, const std::vector<Line2> & lines,
                LineVPScoreFunctorT && scoreFun = LineVPScoreFunctorT()) {
                size_t nlines = lines.size();
                size_t npoints = points.size();
                ImageWithType<double> votes = ImageWithType<double>::zeros(nlines, npoints);
                double maxLineLen = 0;
                for (auto & line : lines) {
                    if (line.length() > maxLineLen)
                        maxLineLen = line.length();
                }
                for (size_t i = 0; i < nlines; i++) {
                    auto & line = lines[i];
                    for (int j = 0; j < npoints; j++) {
                        const HPoint2 & point = points[j];
                        Vec2 mid2vp = (point - HPoint2(line.center())).numerator;
                        double angle = std::min(AngleBetweenDirections(mid2vp, line.direction()),
                            AngleBetweenDirections(-mid2vp, line.direction()));
                        double score = scoreFun(angle, line.length() / maxLineLen);
                        if (std::isinf(score) || std::isnan(score)) {
                            //std::cout << "!" << std::endl;
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


        }


        std::tuple<std::array<HPoint2, 3>, double, std::vector<int>> VanishingPointsDetector::estimateWithProjectionCenterAtOrigin (
            const std::vector<Line2> & lines) const {

            double scoreThreshold = 0.8;
            std::vector<int> lineClasses(lines.size(), -1);

            // get all intersections
            std::vector<std::pair<int, int>> intersectionMakerLineIds;
            auto intersections = ComputeLineIntersections(lines, &intersectionMakerLineIds, true);

            //std::cout << "intersection num: " << intersections.size() << std::endl;        
            RefineIntersections(intersections, intersectionMakerLineIds);
            //std::cout << "intersection num: " << intersections.size() << std::endl;

            // vote all lines for all intersections
            std::vector<double> votesForIntersections(intersections.size(), 0.0);
            // nlines x npoints
            ImageWithType<double> votesPanel = LinesVotesToPoints(intersections, lines);
            for (int i = 0; i < intersections.size(); i++) {
                votesForIntersections[i] = cv::sum(votesPanel.col(i)).val[0];
            }

            // get vp1
            auto intersectionIdWithMaxVotes = std::distance(votesForIntersections.begin(),
                std::max_element(votesForIntersections.begin(), votesForIntersections.end()));
            const HPoint2 & vp1 = intersections[intersectionIdWithMaxVotes];
            //std::cout << "vp1: " << vp1.value() << std::endl;

            // classify lines for vp1 to be 0, and collect remained lines
            std::vector<Line2> remainedLines;
            remainedLines.reserve(lines.size() / 2);
            for (int i = 0; i < lines.size(); i++) {
                if (votesPanel(i, intersectionIdWithMaxVotes) > 0.5) {
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
            lineClasses = ClassifyLines(LinesVotesToPoints(vps, lines), scoreThreshold);
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
                    cv::Point  ijmax;
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



        

        LocalManhattanVanishingPointsDetector::Result LocalManhattanVanishingPointsDetector::estimateWithProjectionCenterAtOrigin(
            const std::vector<Line2> & lines) const {

            Box2 bbox = BoundingBoxOfContainer(lines);
            double scale = norm(bbox.minCorner, bbox.maxCorner);

            Result result;
            result.lineClasses.resize(lines.size(), -1);

            // get all intersections
            std::vector<std::pair<int, int>> intersectionMakerLineIds;
            auto intersections = ComputeLineIntersections(lines, &intersectionMakerLineIds, true);

            //std::cout << "intersection num: " << intersections.size() << std::endl;        
            RefineIntersections(intersections, intersectionMakerLineIds);
            //std::cout << "intersection num: " << intersections.size() << std::endl;

            // vote all lines for all intersections
            std::vector<std::pair<int, double>> intersectionIdsWithVotes(intersections.size());
            
            // nlines x npoints
            ImageWithType<double> votesPanel = LinesVotesToPoints(intersections, lines);
            for (int i = 0; i < intersections.size(); i++) {
                intersectionIdsWithVotes[i].first = i;
                intersectionIdsWithVotes[i].second = cv::sum(votesPanel.col(i)).val[0];
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
                double angle = atan2(abs(direction.numerator[0]), abs(direction.numerator[1]));
                if (angle < _params.verticalVPAngleRange && 
                    norm(direction.numerator) > abs(direction.denominator) * scale * _params.verticalVPMinDistanceRatioToCenter){
                    verticalIntersectionIdWithMaxVotes = idWithVotes;
                    break;
                }
            }

            assert(verticalIntersectionIdWithMaxVotes.first != -1 && 
                "failed to find vertical vp! "
                "try 1) increasing verticalVPAngleRange or 2) decreasing verticalVPMinDistanceRatioToCenter");

            const HPoint2 & vp1 = intersections[verticalIntersectionIdWithMaxVotes.first];
            std::cout << "vp1: " << vp1.value() << std::endl;

            // classify lines for vp1 to be 0, and collect remained lines
            std::vector<Line2> remainedLines;
            remainedLines.reserve(lines.size() / 2);
            for (int i = 0; i < lines.size(); i++) {
                if (votesPanel(i, verticalIntersectionIdWithMaxVotes.first) > 1.2) {
                    result.lineClasses[i] = 0;
                } else {
                    remainedLines.push_back(lines[i]);
                }
            }

            
            // get remained intersections
            std::vector<std::pair<int, int>> remainedIntersectionMakerLineIds;
            auto remainedIntersections = ComputeLineIntersections(remainedLines, &remainedIntersectionMakerLineIds, true);

            //std::cout << "remained intersection num: " << remainedLines.size() << std::endl;
            RefineIntersections(remainedIntersections, remainedIntersectionMakerLineIds, 20.0);
            //std::cout << "remained intersection num: " << remainedLines.size() << std::endl;

            // vote remained lines for remained intersections
            std::vector<std::pair<int, double>> remainedIntersectionIdsWithVotes(remainedIntersections.size());
            // nlines x npoints
            ImageWithType<double> votesRemainedPanel = LinesVotesToPoints(remainedIntersections, lines); // all lines participated!!!!!
            for (int i = 0; i < remainedIntersections.size(); i++) {
                remainedIntersectionIdsWithVotes[i].first = i;
                remainedIntersectionIdsWithVotes[i].second = cv::sum(votesRemainedPanel.col(i)).val[0];
            }

            // sort remained intersection ids
            std::sort(remainedIntersectionIdsWithVotes.begin(), remainedIntersectionIdsWithVotes.end(),
                [](const std::pair<int, double> & a, const std::pair<int, double> & b) -> bool {
                return a.second > b.second;
            });

            // max test count
            int maxTestCount = 5000;

            //// traverse other vp pairs
            struct FocalAndPPData {
                inline FocalAndPPData() {}
                inline FocalAndPPData(double focal, const Point2 & pp, double s) 
                    : focalLength(focal), principlePoint(pp), score(s) {
                    representative = Concat(focal / 3.0, pp);
                }
                Point3 representative;
                double focalLength;
                Point2 principlePoint;
                std::vector<std::pair<int, int>> remainedIntersectionIdPairs;
                double score;
            };
            std::vector<FocalAndPPData> focalAndPPTable;
            focalAndPPTable.reserve(std::min<unsigned long>(
                remainedIntersectionIdsWithVotes.size() * (remainedIntersectionIdsWithVotes.size() - 1), 
                maxTestCount
            ) / 2);
            auto getBBoxOffocalAndPPData = [&focalAndPPTable](int id) -> Box3 {
                return BoundingBox(focalAndPPTable[id].representative);
            };
            RTreeWrapper<int, decltype(getBBoxOffocalAndPPData)> focalAndPPRTree(getBBoxOffocalAndPPData);


            Point2 vp1ccenter = vp1.value() / 2.0;
            double vp1cdist = norm(vp1).value();			
            int testedCount = 0;

            for (int i = 0; i < remainedIntersectionIdsWithVotes.size(); i++) {
                if (testedCount >= maxTestCount)
                    break;

                auto & vp2cand = remainedIntersections[remainedIntersectionIdsWithVotes[i].first];

                if (Distance(vp2cand, vp1) < _params.minFocalLength)
                    continue;
                if (Distance(vp2cand.value(), vp1ccenter) < vp1cdist / 2.0 - _params.minFocalLength)
                    continue;

                Point2 vp12center = (vp1 + vp2cand).value() / 2.0;
                double vp12dist = norm(vp1 - vp2cand).value();

                for (int j = i + 1; j < remainedIntersectionIdsWithVotes.size(); j++) {
                    if (testedCount >= maxTestCount)
                        break;

                    auto & vp3cand = remainedIntersections[remainedIntersectionIdsWithVotes[j].first];

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
                        testedCount++;

                        // get votes from each line
                        double score = 0.0;
                        std::vector<int> vp2lineIds, vp3lineIds;
                        vp2lineIds.reserve(lines.size() / 2);
                        vp3lineIds.reserve(lines.size() / 2);
                        for (int lineId = 0; lineId < lines.size(); lineId++) {
                            double voteForVP1 = votesPanel(lineId, verticalIntersectionIdWithMaxVotes.first);
                            double voteForVP2 = votesRemainedPanel(lineId, i);
                            double voteForVP3 = votesRemainedPanel(lineId, j);
                            score += std::max({ voteForVP1, voteForVP2, voteForVP3 });
                            if (voteForVP2 > 0.8 && voteForVP2 > voteForVP1 && voteForVP2 > voteForVP3){
                                vp2lineIds.push_back(lineId);
                            }
                            else if (voteForVP3 > 0.8 && voteForVP3 > voteForVP1 && voteForVP3 > voteForVP2){
                                vp3lineIds.push_back(lineId);
                            }
                        }

                        // check closeness between vp2lines and vp3lines
                        int closeCount = 0;
                        for (int line2id : vp2lineIds){
                            auto & line2 = lines[line2id];
                            for (auto & line3id : vp3lineIds){
                                auto & line3 = lines[line3id];
                                if (DistanceBetweenTwoLines(line2, line3).first < 10){
                                    closeCount++;
                                }
                            }
                        }

                        score *= Square(double(closeCount) / std::max(vp2lineIds.size(), vp2lineIds.size()));


                        FocalAndPPData currentFocalAndPP(focalLength, principlePoint, score);
                        Box3 bboxOfFocalAndPP = BoundingBox(currentFocalAndPP.representative);
                        static const Vec3 influence(10, 10, 10);
                        bboxOfFocalAndPP.minCorner -= influence;
                        bboxOfFocalAndPP.maxCorner += influence;

                        int idOfExistedSimiliarFocalAndPP = -1;
                        double minDistance = 10;
                        focalAndPPRTree.search(bboxOfFocalAndPP, 
                            [&idOfExistedSimiliarFocalAndPP, &focalAndPPTable, &minDistance, &currentFocalAndPP](int existedId) -> bool{
                            auto & existed = focalAndPPTable[existedId];
                            double distance = Distance(existed.representative, currentFocalAndPP.representative);
                            if (distance < minDistance){
                                minDistance = distance;
                                idOfExistedSimiliarFocalAndPP = existedId;
                            }
                            return true;
                        });

                        if (idOfExistedSimiliarFocalAndPP == -1){ // no similiar focal and pp
                            focalAndPPTable.push_back(currentFocalAndPP);
                            focalAndPPRTree.insert(focalAndPPTable.size() - 1);
                        }
                        else { // found similiar focal and pp
                            focalAndPPTable[idOfExistedSimiliarFocalAndPP].remainedIntersectionIdPairs
                                .emplace_back(remainedIntersectionIdsWithVotes[i].first, remainedIntersectionIdsWithVotes[j].first);
                            focalAndPPTable[idOfExistedSimiliarFocalAndPP].score += score;
                        }
                    }
                }
            }

            assert(!focalAndPPTable.empty() && "failed to find valid vps!");

            size_t bestFocalAndPPId = std::distance(focalAndPPTable.begin(),
                std::max_element(focalAndPPTable.begin(), focalAndPPTable.end(), 
                [](const FocalAndPPData & a, const FocalAndPPData & b){
                return a.score * b.remainedIntersectionIdPairs.size() < 
                    b.score * a.remainedIntersectionIdPairs.size();
            }));

            auto & bestFocalAndPP = focalAndPPTable[bestFocalAndPPId];

            result.vanishingPoints.clear();
            result.vanishingPoints.reserve(bestFocalAndPP.remainedIntersectionIdPairs.size() * 2 + 1);

            result.vanishingPoints.push_back(vp1);
            result.verticalVanishingPointId = 0;

            for (auto & intersectionIdPair : bestFocalAndPP.remainedIntersectionIdPairs){
                result.vanishingPoints.push_back(remainedIntersections[intersectionIdPair.first]);
                int id1 = result.vanishingPoints.size() - 1;
                result.vanishingPoints.push_back(remainedIntersections[intersectionIdPair.second]);
                int id2 = result.vanishingPoints.size() - 1;
                result.horizontalVanishingPointIds.emplace_back(id1, id2);

                std::cout << "vp pair: {" 
                    << result.vanishingPoints[result.vanishingPoints.size() - 2].value() << ", "
                    << result.vanishingPoints[result.vanishingPoints.size() - 1].value() << "}" 
                    << std::endl;
            }


            
            double scoreThreshold = 0.8;
            result.lineClasses = ClassifyLines(LinesVotesToPoints(result.vanishingPoints, lines), scoreThreshold);
            result.focalLength = bestFocalAndPP.focalLength;

            return result;
        }


        LocalManhattanVanishingPointsDetector::Result LocalManhattanVanishingPointsDetector::estimateWithProjectionCenterAtOriginII(
            const std::vector<Line2> & lines) const {

            std::cout << "line num: " << lines.size() << std::endl;

            Box2 bbox = BoundingBoxOfContainer(lines);
            double scale = norm(bbox.minCorner, bbox.maxCorner);

            Result result;
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
                intersectionIdsWithVotes[i].second = cv::sum(votesPanel.col(i)).val[0];
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
                double angle = atan2(abs(direction.numerator[0]), abs(direction.numerator[1]));
                if (angle < _params.verticalVPAngleRange &&
                    norm(direction.numerator) > abs(direction.denominator) * scale * _params.verticalVPMinDistanceRatioToCenter){
                    verticalIntersectionIdWithMaxVotes = idWithVotes;
                    break;
                }
            }

            assert(verticalIntersectionIdWithMaxVotes.first != -1 &&
                "failed to find vertical vp! "
                "try "
                "1) increasing verticalVPAngleRange or "
                "2) decreasing verticalVPMinDistanceRatioToCenter");

            const HPoint2 & vp1 = intersections[verticalIntersectionIdWithMaxVotes.first];
            std::cout << "vp1: " << vp1.value() << std::endl;

            // classify lines for vp1 to be 0, and collect remained lines
            static const double thresholdVotesForVP1 = 0.8;
            std::vector<Line2> remainedLines;
            remainedLines.reserve(lines.size() / 2);
            for (int i = 0; i < lines.size(); i++) {
                if (votesPanel(i, verticalIntersectionIdWithMaxVotes.first) > thresholdVotesForVP1) {
                    result.lineClasses[i] = 0;
                }
                else {
                    remainedLines.push_back(lines[i]);
                }
            }

            std::cout << "remained lines num: " << remainedLines.size() << std::endl;

            static const double closeLinePairDistanceThreshold = 10;
            std::set<std::pair<int, int>> remainedCloseLineIdPairs;
            for (int i = 0; i < lines.size(); i++){
                for (int j = i + 1; j < lines.size(); j++){
                    if (DistanceBetweenTwoLines(lines[i], lines[j]).first < closeLinePairDistanceThreshold){
                        remainedCloseLineIdPairs.emplace(i, j);
                    }
                }
            }
            std::cout << "close line pair num: " << remainedCloseLineIdPairs.size() << std::endl;


            const double fakeFocal = scale * 5;
            const double angleThres = M_PI_4 / 10.0;

            // get remained intersections
            auto remainedIntersections = ComputeLineIntersections(remainedLines, nullptr, true, 
                closeLinePairDistanceThreshold);

            std::cout << "remained intersection num: " << remainedIntersections.size() << std::endl;
            auto remainedIntersectionVs = RefineIntersectionsAndProjectToSpace(remainedIntersections, 
                fakeFocal, angleThres);
            std::cout << "remained intersection vector num after refinement: " 
                << remainedIntersectionVs.size() 
                << std::endl;

            // project onto one panel
            int longitudeDivideNum = 1000, latitudeDivideNum = 500;
            ImageWithType<double> remainedVotesPanel = 
                ImageWithType<double>::zeros(longitudeDivideNum, latitudeDivideNum);
            for (const Vec3& p : remainedIntersectionVs){
                PixelLoc pixel = PixelLocFromGeoCoord(GeoCoord(p), longitudeDivideNum, latitudeDivideNum);
                remainedVotesPanel(pixel.x, pixel.y) += 1.0;
                pixel = PixelLocFromGeoCoord(GeoCoord(-p), longitudeDivideNum, latitudeDivideNum);
                remainedVotesPanel(pixel.x, pixel.y) += 1.0;
            }
            cv::GaussianBlur(remainedVotesPanel, remainedVotesPanel, 
                cv::Size(5, 5),
                4, 4, cv::BORDER_REPLICATE);
            std::cout << "done voting" << std::endl;

            ImageWithType<double> remainedVotesMaxs;
            std::vector<PixelLoc> remainedIntersectionVReps;
            NonMaximaSuppression(remainedVotesPanel, remainedVotesMaxs, 8, &remainedIntersectionVReps, remainedVotesPanel >= 1);

            std::cout << "remained intersection distribution max num: " << remainedIntersectionVReps.size() << std::endl;

            //int radioPanelSize = 2000;
            //ImageWithType<double> radioPanel = ImageWithType<double>::zeros(radioPanelSize, radioPanelSize);
            //for (int x = 0; x < radioPanelSize; x++){
            //    for (int y = 0; y < radioPanelSize; y++){
            //        Vec3 dir(x - radioPanelSize/2, y - radioPanelSize/2, fakeFocal);
            //        PixelLoc pixel = PixelLocFromGeoCoord(GeoCoord(dir), longitudeDivideNum, latitudeDivideNum);
            //        radioPanel(y, x) = remainedVotesPanel(pixel.x, pixel.y);
            //    }
            //}
            //cv::GaussianBlur(radioPanel, radioPanel,
            //    cv::Size(7, 7),
            //    4, 4, cv::BORDER_REPLICATE);

            std::vector<Vec3> remainedHLineCandEqs;
            remainedHLineCandEqs.reserve(Square(remainedIntersectionVReps.size()) / 2);
            for (int i = 0; i < remainedIntersectionVReps.size(); i++){
                HPoint2 interp1 = HPointFromVector(GeoCoordFromPixelLoc(
                    remainedIntersectionVReps[i], longitudeDivideNum, latitudeDivideNum).toVector(), fakeFocal);
                Vec3 interv1 = VectorFromHPoint(interp1);
                for (int j = i + 1; j < remainedIntersectionVReps.size(); j++){
                    HPoint2 interp2 = HPointFromVector(GeoCoordFromPixelLoc(
                        remainedIntersectionVReps[i], longitudeDivideNum, latitudeDivideNum).toVector(), fakeFocal);
                    Vec3 interv2 = VectorFromHPoint(interp2);
                    Vec3 eq = normalize(interv1.cross(interv2));

                    remainedHLineCandEqs.push_back(eq);
                }
            }

            std::cout << "remained hline cand num: " << remainedHLineCandEqs.size() << std::endl;

            // todo

            
            

            return result;
        }


        LocalManhattanVanishingPointsDetector::Result LocalManhattanVanishingPointsDetector::operator()(
            const std::vector<Line2> & lines, const Point2 & projCenter) const {
            std::vector<Line2> offsetedLines = lines;
            for (auto & line : offsetedLines) {
                line.first -= projCenter;
                line.second -= projCenter;
            }
            Result result = estimateWithProjectionCenterAtOriginII(lines);
            for (HPoint2 & vp : result.vanishingPoints) {
                vp = vp + HPoint2(projCenter);
            }
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

    
