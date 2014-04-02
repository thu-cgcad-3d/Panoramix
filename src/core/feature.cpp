#include "feature.hpp"

#include <iostream>
#include <unordered_map>

#include <Eigen/Dense>

#include "utilities.hpp"

namespace panoramix {
    namespace core {

        namespace {

            Mat4 MakeMat4LookAt(const Vec3 & eye, const Vec3 & center,
                const Vec3 & up) {
                Vec3 zaxis = (center - eye); zaxis /= core::norm(zaxis);
                Vec3 xaxis = up.cross(zaxis); xaxis /= core::norm(xaxis);
                Vec3 yaxis = zaxis.cross(xaxis);
                Mat4 m(
                    -xaxis(0), yaxis(0), -zaxis(0), 0,
                    -xaxis(1), yaxis(1), -zaxis(1), 0,
                    -xaxis(2), yaxis(2), -zaxis(2), 0,
                    xaxis.dot(eye), -yaxis.dot(eye), zaxis.dot(eye), 1);
                return m.t();
            }

            Mat4 MakeMat4Perspective(const double & fovyRadians, const double & aspect,
                const double & nearZ, const double & farZ) {
                double cotan = double(1.0) / std::tan(fovyRadians / 2.0);
                Mat4 m(
                    cotan / aspect, 0, 0, 0,
                    0, cotan, 0, 0,
                    0, 0, (farZ + nearZ) / (nearZ - farZ), -1,
                    0, 0, (2 * farZ * nearZ) / (nearZ - farZ), 0);
                return m.t();
            }

        }

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

        Vec3 PerspectiveCamera::spatialDirection(const Vec2 & p2d) const {
            double xratio = (p2d(0) / _screenW - 0.5) * 2;
            double yratio = ((_screenH - p2d(1)) / _screenH - 0.5) * 2;
            Vec4 position(xratio, yratio, 1, 1);
            Vec4 realPosition = _viewProjectionMatrixInv * position;
            return Vec3(realPosition(0) / realPosition(3), realPosition(1) / realPosition(3), realPosition(2) / realPosition(3));
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

           
            void ExtractLinesExperimental(const cv::Mat& im, std::vector<Line2> & lines,
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
                std::set<int> edgePixels; // store all edge pixels in a set
                for (int x = 0; x < imCanny.cols; x++){
                    for (int y = 0; y < imCanny.rows; y++){
                        if (imCanny.at<uchar>(y, x) > 0){
                            double a = atan(dy.at<double>(y, x) / dx.at<double>(y, x));
                            if (isnan(a)){ // NaN
                                continue;
                            }
                            // compute bin id
                            int binId = int((a / M_PI + 0.5) * numDir);
                            imBinIds.at<int32_t>(y, x) = WrapBetween(binId, 0, numDir);
                            edgePixels.insert(y * w + x);
                        }
                        else{
                            imBinIds.at<int32_t>(y, x) = -1;
                        }
                    }
                }

                std::cout << "collecting pixels.." << std::endl;

                cv::Mat imComponentIds(im.size(), CV_32SC1); // component ids int32_t
                cv::Mat imParentsInTraversal(im.size(), CV_32SC2); // 
                std::map<int, bool> componentUsed;
                std::fill(imComponentIds.begin<int32_t>(), imComponentIds.end<int32_t>(), -1); // init to -1

                std::vector<int> xs, ys;
                xs.reserve(2048);
                ys.reserve(2048);

                int componentIdGenerator = 0;
                while (true){
                    bool allPixelsAreUsed = true;
                    PixelLoc unusedPixel;
                    for (int index : edgePixels){
                        int y = index / w;
                        int x = index - y * w;
                        unusedPixel = {x, y};
                        if (imBinIds.at<int32_t>(unusedPixel) == -1) // this is not an edge pixel
                            continue;
                        int32_t componentId = imComponentIds.at<int32_t>(unusedPixel);
                        if (componentUsed.find(componentId)!=componentUsed.end() && componentUsed[componentId]) 
                            // this component is already used to compose a line
                            continue;
                        allPixelsAreUsed = false;
                        break;
                    }
                    if (allPixelsAreUsed)
                        break;

                    const PixelLoc & fistPixel = unusedPixel;
                    int activeComponentId = componentIdGenerator++;
                    int activeBinId = imBinIds.at<int32_t>(fistPixel);
                    int activeSimiliarBinIds[] = { WrapBetween(activeBinId + 1, 0, numDir), WrapBetween(activeBinId - 1, 0, numDir) };
                    int pixelsNum = 0;
                    xs.clear(); 
                    ys.clear();

                    static const int xdirs[] = { 1, 1, 0, -1, -1, -1, 0, 1 };
                    static const int ydirs[] = { 0, 1, 1, 1, 0, -1, -1, -1 };
                    PixelLoc detector = fistPixel;
                    PixelLoc detectorPrev = fistPixel;
                    while (true){ // depth first search
                        imComponentIds.at<int32_t>(detector) = activeComponentId;
                        xs.push_back(detector.x);
                        ys.push_back(detector.y);
                        pixelsNum++;
                        
                        bool allConnectedPixelsAreFound = false;
                        while (true){
                            bool noWayOut = true;
                            for (int k = 0; k < 8; k++){
                                PixelLoc next = detector;
                                next.x += xdirs[k];
                                next.y += ydirs[k];
                                if (!IsBetween(next.x, 0, w - 1) || !IsBetween(next.y, 0, h - 1)) // outsided
                                    continue;

                                int32_t binId = imBinIds.at<int32_t>(next);
                                if (binId == -1)
                                    continue;
                                if (binId != activeBinId && binId != activeSimiliarBinIds[0] && binId != activeSimiliarBinIds[1])
                                    // bin ids are not even similiar
                                    continue;
                                int32_t componentId = imComponentIds.at<int32_t>(next);
                                if (componentId == activeComponentId) // already visited
                                    continue;
                                if (componentUsed.find(componentId) != componentUsed.end() && componentUsed[componentId])
                                    // this component is already used to compose a line
                                    continue;

                                detectorPrev = detector;
                                detector = next;
                                imParentsInTraversal.at<Vec<int32_t, 2>>(detector) = { detectorPrev.x, detectorPrev.y };
                                noWayOut = false;
                                break;
                            }
                            if (!noWayOut){
                                break; // a way is found
                            }else if (detector == fistPixel){ // no way at root pixel
                                allConnectedPixelsAreFound = true;
                                break;
                            }else{ // not root pixel yet, traverse back
                                auto parent = imParentsInTraversal.at<Vec<int32_t, 2>>(detector);
                                detector = { parent[0], parent[1] };
                            }
                        }
                        if (allConnectedPixelsAreFound)
                            break;
                    }

                    std::cout << "pixels num: " << pixelsNum << std::endl;
                    if (pixelsNum < minlen){
                        // declare that this component is abandoned
                        componentUsed[activeComponentId] = true;
                        continue;
                    }

                    // check confidence
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

                    if (confidence < 400){
                        continue;
                    }

                    // declare that this component is used
                    componentUsed[activeComponentId] = true;

                    // install line
                    auto xends = std::minmax_element(xs.begin(), xs.end());
                    auto yends = std::minmax_element(ys.begin(), ys.end());
                    double minx = *xends.first, maxx = *xends.second;
                    double miny = *yends.first, maxy = *yends.second;

                    // ignore boundary line
                    if (maxx <= xborderw || minx >= w - xborderw || maxy <= yborderw || miny >= h - yborderw)
                        continue;

                    double len = sqrt((maxx - minx)*(maxx - minx) + (maxy - miny)*(maxy - miny));
                    double x1 = meanx - cos(theta) * len / 2;
                    double x2 = meanx + cos(theta) * len / 2;
                    double y1 = meany - sin(theta) * len / 2;
                    double y2 = meany + sin(theta) * len / 2;

                    lines.push_back({ Vec2(x1, y1), Vec2(x2, y2) });
                }



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
                std::cout << "done" << std::endl;

            }

        }


        LineSegmentExtractor::Feature LineSegmentExtractor::operator() (const Image & im) const {
            LineSegmentExtractor::Feature lines;
            lines.reserve(1000);
            if (_params.useExperimentalAlgorithm)
                ExtractLinesExperimental(im, lines, _params.minLength, 
                    _params.xBorderWidth, _params.yBorderWidth, _params.numDirs);
            else
                ExtractLines(im, lines, _params.minLength, 
                    _params.xBorderWidth, _params.yBorderWidth, _params.numDirs);
            return lines;
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

    }
}

    
