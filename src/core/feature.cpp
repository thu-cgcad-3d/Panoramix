#include "feature.hpp"

#include <iostream>

#include <Eigen/Dense>

#include "utilities.hpp"

namespace panoramix {
    namespace core {

        PerspectiveCamera::PerspectiveCamera(int w, int h, double focal, const Vec3 & eye,
            const Vec3 & center, const Vec3 & up) 
        : _screenW(w), _screenH(h), _focal(focal), _eye(eye), _center(center), _up(up) {
            Mat4 m;
            m.setIdentity();
            _viewMatrix = Matrix4MakeLookAt(eye, center, up, m);
            
            double verticalViewAngle = atan(_screenH / 2.0 / _focal) * 2;
            double aspect = double(_screenW) / double(_screenH);
            _projectionMatrix = Matrix4MakePerspective(verticalViewAngle, aspect, 0.1, 1e4, m);

            _viewProjectionMatrix = _projectionMatrix * _viewMatrix;
            _viewProjectionMatrixInv = _viewProjectionMatrix.inverse();
        }

        PerspectiveCamera::Vec2 PerspectiveCamera::screenProjection(const Vec3 & p3) const {
            Vec4 p4(p3(0), p3(1), p3(2), 1);
            Vec4 position = _viewProjectionMatrix * p4;
            double xratio = position(0) / position(3) / 2;
            double yratio = position(1) / position(3) / 2;
            double x = (xratio + 0.5) * _screenW;
            double y = _screenH - (yratio + 0.5) * _screenH;
            return Vec2(x, y);
        }

        PerspectiveCamera::Vec3 PerspectiveCamera::spatialDirection(const Vec2 & p2d) const {
            double xratio = (p2d(0) / _screenW - 0.5) * 2;
            double yratio = ((_screenH - p2d(1)) / _screenH - 0.5) * 2;
            Vec4 position(xratio, yratio, 1, 1);
            Vec4 realPosition = _viewProjectionMatrixInv * position;
            return Vec3(realPosition(0) / realPosition(3), realPosition(1) / realPosition(3), realPosition(2) / realPosition(3));
        }


        PanoramicCamera::PanoramicCamera(double focal, const Vec3 & eye,
            const Vec3 & center, const Vec3 & up)
            : _focal(focal), _eye(eye), _center(center), _up(up) {
            _xaxis = (_center - _eye).normalized();
            _yaxis = _up.cross(_xaxis).normalized();
            _zaxis = _xaxis.cross(_yaxis);
        }

        PanoramicCamera::Vec2 PanoramicCamera::screenProjection(const Vec3 & p3) const {
            double xx = p3.dot(_xaxis);
            double yy = p3.dot(_yaxis);
            double zz = p3.dot(_zaxis);
            double longi, lati;
            LongitudeLatitudeFromDirection(Vec3(xx, yy, zz), longi, lati);
            auto sz = screenSize();
            double x = (longi + M_PI) / 2.0 / M_PI * sz.width;
            double y = (lati + M_PI_2) / M_PI * sz.height;
            return Vec2(x, y);
        }

        PanoramicCamera::Vec3 PanoramicCamera::spatialDirection(const Vec2 & p2d) const {
            auto sz = screenSize();
            double longi = p2d(0) / double(sz.width) * 2 * M_PI - M_PI;
            double lati = p2d(1) / double(sz.height) * M_PI - M_PI_2;
            Vec3 dd;
            DirectionFromLongitudeLatitude(longi, lati, dd);
            return dd(0) * _xaxis + dd(1) * _yaxis + dd(2) * _zaxis;
        }


        namespace {

            template <class T, int N>
            void ExtractLines(const cv::Mat& im, std::vector<Line<T, N>> & lines,
                int minlen, int xborderw, int yborderw, int numdir) {

                cv::Mat gim;
                cv::cvtColor(im, gim, CV_BGR2GRAY);
                int h = gim.rows;
                int w = gim.cols;

                cv::Mat dx, dy;
                cv::Mat ggim;
                cv::GaussianBlur(gim, ggim, cv::Size(7, 7), 1.5);
                cv::Sobel(ggim, dx, CV_64F, 1, 0);
                cv::Sobel(ggim, dy, CV_64F, 0, 1);

                cv::Mat imcanny;
                cv::Canny(gim, imcanny, 5, 20);

                std::vector<std::unordered_set<int>> pixelidsset(numdir);

                for (int x = 0; x < imcanny.cols; x++){
                    for (int y = 0; y < imcanny.rows; y++){
                        if (imcanny.at<uchar>(y, x) > 0){
                            double a = atan(dy.at<double>(y, x) / dx.at<double>(y, x));
                            if (a != a){ // NaN
                                continue;
                            }
                            // compute bin id
                            int binid = int((a / M_PI + 0.5) * numdir);
                            if (binid == -1) binid = 0;
                            if (binid == numdir) binid = numdir - 1;

                            int pixelid = y + imcanny.rows * x;

                            pixelidsset[(binid + numdir - 1) % numdir].insert(pixelid);
                            pixelidsset[binid].insert(pixelid);
                            pixelidsset[(binid + 1) % numdir].insert(pixelid);
                        }
                    }
                }

                std::vector<int> xs, ys, ids;
                xs.reserve(512);
                ys.reserve(512);
                ids.reserve(512);

                for (int binid = 0; binid < numdir; binid++){
                    // search in bins
                    auto pixelidnotsearchedyet = pixelidsset[binid];

                    while (true){
                        if (pixelidnotsearchedyet.empty())
                            break;
                        int rootid = *pixelidnotsearchedyet.begin();

                        // BFS
                        xs.clear(); ys.clear(); ids.clear();
                        int x = rootid / imcanny.rows;
                        int y = rootid - imcanny.rows * x;
                        xs.push_back(x);
                        ys.push_back(y);
                        ids.push_back(rootid);

                        // used for search only
                        pixelidnotsearchedyet.erase(rootid);
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
                                int npixelid = ny + imcanny.rows * nx;
                                if (pixelidnotsearchedyet.find(npixelid) != pixelidnotsearchedyet.end()){
                                    xs.push_back(nx);
                                    ys.push_back(ny);
                                    ids.push_back(npixelid);
                                    pixelidnotsearchedyet.erase(npixelid);
                                }
                            }
                            head++;
                        }

                        int edgesize = xs.size();
                        if (edgesize < minlen)
                            continue;

                        cv::Mat xsmat(xs), ysmat(ys);
                        double meanx = cv::mean(xsmat).val[0], meany = cv::mean(ysmat).val[0];
                        cv::Mat zmx = xsmat - meanx, zmy = ysmat - meany;

                        cv::Mat v, lambda;
                        cv::Mat D(2, 2, CV_64FC1);
                        D.at<double>(0, 0) = cv::sum(zmx.mul(zmx)).val[0];
                        D.at<double>(0, 1) = D.at<double>(1, 0) = cv::sum(zmx.mul(zmy)).val[0];
                        D.at<double>(1, 1) = cv::sum(zmy.mul(zmy)).val[0];
                        cv::eigen(D, true, lambda, v);

                        double theta = atan2(v.at<double>(0, 1), v.at<double>(0, 0));
                        double confidence = std::numeric_limits<double>::max();
                        if (lambda.at<double>(1) > 0){
                            confidence = lambda.at<double>(0) / lambda.at<double>(1);
                        }

                        // build line
                        if (confidence >= 400){
                            for (int pid : ids){
                                pixelidsset[binid].erase(pid);
                                pixelidsset[(binid - 1 + numdir) % numdir].erase(pid);
                                pixelidsset[(binid + 1) % numdir].erase(pid);
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

                            lines.push_back({ Vec<T, N>(x1, y1), Vec<T, N>(x2, y2) });
                        }
                    }
                }
            }

        }


        LineSegmentExtractor::Feature LineSegmentExtractor::operator() (const Image & im) const {
            LineSegmentExtractor::Feature lines;
            lines.reserve(200);
            ExtractLines(im, lines, _params.minLength, _params.xBorderWidth, _params.yBorderWidth, _params.numDirs);
            return lines;
        }

    }
}


