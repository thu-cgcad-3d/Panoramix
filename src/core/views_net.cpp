#include "views_net.hpp"

#include "utilities.hpp"

namespace panoramix {
    namespace core {

        ViewsNet::VertHandle ViewsNet::insertPhoto(const Image & im, const PerspectiveCamera & cam) {
            VertData vd;
            vd.camera = vd.originalCamera = cam;
            vd.image = im;
            return insertVertex(vd);
        }

        namespace {

            void LineIntersectons(const std::vector<Line2> & lines,
                std::vector<HPoint2> & hinterps, std::vector<std::pair<int, int>> & lineids,
                bool suppresscross)
            {
                int lnum = lines.size();
                for (int i = 0; i < lnum; i++){
                    auto eqi = cv::Vec3d(lines[i].first[0], lines[i].first[1], 1)
                        .cross(cv::Vec3d(lines[i].second[0], lines[i].second[1], 1));
                    for (int j = i + 1; j < lnum; j++){
                        auto eqj = cv::Vec3d(lines[j].first[0], lines[j].first[1], 1)
                            .cross(cv::Vec3d(lines[j].second[0], lines[j].second[1], 1));
                        auto interp = eqi.cross(eqj);
                        if (interp[0] == 0 && interp[1] == 0 && interp[2] == 0){ // lines overlapped
                            interp[0] = -eqi[1];
                            interp[1] = eqi[0];
                        }
                        interp /= cv::norm(interp);

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
                        lineids.push_back({ i, j });
                    }
                }
            }
        }

        void ViewsNet::computeFeatures(VertHandle h) {
            auto & vd = _views.data(h);
            const Image & im = vd.image;
            vd.featureLineSegment = _params.lineSegmentExtractor(im);
            vd.featureSIFT = _params.siftExtractor(im);
            vd.featureSURF = _params.surfExtractor(im);
            vd.weight = vd.featureLineSegment.size() * _params.lineSegmentWeight +
                vd.featureSIFT.size() * _params.siftWeight +
                vd.featureSURF.size() * _params.surfWeight;

            vd.featureLineIntersections.clear();
            vd.featureLineIntersectionLineIDs.clear();
            LineIntersectons(vd.featureLineSegment, 
                vd.featureLineIntersections, 
                vd.featureLineIntersectionLineIDs, true);
        }

        namespace {
            inline double PerspectiveCameraAngleRadius(const PerspectiveCamera & cam) {
                return atan(sqrt(Square(cam.screenSize().height) + Square(cam.screenSize().width)) /
                    2.0 / cam.focal());
            }
        }

        int ViewsNet::updateConnections(VertHandle h) {
            auto & thisv = _views.data(h);
            const PerspectiveCamera & thisvCam = thisv.originalCamera;
            double thisvCamAngleRadius = PerspectiveCameraAngleRadius(thisvCam);
            thisvCamAngleRadius *= _params.cameraAngleScaler;
            for (auto & v : _views.vertices()){
                if (v.topo.hd == h)
                    continue;
                const PerspectiveCamera & vcam = v.data.originalCamera;
                double vCamAngleRadius = PerspectiveCameraAngleRadius(vcam);
                vCamAngleRadius *= _params.cameraAngleScaler;
                double angleDistance = 0;
                AngleBetweenDirections(thisvCam.center(), vcam.center(), angleDistance);
                if (angleDistance <= thisvCamAngleRadius + vCamAngleRadius){
                    // may overlap
                    HalfData hd;
                    hd.cameraAngleDistance = angleDistance;
                    _views.addEdge(h, v.topo.hd, hd);
                }
            }
            return _views.topo(h).halfedges.size();
        }

        void ViewsNet::computeTransformationOnConnections(VertHandle h) {
            
        }

        void ViewsNet::calibrateCamera(VertHandle h) {
            
        }

        void ViewsNet::calibrateAllCameras() {
            
        }

        namespace {

            void FindVanishingPoints(const std::vector<Vec2>& spinterps, std::vector<Vec2>& vps)
            {
                vps.resize(3);
                int longtdiv = 1000, latdiv = 500;
                cv::Mat votepanel = cv::Mat::zeros(longtdiv, latdiv, CV_32FC1);
                static_assert(sizeof(float) == 4, "invalid double size!");
                int pn = spinterps.size();
                for (const cv::Vec2d& p : spinterps){
                    int longtid = int(p[0] * longtdiv / M_PI / 2);
                    int latid = int((p[1] + M_PI_2) * latdiv / M_PI);
                    longtid = (longtid % longtdiv + longtdiv) % longtdiv;
                    latid = (latid % latdiv + latdiv) % latdiv;
                    votepanel.at<float>(longtid, latid) += 1;
                }
                cv::GaussianBlur(votepanel, votepanel, cv::Size((longtdiv / 50) * 2 + 1, (latdiv / 50) * 2 + 1),
                    4, 4, cv::BORDER_REPLICATE);

                double minval = 0, maxval = 0;
                int maxidx[] = { -1, -1 };
                cv::minMaxIdx(votepanel, &minval, &maxval, 0, maxidx);

                vps[0] = cv::Vec2d(maxidx[0] * M_PI * 2 / longtdiv, maxidx[1] * M_PI / latdiv - M_PI_2);
                cv::Vec3d vec0(cos(vps[0][0])*cos(vps[0][1]), sin(vps[0][0])*cos(vps[0][1]), sin(vps[0][1]));

                // iterate locations orthogonal to vps[0]
                double maxscore = -1;
                for (int x = 0; x < longtdiv; x++){
                    double longt1 = double(x) / longtdiv*M_PI * 2;
                    double lat1 = -atan((vec0[0] * cos(longt1) + vec0[1] * sin(longt1)) / vec0[2]);
                    double longt1rev = longt1 + M_PI;
                    double lat1rev = -lat1;
                    double longt2 = longt1 + M_PI_2;
                    double lat2 = -atan((vec0[0] * cos(longt2) + vec0[1] * sin(longt2)) / vec0[2]);
                    double longt2rev = longt2 + M_PI;
                    double lat2rev = -lat2;

                    double longts[] = { longt1, longt1rev, longt2, longt2rev };
                    double lats[] = { lat1, lat1rev, lat2, lat2rev };

                    double score = 0;
                    for (int k = 0; k < 4; k++){
                        int xx = int(longts[k] / M_PI / 2 * longtdiv) % longtdiv;
                        if (xx == -1) xx = 0;
                        if (xx == longtdiv) xx = longtdiv;
                        int yy = int((lats[k] / M_PI + 0.5)*latdiv);
                        if (yy == -1) yy = 0;
                        if (yy == latdiv) yy = latdiv - 1;
                        score += votepanel.at<float>(xx, yy);

                    }
                    if (score > maxscore){
                        maxscore = score;
                        vps[1] = cv::Vec2d(longt1, lat1);
                        vps[2] = cv::Vec2d(longt2, lat2);
                    }
                }

            }

        }


        void ViewsNet::computeGlobalFeatures() {
            GlobalData newgb;

            // estimate vanishing points


            // classify lines


            // add spatial line segments
            int spatialLineSegmentsNum = 0;
            for (auto & v : _views.vertices())
                spatialLineSegmentsNum += v.data.featureLineSegment.size();
            newgb.spatialLineSegments.resize(spatialLineSegmentsNum);
            auto spatialLineSegmentBegin = newgb.spatialLineSegments.begin();
            for (auto & v : _views.vertices()){
                spatialLineSegmentBegin = std::transform(v.data.featureLineSegment.begin(), v.data.featureLineSegment.end(),
                    spatialLineSegmentBegin, [&v](const Line2 & line) -> Line3{
                    auto & p1 = line.first;
                    auto & p2 = line.second;
                    auto pp1 = v.data.camera.spatialDirection(p1);
                    auto pp2 = v.data.camera.spatialDirection(p2);
                    return Line3{ pp1, pp2 };
                });
            }

        }
 
    }
}