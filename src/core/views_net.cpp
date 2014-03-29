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

            void LineIntersectons(const std::vector<Classified<Line2>> & lines,
                std::vector<HPoint2> & hinterps, std::vector<std::pair<int, int>> & lineids,
                bool suppresscross)
            {
                int lnum = lines.size();
                for (int i = 0; i < lnum; i++){
                    auto eqi = cv::Vec3d(lines[i].component.first[0], lines[i].component.first[1], 1)
                        .cross(cv::Vec3d(lines[i].component.second[0], lines[i].component.second[1], 1));
                    for (int j = i + 1; j < lnum; j++){
                        auto eqj = cv::Vec3d(lines[j].component.first[0], lines[j].component.first[1], 1)
                            .cross(cv::Vec3d(lines[j].component.second[0], lines[j].component.second[1], 1));
                        auto interp = eqi.cross(eqj);
                        if (interp[0] == 0 && interp[1] == 0 && interp[2] == 0){ // lines overlapped
                            interp[0] = -eqi[1];
                            interp[1] = eqi[0];
                        }
                        interp /= norm(interp);

                        if (suppresscross){
                            auto& a1 = lines[i].component.first;
                            auto& a2 = lines[i].component.second;
                            auto& b1 = lines[j].component.first;
                            auto& b2 = lines[j].component.second;
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
            auto lineSegments = _params.lineSegmentExtractor(im);
            vd.lineSegments.resize(lineSegments.size());
            for (size_t i = 0; i < lineSegments.size(); i++){
                vd.lineSegments[i].claz = -1;
                vd.lineSegments[i].component = lineSegments[i];
            }

            vd.SIFTs = _params.siftExtractor(im);
            vd.SURFs = _params.surfExtractor(im);
            vd.weight = vd.lineSegments.size() * _params.lineSegmentWeight +
                vd.SIFTs.size() * _params.siftWeight +
                vd.SURFs.size() * _params.surfWeight;

            vd.lineSegmentIntersections.clear();
            vd.lineSegmentIntersectionLineIDs.clear();
            LineIntersectons(vd.lineSegments, 
                vd.lineSegmentIntersections, 
                vd.lineSegmentIntersectionLineIDs, true);

            // build region net
            vd.regionNet = std::make_shared<RegionsNet>(vd.image);
            vd.regionNet->buildNetAndComputeGeometricFeatures();
            vd.regionNet->computeImageFeatures();
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
                double angleDistance = AngleBetweenDirections(thisvCam.center(), vcam.center());
                if (angleDistance <= thisvCamAngleRadius + vCamAngleRadius){
                    // may overlap
                    HalfData hd;
                    hd.cameraAngleDistance = angleDistance;
                    _views.addEdge(h, v.topo.hd, hd);
                }
            }
            return _views.topo(h).halfedges.size();
        }

        ViewsNet::VertHandle ViewsNet::isTooCloseToAnyExistingView(VertHandle h) const {
            auto & camera = _views.data(h).camera;
            double cameraRadius = PerspectiveCameraAngleRadius(camera);
            auto connections = _views.topo(h).halfedges;
            for (auto & con : connections){
                auto & neighborCamera = _views.data(_views.topo(con).to()).camera;
                double cameraAngle = AngleBetweenDirections(camera.center(), neighborCamera.center());
                double neighborCameraRadius = PerspectiveCameraAngleRadius(camera);
                if (cameraAngle <= (cameraRadius + neighborCameraRadius) * _params.smallCameraAngleScalar){ // too close
                    return _views.topo(con).to();
                }
            }
            return VertHandle();
        }


        void ViewsNet::computeTransformationOnConnections(VertHandle h) {
            
        }

        void ViewsNet::calibrateCamera(VertHandle h) {
            
        }

        void ViewsNet::calibrateAllCameras() {
            
        }

        namespace {

            inline cv::Point PixelIndexFromGeoCoord(const GeoCoord & p, int longidiv, int latidiv) {
                int longtid = static_cast<int>((p.longitude + M_PI) * longidiv / M_PI / 2);
                int latid = static_cast<int>((p.latitude + M_PI_2) * latidiv / M_PI);
                longtid = (longtid % longidiv + longidiv) % longidiv;
                latid = (latid % latidiv + latidiv) % latidiv;
                return cv::Point(longtid, latid);
            }

            inline GeoCoord GeoCoordFromPixelIndex(const cv::Point & pixel, int longidiv, int latidiv) {
                return GeoCoord{ pixel.x * M_PI * 2 / longidiv - M_PI, pixel.y * M_PI / latidiv - M_PI_2 };
            }

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

            std::array<Vec3, 3> FindVanishingPoints(const std::vector<Vec3>& intersections,
                int longitudeDivideNum = 1000, int latitudeDivideNum = 500) {

                std::array<Vec3, 3> vps;

                cv::Mat votePanel = cv::Mat::zeros(longitudeDivideNum, latitudeDivideNum, CV_32FC1);

                int pn = intersections.size();
                for (const Vec3& p : intersections){
                    cv::Point pixel = PixelIndexFromGeoCoord(GeoCoord(p), longitudeDivideNum, latitudeDivideNum);
                    votePanel.at<float>(pixel.x, pixel.y) += 1.0;
                }
                cv::GaussianBlur(votePanel, votePanel, cv::Size((longitudeDivideNum / 50) * 2 + 1, (latitudeDivideNum / 50) * 2 + 1),
                    4, 4, cv::BORDER_REPLICATE);

                double minVal = 0, maxVal = 0;
                int maxIndex[] = { -1, -1 };
                cv::minMaxIdx(votePanel, & minVal, & maxVal, 0, maxIndex);
                cv::Point maxPixel(maxIndex[0], maxIndex[1]);

                vps[0] = GeoCoordFromPixelIndex(maxPixel, longitudeDivideNum, latitudeDivideNum).toVector();
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
                        cv::Point pixel = PixelIndexFromGeoCoord(GeoCoord(v), longitudeDivideNum, latitudeDivideNum);
                        score += votePanel.at<float>(WrapBetween(pixel.x, 0, longitudeDivideNum), WrapBetween(pixel.y, 0, latitudeDivideNum));
                    }
                    if (score > maxScore){
                        maxScore = score;
                        vps[1] = vec1;
                        vps[2] = vec2;
                    }
                }

                if (UnOrthogonality(vps[0], vps[1], vps[2]) < 0.1)
                    return vps;
                
                maxScore = -1;
                for (int y = 0; y < latitudeDivideNum; y++){
                    double lat1 = double(y) / latitudeDivideNum * M_PI - M_PI_2;
                    double longt1s[] = { Longitude1FromLatitudeAndNormalVector(lat1, vec0), Longitude2FromLatitudeAndNormalVector(lat1, vec0) };
                    for (double longt1 : longt1s){
                        Vec3 vec1 = GeoCoord(longt1, lat1).toVector();
                        Vec3 vec1rev = -vec1;
                        Vec3 vec2 = vec0.cross(vec1);
                        Vec3 vec2rev = -vec2;
                        Vec3 vecs[] = { vec1, vec1rev, vec2, vec2rev };

                        double score = 0;
                        for (Vec3 & v : vecs){
                            cv::Point pixel = PixelIndexFromGeoCoord(GeoCoord(v), longitudeDivideNum, latitudeDivideNum);
                            score += votePanel.at<float>(WrapBetween(pixel.x, 0, longitudeDivideNum), WrapBetween(pixel.y, 0, latitudeDivideNum));
                        }
                        if (score > maxScore){
                            maxScore = score;
                            vps[1] = vec1;
                            vps[2] = vec2;
                        }
                    }
                }

                return vps;

            }

            template <class Vec3Container>
            void ClassifyLines(const Vec3Container & points, std::vector<Classified<Line3>> & lines,
                double angleThreshold = M_PI/3, double sigma = 0.1) {
                
                size_t nlines = lines.size();
                size_t npoints = points.size();

                for (size_t i = 0; i < nlines; i++){
                    Vec3 a = lines[i].component.first;
                    Vec3 b = lines[i].component.second;
                    Vec3 normab = a.cross(b);
                    normab /= norm(normab);

                    std::vector<double> lineangles(npoints);
                    std::vector<double> linescores(npoints);

                    for (int j = 0; j < npoints; j++){
                        Vec3 point = points[j];
                        double angle = abs(asin(normab.dot(point)));
                        lineangles[j] = angle;
                    }

                    // get score based on angle
                    for (int j = 0; j < npoints; j++){
                        double angle = lineangles[j];
                        double score = exp(-(angle / angleThreshold) * (angle / angleThreshold) / sigma / sigma / 2);
                        linescores[j] = (angle > angleThreshold) ? 0 : score;
                    }

                    // classify lines
                    lines[i].claz = -1;
                    float curscore = 0.8;
                    for (int j = 0; j < npoints; j++){
                        if (linescores[j] > curscore){
                            lines[i].claz = j;
                            curscore = linescores[j];
                        }
                    }
                }
            }

            inline Vec3 RotateDirectionTo(const Vec3 & from, const Vec3 & toDirection, double angle) {
                Vec3 tovec = from.cross(toDirection).cross(from);
                Vec3 result3 = from + tovec * tan(angle);
                return result3 / norm(result3);
            }
        }


        void ViewsNet::estimateVanishingPointsAndClassifyLines() {

            // pick separated views only
            auto seperatedViewIters = MergeNear(_views.vertices().begin(), _views.vertices().end(), std::false_type(),
                _params.smallCameraAngleScalar, [](const ViewMesh::Vertex & v1, const ViewMesh::Vertex & v2) -> double{
                double angleDistance = AngleBetweenDirections(v1.data.camera.center(), v2.data.camera.center());
                return angleDistance / 
                    (PerspectiveCameraAngleRadius(v1.data.camera) + PerspectiveCameraAngleRadius(v2.data.camera));
            });

            // collect line intersections
            int lineIntersectionsNum = 0;
            for (const auto & vIter : seperatedViewIters) // count intersecton num
                lineIntersectionsNum += vIter->data.lineSegmentIntersections.size();
            std::vector<Vec3> intersections(lineIntersectionsNum);
            auto intersectionsBegin = intersections.begin();
            for (const auto & vIter : seperatedViewIters){ // projection 2d intersections to global GeoCoord
                intersectionsBegin = std::transform(vIter->data.lineSegmentIntersections.begin(), 
                    vIter->data.lineSegmentIntersections.end(),
                    intersectionsBegin, 
                    [&vIter](const HPoint2 & p) -> Vec3 {
                    Vec3 p3 = vIter->data.camera.spatialDirection(p.toPoint());
                    return p3;
                });
            }

            // get merged intersections
            auto mergedIntersectionsIters = MergeNear(intersections.begin(), intersections.end(), std::false_type(), 
                M_PI / 150.0, [](const Vec3 & p1, const Vec3 & p2) -> double {
                return AngleBetweenDirections(p1, p2);
            });
            _globalData.mergedSpatialLineSegmentIntersections.clear();
            _globalData.mergedSpatialLineSegmentIntersections.reserve(mergedIntersectionsIters.size());
            for (auto intersectionIter : mergedIntersectionsIters) {
                _globalData.mergedSpatialLineSegmentIntersections.push_back(*intersectionIter);
            }

            // find vanishing points;
            _globalData.vanishingPoints = FindVanishingPoints(intersections);                        

            // add spatial line segments from line segments of all views
            int spatialLineSegmentsNum = 0;
            for (auto & v : _views.vertices())
                spatialLineSegmentsNum += v.data.lineSegments.size();
            _globalData.spatialLineSegments.resize(spatialLineSegmentsNum);
            auto spatialLineSegmentBegin = _globalData.spatialLineSegments.begin();
            for (auto & v : _views.vertices()){
                spatialLineSegmentBegin = std::transform(v.data.lineSegments.begin(), v.data.lineSegments.end(),
                    spatialLineSegmentBegin, [&v](const Classified<Line2> & line) -> Classified<Line3>{
                    auto & p1 = line.component.first;
                    auto & p2 = line.component.second;
                    auto pp1 = v.data.camera.spatialDirection(p1);
                    auto pp2 = v.data.camera.spatialDirection(p2);
                    Classified<Line3> cline3;
                    cline3.claz = line.claz;
                    cline3.component = Line3{ pp1, pp2 };
                    return cline3;
                });
            }

            // classify lines
            ClassifyLines(_globalData.vanishingPoints, _globalData.spatialLineSegments);
     
            // project line classes back to perspective views
            spatialLineSegmentBegin = _globalData.spatialLineSegments.begin();
            for (auto & v : _views.vertices()){
                for (auto & line : v.data.lineSegments){
                    line.claz = spatialLineSegmentBegin->claz;
                    ++ spatialLineSegmentBegin;
                }
            }


        }


        namespace {

            std::vector<Classified<Line3>> MergeColinearSpatialLines(const std::vector<Classified<Line3>> & oldLines) {
                
                std::vector<Classified<Line3>> lines(oldLines);

                // normalize line point directions
                for (auto & line : lines) {
                    line.component.first /= norm(line.component.first);
                    line.component.second /= norm(line.component.second);
                }

                // group all colinear spatial lines
                auto colinearLineIters = MergeNear(lines.begin(), lines.end(), std::true_type(), 
                    0.01, [](const Classified<Line3> & line1, const Classified<Line3> & line2) -> double{
                    if (line1.claz != line2.claz)
                        return 100.0;
                    auto normal1 = line1.component.first.cross(line1.component.second);
                    normal1 /= norm(normal1);
                    auto normal2 = line2.component.first.cross(line2.component.second);
                    normal2 /= norm(normal2);
                    double angle = abs(asin(norm(normal1.cross(normal2))));
                    return angle;
                });


                std::vector<Classified<Line3>> mergedLines;
                mergedLines.reserve(oldLines.size());
                
                assert(!colinearLineIters.empty());
                colinearLineIters.push_back(lines.end());
                auto colinearedBeginIter = colinearLineIters.begin();                
                for (; *colinearedBeginIter != lines.end(); ++colinearedBeginIter) {
                    auto colinearedBegin = *colinearedBeginIter;
                    auto colinearedEnd = *std::next(colinearedBeginIter);
                    assert(colinearedBegin != colinearedEnd); // not empty
                    int claz = colinearedBegin->claz;

                    size_t lineNum = std::distance(colinearedBegin, colinearedEnd);
                    if (lineNum == 1){ // only one line, no need to merge
                        mergedLines.push_back(*colinearedBegin);
                        continue;
                    }
                    
                    auto & firstLine = *colinearedBegin; // get the normal direction of first line
                    /*
                        second point
                        ^     \
                        |      \
                        x ----- >first point
                        normal
                    */
                    auto firstNormal = firstLine.component.first.cross(firstLine.component.second);
                    firstNormal /= norm(firstNormal);
                    
                    for (auto lineIter = colinearedBegin; lineIter != colinearedEnd; ++lineIter) {
                        assert(lineIter->claz == claz); // all lines in the same group should be in the same class
                        auto & line = *lineIter;
                        auto normal = line.component.first.cross(line.component.second);
                        normal /= norm(normal);
                        if (normal.dot(firstNormal) < 0) { // not same direction, then must be opposite
                            // swap order of the line to make all line directions in this group consistent
                            std::swap(line.component.first, line.component.second); 
                        }
                    }

                    auto firstPointDirection = firstLine.component.first;
                    auto firstNormalCrossPoint = firstNormal.cross(firstPointDirection);

                    // compute line projection angles
                    std::vector<std::pair<double, double>> lineAngleSegments(lineNum);
                    size_t count = 0;
                    for (auto lineIter = colinearedBegin; lineIter != colinearedEnd; ++lineIter){
                        const auto & line = *lineIter;
                        auto & pdir1 = line.component.first;
                        Vec2 pdv1(pdir1.dot(firstPointDirection), pdir1.dot(firstNormalCrossPoint));
                        double angle1 = SignedAngleBetweenDirections(pdv1, Vec2(1, 0), true);
                        double angle2 = angle1 + AngleBetweenDirections(line.component.first, line.component.second); // may be higher than + M_PI
                        lineAngleSegments[count++] = std::make_pair(angle1, angle2);
                    }

                    // sort the angles
                    std::sort(lineAngleSegments.begin(), lineAngleSegments.end(), 
                        [](const std::pair<double, double> & a1, const std::pair<double, double> & a2){
                        return a1.first < a2.first;
                    });

                    // now merge line angles
                    std::vector<std::pair<double, double>> mergedLineAngleSegments;
                    double curFrom = lineAngleSegments.front().first;
                    double curTo = lineAngleSegments.front().second;
                    for (auto & lineAngle : lineAngleSegments) {
                        if (lineAngle.first <= curTo){
                            curTo = lineAngle.second;
                        } else {
                            if (curTo - curFrom >= M_PI){
                                mergedLineAngleSegments.push_back(std::make_pair(curFrom, (curFrom + curTo) / 2.0));
                                mergedLineAngleSegments.push_back(std::make_pair((curFrom + curTo) / 2.0, curTo));
                            }else
                                mergedLineAngleSegments.push_back(std::make_pair(curFrom, curTo));
                            curFrom = lineAngle.first;
                            curTo = lineAngle.second;
                        }
                    }
                    if (curTo - curFrom >= M_PI){
                        mergedLineAngleSegments.push_back(std::make_pair(curFrom, (curFrom + curTo) / 2.0));
                        mergedLineAngleSegments.push_back(std::make_pair((curFrom + curTo) / 2.0, curTo));
                    }
                    else
                        mergedLineAngleSegments.push_back(std::make_pair(curFrom, curTo));


                    // recover lines from angle pairs
                    for (auto & lineAngle : mergedLineAngleSegments) {
                        Vec3 d1 = Vec3(1, 0, 0) * cos(lineAngle.first) + Vec3(0, 1, 0) * sin(lineAngle.first);
                        Vec3 dd1 = d1(0) * firstPointDirection + d1(1) * firstNormalCrossPoint + d1(2) * firstNormal;
                        Vec3 d2 = Vec3(1, 0, 0) * cos(lineAngle.second) + Vec3(0, 1, 0) * sin(lineAngle.second);
                        Vec3 dd2 = d2(0) * firstPointDirection + d2(1) * firstNormalCrossPoint + d2(2) * firstNormal;
                        Classified<Line3> line;
                        line.claz = claz;
                        line.component.first = dd1;
                        line.component.second = dd2;
                        mergedLines.push_back(line);
                    }
                }

                return mergedLines;
            }

        }


        void ViewsNet::rectifySpatialLines() {
            // merge lines
            _globalData.mergedSpatialLineSegments = MergeColinearSpatialLines(_globalData.spatialLineSegments);

            // build constraint graph
            struct ConstriantData {
                int lineids[2];
                Vec3 position;
                double weight;
                enum {
                    Intersection,
                    Incidence
                } type;
            };
            std::vector<ConstriantData> constraints;
            constraints.reserve(Square(_globalData.spatialLineSegments.size()));
        }
 
    }
}