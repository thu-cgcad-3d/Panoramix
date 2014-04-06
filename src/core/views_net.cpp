#include "views_net.hpp"

#include "utilities.hpp"
#include "optimization.hpp"

namespace panoramix {
    namespace core {

        ViewsNet::Params::Params() 
            : camera(250.0), lineSegmentWeight(1.0), siftWeight(1.0),
            surfWeight(1.0), cameraAngleScaler(1.8), smallCameraAngleScalar(0.05),
            intersectionConstraintLineDistanceAngleThreshold(0.06),
            incidenceConstraintLineDistanceAngleThreshold(0.2),
            mergeLineDistanceAngleThreshold(0.05),
            mjWeightTriplet(5.0), mjWeightX(5.0), mjWeightT(2.0), mjWeightL(1.0), mjWeightI(2.0) {}

        ViewsNet::VertHandle ViewsNet::insertPhoto(const Image & im, const PerspectiveCamera & cam,
            double cameraDirectionErrorScale,
            double cameraPositionErrorScale) {
            VertData vd;
            vd.camera = vd.originalCamera = cam;
            vd.cameraDirectionErrorScale = cameraDirectionErrorScale;
            vd.cameraPositionErrorScale = cameraPositionErrorScale;
            vd.image = im;
            return insertVertex(vd);
        }

        namespace {

            void LineIntersectons(const std::vector<Classified<Line2>> & lines,
                std::vector<HPoint2> & hinterps, std::vector<std::pair<int, int>> & lineids,
                bool suppresscross)
            {
                size_t lnum = lines.size();
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
            return static_cast<int>(_views.topo(h).halfedges.size());
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




        namespace {

            

        }
        


        void ViewsNet::calibrateCamera(VertHandle h) {
            auto & vd = _views.data(h);
            if (vd.cameraDirectionErrorScale == 0 && vd.cameraPositionErrorScale == 0)
                return;
            // TODO
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
                return GeoCoord{ pixel .x * M_PI * 2 / longidiv - M_PI, pixel.y * M_PI / latidiv - M_PI_2 };
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

                std::cout << "begin voting ..." << std::endl;
                size_t pn = intersections.size();
                for (const Vec3& p : intersections){
                    cv::Point pixel = PixelIndexFromGeoCoord(GeoCoord(p), longitudeDivideNum, latitudeDivideNum);
                    votePanel.at<float>(pixel.x, pixel.y) += 1.0;
                }
                std::cout << "begin gaussian bluring ..." << std::endl;
                cv::GaussianBlur(votePanel, votePanel, cv::Size((longitudeDivideNum / 50) * 2 + 1, (latitudeDivideNum / 50) * 2 + 1),
                    4, 4, cv::BORDER_REPLICATE);
                std::cout << "done voting" << std::endl;

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
                
                // failed, then use y instead of x
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


















        ViewsNet::ConstraintData::ConstraintData(){
            mergedSpatialLineSegmentIds[0] = mergedSpatialLineSegmentIds[1] = 0;
            weight = 0.0;
            std::memset(lineVotings, 0, sizeof(lineVotings));
            junctionWeights.I = junctionWeights.L = junctionWeights.T = junctionWeights.Triplet = junctionWeights.X = 0.0;
        }

        namespace {

            // merge colinear spatial lines and find some incidence constraints
            std::vector<Classified<Line3>> MergeColinearSpatialLinesAndAppendIncidenceConstraints(
                const std::vector<Classified<Line3>> & oldLines, 
                std::vector<int> & chainIds, std::vector<ViewsNet::ConstraintData> & constraints, 
                double mergeAngleThres, double incidenceAngleThres) {
                
                std::vector<Classified<Line3>> lines(oldLines);

                // normalize line point directions
                for (auto & line : lines) {
                    line.component.first /= norm(line.component.first);
                    line.component.second /= norm(line.component.second);
                }

                // group all colinear spatial lines
                auto colinearLineIters = MergeNear(lines.begin(), lines.end(), std::true_type(), 
                    mergeAngleThres,
                    [](const Classified<Line3> & line1, const Classified<Line3> & line2) -> double{
                    if (line1.claz != line2.claz)
                        return 100.0; // never merge lines with different classes
                    auto normal1 = line1.component.first.cross(line1.component.second);
                    auto normal2 = line2.component.first.cross(line2.component.second);
                    return std::min(AngleBetweenDirections(normal1, normal2), 
                        AngleBetweenDirections(normal1, -normal2));
                });


                std::vector<Classified<Line3>> mergedLines;
                mergedLines.reserve(oldLines.size());
                chainIds.clear();
                chainIds.reserve(mergedLines.capacity());
                int chainId = 0;
                
                assert(!colinearLineIters.empty());
                colinearLineIters.push_back(lines.end());
                auto colinearedBeginIter = colinearLineIters.begin();                
                for (; *colinearedBeginIter != lines.end(); ++colinearedBeginIter, ++chainId) {
                    auto colinearedBegin = *colinearedBeginIter;
                    auto colinearedEnd = *std::next(colinearedBeginIter);
                    assert(colinearedBegin != colinearedEnd); // not empty
                    int claz = colinearedBegin->claz;

                    size_t lineNum = std::distance(colinearedBegin, colinearedEnd);
                    if (lineNum == 1){ // only one line, no need to merge
                        mergedLines.push_back(*colinearedBegin);
                        chainIds.push_back(chainId);
                        continue;
                    }
                    
                    auto & firstLine = *colinearedBegin; // get the normal direction of first line
                    /*
                        second point
                        ^     \
                        |      \ the line
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

                    // compute line projection angle spans
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

                    // sort lines in start-angle's ascending order
                    std::sort(lineAngleSegments.begin(), lineAngleSegments.end(), 
                        [](const std::pair<double, double> & a1, const std::pair<double, double> & a2){
                        return a1.first < a2.first;
                    });

                    // now merge line angle spans
                    std::vector<std::pair<double, double>> mergedLineAngleSegments;
                    double curFrom = lineAngleSegments.front().first;
                    double curTo = lineAngleSegments.front().second;
                    for (auto & lineAngle : lineAngleSegments) {
                        if (lineAngle.first <= curTo){
                            curTo = lineAngle.second;
                        } else {
                            if (curTo - curFrom >= M_PI){ // break major arcs
                                mergedLineAngleSegments.push_back(std::make_pair(curFrom, (curFrom + curTo) / 2.0));
                                mergedLineAngleSegments.push_back(std::make_pair((curFrom + curTo) / 2.0, curTo));
                            }else
                                mergedLineAngleSegments.push_back(std::make_pair(curFrom, curTo));
                            curFrom = lineAngle.first;
                            curTo = lineAngle.second;
                        }
                    }
                    if (curTo - curFrom >= M_PI){ // break major arcs
                        mergedLineAngleSegments.push_back(std::make_pair(curFrom, (curFrom + curTo) / 2.0));
                        mergedLineAngleSegments.push_back(std::make_pair((curFrom + curTo) / 2.0, curTo));
                    }
                    else
                        mergedLineAngleSegments.push_back(std::make_pair(curFrom, curTo));


                    // recover lines from angle pairs
                    size_t firstLineIdInThisChain = chainIds.size();
                    for (auto & lineAngle : mergedLineAngleSegments) {
                        Vec3 d1 = Vec3(1, 0, 0) * cos(lineAngle.first) + Vec3(0, 1, 0) * sin(lineAngle.first);
                        Vec3 dd1 = d1(0) * firstPointDirection + d1(1) * firstNormalCrossPoint + d1(2) * firstNormal;
                        Vec3 d2 = Vec3(1, 0, 0) * cos(lineAngle.second) + Vec3(0, 1, 0) * sin(lineAngle.second);
                        Vec3 dd2 = d2(0) * firstPointDirection + d2(1) * firstNormalCrossPoint + d2(2) * firstNormal;
                        Classified<Line3> line;
                        line.claz = claz;
                        line.component.first = dd1;
                        line.component.second = dd2;
                        
                        if (!chainIds.empty() && chainIds.back() == chainId) { // insert a new incidence constraint
                            Vec3 lastEndPoint = mergedLines.back().component.second;
                            if (AngleBetweenDirections(lastEndPoint, line.component.first) <= incidenceAngleThres){
                                Vec3 mid = (lastEndPoint + line.component.first) / 2.0;
                                ViewsNet::ConstraintData cons;
                                cons.mergedSpatialLineSegmentIds[0] = mergedLines.size() - 1;
                                cons.mergedSpatialLineSegmentIds[1] = mergedLines.size();
                                cons.position = mid / norm(mid);
                                cons.type = ViewsNet::ConstraintData::Incidence;
                                cons.weight = 0.0;
                                constraints.push_back(cons);
                            }
                        }                        
                        mergedLines.push_back(line); // insert new merged line
                        chainIds.push_back(chainId); // insert the chainId of this new merged line  
                    }
                    // check whether the head and tail of this chain can form an incidence constraint
                    assert(chainIds[firstLineIdInThisChain] == chainIds.back());
                    Vec3 tailPoint = mergedLines.back().component.second;
                    Vec3 headPoint = mergedLines[firstLineIdInThisChain].component.first;
                    if (AngleBetweenDirections(tailPoint, headPoint) <= incidenceAngleThres){
                        Vec3 mid = (tailPoint + headPoint) / 2.0;
                        ViewsNet::ConstraintData cons;
                        cons.mergedSpatialLineSegmentIds[0] = mergedLines.size() - 1;
                        cons.mergedSpatialLineSegmentIds[1] = firstLineIdInThisChain;
                        cons.position = mid / norm(mid);
                        cons.type = ViewsNet::ConstraintData::Incidence;
                        cons.weight = 0.0;
                        constraints.push_back(cons);
                    }
                }

                return mergedLines;
            }

            // find intersection and incidence constraints
            void AppendIntersectionAndOptionallyIncidenceConstraints(const std::vector<Classified<Line3>> & mergedLines,
                std::vector<ViewsNet::ConstraintData> & constraints, 
                double intersectionAngleThres, double incidenceAngleThres,
                bool appendIncidenceCons = true){
                for (size_t i = 0; i < mergedLines.size(); i++){
                    auto & linei = mergedLines[i];
                    if (linei.claz == -1)
                        continue;
                    auto normali = linei.component.first.cross(linei.component.second);
                    for (size_t j = i + 1; j < mergedLines.size(); j++){
                        auto & linej = mergedLines[j];
                        if (linej.claz == -1)
                            continue;
                        if (linei.claz == linej.claz && !appendIncidenceCons)
                            continue;
                        auto normalj = linej.component.first.cross(linej.component.second);

                        // we MUST ENSURE that all point direction vectors are NORMALIED!!!!!!!!!!!

                        auto nearestPoss = DistanceBetweenTwoLines(linei.component, linej.component).second;
                        double angleDist = AngleBetweenDirections(nearestPoss.first.position, nearestPoss.second.position);

                        if (linei.claz == linej.claz && angleDist <= incidenceAngleThres){ // incidence
                            // the angle distance vertical to lines direction must be lower than the intersectionAngleThres
                            if (std::min(AngleBetweenDirections(normali, normalj), 
                                    AngleBetweenDirections(normali, -normalj)) > 
                                intersectionAngleThres)
                                continue;
                            ViewsNet::ConstraintData cons;
                            cons.mergedSpatialLineSegmentIds[0] = i;
                            cons.mergedSpatialLineSegmentIds[1] = j;
                            cons.position = (nearestPoss.first.position + nearestPoss.second.position) / 2.0;
                            cons.weight = 0.0;
                            cons.type = ViewsNet::ConstraintData::Incidence;
                            constraints.push_back(cons);
                            continue;
                        }
                        else if (angleDist <= intersectionAngleThres){
                            // get intersection point direction
                            Vec3 intersection = normali.cross(normalj); // the direction may be opposited
                            intersection /= norm(intersection);
                            auto nearestPosMid = (nearestPoss.first.position + nearestPoss.second.position) / 2.0;
                            if (AngleBetweenDirections(nearestPosMid, intersection) >
                                AngleBetweenDirections(nearestPosMid, -intersection))
                                intersection = -intersection; // fix the direction

                            ViewsNet::ConstraintData cons;
                            cons.mergedSpatialLineSegmentIds[0] = i;
                            cons.mergedSpatialLineSegmentIds[1] = j;
                            cons.position = intersection;
                            cons.weight = 0.0;
                            cons.type = ViewsNet::ConstraintData::Intersection;
                            constraints.push_back(cons);
                        }
                    }
                }
            }

            // compute manhattan junction weights
            void VoteManhattanJunctionWeightsOnConstraints(const std::vector<Classified<Line3>> & mergedLines,
                const std::array<Vec3, 3> & vanishingPoints,
                std::vector<ViewsNet::ConstraintData> & constraints) {
                // compute line votings
                for (auto & cons : constraints){
                    Vec3 position = cons.position;
                    for (auto & line : mergedLines) {
                        if (line.claz == -1)
                            continue;
                        int vpid = line.claz;
                        const Vec3 & vp = vanishingPoints[vpid];

                        Vec3 p1 = line.component.first;
                        Vec3 p2 = line.component.second;

                        if (vp.cross(position).dot(p1.cross(p2)) < 0)
                            std::swap(p1, p2); // make vp->position and p1->p2 share the same order

                        Vec3 vpPositionPlaneNormal = vp.cross(position);
                        vpPositionPlaneNormal /= norm(vpPositionPlaneNormal);
                        Vec3 p12PlaneNormal = p1.cross(p2);
                        p12PlaneNormal /= norm(p12PlaneNormal);
                        double angle = AngleBetweenDirections(vpPositionPlaneNormal, p12PlaneNormal);
                        if (angle > M_PI / 15.0)
                            continue;  // this position is not on the line (nor its extension), continue

                        double weight = exp(-Square(angle / (M_PI / 30.0)));

                        Line3 spans[2] = { { vp, position }, { position, -vp } };
                        for (size_t i = 0; i < 2; i++){
                            auto proj1 = ProjectionOfPointOnLine(p1, spans[i]);
                            if (proj1.ratio > 1.0)
                                continue;
                            Vec3 lowb = proj1.ratio < 0 ? spans[i].first : p1; // lower bound of the overlapped arc
                            auto proj2 = ProjectionOfPointOnLine(p2, spans[i]);
                            Vec3 highb = proj2.ratio > 1.0 ? spans[i].second : p2; // higher bound of the overlapped arc
                            cons.lineVotings[vpid][i] += AngleBetweenDirections(lowb, highb) * weight;
                        }
                    }
                }

                // compute junction weights
                for (auto & cons : constraints){
                    cons.weight = 0.0;
                    auto & v = cons.lineVotings;                   

                    // compute Triplet junctions
                    double Tp = 0.0;
                    for (int i = 0; i < 3; i++){
                        for (int j = i + 1; j < 3; j++){
                            int k = 3 - i - j;
                            Tp += (v[i][0] + v[i][1]) * (v[j][0] + v[j][1]) * (v[k][0] + v[k][1]);
                        }
                    }
                    cons.junctionWeights.Triplet = Tp;
                    // compute X junction
                    double X = 0.0;
                    for (int i = 0; i < 3; i++){
                        for (int j = i + 1; j < 3; j++){
                            int k = 3 - i - j;
                            X += v[i][0] * v[i][1] * v[j][0] * v[j][1] * DiracDelta(v[k][0] + v[k][1]);
                        }
                    }
                    cons.junctionWeights.X = X;
                    // compute T junction
                    double T = 0.0;
                    for (int i = 0; i < 3; i++){
                        for (int j = 0; j < 3; j++){
                            if (i == j)
                                continue;
                            int k = 3 - i - j;
                            T += v[i][0] * v[i][1] * v[j][0] * DiracDelta(v[j][1] + v[k][0] + v[k][1]);
                            T += v[i][0] * v[i][1] * v[j][1] * DiracDelta(v[j][0] + v[k][0] + v[k][1]);
                        }
                    }
                    cons.junctionWeights.T = T;
                    // compute L junction
                    double L = 0.0;
                    for (int i = 0; i < 3; i++){
                        for (int j = i + 1; j < 3; j++){
                            int k = 3 - i - j;
                            for (int a = 0; a < 2; a++){
                                int nota = 1 - a;
                                for (int b = 0; b < 2; b++){
                                    int notb = 1 - b;
                                    L += v[i][a] * v[j][b] * DiracDelta(v[i][nota] + v[j][notb] + v[k][0] + v[k][1]);
                                }
                            }
                        }
                    }
                    cons.junctionWeights.L = L;
                    // compute I junction
                    double I = 0.0;
                    for (int i = 0; i < 3; i++){
                        for (int j = i + 1; j < 3; j++){
                            int k = 3 - i - j;
                            I += (v[i][0] + v[i][1]) * DiracDelta(v[j][0] + v[j][1] + v[k][0] + v[k][1]);
                        }
                    }
                    cons.junctionWeights.I = I;
                }
            }

            inline bool MaybeVanishingPoint(const Vec3 & p, const std::array<Vec3, 3> & vps, double thres) {
                return AngleBetweenDirections(vps[0], p) < thres ||
                    AngleBetweenDirections(vps[1], p) < thres ||
                    AngleBetweenDirections(vps[2], p) < thres ||
                    AngleBetweenDirections(-vps[0], p) < thres ||
                    AngleBetweenDirections(-vps[1], p) < thres ||
                    AngleBetweenDirections(-vps[2], p) < thres;
            }

            // optimize lines using constraints
            void OptimizeLines(std::vector<Classified<Line3>> & lines, std::vector<ViewsNet::ConstraintData> & constraints, 
                const std::array<Vec3, 3> & vps ) {
                // build equations
                // get all line factors
                struct LineDeterminer {
                    Vec3 firstPointFactor; // \lambda * firstPointFactor = line.first
                    Vec3 secondPointFactor; // \lambda * secondPointFactor = line.second
                    inline Line3 operator()(double lambda) const {
                        return Line3(firstPointFactor * lambda, secondPointFactor * lambda);
                    }
                };
                std::vector<LineDeterminer> lineDeterminers(lines.size());
                for (size_t i = 0; i < lineDeterminers.size(); i++){
                    auto & line = lines[i];
                    if (line.claz == -1)
                        continue;
                    auto vp = vps[line.claz];
                    auto p1 = line.component.first;
                    auto p2 = line.component.second;
                    if (vp.dot(p2 - p1) < 0){
                        vp = -vp; // make p1->p2 and vp direction consistent
                    }
                    lineDeterminers[i].firstPointFactor = p1 / norm(p1);
                    // use the Law of sines
                    double innerAngleOfLine = AngleBetweenDirections(p1, p2);
                    double innerAngleAtFirstPoint = AngleBetweenDirections(vp, -p1); // corresponding to second point factor
                    double innerAngleAtSecondPoint = AngleBetweenDirections(-vp, -p2); // corresponding to first point factor
                    double innerAngleSum = innerAngleOfLine + innerAngleAtFirstPoint + innerAngleAtSecondPoint;
                    assert(FuzzyEquals(innerAngleSum, M_PI, 1e-1));
                    lineDeterminers[i].secondPointFactor = p2 / norm(p2) *
                        sin(innerAngleAtFirstPoint) / sin(innerAngleAtSecondPoint);
                    assert(FuzzyEquals(AngleBetweenDirections(lineDeterminers[i](1.0).direction(), vp), 0.0, 1e-1));
                }

                // build constraint equations
                struct VariableInfo {
                    int varId;
                    bool isSlackVariable;
                    int constraintId; // if is slack variable
                    int lineId; // if is not slack variable
                    double weightInObjectiveFunction;
                };

                /*
                    min 0.5 * x G x + g0 x //G=0
                    s.t.
                        CE^T x + ce0 = 0
                        CI^T x + ci0 >= 0
                */
                /***
                    min \Sum w_ij * [s_ij]
                    s.t.
                        [\lambda_i] * d_ia - [\lambda_j] * d_ja <= s_ij
                        [\lambda_j] * d_ja - [\lambda_i] * d_ia <= s_ij \foreach constraint (i, j)
                        [\lambda_i] >= 1 \foreach line (i)
                ***/
                struct ConstraintInequationInfo { // [\lambda_i] * d_ia - [\lambda_j] * d_ja <= s_ij
                    int equationId;
                    int lambda_iId;
                    double d_ia;
                    int lambda_jId;
                    double d_ja;
                    int s_ijId;
                };
                struct ScaleInequationInfo { // [\lambda_i] >= 1 \foreach line(i)
                    int equationId;
                    int lambdaId; 
                };

                int varIdGenerator = 0;
                int equationIdGenerator = 0;

                std::vector<VariableInfo> lambdas(lines.size()); // \lambda depth variables
                std::vector<ScaleInequationInfo> scaleInequations(lines.size());
                for (size_t i = 0; i < lines.size(); i++){
                    // \lambda_i
                    lambdas[i].varId = varIdGenerator++;
                    lambdas[i].isSlackVariable = false;
                    lambdas[i].lineId = i;
                    lambdas[i].weightInObjectiveFunction = 0.0; // not appeared in objective function
                    scaleInequations[i].equationId = equationIdGenerator++;
                    scaleInequations[i].lambdaId = lambdas[i].varId;
                }

                std::vector<VariableInfo> slacks; // slack variables
                std::vector<ConstraintInequationInfo> constraintInequations;
                slacks.reserve(constraints.size());
                constraintInequations.reserve(constraints.size() * 3);
                for (size_t i = 0; i < constraints.size(); i++){
                    auto & cons = constraints[i];
                    bool maybeVanishingPoint =
                        MaybeVanishingPoint(cons.position, vps, M_PI / 100.0);
                    if (maybeVanishingPoint)
                        continue;

                    const auto & line1 = lines[cons.mergedSpatialLineSegmentIds[0]];
                    const auto & line2 = lines[cons.mergedSpatialLineSegmentIds[1]];
                    static const int nearestPointPairsTable[4][2] = {
                        {0, 0},
                        {0, 1},
                        {1, 0},
                        {1, 1}
                    };
                    int minId = -1;
                    double minAngle = std::numeric_limits<double>::max();
                    for (int j = 0; j < 4; j++){
                        double angle = AngleBetweenDirections(nearestPointPairsTable[j][0] == 0 ?
                            line1.component.first : line1.component.second,
                            nearestPointPairsTable[j][1] == 0 ?
                            line2.component.first : line2.component.second);
                        if (angle < minAngle){
                            minAngle = angle;
                            minId = j;
                        }
                    }
                    assert(minId != -1);
                    auto nearestPointPair = nearestPointPairsTable[minId];
                    const auto & lineDet1 = lineDeterminers[cons.mergedSpatialLineSegmentIds[0]];
                    const auto & lineDet2 = lineDeterminers[cons.mergedSpatialLineSegmentIds[1]];
                    Vec3 di = nearestPointPair[0] == 0 ? lineDet1.firstPointFactor : lineDet1.secondPointFactor;
                    Vec3 dj = nearestPointPair[1] == 0 ? lineDet2.firstPointFactor : lineDet2.secondPointFactor;
                    Vec3 ddi(di.dot(vps[0]), di.dot(vps[1]), di.dot(vps[2]));
                    Vec3 ddj(dj.dot(vps[0]), dj.dot(vps[1]), dj.dot(vps[2]));

                    // slack variable
                    VariableInfo slackVar;
                    slackVar.varId = varIdGenerator++;
                    slackVar.isSlackVariable = true;
                    slackVar.constraintId = i;
                    slackVar.weightInObjectiveFunction = cons.weight; // weight of constraint
                    slacks.push_back(slackVar);                    
                    
                    if (cons.type == ViewsNet::ConstraintData::Intersection){
                        // [\lambda_i] * d_ia - [\lambda_j] * d_ja <= s_ij
                        for (int k = 0; k < 3; k++){
                            ConstraintInequationInfo eq;
                            eq.equationId = equationIdGenerator++;
                            eq.lambda_iId = lambdas[cons.mergedSpatialLineSegmentIds[0]].varId;
                            eq.d_ia = ddi[k];
                            eq.lambda_jId = lambdas[cons.mergedSpatialLineSegmentIds[1]].varId;
                            eq.d_ja = ddj[k];
                            eq.s_ijId = slacks.back().varId;
                            constraintInequations.push_back(eq);
                        }
                        // [\lambda_j] * d_ja - [\lambda_i] * d_ia <= s_ij
                        for (int k = 0; k < 3; k++){
                            ConstraintInequationInfo eq;
                            eq.equationId = equationIdGenerator++;
                            eq.lambda_jId = lambdas[cons.mergedSpatialLineSegmentIds[0]].varId;
                            eq.d_ja = ddi[k];
                            eq.lambda_iId = lambdas[cons.mergedSpatialLineSegmentIds[1]].varId;
                            eq.d_ia = ddj[k];
                            eq.s_ijId = slacks.back().varId;
                            constraintInequations.push_back(eq);
                        }
                    }
                    else {
                        // [\lambda_i] * d_ia - [\lambda_j] * d_ja <= s_ij
                        for (int k = 0; k < 3; k++){
                            if (k == line1.claz)
                                continue;
                            ConstraintInequationInfo eq;
                            eq.equationId = equationIdGenerator++;
                            eq.lambda_iId = lambdas[cons.mergedSpatialLineSegmentIds[0]].varId;
                            eq.d_ia = ddi[k];
                            eq.lambda_jId = lambdas[cons.mergedSpatialLineSegmentIds[1]].varId;
                            eq.d_ja = ddj[k];
                            eq.s_ijId = slacks.back().varId;
                            constraintInequations.push_back(eq);
                        }
                        // [\lambda_j] * d_ja - [\lambda_i] * d_ia <= s_ij
                        for (int k = 0; k < 3; k++){
                            if (k == line1.claz)
                                continue;
                            ConstraintInequationInfo eq;
                            eq.equationId = equationIdGenerator++;
                            eq.lambda_jId = lambdas[cons.mergedSpatialLineSegmentIds[0]].varId;
                            eq.d_ja = ddi[k];
                            eq.lambda_iId = lambdas[cons.mergedSpatialLineSegmentIds[1]].varId;
                            eq.d_ia = ddj[k];
                            eq.s_ijId = slacks.back().varId;
                            constraintInequations.push_back(eq);
                        }
                    }
                }

                /***
                    min \Sum w_ij * [s_ij]
                    s.t.
                        [\lambda_i] * d_ia - [\lambda_j] * d_ja <= s_ij
                        [\lambda_j] * d_ja - [\lambda_i] * d_ia <= s_ij \foreach constraint (i, j)
                        [\lambda_i] >= 1 \foreach line (i)
                ***/

                int varNum = varIdGenerator;
                int equationNum = equationIdGenerator;

                sinfo *info;
                info = (sinfo*)malloc(sizeof(sinfo));
                info->env = (jmp_buf *)malloc(sizeof(jmp_buf));
                info->text = "This information was passed to the hook function.";
                setjmp(*(info->env));
                glp_error_hook(glpErrorHook, info);


                glp_prob * problem = glp_create_prob();
                glp_set_prob_name(problem, "Optimize Lines");
                glp_set_obj_name(problem, "Optimize Lines: Objective");
                glp_set_obj_dir(problem, GLP_MIN);

                glp_add_rows(problem, equationNum+1); // add rows
                glp_add_cols(problem, varNum+1); // add cols

                // add scale inequations
                for (auto & var : lambdas){
                    glp_set_col_bnds(problem, var.varId + 1, GLP_LO, 3.0, 1e5);
                }
                for (auto & var : slacks){
                    glp_set_col_bnds(problem, var.varId + 1, GLP_LO, 0.0, 1e5);
                }

                // add constrant iequations
                for (auto & ie : constraintInequations){
                    glp_set_row_bnds(problem, ie.equationId+1, GLP_UP, -1e5, 0.0); // ... <= 0.0
                }
                for (auto & ie : constraintInequations){
                    int varIds[] = { -1, ie.lambda_iId+1, ie.lambda_jId+1, ie.s_ijId+1 };
                    //std::cout << "constraint: lambda_i-" << (ie.lambda_iId + 1) << 
                    //    "  lambda_j-" << (ie.lambda_jId + 1) << 
                    //    " s_ij-" << (ie.s_ijId + 1) << std::endl;
                    double coefs[] = { 0.0, ie.d_ia, -ie.d_ja, -1.0 };
                    if (ie.lambda_iId == ie.lambda_jId || ie.lambda_iId == ie.s_ijId || ie.lambda_jId == ie.s_ijId){
                        std::cout << "f" << std::endl;
                    }
                    glp_set_mat_row(problem, ie.equationId + 1, 3, varIds, coefs);
                }

                // set objective coefs
                // only slacks appear in objective function
                for (auto & var : slacks){
                    glp_set_obj_coef(problem, var.varId+1, var.weightInObjectiveFunction);
                }
                
                glp_adv_basis(problem, 0);

                bool useSimplex = true;
                if (useSimplex){
                    glp_smcp params;
                    glp_init_smcp(&params);
                    params.msg_lev = GLP_MSG_ON;
                    int stat = glp_simplex(problem, &params);

                    // retrieve results
                    for (auto & var : lambdas){
                        double lambda = glp_get_col_prim(problem, var.varId+1);
                        // update line coordinates
                        lines[var.lineId].component = lineDeterminers[var.lineId](lambda);
                    }
                    for (auto & var : slacks){
                        double s = glp_get_col_prim(problem, var.varId+1);
                        // set the slack value for constraint
                        constraints[var.constraintId].slackValue = s;
                    }
                }
                else{
                    glp_iptcp params;
                    glp_init_iptcp(&params);
                    params.msg_lev = GLP_MSG_ON;
                    int stat = glp_interior(problem, &params);

                    // retrieve results
                    for (auto & var : lambdas){
                        double lambda = glp_ipt_col_prim(problem, var.varId+1);
                        // update line coordinates
                        lines[var.lineId].component = lineDeterminers[var.lineId](lambda);
                    }
                    for (auto & var : slacks){
                        double s = glp_ipt_col_prim(problem, var.varId+1);
                        // set the slack value for constraint
                        if (isnan(s))
                            std::cout << "slack value of constraint " << var.constraintId << "is NaN" << std::endl;
                        constraints[var.constraintId].slackValue = isnan(s) ? 1e10 : s;
                    }
                }

                glp_delete_prob(problem);


                glp_error_hook(NULL, NULL);
                free(info->env);
                free(info);
                
            }


        }


        void ViewsNet::rectifySpatialLines() {
            _globalData.constraints.reserve(Square(_globalData.mergedSpatialLineSegments.size()) / 4);
            _globalData.constraints.clear();

            // merge lines, and get (some) incidence constraints
            _globalData.mergedSpatialLineSegments = 
                MergeColinearSpatialLinesAndAppendIncidenceConstraints(
                _globalData.spatialLineSegments,
                _globalData.mergedSpatialLineSegmentChainIds, 
                _globalData.constraints, 
                _params.mergeLineDistanceAngleThreshold,
                _params.incidenceConstraintLineDistanceAngleThreshold);

            // make all point directions of lines normalized
            for (auto & line : _globalData.mergedSpatialLineSegments){
                line.component.first /= norm(line.component.first);
                line.component.second /= norm(line.component.second);
            }

            // get all intersection and incidence constraints
            AppendIntersectionAndOptionallyIncidenceConstraints(_globalData.mergedSpatialLineSegments,
                _globalData.constraints, 
                _params.intersectionConstraintLineDistanceAngleThreshold, 
                _params.incidenceConstraintLineDistanceAngleThreshold, false);

            // remove all duplicated constraints
            // remove all self connected lines
            // remove all constraints close to any vanishing points
            auto uniqueConsIters = MergeNear(_globalData.constraints.begin(), _globalData.constraints.end(), std::false_type(),
                1, [](const ConstraintData & cons1, const ConstraintData & cons2){
                return (cons1.type == cons2.type && 
                    std::is_permutation(std::begin(cons1.mergedSpatialLineSegmentIds), std::end(cons1.mergedSpatialLineSegmentIds), 
                    std::begin(cons2.mergedSpatialLineSegmentIds))) ? 0 : 2;
            });
            std::vector<ConstraintData> uniqueCons;
            for (auto consIter : uniqueConsIters){
                auto & cons = *consIter;
                if (cons.mergedSpatialLineSegmentIds[0] == cons.mergedSpatialLineSegmentIds[1] ||
                    MaybeVanishingPoint(cons.position, _globalData.vanishingPoints, 
                    _params.intersectionConstraintLineDistanceAngleThreshold))
                    continue;
                uniqueCons.push_back(cons);
            }
            _globalData.constraints = uniqueCons;
            
            // vote for the position of each constraint
            VoteManhattanJunctionWeightsOnConstraints(_globalData.mergedSpatialLineSegments,
                _globalData.vanishingPoints, _globalData.constraints);

            // compute weights of contraints 
            for (auto & cons : _globalData.constraints){
                cons.weight = cons.junctionWeights.Triplet * _params.mjWeightTriplet + 
                    cons.junctionWeights.T * _params.mjWeightT + 
                    cons.junctionWeights.X * _params.mjWeightX + 
                    cons.junctionWeights.L * _params.mjWeightL +
                    cons.junctionWeights.I * _params.mjWeightI;
                //std::cout << "cons weight: " << cons.weight << std::endl;
                assert(cons.weight >= 0);
            }           
            std::sort(_globalData.constraints.begin(), _globalData.constraints.end(),
                [](const ConstraintData & cons1, const ConstraintData & cons2){
                return cons1.weight > cons2.weight;
            });

            std::cout << "line num: " << _globalData.mergedSpatialLineSegments.size() << std::endl;
            std::cout << "constraint num: " << _globalData.constraints.size() << std::endl;

            // optimize lines
            OptimizeLines(_globalData.mergedSpatialLineSegments, 
                _globalData.constraints, _globalData.vanishingPoints);

            // find necessary constraints using MST with slackValues
            std::vector<size_t> lineIds(_globalData.mergedSpatialLineSegments.size()), 
                consIds(_globalData.constraints.size());
            for (size_t i = 0; i < lineIds.size(); i++)
                lineIds[i] = i;
            for (size_t i = 0; i < consIds.size(); i++)
                consIds[i] = i;
            std::vector<size_t> MSTconsIds;
            MSTconsIds.reserve(_globalData.constraints.size());
            MinimumSpanningTree(lineIds.begin(), lineIds.end(), consIds.begin(), consIds.end(),
                std::back_inserter(MSTconsIds),
                [this](size_t e){return std::make_pair(_globalData.constraints[e].mergedSpatialLineSegmentIds[0], 
                    _globalData.constraints[e].mergedSpatialLineSegmentIds[1]); },
                [this](size_t e){return _globalData.constraints[e].slackValue; }
            );
            std::vector<ConstraintData> MSTconstraints(MSTconsIds.size());
            for (size_t i = 0; i < MSTconsIds.size(); i++){
                MSTconstraints[i] = _globalData.constraints[MSTconsIds[i]];
            }
            std::cout << "line num: " << _globalData.mergedSpatialLineSegments.size() << std::endl;
            std::cout << "mst constraint num: " << MSTconstraints.size() << std::endl;
            _globalData.refinedConstraints = MSTconstraints;

            // optimize lines again
            OptimizeLines(_globalData.mergedSpatialLineSegments,
                _globalData.refinedConstraints, _globalData.vanishingPoints);
        }
 


    }
}