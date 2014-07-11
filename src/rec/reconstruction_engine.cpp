#include <filesystem>

extern "C" {
    #include <gpc.h>
}

#include <Eigen/StdVector>
#include <dlib/matrix.h>
#include <dlib/optimization.h>

#include "../vis/visualize2d.hpp"
#include "../vis/visualize3d.hpp"

#include "optimization.hpp"
#include "reconstruction_engine.hpp"

namespace panoramix {
    namespace rec {

        ReconstructionEngine::Params::Params() 
            : camera(250.0), lineSegmentWeight(1.0), siftWeight(1.0),
            surfWeight(1.0), cameraAngleScaler(1.8), smallCameraAngleScalar(0.05),
            intersectionConstraintLineDistanceAngleThreshold(0.06),
            incidenceConstraintLineDistanceAngleThreshold(0.2),
            mergeLineDistanceAngleThreshold(0.05),
            mjWeightTriplet(5.0), mjWeightX(5.0), mjWeightT(2.0), mjWeightL(1.0), mjWeightI(2.0) {
        }

        ReconstructionEngine::ComponentData::ComponentData(ReconstructionEngine::ComponentData::Type t) : type(t) {}

        ReconstructionEngine::LineStructureConnectivityConstraintData::LineStructureConnectivityConstraintData() {
            mergedSpatialLineSegmentIds[0] =
                mergedSpatialLineSegmentIds[1] = 0;
            weight = 0.0;
            std::memset(lineVotings, 0, sizeof(lineVotings));
            junctionWeights.I =
                junctionWeights.L =
                junctionWeights.T =
                junctionWeights.Triplet =
                junctionWeights.X = 0.0;
        }

        ReconstructionEngine::ConstraintData::ConstraintData(ReconstructionEngine::ConstraintData::Type t) : type(t) {}


        ReconstructionEngine::ViewHandle ReconstructionEngine::insertPhoto(const Image & im, const PerspectiveCamera & cam,
            double cameraDirectionErrorScale) {
            ViewData vd;
            vd.camera = vd.originalCamera = cam;
            vd.cameraDirectionErrorScale = cameraDirectionErrorScale;
            vd.image = im;
            return insertView(vd);
        }

        void ReconstructionEngine::insertPanorama(const Image & panorama, const std::vector<PerspectiveCamera> & viewCams,
            const PanoramicCamera & panCam) {
            for (int i = 0; i < viewCams.size(); i++) {
                auto & camera = viewCams[i];
                const auto im =
                    core::CameraSampler<core::PerspectiveCamera, core::PanoramicCamera>(camera, panCam)(panorama);
                auto viewHandle = insertPhoto(im, camera);
                updateConnections(viewHandle);
            }
            _globalData.panorama = panorama;
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

        void ReconstructionEngine::computeFeatures(ViewHandle h) {
            auto & vd = _views.data(h);
            const Image & im = vd.image;
            auto lineSegments = _params.lineSegmentExtractor(im);
            vd.lineSegments.resize(lineSegments.size());
            for (size_t i = 0; i < lineSegments.size(); i++){
                vd.lineSegments[i].claz = -1;
                vd.lineSegments[i].component = lineSegments[i];
            }

            // compute line intersections
            vd.lineSegmentIntersections.clear();
            vd.lineSegmentIntersectionLineIDs.clear();
            LineIntersectons(vd.lineSegments, 
                vd.lineSegmentIntersections, 
                vd.lineSegmentIntersectionLineIDs, true);
            
            // collect keypoints/decriptors for matching
            vd.keypointsForMatching = _params.surfExtractor(im, cv::Mat(), vd.descriptorsForMatching);
            // and build RTree for the keypoints
            /*vd.lineSegmentIntersectionsRTree =
                RTreeWrapper<HPoint2>(vd.lineSegmentIntersections.begin(), vd.lineSegmentIntersections.end());
            vd.keypointsForMatchingRTree = 
                RTreeWrapper<KeyPoint>(vd.keypointsForMatching.begin(), vd.keypointsForMatching.end());*/
        }

        void ReconstructionEngine::buildRegionNet(ViewHandle h) {
            auto & vd = _views.data(h);
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

        size_t ReconstructionEngine::updateConnections(ViewHandle h) {
            auto & thisv = _views.data(h);
            const PerspectiveCamera & thisvCam = thisv.originalCamera;
            double thisvCamAngleRadius = PerspectiveCameraAngleRadius(thisvCam);
            thisvCamAngleRadius *= _params.cameraAngleScaler;
            for (auto & v : _views.elements<0>()){
                if (v.topo.hd == h)
                    continue;
                const PerspectiveCamera & vcam = v.data.originalCamera;
                double vCamAngleRadius = PerspectiveCameraAngleRadius(vcam);
                vCamAngleRadius *= _params.cameraAngleScaler;
                double angleDistance = AngleBetweenDirections(thisvCam.center(), vcam.center());
                if (angleDistance <= thisvCamAngleRadius + vCamAngleRadius){
                    // may overlap
                    ViewConnectionData hd;
                    //hd.cameraAngleDistance = angleDistance;
                    _views.add<1>({ h, v.topo.hd }, hd);
                }
            }
            return _views.topo(h).uppers.size();
        }

        ReconstructionEngine::ViewHandle ReconstructionEngine::isTooCloseToAnyExistingView(ViewHandle h) const {
            auto & camera = _views.data(h).camera;
            double cameraRadius = PerspectiveCameraAngleRadius(camera);
            auto connections = _views.topo(h).uppers;
            for (auto & con : connections){
                auto to = _views.topo(con).lowers[0];
                if (to == h)
                    to = _views.topo(con).lowers[1];
                auto & neighborCamera = _views.data(to).camera;
                double cameraAngle = AngleBetweenDirections(camera.center(), neighborCamera.center());
                double neighborCameraRadius = PerspectiveCameraAngleRadius(camera);
                if (cameraAngle <= (cameraRadius + neighborCameraRadius) * _params.smallCameraAngleScalar){ // too close
                    return to;
                }
            }
            return ViewHandle();
        }

        void ReconstructionEngine::findMatchesToConnectedViews(ViewHandle h) {
            cv::detail::BestOf2NearestMatcher matcher(false, 0.3f);
            auto & connections = _views.topo(h).uppers;
            for (auto & con : connections) {
                auto & conData = _views.data(con);
                auto & thisVD = _views.data(_views.topo(con).lowers[0]);
                auto & neighborVD = _views.data(_views.topo(con).lowers[1]);
                
                // find matches and estimate homography using opencv
                cv::detail::ImageFeatures thisFea;
                thisFea.descriptors = thisVD.descriptorsForMatching;
                thisFea.keypoints = thisVD.keypointsForMatching;

                cv::detail::ImageFeatures neighborFea;
                neighborFea.descriptors = neighborVD.descriptorsForMatching;
                neighborFea.keypoints = neighborVD.keypointsForMatching;

                matcher(thisFea, neighborFea, conData.matchInfo);
                conData.matchInfo.src_img_idx = _views.topo(con).lowers[0].id;
                conData.matchInfo.dst_img_idx = _views.topo(con).lowers[1].id;
            }
        }

        namespace {

            void DecodeLookAtMatGradient(const Eigen::MatrixXd & m4, Vec3 & deye, Vec3 & dcenter, Vec3 & dup) {

            }

        }

        void ReconstructionEngine::calibrateAllCameras() {
            
            //using namespace deriv;
            //
            //ExpressionGraph graph;
            //std::vector<Expression<Eigen::MatrixXd>> cameraViewMats(_views.internalElements<0>().size());
            ////std::vector<Eigen::MatrixXd> 

            //for (auto & v : _views.elements<0>()){
            //    cameraViewMats[v.topo.hd.id] = deriv::composeFunction(graph, 
            //        [&v](){
            //        return deriv::CVMatToEigenMatX(v.data.camera.viewMatrix()); 
            //    });
            //}

            //for (auto & c : _views.elements<1>()){

            //}

            NOT_IMPLEMENTED_YET();

        }


        void ReconstructionEngine::stitchPanorama() {
            NOT_IMPLEMENTED_YET();
        }











        namespace {

            inline PixelLoc PixelIndexFromGeoCoord(const GeoCoord & p, int longidiv, int latidiv) {
                int longtid = static_cast<int>((p.longitude + M_PI) * longidiv / M_PI / 2);
                int latid = static_cast<int>((p.latitude + M_PI_2) * latidiv / M_PI);
                longtid = (longtid % longidiv + longidiv) % longidiv;
                latid = (latid % latidiv + latidiv) % latidiv;
                return PixelLoc(longtid, latid);
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
                    PixelLoc pixel = PixelIndexFromGeoCoord(GeoCoord(p), longitudeDivideNum, latitudeDivideNum);
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
                        PixelLoc pixel = PixelIndexFromGeoCoord(GeoCoord(v), longitudeDivideNum, latitudeDivideNum);
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
                            PixelLoc pixel = PixelIndexFromGeoCoord(GeoCoord(v), longitudeDivideNum, latitudeDivideNum);
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
                    double curscore = 0.8;
                    for (int j = 0; j < npoints; j++){
                        if (linescores[j] > curscore){
                            lines[i].claz = j;
                            curscore = linescores[j];
                        }
                    }
                }
            }

            inline Vec3 RotateDirectionTo(const Vec3 & originalDirection, const Vec3 & toDirection, double angle) {
                Vec3 tovec = originalDirection.cross(toDirection).cross(originalDirection);
                Vec3 result3 = originalDirection + tovec * tan(angle);
                return result3 / norm(result3);
            }
            
        }


        void ReconstructionEngine::estimateVanishingPointsAndClassifyLines() {

            // pick separated views only
            std::vector<decltype(_views.elements<0>().begin())> seperatedViewIters;
            MergeNearNaive(_views.elements<0>().begin(), _views.elements<0>().end(), std::back_inserter(seperatedViewIters), std::false_type(),
                _params.smallCameraAngleScalar,
                [](const ViewsGraph::TripletType<0> & v1, const ViewsGraph::TripletType<0> & v2) -> double{
                double angleDistance = AngleBetweenDirections(v1.data.camera.center(), v2.data.camera.center());
                return angleDistance / 
                    (PerspectiveCameraAngleRadius(v1.data.camera) + PerspectiveCameraAngleRadius(v2.data.camera));
            });

            // collect line intersections
            size_t lineIntersectionsNum = 0;
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
                    return p3 / norm(p3); // normalized
                });
            }

            // get merged intersections
            // normalize spatial intersections
            auto intersectionValidEnd = std::remove_if(intersections.begin(), intersections.end(),
                [](const Vec3 & v){return std::isnan(v[0]) || std::isnan(v[1]) || std::isnan(v[2]); });
            std::vector<decltype(intersections.begin())> mergedIntersectionsIters;
            mergedIntersectionsIters.reserve(intersections.size() / 2);
            MergeNearRTree(intersections.begin(), intersectionValidEnd, std::back_inserter(mergedIntersectionsIters), std::false_type(),
                2 * sin(M_PI / 150.0 / 2.0));
            _globalData.mergedSpatialLineSegmentIntersections.clear();
            _globalData.mergedSpatialLineSegmentIntersections.reserve(mergedIntersectionsIters.size());
            for (auto intersectionIter : mergedIntersectionsIters) {
                _globalData.mergedSpatialLineSegmentIntersections.push_back(*intersectionIter);
            }


            // find vanishing points;
            _globalData.vanishingPoints = FindVanishingPoints(intersections);                        

            // add spatial line segments from line segments of all views
            size_t spatialLineSegmentsNum = 0;
            for (auto & v : _views.elements<0>())
                spatialLineSegmentsNum += v.data.lineSegments.size();
            _globalData.spatialLineSegments.resize(spatialLineSegmentsNum);
            auto spatialLineSegmentBegin = _globalData.spatialLineSegments.begin();
            for (auto & v : _views.elements<0>()){
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
            for (auto & v : _views.elements<0>()){
                for (auto & line : v.data.lineSegments){
                    line.claz = spatialLineSegmentBegin->claz;
                    ++ spatialLineSegmentBegin;
                }
            }
        }

















        namespace {

            using LineStructureConnectivityConstraintData = ReconstructionEngine::LineStructureConnectivityConstraintData;


            // merge colinear spatial lines and find some incidence constraints
            std::vector<Classified<Line3>> MergeColinearSpatialLinesAndAppendIncidenceConstraints(
                const std::vector<Classified<Line3>> & oldLines, 
                std::vector<int> & chainIds, std::vector<LineStructureConnectivityConstraintData> & constraints, 
                double mergeAngleThres, double incidenceAngleThres) {
                
                std::vector<Classified<Line3>> lines(oldLines);

                // normalize line point directions
                for (auto & line : lines) {
                    line.component.first /= norm(line.component.first);
                    line.component.second /= norm(line.component.second);
                }

                // group all colinear spatial lines
                std::vector<decltype(lines.begin())> colinearLineIters;
                MergeNearNaive(lines.begin(), lines.end(), std::back_inserter(colinearLineIters), std::true_type(), 
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
                                LineStructureConnectivityConstraintData cons;
                                cons.mergedSpatialLineSegmentIds[0] = mergedLines.size() - 1;
                                cons.mergedSpatialLineSegmentIds[1] = mergedLines.size();
                                cons.position = mid / norm(mid);
                                cons.type = LineStructureConnectivityConstraintData::Incidence;
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
                        LineStructureConnectivityConstraintData cons;
                        cons.mergedSpatialLineSegmentIds[0] = mergedLines.size() - 1;
                        cons.mergedSpatialLineSegmentIds[1] = firstLineIdInThisChain;
                        cons.position = mid / norm(mid);
                        cons.type = LineStructureConnectivityConstraintData::Incidence;
                        cons.weight = 0.0;
                        constraints.push_back(cons);
                    }
                }

                return mergedLines;
            }

            // find intersection and incidence constraints
            void AppendIntersectionAndOptionallyIncidenceConstraints(const std::vector<Classified<Line3>> & mergedLines,
                std::vector<LineStructureConnectivityConstraintData> & constraints, 
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
                            LineStructureConnectivityConstraintData cons;
                            cons.mergedSpatialLineSegmentIds[0] = i;
                            cons.mergedSpatialLineSegmentIds[1] = j;
                            cons.position = (nearestPoss.first.position + nearestPoss.second.position) / 2.0;
                            cons.weight = 0.0;
                            cons.type = LineStructureConnectivityConstraintData::Incidence;
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

                            LineStructureConnectivityConstraintData cons;
                            cons.mergedSpatialLineSegmentIds[0] = i;
                            cons.mergedSpatialLineSegmentIds[1] = j;
                            cons.position = intersection;
                            cons.weight = 0.0;
                            cons.type = LineStructureConnectivityConstraintData::Intersection;
                            constraints.push_back(cons);
                        }
                    }
                }
            }

            // compute manhattan junction weights
            void VoteManhattanJunctionWeightsOnConstraints(const std::vector<Classified<Line3>> & mergedLines,
                const std::array<Vec3, 3> & vanishingPoints,
                std::vector<LineStructureConnectivityConstraintData> & constraints) {
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
            void OptimizeLines(std::vector<Classified<Line3>> & lines, std::vector<LineStructureConnectivityConstraintData> & constraints, 
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
                    int lineId; // if is not slack variable (the lambda variable)
                    double weightInObjectiveFunction;
                };


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
                    lambdas[i].lineId = static_cast<int>(i);
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
                    Vec3 ddi(di.dot(vps[0] / norm(vps[0])), di.dot(vps[1] / norm(vps[1])), di.dot(vps[2] / norm(vps[2])));
                    Vec3 ddj(dj.dot(vps[0] / norm(vps[0])), dj.dot(vps[1] / norm(vps[1])), dj.dot(vps[2] / norm(vps[2])));

                    // slack variable
                    VariableInfo slackVar;
                    slackVar.varId = varIdGenerator++;
                    slackVar.isSlackVariable = true;
                    slackVar.constraintId = static_cast<int>(i);
                    slackVar.weightInObjectiveFunction = cons.weight; // weight of constraint
                    slacks.push_back(slackVar);                    
                    
                    if (cons.type == LineStructureConnectivityConstraintData::Intersection){
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

                glp_add_rows(problem, equationNum); // add rows
                glp_add_cols(problem, varNum); // add cols

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

                enum Method {Simplex, Exact, Interior};
                Method method = Simplex;
                if (method == Simplex){
                    glp_smcp params;
                    glp_init_smcp(&params);
                    params.tol_bnd = 1e-10;
                    params.tol_dj = 1e-10;
                    params.tol_piv = 1e-13;
                    params.presolve = GLP_ON;
                    params.meth = GLP_DUALP;
                    //params.r_test = GLP_RT_STD;

                    params.msg_lev = GLP_MSG_ON;
                    int stat = glp_simplex(problem, &params);
                    switch (stat) {
                    case 0: std::cout << "solved" << std::endl; break;
                    case GLP_EBADB: std::cout << "GLP_EBADB" << std::endl; break;
                    case GLP_ESING: std::cout << "GLP_ESING" << std::endl; break;
                    case GLP_ECOND: std::cout << "GLP_ECOND" << std::endl; break;
                    case GLP_EBOUND: std::cout << "GLP_EBOUND" << std::endl; break;
                    case GLP_EFAIL: std::cout << "GLP_EFAIL" << std::endl; break;
                    case GLP_EOBJLL: std::cout << "GLP_EOBJLL" << std::endl; break;
                    case GLP_EOBJUL: std::cout << "GLP_EOBJUL" << std::endl; break;
                    case GLP_EITLIM: std::cout << "GLP_EITLIM" << std::endl; break;
                    case GLP_ETMLIM: std::cout << "GLP_ETMLIM" << std::endl; break;
                    case GLP_ENOPFS: std::cout << "GLP_ENOPFS" << std::endl; break;
                    case GLP_ENODFS: std::cout << "GLP_ENODFS" << std::endl; break;
                    }

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
                else if (method == Exact) {
                    glp_smcp params;
                    glp_init_smcp(&params);

                    params.msg_lev = GLP_MSG_ON;
                    int stat = glp_exact(problem, &params);
                    switch (stat) {
                    case 0: std::cout << "solved" << std::endl; break;
                    case GLP_EBADB: std::cout << "GLP_EBADB" << std::endl; break;
                    case GLP_ESING: std::cout << "GLP_ESING" << std::endl; break;
                    case GLP_ECOND: std::cout << "GLP_ECOND" << std::endl; break;
                    case GLP_EBOUND: std::cout << "GLP_EBOUND" << std::endl; break;
                    case GLP_EFAIL: std::cout << "GLP_EFAIL" << std::endl; break;
                    case GLP_EOBJLL: std::cout << "GLP_EOBJLL" << std::endl; break;
                    case GLP_EOBJUL: std::cout << "GLP_EOBJUL" << std::endl; break;
                    case GLP_EITLIM: std::cout << "GLP_EITLIM" << std::endl; break;
                    case GLP_ETMLIM: std::cout << "GLP_ETMLIM" << std::endl; break;
                    case GLP_ENOPFS: std::cout << "GLP_ENOPFS" << std::endl; break;
                    case GLP_ENODFS: std::cout << "GLP_ENODFS" << std::endl; break;
                    }

                    // retrieve results
                    for (auto & var : lambdas) {
                        double lambda = glp_get_col_prim(problem, var.varId + 1);
                        // update line coordinates
                        lines[var.lineId].component = lineDeterminers[var.lineId](lambda);
                    }
                    for (auto & var : slacks) {
                        double s = glp_get_col_prim(problem, var.varId + 1);
                        // set the slack value for constraint
                        constraints[var.constraintId].slackValue = s;
                    }
                }
                else if(method == Interior){
                    glp_iptcp params;
                    
                    glp_init_iptcp(&params);
                    params.msg_lev = GLP_MSG_ON;
                    int stat = glp_interior(problem, &params);
                    switch (stat) {
                    case 0: std::cout << "solved" << std::endl; break;
                    case GLP_EFAIL: std::cout << "GLP_EFAIL" << std::endl; break;
                    case GLP_ENOCVG: std::cout << "GLP_ENOCVG" << std::endl; break;
                    case GLP_EITLIM: std::cout << "GLP_EITLIM" << std::endl; break;
                    case GLP_EINSTAB: std::cout << "GLP_EINSTAB" << std::endl; break;
                    }

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


            // optimize lines using constraints based on expressions
            void OptimizeLinesUsingConstraintGraph(std::vector<Classified<Line3>> & lines, 
                std::vector<LineStructureConnectivityConstraintData> & constraints,
                const std::array<Vec3, 3> & vps) {
                using namespace deriv;
                using namespace Eigen;

                using ConsGraph = ConstraintGraph<int, int>;
                ConsGraph consGraph;
                for (int i = 0; i < lines.size(); i++)
                    consGraph.addComponent(i);
                for (int i = 0; i < constraints.size(); i++)
                    consGraph.addConstraint({ 
                    ConsGraph::ComponentHandle(constraints[i].mergedSpatialLineSegmentIds[0]), 
                    ConsGraph::ComponentHandle(constraints[i].mergedSpatialLineSegmentIds[1]) 
                }, i);
                
                ExpressionGraph graph;

                // build equations
                // get all line factors
                MatrixX3d firstPointFactorValues = MatrixX3d::Zero(lines.size(), 3);
                MatrixX3d secondPointFactorValues = MatrixX3d::Zero(lines.size(), 3);
                auto firstPointFactors = graph.addRef(firstPointFactorValues, "firstPointFactors");
                auto secondPointFactors = graph.addRef(secondPointFactorValues, "secondPointFactors");

                for (size_t i = 0; i < lines.size(); i++){
                    auto & line = lines[i];
                    if (line.claz == -1)
                        continue;
                    auto vp = vps[line.claz];
                    auto p1 = line.component.first;
                    auto p2 = line.component.second;
                    if (vp.dot(p2 - p1) < 0){
                        vp = -vp; // make p1->p2 and vp direction consistent
                    }
                    firstPointFactorValues.row(i) = CVMatToEigenMat(p1 / norm(p1));
                    // use the Law of sines
                    double innerAngleOfLine = AngleBetweenDirections(p1, p2);
                    double innerAngleAtFirstPoint = AngleBetweenDirections(vp, -p1); // corresponding to second point factor
                    double innerAngleAtSecondPoint = AngleBetweenDirections(-vp, -p2); // corresponding to first point factor
                    double innerAngleSum = innerAngleOfLine + innerAngleAtFirstPoint + innerAngleAtSecondPoint;
                    assert(FuzzyEquals(innerAngleSum, M_PI, 1e-1));
                    secondPointFactorValues.row(i) = CVMatToEigenMat(p2 / norm(p2) *
                        sin(innerAngleAtFirstPoint) / sin(innerAngleAtSecondPoint));
                }


                VectorXd lambdaValues = VectorXd::Ones(lines.size());               
                auto lambdas = graph.addRef(lambdaValues, "lambdas").assign<VectorXd>();
                
                // broad cast for multiplication
                auto broadcast3 = [](const VectorXd & v) {
                    MatrixX3d mat;
                    mat << v, v, v;
                    return mat;
                };
                auto broadcast3Back = [](const MatrixX3d & mat) {
                    VectorXd v = mat.rowwise().sum();
                    return v;
                };

                auto firstPoints = cwiseProd(
                    composeMappingFunction<MatrixX3d, VectorXd>(lambdas.assign<VectorXd>(), 
                        broadcast3, broadcast3Back, "broadcast3", "broadcast3back"), 
                    firstPointFactors).eval();
                auto secondPoints = cwiseProd(
                    composeMappingFunction<MatrixX3d, VectorXd>(lambdas.assign<VectorXd>(),
                        broadcast3, broadcast3Back, "broadcast3", "broadcast3back"),
                    secondPointFactors).eval();

                // TODO

                // select lines for constrained lines
                // lines1
                //auto selectLambdasOfConstrainedLines1 = [&constraints](const MatrixX3d & lmds) {
                //    MatrixX3d slmds = MatrixX3d::Zero(constraints.size(), 3);
                //    for (int i = 0; i < constraints.size(); i++) {
                //        slmds.row(i) = lmds.row(constraints[i].mergedSpatialLineSegmentIds[0]);
                //    }
                //    return slmds;
                //};
                //auto selectLambdasOfConstrainedLines1Back = [&constraints, &lines](const MatrixX3d & slmds) {
                //    VectorXd lmds = VectorXd::Zero(lines.size());
                //    for (int i = 0; i < constraints.size(); i++) {
                //        lmds(constraints[i].mergedSpatialLineSegmentIds[0]) += slmds(i);
                //    }
                //    return lmds;
                //};

                //// lines2
                //auto selectLambdasOfConstrainedLines2 = [&constraints](const MatrixX3d & lmds) {
                //    VectorXd slmds = VectorXd::Zero(constraints.size());
                //    for (int i = 0; i < constraints.size(); i++) {
                //        slmds(i) = lmds(constraints[i].mergedSpatialLineSegmentIds[1]);
                //    }
                //    return slmds;
                //};
                //auto selectLambdasOfConstrainedLines2Back = [&constraints, &lines](const MatrixX3d & slmds) {
                //    VectorXd lmds = VectorXd::Zero(lines.size());
                //    for (int i = 0; i < constraints.size(); i++) {
                //        lmds(constraints[i].mergedSpatialLineSegmentIds[1]) += slmds(i);
                //    }
                //    return lmds;
                //};

                //auto lambdas1 = composeMappingFunction<VectorXd, VectorXd>(lambdas,
                //    selectLambdasOfConstrainedLines1, selectLambdasOfConstrainedLines1Back);
                //auto lambdas2 = composeMappingFunction<VectorXd, VectorXd>(lambdas,
                //    selectLambdasOfConstrainedLines2, selectLambdasOfConstrainedLines2Back);
            }
        }


        void ReconstructionEngine::rectifySpatialLines() {
            // constraints for 3D lines reconstruction
            std::vector<LineStructureConnectivityConstraintData> constraints, refinedConstraints;
            constraints.reserve(Square(_globalData.mergedSpatialLineSegments.size()) / 4);
            constraints.clear();

            // merge lines, and get (some) incidence constraints
            _globalData.mergedSpatialLineSegments = 
                MergeColinearSpatialLinesAndAppendIncidenceConstraints(
                _globalData.spatialLineSegments,
                _globalData.mergedSpatialLineSegmentChainIds, 
                constraints, 
                _params.mergeLineDistanceAngleThreshold,
                _params.incidenceConstraintLineDistanceAngleThreshold);

            // make all point directions of lines normalized
            for (auto & line : _globalData.mergedSpatialLineSegments){
                line.component.first /= norm(line.component.first);
                line.component.second /= norm(line.component.second);
            }

            // get all intersection and incidence constraints
            AppendIntersectionAndOptionallyIncidenceConstraints(_globalData.mergedSpatialLineSegments,
                constraints, 
                _params.intersectionConstraintLineDistanceAngleThreshold, 
                _params.incidenceConstraintLineDistanceAngleThreshold, false);

            // remove all duplicated constraints
            // remove all self connected lines
            // remove all constraints close to any vanishing points
            std::vector<decltype(constraints.begin())> uniqueConsIters;
            MergeNearNaive(constraints.begin(), constraints.end(), std::back_inserter(uniqueConsIters), std::false_type(),
                1, [](const LineStructureConnectivityConstraintData & cons1, const LineStructureConnectivityConstraintData & cons2){
                return (cons1.type == cons2.type && 
                    std::is_permutation(std::begin(cons1.mergedSpatialLineSegmentIds), std::end(cons1.mergedSpatialLineSegmentIds), 
                    std::begin(cons2.mergedSpatialLineSegmentIds))) ? 0 : 2;
            });
            std::vector<LineStructureConnectivityConstraintData> uniqueCons;
            for (auto consIter : uniqueConsIters){
                auto & cons = *consIter;
                if (cons.mergedSpatialLineSegmentIds[0] == cons.mergedSpatialLineSegmentIds[1] ||
                    MaybeVanishingPoint(cons.position, _globalData.vanishingPoints, 
                    _params.intersectionConstraintLineDistanceAngleThreshold))
                    continue;
                uniqueCons.push_back(cons);
            }
            constraints = uniqueCons;

            IF_DEBUG_USING_VISUALIZERS {
                // debug constraints
                std::vector<Line3> consLines;
                std::vector<Point3> consPoints;
                consLines.reserve(constraints.size());
                consPoints.reserve(constraints.size());
                for (auto & cons : constraints){
                    auto & line1 = _globalData.mergedSpatialLineSegments[cons.mergedSpatialLineSegmentIds[0]].component;
                    auto & line2 = _globalData.mergedSpatialLineSegments[cons.mergedSpatialLineSegmentIds[1]].component;
                    auto pp = DistanceBetweenTwoLines(line1, line2);
                    consLines.push_back(Line3(pp.second.first.position, pp.second.second.position));
                    consPoints.push_back(cons.position);
                }
                vis::Visualizer3D viz;
                viz << vis::manip3d::SetWindowName("show constraints recognized")
                    << vis::manip3d::SetDefaultColor(vis::ColorTag::Yellow)
                    << vis::manip3d::SetColorTableDescriptor(vis::ColorTableDescriptor::RGB)
                    << _globalData.mergedSpatialLineSegments;
                viz << vis::manip3d::SetDefaultColor(vis::ColorTag::DimGray)
                    << consLines
                    << vis::manip3d::AutoSetCamera
                    << vis::manip3d::Show(false);
            }
            
            // vote for the position of each constraint
            VoteManhattanJunctionWeightsOnConstraints(_globalData.mergedSpatialLineSegments,
                _globalData.vanishingPoints, constraints);

            // compute weights of contraints 
            for (auto & cons : constraints){
                cons.weight = cons.junctionWeights.Triplet * _params.mjWeightTriplet + 
                    cons.junctionWeights.T * _params.mjWeightT + 
                    cons.junctionWeights.X * _params.mjWeightX + 
                    cons.junctionWeights.L * _params.mjWeightL +
                    cons.junctionWeights.I * _params.mjWeightI;
                //std::cout << "cons weight: " << cons.weight << std::endl;
                assert(cons.weight >= 0);
            }           
            std::sort(constraints.begin(), constraints.end(),
                [](const LineStructureConnectivityConstraintData & cons1, const LineStructureConnectivityConstraintData & cons2){
                return cons1.weight > cons2.weight;
            });

            std::cout << "line num: " << _globalData.mergedSpatialLineSegments.size() << std::endl;
            std::cout << "constraint num: " << constraints.size() << std::endl;

            // optimize lines
            OptimizeLines(_globalData.mergedSpatialLineSegments, 
                constraints, _globalData.vanishingPoints);

            IF_DEBUG_USING_VISUALIZERS {
                std::vector<Line3> consLines;
                consLines.reserve(constraints.size());
                double maxAngle = 0;
                LineStructureConnectivityConstraintData * maxCons = nullptr;
                for (auto & cons : constraints){
                    auto & line1 = _globalData.mergedSpatialLineSegments[cons.mergedSpatialLineSegmentIds[0]].component;
                    auto & line2 = _globalData.mergedSpatialLineSegments[cons.mergedSpatialLineSegmentIds[1]].component;
                    auto pp = DistanceBetweenTwoLines(line1, line2);
                    if (pp.first > maxAngle){
                        maxAngle = pp.first;
                        maxCons = &cons;
                    }
                    consLines.push_back(Line3(pp.second.first.position, pp.second.second.position));
                }
                std::cout << "max distance between constrained lines: " << maxAngle << std::endl;
                vis::Visualizer3D viz;
                viz << vis::manip3d::SetWindowName("constraints and optimized lines")
                    << vis::manip3d::SetDefaultColor(vis::ColorTag::Yellow)
                    << vis::manip3d::SetColorTableDescriptor(vis::ColorTableDescriptor::RGB)
                    << _globalData.mergedSpatialLineSegments;
                viz << vis::manip3d::SetDefaultColor(vis::ColorTag::DimGray)
                    << consLines
                    << vis::manip3d::AutoSetCamera
                    << vis::manip3d::Show(false);
            }

            // find necessary constraints using MST with slackValues
            std::vector<size_t> lineIds(_globalData.mergedSpatialLineSegments.size()), 
                consIds(constraints.size());
            for (size_t i = 0; i < lineIds.size(); i++)
                lineIds[i] = i;
            for (size_t i = 0; i < consIds.size(); i++)
                consIds[i] = i;
            std::vector<size_t> MSTconsIds;
            MSTconsIds.reserve(constraints.size());
            MinimumSpanningTree(lineIds.begin(), lineIds.end(), consIds.begin(), consIds.end(),
                std::back_inserter(MSTconsIds),
                [&constraints](size_t e){ return std::make_pair(constraints[e].mergedSpatialLineSegmentIds[0], constraints[e].mergedSpatialLineSegmentIds[1]); },
                [&constraints](size_t e1, size_t e2){ return abs(constraints[e1].slackValue) < abs(constraints[e2].slackValue);}
            );
            std::vector<LineStructureConnectivityConstraintData> MSTconstraints(MSTconsIds.size());
            for (size_t i = 0; i < MSTconsIds.size(); i++){
                MSTconstraints[i] = constraints[MSTconsIds[i]];
            }
            std::cout << "line num: " << _globalData.mergedSpatialLineSegments.size() << std::endl;
            std::cout << "mst constraint num: " << MSTconstraints.size() << std::endl;
            refinedConstraints = MSTconstraints;

            // optimize lines again
            OptimizeLines(_globalData.mergedSpatialLineSegments,
                refinedConstraints, _globalData.vanishingPoints);


            IF_DEBUG_USING_VISUALIZERS {
                std::vector<Line3> consLines;
                consLines.reserve(refinedConstraints.size());
                double maxAngle = 0;
                LineStructureConnectivityConstraintData * maxCons = nullptr;
                for (auto & cons : refinedConstraints){
                    auto & line1 = _globalData.mergedSpatialLineSegments[cons.mergedSpatialLineSegmentIds[0]].component;
                    auto & line2 = _globalData.mergedSpatialLineSegments[cons.mergedSpatialLineSegmentIds[1]].component;
                    auto pp = DistanceBetweenTwoLines(line1, line2);
                    if (pp.first > maxAngle){
                        maxAngle = pp.first;
                        maxCons = &cons;
                    }
                    consLines.push_back(Line3(pp.second.first.position, pp.second.second.position));
                }
                std::cout << "max distance between constrained lines (after refinement): " << maxAngle << std::endl;
                vis::Visualizer3D viz;
                viz << vis::manip3d::SetWindowName("refined constraints and again-optimized lines")
                    << vis::manip3d::SetDefaultColor(vis::ColorTag::Yellow)
                    << vis::manip3d::SetColorTableDescriptor(vis::ColorTableDescriptor::RGB)
                    << _globalData.mergedSpatialLineSegments;
                viz << vis::manip3d::SetDefaultColor(vis::ColorTag::DimGray)
                    << consLines
                    << vis::manip3d::AutoSetCamera
                    << vis::manip3d::Show(true);
            }


            // compute connected components of mergedSpatialLineSegments using constraint connections
            std::vector<int> mergedSpatialLineSegmentIds(_globalData.mergedSpatialLineSegments.size());
            for (int i = 0; i < mergedSpatialLineSegmentIds.size(); i++){
                mergedSpatialLineSegmentIds[i] = i;
            }
            _globalData.mergedSpatialLineSegmentsClassifiedWithStructureIds = _globalData.mergedSpatialLineSegments;
            ConnectedComponents(mergedSpatialLineSegmentIds.begin(), mergedSpatialLineSegmentIds.end(), 
                [&refinedConstraints](int lineid){
                std::vector<int> lineidsInRelation;
                for (auto & con : refinedConstraints) {
                    if (con.mergedSpatialLineSegmentIds[0] == lineid)
                        lineidsInRelation.push_back(con.mergedSpatialLineSegmentIds[1]);
                    if (con.mergedSpatialLineSegmentIds[1] == lineid)
                        lineidsInRelation.push_back(con.mergedSpatialLineSegmentIds[0]);
                }
                return lineidsInRelation;
            }, 
                [this](int lineid, int ccid){
                _globalData.spatialStructuresOfMergedSpatialLineIds[ccid].push_back(lineid);
                _globalData.mergedSpatialLineSegmentsClassifiedWithStructureIds[lineid].claz = ccid;
            });

            // install to constraint graph
            // ignore the connected components, insert each line as individual

            // reset expression graph
            _exprGraph = deriv::ExpressionGraph();
            auto & graph = _exprGraph;

            _constraints.clear();

            // add line structure components
            _constraints.internalElements<0>().reserve(_globalData.mergedSpatialLineSegments.size());
            _constraints.internalElements<1>().reserve(refinedConstraints.size());
            for (int i = 0; i < _globalData.mergedSpatialLineSegments.size(); i++) {
                ComponentData cd(ComponentData::Type::LineStructure);
                auto h = _constraints.add(cd);
                LineStructureComponentData & lineStruct = _constraints.data(h).asLineStructure;
                lineStruct.mergedSpatialLineSegmentId = i;
                lineStruct.eta = norm(_globalData.mergedSpatialLineSegments[i].component.first);
                lineStruct.etaExpr = graph.addRef(lineStruct.eta);
            }

            // add line structure connectivity constraints
            for (auto & con : refinedConstraints) {
                ConstraintData cd(ConstraintData::Type::LineStructureConnectivity);
                cd.asLineStructureConnectivity = con;
                auto h = _constraints.add<1>({ ComponentHandle(con.mergedSpatialLineSegmentIds[0]), 
                    ComponentHandle(con.mergedSpatialLineSegmentIds[1]) }, cd);
                // expression
                //_constraints.data(h)
            }
        }


        namespace {

            // index
            struct RegionIndex {
                ReconstructionEngine::ViewHandle viewHandle;
                RegionsNet::RegionHandle regionHandle;
            };

            inline bool operator == (const RegionIndex & a, const RegionIndex & b) {
                return a.viewHandle == b.viewHandle && a.regionHandle == b.regionHandle;
            }

            inline bool operator < (const RegionIndex & a, const RegionIndex & b) {
                if (a.viewHandle.id != b.viewHandle.id)
                    return a.viewHandle.id < b.viewHandle.id;
                return a.regionHandle.id < b.regionHandle.id;
            }

            struct RegionBoundaryIndex {
                ReconstructionEngine::ViewHandle viewHandle;
                RegionsNet::BoundaryHandle boundaryHandle;
            };    

            inline bool operator == (const RegionBoundaryIndex & a, const RegionBoundaryIndex & b) {
                return a.viewHandle == b.viewHandle && a.boundaryHandle == b.boundaryHandle;
            }


            // region map
            struct RegionMapVertex {
                RegionIndex regionIndex;

                // (a,b,c) representing the plane: ax + by + cz = 1 
                // => a(kA) + b(kB) + c(kC) = 1 
                // => k = 1/(aA + bB + cC), 
                // (A, B, C) is the spatial direction of a given point, norm(A, B, C) = 1, k is the depth
                Vec3 theta; 

                bool isVoid;
            };

            struct RegionMapEdge {
                bool isOverlap;
                inline bool isBoundary() const { return !isOverlap; }
                struct {
                    RegionBoundaryIndex boundaryIndex;
                    bool isOccludingBoundary; // true/false -> disconnected/connected spatially
                } asBoundary; // isOverlap
                struct {
                    double overlapRatio;
                } asOverlap; // !isOverlap
            };

            using HolisticRegionMap = GraphicalModel02<RegionMapVertex, RegionMapEdge>;


            // polygon conversion
            void ConvertToGPCPolygon(const std::vector<PixelLoc> & pts, gpc_polygon & poly) {
                poly.num_contours = 1;
                poly.contour = new gpc_vertex_list[1];
                poly.contour[0].num_vertices = pts.size();
                poly.contour[0].vertex = new gpc_vertex[pts.size()];
                for (int i = 0; i < pts.size(); i++) {
                    poly.contour[0].vertex[i].x = pts[i].x;
                    poly.contour[0].vertex[i].y = pts[i].y;
                }
                poly.hole = new int[1];
                poly.hole[0] = 0;
            }

            void ConvertToPixelVector(const gpc_polygon & poly, std::vector<PixelLoc> & pts) {
                pts.clear();
                pts.resize(poly.contour[0].num_vertices);
                for (int i = 0; i < pts.size(); i++) {
                    pts[i].x = static_cast<int>(poly.contour[0].vertex[i].x);
                    pts[i].y = static_cast<int>(poly.contour[0].vertex[i].y);
                }
            }

            /*class Polygon2 {
            public:
                inline Polygon2() : _p(nullptr) {}
                inline Polygon2(const std::vector<PixelLoc> & pts) {
                    _p = new gpc_polygon;
                    ConvertToGPCPolygon(pts, *_p);
                }
                inline Polygon2(const Image & mask) {
                    NOT_IMPLEMENTED_YET();
                }

                inline ~Polygon2() {
                    gpc_free_polygon(_p);
                    delete _p;
                }
                inline Polygon2(const Polygon2 & p) {
                    NOT_IMPLEMENTED_YET();
                }
                inline Polygon2(Polygon2 && p) {
                    std::swap(_p, p._p);
                }

                inline std::vector<PixelLoc> toPixels() const {
                    std::vector<PixelLoc> pixels;
                    ConvertToPixelVector(*_p, pixels);
                    return pixels;
                }
                inline double area() const {
                    return cv::contourArea(toPixels());
                }

                inline Polygon2 & operator &= (const Polygon2 & p) {
                    NOT_IMPLEMENTED_YET();
                }
                inline Polygon2 & operator |= (const Polygon2 & p) {
                    NOT_IMPLEMENTED_YET();
                }
                inline Polygon2 & operator -= (const Polygon2 & p) {
                    NOT_IMPLEMENTED_YET();
                }

            private:
                gpc_polygon * _p;
            };

            inline Polygon2 operator & (const Polygon2 & a, const Polygon2 & b) {
                Polygon2 r = a;
                r &= b;
                return r;
            }
            inline Polygon2 operator | (const Polygon2 & a, const Polygon2 & b) {
                Polygon2 r = a;
                r |= b;
                return r;
            }
            inline Polygon2 operator - (const Polygon2 & a, const Polygon2 & b) {
                Polygon2 r = a;
                r -= b;
                return r;
            }*/


            struct MinOf6Traits : public deriv::OpTraitsBase<double, double, double, double, double, double, double> {
                inline double value(const double & e1, const double & e2, const double & e3, const double & e4, const double & e5, const double & e6) const {
                    return std::min({ e1, e2, e3, e4, e5, e5 });
                }
                inline void derivatives(
                    deriv::Expression<double> output,
                    deriv::DerivativeExpression<double> sumOfDOutputs,
                    deriv::OriginalAndDerivativeExpression<double> input1,
                    deriv::OriginalAndDerivativeExpression<double> input2,
                    deriv::OriginalAndDerivativeExpression<double> input3,
                    deriv::OriginalAndDerivativeExpression<double> input4,
                    deriv::OriginalAndDerivativeExpression<double> input5,
                    deriv::OriginalAndDerivativeExpression<double> input6) const {

                    std::array<deriv::OriginalAndDerivativeExpression<double>, 6> inputs = { {
                        input1,
                        input2,
                        input3,
                        input4,
                        input5,
                        input6
                    } };
                    for (int i = 0; i < 6; i++) {
                        auto selectMinumunToBeSumOfDOutputs = [output, inputs, i](double sumOfDOutputsVal) -> double {
                            auto maxPos = std::min_element(inputs.begin(), inputs.end(), 
                                [](const deriv::OriginalAndDerivativeExpression<double> & a,
                                const deriv::OriginalAndDerivativeExpression<double> & b) {
                                return a.first.result() < b.first.result();
                            });
                            auto maxId = std::distance(inputs.begin(), maxPos);
                            if (maxId == i) {
                                return sumOfDOutputsVal;
                            } else {
                                return 0.0;
                            }
                        };
                        inputs[i].second = deriv::ComposeExpressionWithoutDerivativeDefinition(
                            selectMinumunToBeSumOfDOutputs, sumOfDOutputs);
                    }                   
                }
                virtual std::ostream & toString(std::ostream & os) const { os << "min6"; return os; }
            };

            inline deriv::Expression<double> minOf6(const std::array<deriv::Expression<double>, 6> & inputs) {
                return deriv::ComposeExpression(MinOf6Traits(), 
                    inputs[0], inputs[1], inputs[2], inputs[3], inputs[4], inputs[5]);
            }

            inline deriv::Expression<double> gaussianFunc(const deriv::Expression<double> & input, double sigma) {
                return exp(-input * input / 2.0 / sigma / sigma);
            }

            struct HFuncTraits : public deriv::OpTraitsBase<double, double> {
                explicit HFuncTraits(double s) : sigma(s) {}
                inline double value(const double & e) const {
                    return abs(e) <= sigma ? Square(e / sigma) : 1.0;
                }
                inline void derivatives(deriv::Expression<double> output,
                    deriv::DerivativeExpression<double> sumOfDOutputs,
                    deriv::OriginalAndDerivativeExpression<double> input) const {
                    input.second = deriv::ComposeExpressionWithoutDerivativeDefinition(
                        [input, this](double sumOfDOutputsVal) -> double{
                        double inputv = input.first.result();
                        if (abs(inputv) < sigma) {
                            return sumOfDOutputsVal * 2.0 * inputv / sigma / sigma;
                        } else {
                            return 0.0;
                        }
                    }, sumOfDOutputs);
                }
                double sigma;
            };

            inline deriv::Expression<double> hFunc(const deriv::Expression<double> & input, double sigma) {
                return deriv::ComposeExpression(HFuncTraits(sigma), input);
            }

            inline deriv::Expression<double> angleFunc(const deriv::Expression<Eigen::Vector3d> & a, const deriv::Expression<Eigen::Vector3d> & b) {
                return deriv::acos(deriv::dotProd(a, b) / deriv::norm(a) / deriv::norm(b));
            }
        }
 

        void ReconstructionEngine::reconstructFaces() {

            // compute spatial positions of each region
            struct {
                inline uint64_t operator()(const RegionIndex & ri) const {
                    return static_cast<uint64_t>((ri.viewHandle.id << 4) + ri.regionHandle.id);
                }
            } hashRegionIndex;
            std::unordered_map<RegionIndex, std::vector<Vec3>, decltype(hashRegionIndex)> 
                regionSpatialContours(10000, hashRegionIndex);          
            for (auto & view : _views.elements<0>()) {
                const auto & regions = * view.data.regionNet;
                for (auto & region : regions.regions().elements<0>()) {
                    RegionIndex ri = { view.topo.hd, region.topo.hd };
                    const ReconstructionEngine::ViewData & vd = view.data;
                    const RegionsNet::RegionData & rd = region.data;

                    assert(!rd.contours.empty() && "Region contour not initialized yet?");
                    std::vector<Vec3> spatialContour;
                    spatialContour.reserve(rd.contours.front().size());
                    std::transform(rd.contours.front().begin(), rd.contours.front().end(), 
                        std::back_inserter(spatialContour), 
                        [&vd](const PixelLoc & p) {
                        auto direction = vd.camera.spatialDirection(p);
                        return direction / norm(direction);
                    });
                    regionSpatialContours[ri] = spatialContour;
                }
            }

            // build spatial rtree for regions
            auto lookupRegionBB = [&regionSpatialContours](const RegionIndex& ri) {
                return BoundingBoxOfContainer(regionSpatialContours[ri]);
            };
            RTreeWrapper<RegionIndex, decltype(lookupRegionBB)> regionsRTree(lookupRegionBB);
            for (auto & region : regionSpatialContours) {
                regionsRTree.insert(region.first);
            }
            
            // store overlapping ratios between overlapped regions
            struct {
                inline uint64_t operator()(const std::pair<RegionIndex, RegionIndex> & vhp) const {
                    return static_cast<uint64_t>(
                        (((vhp.first.regionHandle.id << 3) + vhp.first.viewHandle.id) << 10) + 
                        ((vhp.second.regionHandle.id << 3) + vhp.second.viewHandle.id));
                }
            } hashRegionIndexPair;
            std::unordered_map<std::pair<RegionIndex, RegionIndex>, double, decltype(hashRegionIndexPair)> 
                overlappedRegionIndexPairs(1000, hashRegionIndexPair);
            
            for (auto & rip : regionSpatialContours) {
                auto & ri = rip.first;
                auto & riContour2d = _views.data(ri.viewHandle).regionNet->regions().data(ri.regionHandle).contours.front();
                auto & riCamera = _views.data(ri.viewHandle).camera;
                double riArea = _views.data(ri.viewHandle).regionNet->regions().data(ri.regionHandle).area;

                gpc_polygon riPoly;
                ConvertToGPCPolygon(riContour2d, riPoly);

                regionsRTree.search(lookupRegionBB(ri), 
                    [&ri, &riContour2d, &riPoly, &riCamera, riArea, &overlappedRegionIndexPairs, &regionSpatialContours](
                    const RegionIndex & relatedRi) {

                    if (ri.viewHandle == relatedRi.viewHandle) {
                        return true;
                    }

                    // project relatedRi contour to ri's camera plane
                    auto & relatedRiContour3d = regionSpatialContours[relatedRi];
                    std::vector<PixelLoc> relatedRiContour2d(relatedRiContour3d.size());
                    for (int i = 0; i < relatedRiContour3d.size(); i++) {
                        auto p = riCamera.screenProjection(relatedRiContour3d[i]);
                        relatedRiContour2d[i] = PixelLoc(p);
                    }
                    gpc_polygon relatedRiPoly;
                    ConvertToGPCPolygon(relatedRiContour2d, relatedRiPoly);

                    // compute overlapping area ratio
                    gpc_polygon intersectedPoly;
                    gpc_polygon_clip(GPC_INT, &relatedRiPoly, &riPoly, &intersectedPoly);

                    if (intersectedPoly.num_contours > 0 && intersectedPoly.contour[0].num_vertices > 0) {
                        std::vector<PixelLoc> intersected;
                        ConvertToPixelVector(intersectedPoly, intersected);
                        double intersectedArea = cv::contourArea(intersected);

                        double overlapRatio = intersectedArea / riArea;
                        assert(overlapRatio <= 1.0 && "Invalid overlap ratio!");

                        if (overlapRatio > 0.2)
                            overlappedRegionIndexPairs[std::make_pair(relatedRi, ri)] = overlapRatio;
                    }

                    gpc_free_polygon(&relatedRiPoly);
                    gpc_free_polygon(&intersectedPoly);

                    return true;
                });

                gpc_free_polygon(&riPoly);
            }

            for (auto & riPair : overlappedRegionIndexPairs) {
                auto revRiPair = std::make_pair(riPair.first.second, riPair.first.first);
                std::cout << "a-b: " << riPair.second;
                if (overlappedRegionIndexPairs.find(revRiPair) != overlappedRegionIndexPairs.end())
                    std::cout << "   b-a: " << overlappedRegionIndexPairs[revRiPair];
                std::cout << std::endl;
            }


            // build spatial rtree for mergedlines
            // find line-region relations
            std::map<std::pair<int, RegionIndex>, std::vector<Point2>> lineRegionSharedPoints; // (lineid, regionidx) -> points
            std::map<std::pair<int, RegionIndex>, std::vector<Vec3>> lineRegionSharedDirections;
            for (int i = 0; i < _globalData.mergedSpatialLineSegments.size(); i++) {
                if (_globalData.mergedSpatialLineSegments[i].claz == -1) // ignore non manhattan lines
                    continue;
                auto & line = _globalData.mergedSpatialLineSegments[i].component;
                static const double angleStep = M_PI / 128.0;
                double spanAngle = AngleBetweenDirections(line.first, line.second);
                std::vector<Vec3> sampledDirectionsInRegions(static_cast<size_t>(std::ceil(spanAngle / angleStep)));
                for (int k = 0; k < sampledDirectionsInRegions.size(); k++) {
                    sampledDirectionsInRegions[k] = RotateDirectionTo(line.first, line.second, k * angleStep / spanAngle);
                    sampledDirectionsInRegions[k] /= norm(sampledDirectionsInRegions[k]);
                }
                for (auto & sampledDir : sampledDirectionsInRegions) {
                    auto sampledBB = DefaultInfluenceBoxFunctor<Vec3>(0.5)(sampledDir / norm(sampledDir));                    
                    regionsRTree.search(sampledBB, 
                        [&line, this, i, &sampledDir, &lineRegionSharedPoints, &lineRegionSharedDirections](const RegionIndex & ri) -> bool {
                        // project line to ri view
                        auto & cam = _views.data(ri.viewHandle).camera;
                        auto p = cam.screenProjection(sampledDir);
                        auto d = cv::pointPolygonTest(_views.data(ri.viewHandle).regionNet->regions().data(ri.regionHandle).contours.front(),
                            PixelLoc(p[0], p[1]), false);
                        if (d > 0) { // inside
                            lineRegionSharedPoints[std::make_pair(i, ri)].push_back(p);
                            lineRegionSharedDirections[std::make_pair(i, ri)].push_back(sampledDir);
                        }
                        return true;
                    });
                }
            }



            // append to constraint graph
            _constraints.internalElements<0>().reserve(_constraints.internalElements<0>().capacity() + 
                regionsRTree.size());
            _constraints.internalElements<1>().reserve(_constraints.internalElements<1>().capacity() + 
                overlappedRegionIndexPairs.size() + regionsRTree.size());

            // add region components
            std::unordered_map<RegionIndex, HandleAtLevel<0>, decltype(hashRegionIndex)> ri2Handle(50000, hashRegionIndex);
            for (auto & r : regionSpatialContours) {
                ComponentData cd(ComponentData::Type::Region);
                //cd.asRegion.isVoid = false;
                cd.asRegion.regionHandle = r.first.regionHandle;
                cd.asRegion.viewHandle = r.first.viewHandle;
                auto & cam = _views.data(r.first.viewHandle).camera;
                cd.asRegion.theta = (cam.center() - cam.eye()) / (cam.center() - cam.eye()).dot(cam.center());
                ri2Handle[r.first] = _constraints.add(cd);
            }
            //// add line structure components
            //std::vector<HandleAtLevel<0>> mergedLineId2Handle(_globalData.mergedSpatialLineSegments.size());
            //for (int i = 0; i < _globalData.mergedSpatialLineSegments.size(); i++){
            //    ComponentData cd(ComponentData::Type::LineStructure);
            //    cd.asLineStructure.eta = 1.0;
            //    cd.asLineStructure.mergedSpatialLineSegmentId = i;
            //    mergedLineId2Handle[i] = _constraints.add(cd);
            //}
            // add region overlap constraints
            for (auto & o : overlappedRegionIndexPairs) {
                auto h1 = ri2Handle[o.first.first];
                auto h2 = ri2Handle[o.first.second];

                ConstraintData cd(ConstraintData::Type::RegionOverlap);
                cd.asRegionOverlap.overlapRatio = o.second;
                _constraints.add<1>({ h1, h2 }, cd);
            }
            // add region connectivity constraints
            for (auto & vd : _views.elements<0>()) {
                auto & regions = * vd.data.regionNet;
                for (auto & regionBoundaryData : regions.regions().elements<1>()) {
                    auto & rids = regionBoundaryData.topo.lowers;
                    RegionIndex ri1 = { vd.topo.hd, rids[0] };
                    RegionIndex ri2 = { vd.topo.hd, rids[1] };
                    auto rh1 = ri2Handle[ri1];
                    auto rh2 = ri2Handle[ri2];

                    ConstraintData cd(ConstraintData::Type::RegionConnectivity);
                    cd.asRegionPairConsistency.viewHandle = vd.topo.hd;
                    cd.asRegionPairConsistency.boundaryHandle = regionBoundaryData.topo.hd;
                    _constraints.add<1>({ rh1, rh2 }, cd);
                }
            }
            // add region line structure connectivity constraints
            for (auto & lineRegionData : lineRegionSharedDirections) {
                if (lineRegionData.second.empty())
                    continue;
                ConstraintData cd(ConstraintData::Type::RegionLineStructureConnectivity);
                cd.asRegionLineStructureConnectivity.sampledPoints = lineRegionData.second;
                
                auto lineHandle = HandleAtLevel<0>(lineRegionData.first.first);
                auto regionHandle = ri2Handle[lineRegionData.first.second];
                _constraints.add<1>({ lineHandle, regionHandle }, cd);
            }

            std::cout << "component num: " << _constraints.internalElements<0>().size() << std::endl;
            std::cout << "constraint num: " << _constraints.internalElements<1>().size() << std::endl;

            std::cout << "begin building constraints" << std::endl;
            
            
            // inference region orientations and spatial connectivity of boundaries
            using namespace deriv;
            ExpressionGraph & graph = _exprGraph;

            std::array<Expression<Eigen::Vector3d>, 6> principleDirectionExprs;
            for (int i = 0; i < 3; i++) {
                const Vec3 & vp = _globalData.vanishingPoints[i];
                principleDirectionExprs[i * 2] = composeFunction(graph, [&vp]() {
                    return CVMatToEigenMat(vp / norm(vp)); 
                });
                principleDirectionExprs[i * 2 + 1] = composeFunction(graph, [&vp]() {
                    return CVMatToEigenMat(- vp / norm(vp));
                });
            }
            
            // component expressions
            for (auto & r : _constraints.elements<0>()) {
                if (r.data.type == ComponentData::Type::Region) {
                    // theta
                    r.data.asRegion.thetaExpr = composeFunction(graph, [&r]() {
                        return CVMatToEigenMat(r.data.asRegion.theta); 
                    });
                    // E_manh
                    std::array<Expression<double>, 6> angles;
                    for (int i = 0; i < 6; i++) {
                        angles[i] = angleFunc(r.data.asRegion.thetaExpr, principleDirectionExprs[i]);
                    }
                    Expression<double> minAngle = minOf6(angles);
                    static const double wManh = 1.0;
                    r.data.asRegion.manhattanEnergyExpr = 
                        DisableableExpression<double>(hFunc(minAngle, M_PI / 6) * wManh, graph);
                } else {
                    // eta
                    r.data.asLineStructure.etaExpr = graph.addRef(r.data.asLineStructure.eta);
                }
            }

            // constraint energy expressions
            for (auto & con : _constraints.elements<1>()) {
                auto & rd1 = _constraints.data(con.topo.lowers[0]);
                auto & rd2 = _constraints.data(con.topo.lowers[1]);

                // region region overlap
                if (con.data.type == ConstraintData::Type::RegionOverlap) {
                    assert(rd1.type == ComponentData::Type::Region &&
                        rd2.type == ComponentData::Type::Region &&
                        "invalid component type!");

                    auto thetaDiffExpr = rd1.asRegion.thetaExpr - rd2.asRegion.thetaExpr;
                    auto thetaDiffSquaredSum = abs(dotProd(thetaDiffExpr, thetaDiffExpr));
                    
                    static const double wOverlap = 5.0;

                    // E_overlap
                    auto overlapRatio = con.data.asRegionOverlap.overlapRatio;
                    assert(overlapRatio > 0);
                    con.data.constraintEnergyExpr =
                        DisableableExpression<double>(wOverlap * overlapRatio * thetaDiffSquaredSum, graph);

                }
                // region region connectivity
                else if (con.data.type == ConstraintData::Type::RegionConnectivity) {
                    assert(rd1.type == ComponentData::Type::Region &&
                        rd2.type == ComponentData::Type::Region &&
                        "invalid component type!");

                    auto & camera = _views.data(con.data.asRegionPairConsistency.viewHandle).camera;
                    auto & boundaryData = _views.data(con.data.asRegionPairConsistency.viewHandle).regionNet
                        ->regions().data(con.data.asRegionPairConsistency.boundaryHandle);
                    auto & sampledPoints = boundaryData.sampledPoints;
                    auto straightness = boundaryData.straightness;

                    // get sampled direction expressions
                    std::vector<Expression<Eigen::Vector3d>> sampledDirectionExprs;
                    sampledDirectionExprs.reserve(boundaryData.length / 5);
                    for (auto & ps : sampledPoints) {
                        for (auto & p : ps) {
                            auto sampledDirectionExpr = composeFunction(graph, [&p, &camera]() {
                                auto direction = camera.spatialDirection(p);
                                return CVMatToEigenMat(direction / norm(direction));
                            });
                            sampledDirectionExprs.push_back(sampledDirectionExpr);
                        }
                    }
                    // get depth expressions of sampled direction
                    std::vector<Expression<double>> lambdasR1(sampledDirectionExprs.size()),
                        lambdasR2(sampledDirectionExprs.size());
                    std::vector<EHandle> lambdaSquaredDiffHandles(sampledDirectionExprs.size());
                    for (int i = 0; i < sampledDirectionExprs.size(); i++) {
                        lambdasR1[i] = dotProd(rd1.asRegion.thetaExpr, rd1.asRegion.thetaExpr) /
                            dotProd(rd1.asRegion.thetaExpr, sampledDirectionExprs[i]);
                        lambdasR2[i] = dotProd(rd2.asRegion.thetaExpr, rd2.asRegion.thetaExpr) /
                            dotProd(rd2.asRegion.thetaExpr, sampledDirectionExprs[i]);
                        lambdaSquaredDiffHandles[i] = Square(lambdasR1[i] - lambdasR2[i]).handle();
                    }
                    EHandle lambdaSquaredDiffSumHandle = HSum<double>(&graph, lambdaSquaredDiffHandles);
                    auto lambdaSquaredDiffSum = graph.as<double>(lambdaSquaredDiffSumHandle);

                    // E_connect
                    double occludedness = straightness;
                    assert(occludedness > 0);
                    static const double wConnect = 1.0;
                    con.data.constraintEnergyExpr =
                        DisableableExpression<double>(wConnect * lambdaSquaredDiffSum * occludedness, graph);

                }
                // line line connectivity
                else if (con.data.type == ConstraintData::Type::LineStructureConnectivity) {

                    assert(rd1.type == ComponentData::Type::LineStructure &&
                        rd2.type == ComponentData::Type::LineStructure &&
                        "invalid component type!");

                    auto & lineData1 = rd1.asLineStructure;
                    auto & lineData2 = rd2.asLineStructure;
                    auto & line1 = _globalData.mergedSpatialLineSegments[lineData1.mergedSpatialLineSegmentId];
                    auto & line2 = _globalData.mergedSpatialLineSegments[lineData2.mergedSpatialLineSegmentId];

                    auto nearestData = DistanceBetweenTwoLines(line1.component, line2.component).second;
                    auto sampledPoint = normalize((nearestData.first.position + nearestData.second.position) / 2.0);
                    double depthOfSampledPointOnLine1 =
                        norm(DistanceBetweenTwoLines(InfiniteLine3(Point3(0, 0, 0), sampledPoint),
                        line1.component.infinieLine()).second.second);
                    double depthOfSampledPointOnLine2 =
                        norm(DistanceBetweenTwoLines(InfiniteLine3(Point3(0, 0, 0), sampledPoint),
                        line2.component.infinieLine()).second.second);
                    double ratio1 = depthOfSampledPointOnLine1 / norm(line1.component.first);
                    assert(ratio1 > 0);
                    double ratio2 = depthOfSampledPointOnLine2 / norm(line2.component.first);
                    assert(ratio2 > 0);
                    auto lambda1 = lineData1.etaExpr * ratio1;
                    auto lambda2 = lineData2.etaExpr * ratio2;

                    static const double wLineConnect = 2.0;
                    auto lambdaSquaredDiff = Square(lambda1 - lambda2);
                    con.data.constraintEnergyExpr =
                        DisableableExpression<double>((wLineConnect * lambdaSquaredDiff).cast<double>(), graph);

                }
                // region line connectivity
                else if (con.data.type == ConstraintData::Type::RegionLineStructureConnectivity) {

                    assert((rd1.type == ComponentData::Type::LineStructure && rd2.type == ComponentData::Type::Region ||
                        rd1.type == ComponentData::Type::Region && rd2.type == ComponentData::Type::LineStructure) &&
                        "invalid component type!");

                    auto prrd1 = &rd1, prrd2 = &rd2;
                    if (rd1.type == ComponentData::Type::LineStructure){
                        std::swap(prrd1, prrd2);
                    }
                    auto & rrd1 = *prrd1; // region
                    auto & rrd2 = *prrd2; // linestruct

                    assert(rrd1.type == ComponentData::Type::Region && 
                        rrd2.type == ComponentData::Type::LineStructure &&
                        "invalid component type!");

                    auto & sampledPoints = con.data.asRegionLineStructureConnectivity.sampledPoints;

                    // get sampled direction expressions
                    std::vector<Expression<Eigen::Vector3d>> sampledDirectionExprs;
                    sampledDirectionExprs.reserve(sampledPoints.size());
                    for (auto & p : sampledPoints) {
                        auto sampledDirectionExpr = composeFunction(graph, [&p]() {
                            return CVMatToEigenMat(normalize(p));
                        });
                        sampledDirectionExprs.push_back(sampledDirectionExpr);
                    }
                    
                    // get depth expressions of sampled direction
                    std::vector<Expression<double>> lambdasR1(sampledDirectionExprs.size()),
                        lambdasR2(sampledDirectionExprs.size());
                    std::vector<EHandle> lambdaSquaredDiffHandles(sampledDirectionExprs.size());
                    for (int i = 0; i < sampledDirectionExprs.size(); i++) {
                        lambdasR1[i] = dotProd(rrd1.asRegion.thetaExpr, rrd1.asRegion.thetaExpr) /
                            dotProd(rrd1.asRegion.thetaExpr, sampledDirectionExprs[i]);
                        /*lambdasR2[i] = dotProd(rd2.asRegion.thetaExpr, rd2.asRegion.thetaExpr) /
                            dotProd(rd2.asRegion.thetaExpr, sampledDirectionExprs[i]);*/
                        auto & dir = sampledPoints[i]; // we don't need expression here
                        auto & line = _globalData.mergedSpatialLineSegments[rrd2.asLineStructure.mergedSpatialLineSegmentId].component;
                        auto depth = norm(DistanceBetweenTwoLines(InfiniteLine3({ 0, 0, 0 }, dir),
                            line.infinieLine()).second.second);
                        double ratio = depth / norm(line.first); // we only need the ratio of the depths
                        lambdasR2[i] = (rrd2.asLineStructure.etaExpr * ratio).cast<double>();

                        lambdaSquaredDiffHandles[i] = Square(lambdasR1[i] - lambdasR2[i]).handle();
                    }
                    EHandle lambdaSquaredDiffSumHandle = HSum<double>(&graph, lambdaSquaredDiffHandles);
                    auto lambdaSquaredDiffSum = graph.as<double>(lambdaSquaredDiffSumHandle);

                    // E_connect
                    static const double wConnect = 1.0;
                    con.data.constraintEnergyExpr =
                        DisableableExpression<double>(wConnect * lambdaSquaredDiffSum, graph);
                    
                }


            }

            // sum all energy
            std::vector<EHandle> allEnergyHandles;
            allEnergyHandles.reserve(_constraints.internalElements<1>().size() + _constraints.internalElements<0>().size());
            // add all E_manh
            for (auto & rd : _constraints.elements<0>()){
                if (rd.data.type == ComponentData::Type::Region){
                    allEnergyHandles.push_back(rd.data.asRegion.manhattanEnergyExpr.toExpression().handle());
                }
            }
            // add all connectivity energies
            for (auto & cd : _constraints.elements<1>()){
                allEnergyHandles.push_back(cd.data.constraintEnergyExpr.toExpression().handle());
            }
            // sum all
            auto completeEngergyExpr = graph.as<double>(HSum<double>(&graph, allEnergyHandles));
            std::cout << "current energy is: " << completeEngergyExpr.execute() << std::endl;

            // compute derivative


            std::cout << "done building constraints" << std::endl;

            // reconstruct faces
            // Eigen::ArrayXXd

        }


    }
}