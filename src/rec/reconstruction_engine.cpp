#include <filesystem>

extern "C" {
    #include <gpc.h>
}

#include <glpk.h>
#include <setjmp.h>

#include <Eigen/StdVector>
#include <dlib/matrix.h>
#include <dlib/optimization.h>

#include "../vis/visualize2d.hpp"
#include "../vis/visualize3d.hpp"

#include "optimization.hpp"
#include "reconstruction_engine.hpp"
#include "reconstruction_engine_visualize.hpp"

#include "../core/debug.hpp"

namespace panoramix {
    namespace rec {

        // for glpk
        struct sinfo {
            char * text;
            jmp_buf * env;
        };

        void glpErrorHook(void * in){
            sinfo * info = (sinfo*)in;
            glp_free_env();
            longjmp(*(info->env), 1);
        }


        ReconstructionEngine::Params::Params() 
            : camera(250.0), cameraAngleScaler(1.8), smallCameraAngleScalar(0.05),
            samplingStepLengthOnRegionBoundaries(16.0),
            samplingStepLengthOnLines(8.0),
            intersectionDistanceThreshold(30),
            incidenceDistanceAlongDirectionThreshold(50),
            incidenceDistanceVerticalDirectionThreshold(8)
        {}

        ReconstructionEngine::ComponentData::ComponentData(ReconstructionEngine::ComponentData::Type t) : type(t) {}

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

            RegionsNet::Params regionsNetParams;
            regionsNetParams.samplingStepLengthOnBoundary = _params.samplingStepLengthOnRegionBoundaries;
            vd.regionNet = std::make_shared<RegionsNet>(vd.image, regionsNetParams);
            vd.regionNet->buildNetAndComputeGeometricFeatures();
            vd.regionNet->computeImageFeatures();

            LinesNet::Params linesNetParams;
            linesNetParams.intersectionDistanceThreshold = _params.intersectionDistanceThreshold;
            linesNetParams.incidenceDistanceVerticalDirectionThreshold = _params.incidenceDistanceVerticalDirectionThreshold;
            linesNetParams.incidenceDistanceAlongDirectionThreshold = _params.incidenceDistanceAlongDirectionThreshold;
            vd.lineNet = std::make_shared<LinesNet>(vd.image, linesNetParams);
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

        void ReconstructionEngine::updateConnections() {
            NOT_IMPLEMENTED_YET();
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
                        score += votePanel.at<float>(WrapBetween(pixel.x, 0, longitudeDivideNum), 
                            WrapBetween(pixel.y, 0, latitudeDivideNum));
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
                            PixelLoc pixel = PixelIndexFromGeoCoord(GeoCoord(v), longitudeDivideNum, latitudeDivideNum);
                            score += votePanel.at<float>(WrapBetween(pixel.x, 0, longitudeDivideNum), 
                                WrapBetween(pixel.y, 0, latitudeDivideNum));
                        }
                        if (score > maxScore){
                            maxScore = score;
                            vps[1] = vec1;
                            vps[2] = vec2;
                        }
                    }
                }

                assert(UnOrthogonality(vps[0], vps[1], vps[2]) < 0.1);

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
                lineIntersectionsNum += vIter->data.lineNet->lineSegmentIntersections().size();
            std::vector<Vec3> intersections(lineIntersectionsNum);
            auto intersectionsBegin = intersections.begin();
            for (const auto & vIter : seperatedViewIters){ // projection 2d intersections to global GeoCoord
                intersectionsBegin = std::transform(vIter->data.lineNet->lineSegmentIntersections().begin(), 
                    vIter->data.lineNet->lineSegmentIntersections().end(),
                    intersectionsBegin, 
                    [&vIter](const HPoint2 & p) -> Vec3 {
                    Vec3 p3 = vIter->data.camera.spatialDirection(p.toPoint());
                    return p3 / norm(p3); // normalized
                });
            }

            // find vanishing points;
            _globalData.vanishingPoints = FindVanishingPoints(intersections);           


            // add spatial line segments from line segments of all views
            size_t spatialLineSegmentsNum = 0;
            for (auto & v : _views.elements<0>())
                spatialLineSegmentsNum += v.data.lineNet->lineSegments().size();
            std::vector<Classified<Line3>> spatialLineSegments(spatialLineSegmentsNum);
            auto spatialLineSegmentBegin = spatialLineSegments.begin();
            for (auto & v : _views.elements<0>()){
                spatialLineSegmentBegin = std::transform(v.data.lineNet->lineSegments().begin(),
                    v.data.lineNet->lineSegments().end(),
                    spatialLineSegmentBegin, [&v](const Line2 & line) -> Classified<Line3>{
                    auto & p1 = line.first;
                    auto & p2 = line.second;
                    auto pp1 = v.data.camera.spatialDirection(p1);
                    auto pp2 = v.data.camera.spatialDirection(p2);
                    Classified<Line3> cline3;
                    cline3.claz = -1;
                    cline3.component = Line3{ pp1, pp2 };
                    return cline3;
                });
            }

            // classify lines
            ClassifyLines(_globalData.vanishingPoints, spatialLineSegments);

            // build lines net and compute features
            spatialLineSegmentBegin = spatialLineSegments.begin();
            for (auto & v : _views.elements<0>()){
                std::array<HPoint2, 3> projectedVPs;
                for (int i = 0; i < 3; i++){
                    projectedVPs[i] = v.data.camera.screenProjectionInHPoint(_globalData.vanishingPoints[i]);
                }
                std::vector<int> lineClasses(v.data.lineNet->lineSegments().size());
                for (auto & lineClass : lineClasses){
                    lineClass = spatialLineSegmentBegin->claz;
                    ++spatialLineSegmentBegin;
                }
                v.data.lineNet->buildNetAndComputeFeaturesUsingVanishingPoints(projectedVPs, lineClasses);
            }


        }

















        namespace {

            struct LineStructureConnectivityConstraintData {

                LineStructureConnectivityConstraintData(){
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

                size_t mergedSpatialLineSegmentIds[2]; // corresponded mergedSpatialLineSegments ids
                PositionOnLine3 positionOnLines[2];
                Vec3 position; // location of intersecion

                // [i][0] -> line lengths with class i lying between vp[i] and position
                // [i][1] -> line lengths with class i lying between position and anti-vp[i]
                double lineVotings[3][2];
                double weight;
                struct { double I, L, X, T, Triplet; } junctionWeights;
                enum { Intersection, Incidence } type;
                double slackValue; // retreived after optimization
            };


            // merge colinear spatial lines and find some incidence constraints
            std::vector<Classified<Line3>> MergeColinearSpatialLinesAndAppendIncidenceConstraints(
                const std::vector<Classified<Line3>> & oldLines,
                std::vector<int> & chainIds,
                std::vector<LineStructureConnectivityConstraintData> & constraints,
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
                        if (normal.dot(firstNormal) < 0) { // not same direction, reverse it
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
                        Vec2 pdv1(pdir1.dot(firstPointDirection), pdir1.dot(firstNormalCrossPoint)); // transform to a coordinate system defined by the first line
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
                        if (lineAngle.first <= curTo){ // can extend
                            curTo = lineAngle.second;
                        }
                        else { // cannot extend
                            if (curTo - curFrom >= M_PI){ // break major arcs
                                mergedLineAngleSegments.push_back(std::make_pair(curFrom, (curFrom + curTo) / 2.0));
                                mergedLineAngleSegments.push_back(std::make_pair((curFrom + curTo) / 2.0, curTo));
                            }
                            else
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
            void OptimizeLines(std::vector<Classified<Line3>> & lines, 
                std::vector<LineStructureConnectivityConstraintData> & constraints,
                const std::array<Vec3, 3> & vps) {
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
                            { 0, 0 },
                            { 0, 1 },
                            { 1, 0 },
                            { 1, 1 }
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
                    glp_set_row_bnds(problem, ie.equationId + 1, GLP_UP, -1e5, 0.0); // ... <= 0.0
                }
                for (auto & ie : constraintInequations){
                    int varIds[] = { -1, ie.lambda_iId + 1, ie.lambda_jId + 1, ie.s_ijId + 1 };
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
                    glp_set_obj_coef(problem, var.varId + 1, var.weightInObjectiveFunction);
                }

                glp_adv_basis(problem, 0);

                enum Method { Simplex, Exact, Interior };
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
                        double lambda = glp_get_col_prim(problem, var.varId + 1);
                        // update line coordinates
                        lines[var.lineId].component = lineDeterminers[var.lineId](lambda);
                    }
                    for (auto & var : slacks){
                        double s = glp_get_col_prim(problem, var.varId + 1);
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
                else if (method == Interior){
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
                        double lambda = glp_ipt_col_prim(problem, var.varId + 1);
                        // update line coordinates
                        lines[var.lineId].component = lineDeterminers[var.lineId](lambda);
                    }
                    for (auto & var : slacks){
                        double s = glp_ipt_col_prim(problem, var.varId + 1);
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

        /*
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

        }
        */


        namespace {

            template <class HandleT>
            struct IndexOfSubStructureInView {
                ReconstructionEngine::ViewHandle viewHandle;
                HandleT handle;
            };

            template <class HandleT>
            inline bool operator == (const IndexOfSubStructureInView<HandleT> & a, 
                const IndexOfSubStructureInView<HandleT> & b) {
                return a.viewHandle == b.viewHandle && a.handle == b.handle;
            }

            template <class HandleT>
            inline bool operator < (const IndexOfSubStructureInView<HandleT> & a, 
                const IndexOfSubStructureInView<HandleT> & b) {
                if (a.viewHandle.id != b.viewHandle.id)
                    return a.viewHandle.id < b.viewHandle.id;
                return a.handle.id < b.handle.id;
            }

            template <class IndexT>
            struct IndexOfSubStructureInViewHasher {
                inline size_t operator()(const IndexT & idx) const {
                    return ((idx.viewHandle.id) << 4) + (idx.handle.id);
                }
            };

            template <class IndexT1, class IndexT2>
            struct IndexOfSubStructureInViewHasher < std::pair<IndexT1, IndexT2> > {
                inline size_t operator()(const std::pair<IndexT1, IndexT2> & idx) const {
                    return (((idx.first.viewHandle.id << 3) + idx.first.handle.id) << 10) +
                        (idx.second.viewHandle.id << 3) + idx.second.handle.id;
                }
            };

            using RegionIndex = IndexOfSubStructureInView <RegionsNet::RegionHandle> ;
            using RegionBoundaryIndex = IndexOfSubStructureInView < RegionsNet::BoundaryHandle > ;
            using LineIndex = IndexOfSubStructureInView <LinesNet::LineHandle> ;
            using LineRelationIndex = IndexOfSubStructureInView < LinesNet::LineRelationHandle > ;

            template <class T>
            using IndexHashSet = std::unordered_set<T, IndexOfSubStructureInViewHasher<T>> ;

            template <class KeyT, class ValueT>
            using IndexHashMap = std::unordered_map<KeyT, ValueT, IndexOfSubStructureInViewHasher<KeyT>> ;



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

            // line depth ratio
            double ComputeDepthRatioOfPointOnSpatialLine(Vec3 lineFirstPointDir,
                Vec3 p, Vec3 vp) {
                // firstp -> p vp
                //  \      /
                //   \    /
                //    center
                lineFirstPointDir /= norm(lineFirstPointDir);
                p /= norm(p);
                vp /= norm(vp);

                if ((p - lineFirstPointDir).dot(vp) < 0)
                    vp = -vp;
                double angleCenter = AngleBetweenDirections(lineFirstPointDir, p);
                double angleFirstP = AngleBetweenDirections(-lineFirstPointDir, vp);
                double angleP = AngleBetweenDirections(-p, -vp);
                //assert(FuzzyEquals(angleCenter + angleFirstP + angleP, M_PI, 0.1));
                return sin(angleFirstP) / sin(angleP);
            }


            inline deriv::Expression<double> GaussianFunc(const deriv::Expression<double> & input, double sigma) {
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

            inline deriv::Expression<double> HFunc(const deriv::Expression<double> & input, double sigma) {
                return deriv::ComposeExpression(HFuncTraits(sigma), input);
            }

            inline deriv::Expression<double> AngleFunc(const deriv::Expression<Eigen::Vector3d> & a, const deriv::Expression<Eigen::Vector3d> & b) {
                return deriv::acos(deriv::dotProd(a, b) / deriv::norm(a) / deriv::norm(b));
            }


            template <class CallbackFunctorT>
            void OptimizeConstraintGraphUsingEnergyExpression (
                const deriv::Expression<double> & energy, 
                ReconstructionEngine::ConstraintGraph & constraints,
                double delta, double momentum, int nepoches,
                CallbackFunctorT callback
                ){

                // compute derivative
                EHandleTable handleTable;
                handleTable.reserve(constraints.internalElements<0>().size());
                for (auto & vd : constraints.elements<0>()){
                    if (vd.data.type == ReconstructionEngine::ComponentData::Type::Region){
                        vd.data.asRegion.thetaExpr.registerHandleTable(handleTable);
                        for (int i = 0; i < 3; i++)
                            assert(!isnan(vd.data.asRegion.thetaExpr.expression().execute()[i]));
                    }
                    else if (vd.data.type == ReconstructionEngine::ComponentData::Type::Line){
                        vd.data.asLine.etaExpr.registerHandleTable(handleTable);
                        assert(!isnan(vd.data.asLine.etaExpr.expression().execute()));
                    }
                    else{
                        assert(false && "invalid component type!");
                    }
                }

                EHandleTable derivHandleTable;
                derivHandleTable.reserve(handleTable.size());
                energy.derivativesHandlesRange(handleTable.begin(), handleTable.end(),
                    std::back_inserter(derivHandleTable));

                for (auto & vd : constraints.elements<0>()){
                    if (vd.data.type == ReconstructionEngine::ComponentData::Type::Region){
                        vd.data.asRegion.thetaExpr.getDerivative(derivHandleTable);
                    }
                    else if (vd.data.type == ReconstructionEngine::ComponentData::Type::Line){
                        vd.data.asLine.etaExpr.getDerivative(derivHandleTable);
                    }
                    else{
                        assert(false && "invalid component type!");
                    }
                }

                // reconstruct faces
                auto startTime = std::chrono::high_resolution_clock::now();
                double energyVal = std::numeric_limits<double>::max();
                double energyValChange = energyVal;
                for (int i = 0; i < nepoches; i++){
                    double curEnergy = energy.execute();
                    if (curEnergy > energyVal || isnan(curEnergy)){
                        delta *= 0.6;
                        std::cout << "delta set to " << delta << std::endl;
                        for (auto & vd : constraints.elements<0>()){
                            if (vd.data.type == ReconstructionEngine::ComponentData::Type::Region){
                                vd.data.asRegion.thetaExpr.deOptimizeData();
                            }
                            else if (vd.data.type == ReconstructionEngine::ComponentData::Type::Line){
                                vd.data.asLine.etaExpr.deOptimizeData();
                            }
                            else{
                                assert(false && "invalid component type!");
                            }
                        }
                    }
                    else{
                        energyValChange = curEnergy - energyVal;
                        energyVal = curEnergy;
                    }
                    std::cout << "[" << i << "] current energy: " << energyVal << std::endl;
                    for (auto & vd : constraints.elements<0>()){
                        if (vd.data.type == ReconstructionEngine::ComponentData::Type::Region){
                            vd.data.asRegion.thetaExpr.optimizeData(delta, momentum, handleTable);
                        }
                        else if (vd.data.type == ReconstructionEngine::ComponentData::Type::Line){
                            vd.data.asLine.etaExpr.optimizeData(delta, momentum, handleTable);
                        }
                        else{
                            assert(false && "invalid component type!");
                        }
                    }

                    auto timeCost = std::chrono::high_resolution_clock::now() - startTime;
                    std::cout << "time cost: "
                        << std::chrono::duration_cast<std::chrono::seconds>(timeCost).count()
                        << std::endl;

                    if (!callback(i, energyVal, energyValChange))
                        break;
                }

            }
        
        }
 

        void ReconstructionEngine::reconstructLinesAndFaces() {

            //// REGIONS ////

            // compute spatial positions of each region
            IndexHashMap<RegionIndex, std::vector<Vec3>> 
                regionSpatialContours;
            for (auto & view : _views.elements<0>()) {
                const auto & regions = * view.data.regionNet;
                for (auto & region : regions.regions().elements<0>()) {
                    RegionIndex ri = { view.topo.hd, region.topo.hd };
                    const ReconstructionEngine::ViewData & vd = view.data;
                    const RegionsNet::RegionData & rd = region.data;

                    assert(!rd.contours.empty() && "Region contour not initialized yet?");
                    std::vector<Vec3> spatialContour;
                    for (auto & p : rd.dilatedContours.back()){
                        auto direction = vd.camera.spatialDirection(p);
                        spatialContour.push_back(direction / norm(direction));
                    }
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
            IndexHashMap<std::pair<RegionIndex, RegionIndex>, double> 
                overlappedRegionIndexPairs;
            
            for (auto & rip : regionSpatialContours) {
                auto & ri = rip.first;
                auto & riContour2d = _views.data(ri.viewHandle).regionNet->regions().data(ri.handle).contours.front();
                auto & riCamera = _views.data(ri.viewHandle).camera;
                double riArea = _views.data(ri.viewHandle).regionNet->regions().data(ri.handle).area;

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
                    std::vector<core::PixelLoc> relatedRiContour2d(relatedRiContour3d.size());
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
                        std::vector<core::PixelLoc> intersected;
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


            
            //// LINES ////
            // compute spatial normal directions for each line
            IndexHashMap<LineIndex, Classified<Line3>>
                lineSpatialAvatars;
            for (auto & vd : _views.elements<0>()){
                auto & lines = vd.data.lineNet->lines();
                LineIndex li;
                li.viewHandle = vd.topo.hd;
                auto & cam = vd.data.camera;
                for (auto & ld : lines.elements<0>()){
                    li.handle = ld.topo.hd;
                    auto & line = ld.data.line;
                    Classified<Line3> avatar;
                    avatar.claz = line.claz;
                    avatar.component = Line3(
                        cam.spatialDirection(line.component.first),
                        cam.spatialDirection(line.component.second)
                    );
                    lineSpatialAvatars[li] = avatar;
                }
            }

            // build rtree for lines
            auto lookupLineNormal = [&lineSpatialAvatars](const LineIndex & li) -> Box3 {
                auto normal = lineSpatialAvatars[li].component.first.cross(lineSpatialAvatars[li].component.second);
                Box3 b = BoundingBox(normalize(normal));
                static const double s = 0.2;
                b.minCorner = b.minCorner - Vec3(s, s, s);
                b.maxCorner = b.maxCorner + Vec3(s, s, s);
                return b;
            };

            RTreeWrapper<LineIndex, decltype(lookupLineNormal)> linesRTree(lookupLineNormal);
            for (auto & i : lineSpatialAvatars){
                linesRTree.insert(i.first);
            }

            // recognize incidence constraints between lines of different views
            IndexHashMap<std::pair<LineIndex, LineIndex>, Vec3> lineIncidenceRelations;

            for (auto & i : lineSpatialAvatars){
                auto li = i.first;
                auto & lineData = i.second;
                auto & views = _views;
                linesRTree.search(lookupLineNormal(li), [&li, &lineSpatialAvatars, &views, &lineIncidenceRelations](const LineIndex & relatedLi) -> bool {
                    if (li.viewHandle == relatedLi.viewHandle)
                        return true;
                    if (relatedLi < li) // make sure one relation is stored only once, avoid storing both a-b and b-a
                        return true;
                    auto & line1 = lineSpatialAvatars[li];
                    auto & line2 = lineSpatialAvatars[relatedLi];
                    if (line1.claz != line2.claz)
                        return true;
                    auto normal1 = normalize(line1.component.first.cross(line1.component.second));
                    auto normal2 = normalize(line2.component.first.cross(line2.component.second));
                    auto & vd1 = views.data(li.viewHandle);
                    auto & vd2 = views.data(relatedLi.viewHandle);
                    if (std::min(AngleBetweenDirections(normal1, normal2), 
                            AngleBetweenDirections(normal1, -normal2)) <
                        vd1.lineNet->params().incidenceDistanceVerticalDirectionThreshold / vd1.camera.focal() +
                        vd2.lineNet->params().incidenceDistanceVerticalDirectionThreshold / vd2.camera.focal()){
                        
                        auto nearest = DistanceBetweenTwoLines(line1.component, line2.component);
                        auto relationCenter = (nearest.second.first.position + nearest.second.second.position) / 2.0;
                        relationCenter /= norm(relationCenter);
                        
                        lineIncidenceRelations[std::make_pair(li, relatedLi)] = relationCenter;
                    }
                    return true;
                });
            }
            

            // generate sampled points for line-region connections
            IndexHashMap<std::pair<RegionIndex, LineIndex>, std::vector<Vec3>> 
                regionLineIntersectionSampledPoints;
            
            static const int extendSize = 3;
            std::vector<int> dx, dy;
            dx.reserve(2 * extendSize + 1);
            dy.reserve(2 * extendSize + 1);
            for (int a = -extendSize; a <= extendSize; a++){
                for (int b = -extendSize; b <= extendSize; b++){
                    dx.push_back(a);
                    dy.push_back(b);
                }
            }

            for (auto & vd : _views.elements<0>()){                
                RegionIndex ri;
                ri.viewHandle = vd.topo.hd;
                
                LineIndex li;
                li.viewHandle = vd.topo.hd;

                const Image & segmentedRegions = vd.data.regionNet->segmentedRegions();
                auto & cam = vd.data.camera;
                
                for (auto & ld : vd.data.lineNet->lines().elements<0>()){
                    li.handle = ld.topo.hd;
                    
                    auto & line = ld.data.line.component;
                    auto lineDir = normalize(line.direction());
                    double sampleStep = _params.samplingStepLengthOnLines;
                    int sampledNum = static_cast<int>(std::floor(line.length() / sampleStep));
                    
                    for (int i = 0; i < sampledNum; i++){
                        auto sampledPoint = line.first + lineDir * i * sampleStep;                        
                        
                        std::set<int32_t> rhids;
                        for (int k = 0; k < dx.size(); k++){
                            int x = BoundBetween(static_cast<int>(std::round(sampledPoint[0] + dx[k])), 0, segmentedRegions.cols - 1);
                            int y = BoundBetween(static_cast<int>(std::round(sampledPoint[1] + dy[k])), 0, segmentedRegions.rows - 1);
                            PixelLoc p(x, y);
                            rhids.insert(segmentedRegions.at<int32_t>(p));
                        }
                        
                        for (int32_t rhid : rhids){
                            ri.handle = RegionsNet::RegionHandle(rhid);
                            regionLineIntersectionSampledPoints[std::make_pair(ri, li)]
                                .push_back(normalize(cam.spatialDirection(sampledPoint)));
                        }
                    }
                }
            }


            IF_DEBUG_USING_VISUALIZERS {
                // visualize connections between regions and lines
                std::unordered_map<ViewHandle, vis::Visualizer2D, HandleHasher<AtLevel<0>>> vizs;
                for (auto & vd : _views.elements<0>()){
                    //vis::Visualizer2D viz(vd.data.regionNet->image);
                    int height = vd.data.regionNet->image().rows;
                    int width = vd.data.regionNet->image().cols;

                    Image coloredOutput(vd.data.regionNet->segmentedRegions().size(), CV_8UC3);
                    std::vector<cv::Vec<uint8_t, 3>> colors(vd.data.regionNet->regions().internalElements<0>().size());
                    std::generate(colors.begin(), colors.end(), [](){
                        return cv::Vec<uint8_t, 3>(uint8_t(std::rand() % 256),
                            uint8_t(std::rand() % 256),
                            uint8_t(std::rand() % 256));
                    });
                    for (int y = 0; y < height; y++) {
                        for (int x = 0; x < width; x++) {
                            coloredOutput.at<cv::Vec<uint8_t, 3>>(cv::Point(x, y)) =
                                colors[vd.data.regionNet->segmentedRegions().at<int32_t>(cv::Point(x, y))];
                        }
                    }
                    vizs[vd.topo.hd].setImage(vd.data.regionNet->image());
                    vizs[vd.topo.hd].params.alphaForNewImage = 0.5;
                    vizs[vd.topo.hd] << coloredOutput;
                }

                    //viz.params.alphaForNewImage = 0.8f;
                    //viz << coloredOutput;

                    //viz.params.thickness = 1;
                    //viz.params.color = vis::ColorFromTag(vis::ColorTag::Black);
                    
                for (auto & lineIdRi : regionLineIntersectionSampledPoints){
                    auto & ri = lineIdRi.first.first;
                    auto & li = lineIdRi.first.second;
                    auto & cline2 = _views.data(li.viewHandle).lineNet->lines().data(li.handle).line;
                    auto & cam = _views.data(ri.viewHandle).camera;
                    auto & viz = vizs[ri.viewHandle];

                    /*viz << vis::manip2d::SetColor(vis::ColorTag::White)
                        << vis::manip2d::SetThickness(2);*/
                    //viz << cline2;

                    viz << vis::manip2d::SetColor(vis::ColorTag::Black)
                        << vis::manip2d::SetThickness(1);
                    auto & regionCenter = _views.data(ri.viewHandle).regionNet->regions().data(ri.handle).center;
                    for (auto & d : lineIdRi.second){
                        auto p = cam.screenProjection(d);
                        viz << Line2(regionCenter, p);
                    }
                }

                for (auto & viz : vizs){
                    viz.second << vis::manip2d::Show();
                }
                
            }



            using namespace deriv;
            ExpressionGraph graph;
            
            // append to constraint graph
            ConstraintGraph constraints;
            auto & _constraints = constraints;
            _constraints.internalElements<0>().reserve(_constraints.internalElements<0>().capacity() + 
                regionsRTree.size());
            _constraints.internalElements<1>().reserve(_constraints.internalElements<1>().capacity() + 
                overlappedRegionIndexPairs.size() + regionsRTree.size());

            // add region components
            IndexHashMap<RegionIndex, HandleAtLevel<0>> ri2Handle;
            for (auto & r : regionSpatialContours) {
                ComponentData cd(ComponentData::Type::Region);
                //cd.asRegion.isVoid = false;
                cd.asRegion.regionHandle = r.first.handle;
                cd.asRegion.viewHandle = r.first.viewHandle;
                auto & cam = _views.data(r.first.viewHandle).camera;
                auto direction = 
                    cam.spatialDirection(_views.data(r.first.viewHandle).regionNet->regions().data(r.first.handle).center);
                cd.asRegion.thetaExpr = OptimizibleExpression<Eigen::Vector3d>(
                    deriv::MakeEigenMat(normalize(direction)), 
                    graph);
                ri2Handle[r.first] = _constraints.add(cd);
            }
            // add line components
            IndexHashMap<LineIndex, ComponentHandle> li2Handle;
            for (auto & l : lineSpatialAvatars) {
                ComponentData cd(ComponentData::Type::Line);
                cd.asLine.lineHandle = l.first.handle;
                cd.asLine.viewHandle = l.first.viewHandle;
                cd.asLine.etaExpr = OptimizibleExpression<double>(1.0, graph);
                li2Handle[l.first] = _constraints.add(cd);
            }
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
                    cd.asRegionConnectivity.viewHandle = vd.topo.hd;
                    cd.asRegionConnectivity.boundaryHandle = regionBoundaryData.topo.hd;
                    _constraints.add<1>({ rh1, rh2 }, cd);
                }
            }
            // add line connectivity
            for (auto & vd : _views.elements<0>()){
                auto & lines = *vd.data.lineNet;
                for (auto & lineRelationData : lines.lines().elements<1>()){
                    auto & lids = lineRelationData.topo.lowers;
                    LineIndex li1 = { vd.topo.hd, lids[0] };
                    LineIndex li2 = { vd.topo.hd, lids[1] };
                    auto lh1 = li2Handle[li1];
                    auto lh2 = li2Handle[li2];

                    ConstraintData cd(ConstraintData::Type::LineConnectivity);
                    cd.asLineConnectivity.viewHandle = vd.topo.hd;
                    cd.asLineConnectivity.lineRelationHandle = lineRelationData.topo.hd;
                    _constraints.add<1>({ lh1, lh2 }, cd);
                }
            }
            // add line inter view incidence
            for (auto & iv : lineIncidenceRelations){
                auto h1 = li2Handle[iv.first.first];
                auto h2 = li2Handle[iv.first.second];
                ConstraintData cd(ConstraintData::Type::LineInterViewIncidence);
                cd.asLineInterViewIncidence.relationCenter = iv.second;
                _constraints.add<1>({ h1, h2 }, cd);
            }
            // add region line connectivity constraints
            for (auto & regionLineData : regionLineIntersectionSampledPoints) {
                if (regionLineData.second.empty())
                    continue;
                ConstraintData cd(ConstraintData::Type::RegionLineConnectivity);
                cd.asRegionLineConnectivity.sampledPoints = regionLineData.second;
                
                auto lineHandle = li2Handle[regionLineData.first.second];
                auto regionHandle = ri2Handle[regionLineData.first.first];
                _constraints.add<1>({ lineHandle, regionHandle }, cd);
            }

            std::cout << "component num: " << _constraints.internalElements<0>().size() << std::endl;
            std::cout << "constraint num: " << _constraints.internalElements<1>().size() << std::endl;


            IF_DEBUG_USING_VISUALIZERS{
                // visualize constraint graph
                auto getRegionCenter = [this](const RegionComponentData & rd) -> Point2 {
                    auto & regionData = _views.data(rd.viewHandle).regionNet->regions().data(rd.regionHandle);
                    return regionData.center;
                };
                std::map<ViewHandle, vis::Visualizer2D> viewVizs;
                for (auto & vd : _views.elements<0>()){
                    viewVizs[vd.topo.hd] = viewVizs[vd.topo.hd] << vd.data;
                    assert(!viewVizs[vd.topo.hd].image().empty());
                }
                std::map<std::pair<ViewHandle, ViewHandle>, Image> viewPairVizs;
                for (auto & cd : _constraints.elements<1>()){
                    if (cd.data.type == ConstraintData::Type::RegionOverlap){
                        auto & rd1 = _constraints.data(cd.topo.lowers[0]).asRegion;
                        auto & rd2 = _constraints.data(cd.topo.lowers[1]).asRegion;
                        ViewHandle v1 = rd1.viewHandle;
                        ViewHandle v2 = rd2.viewHandle;
                        
                        Point2 c1 = getRegionCenter(rd1);
                        Point2 c2 = getRegionCenter(rd2);
                        if (v2 < v1){
                            std::swap(v1, v2);
                            std::swap(c1, c2);
                        }
                        Image im1 = viewVizs[v1].image();
                        Image im2 = viewVizs[v2].image();
                        if (viewPairVizs.find(std::make_pair(v1, v2)) == viewPairVizs.end()){                           
                            assert(im1.size == im2.size);
                            assert(im1.type() == im2.type());
                            Image im12 = Image::zeros(im1.rows, im1.cols + im2.cols, im1.type());
                            im12(cv::Range(0, im1.rows), cv::Range(0, im1.cols)) += im1;
                            im12(cv::Range(0, im2.rows), cv::Range(im1.cols, im1.cols + im2.cols)) += im2;
                            viewPairVizs[std::make_pair(v1, v2)] = im12;
                        }

                        vis::Color color(rand() % 255, rand() % 255, rand() % 255, 1);
                        c2[0] += im1.cols;
                        cv::line(viewPairVizs[std::make_pair(v1, v2)], PixelLoc(c1[0], c1[1]), PixelLoc(c2[0], c2[1]),
                            color, 2);
                    }
                }

                for (auto & vp : viewPairVizs){
                    cv::imshow("overlapped regions", vp.second);
                    cv::waitKey();
                }
               
            }




            std::cout << "begin building constraints" << std::endl;
            
            
            // inference region orientations and spatial connectivity of boundaries


            std::array<Expression<Eigen::Vector3d>, 6> principleDirectionExprs;
            for (int i = 0; i < 3; i++) {
                const Vec3 & vp = _globalData.vanishingPoints[i];
                principleDirectionExprs[i * 2] = composeFunction(graph, [&vp]() {
                    return MakeEigenMat(vp / norm(vp)); 
                });
                principleDirectionExprs[i * 2 + 1] = composeFunction(graph, [&vp]() {
                    return MakeEigenMat(- vp / norm(vp));
                });
            }
            
            // component expressions
            for (auto & r : _constraints.elements<0>()) {
                if (r.data.type == ComponentData::Type::Region) {
                    // E_manh
                    std::array<Expression<double>, 6> angles;
                    for (int i = 0; i < 6; i++) {
                        angles[i] = AngleFunc(r.data.asRegion.thetaExpr.expression(), principleDirectionExprs[i]);
                    }
                    Expression<double> minAngle = deriv::minInRange(angles.begin(), angles.end());
                    double wManh = 1.0;
                    r.data.asRegion.manhattanEnergyExpr = HFunc(minAngle, M_PI / 6) * wManh;
                    auto thetaNorm = norm(r.data.asRegion.thetaExpr.expression());
                    r.data.reserveScaleEnergyExpr = cwiseSelect(thetaNorm - 1.0, 1.0, 1.0 / (thetaNorm + 1e-6));
                }
                else if (r.data.type == ComponentData::Type::Line) {
                    auto etaNorm = abs(r.data.asLine.etaExpr.expression());
                    r.data.reserveScaleEnergyExpr = cwiseSelect(etaNorm - 1.0, 1.0, 1.0 / (etaNorm + 1e-6));
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

                    auto thetaDiffExpr = rd1.asRegion.thetaExpr.expression() - rd2.asRegion.thetaExpr.expression();
                    auto thetaDiffSquaredSum = abs(dotProd(thetaDiffExpr, thetaDiffExpr));
                    
                    double wOverlap = 20.0;

                    // E_overlap
                    auto overlapRatio = con.data.asRegionOverlap.overlapRatio;
                    assert(overlapRatio > 0);
                    con.data.invalidityExpr = graph.addConst(0.0);
                    con.data.constraintEnergyExpr = wOverlap * thetaDiffSquaredSum;
                    if (overlapRatio < 0.15)
                        con.data.constraintEnergyExpr.disable();                    
                }
                // region connectivity
                else if (con.data.type == ConstraintData::Type::RegionConnectivity) {
                    assert(rd1.type == ComponentData::Type::Region &&
                        rd2.type == ComponentData::Type::Region &&
                        "invalid component type!");

                    auto & camera = _views.data(con.data.asRegionConnectivity.viewHandle).camera;
                    auto & boundaryData = _views.data(con.data.asRegionConnectivity.viewHandle).regionNet
                        ->regions().data(con.data.asRegionConnectivity.boundaryHandle);
                    auto & sampledPoints = boundaryData.sampledPoints;
                    auto straightness = boundaryData.straightness;
                    auto & fittedLine = boundaryData.fittedLine;

                    // get sampled direction expressions
                    std::vector<Expression<Eigen::Vector3d>> sampledDirectionExprs; // normalized directions
                    sampledDirectionExprs.reserve(boundaryData.length / 5);
                    for (auto & ps : sampledPoints) {
                        for (auto & p : ps) {
                            auto sampledDirectionExpr = composeFunction(graph, [&p, &camera]() {
                                auto direction = camera.spatialDirection(p);
                                assert(norm(direction) != 0);
                                return MakeEigenMat(direction / norm(direction));
                            });
                            sampledDirectionExprs.push_back(sampledDirectionExpr);
                        }
                    }
                    // get depth expressions of sampled direction
                    std::vector<Expression<double>> lambdasR1(sampledDirectionExprs.size()),
                        lambdasR2(sampledDirectionExprs.size());
                    std::vector<EHandle> lambdaSquaredDiffHandles(sampledDirectionExprs.size());
                    for (int i = 0; i < sampledDirectionExprs.size(); i++) {
                        lambdasR1[i] = dotProd(rd1.asRegion.thetaExpr.expression(), rd1.asRegion.thetaExpr.expression()) /
                            dotProd(rd1.asRegion.thetaExpr.expression(), sampledDirectionExprs[i]);
                        lambdasR2[i] = dotProd(rd2.asRegion.thetaExpr.expression(), rd2.asRegion.thetaExpr.expression()) /
                            dotProd(rd2.asRegion.thetaExpr.expression(), sampledDirectionExprs[i]);
                        lambdaSquaredDiffHandles[i] = Square(lambdasR1[i] - lambdasR2[i]).handle();
                    }
                    EHandle lambdaSquaredDiffSumHandle = HSum<double>(&graph, lambdaSquaredDiffHandles);
                    auto lambdaSquaredDiffSum = graph.as<double>(lambdaSquaredDiffSumHandle);


                    // E_connect
                    double wConnect = 10.0;
                    assert(straightness >= 0.0 && straightness <= 1.0);
                    con.data.invalidityExpr = lambdaSquaredDiffSum / sampledPoints.size();
                    con.data.constraintEnergyExpr = wConnect * con.data.invalidityExpr * sampledPoints.size();


                    IF_DEBUG_USING_VISUALIZERS{
                        double v = con.data.constraintEnergyExpr.expression().execute();
                        assert(!isnan(v) && v >= 0);
                    }

                }
                // line line connectivity
                else if (con.data.type == ConstraintData::Type::LineConnectivity) {

                    assert(rd1.type == ComponentData::Type::Line &&
                        rd2.type == ComponentData::Type::Line &&
                        "invalid component type!");

                    auto & lineConnectivityData = con.data.asLineConnectivity;
                    assert(lineConnectivityData.viewHandle == rd1.asLine.viewHandle);
                    assert(lineConnectivityData.viewHandle == rd2.asLine.viewHandle);

                    auto & cam = _views.data(con.data.asLineConnectivity.viewHandle).camera;
                    auto & linesNet = *_views.data(con.data.asLineConnectivity.viewHandle).lineNet;
                    auto & lineRelationData = linesNet.lines().data(lineConnectivityData.lineRelationHandle);
                    auto & lineVotingDistribution = linesNet.lineVotingDistribution()(
                        PixelLoc(lineRelationData.relationCenter[0], lineRelationData.relationCenter[1]));

                    // compute junction weights
                    auto & v = lineVotingDistribution;
                    double junctionWeight = 0;
                    {
                        // Y
                        double Y = 0.0;
                        for (int s = 0; s < 2; s++){
                            Y += v(0, s) * v(1, s) * v(2, s) * DiracDelta(v(0, 1 - s) + v(1, 1 - s) + v(2, 1 - s));
                        }

                        // W
                        double W = 0.0;
                        for (int i = 0; i < 3; i++){
                            for (int j = 0; j < 3; j++){
                                if (i == j)
                                    continue;
                                int k = 3 - i - j;
                                for (int s = 0; s < 2; s++){
                                    W += v(i, s) * v(j, 1 - s) * v(k, 1 - s) * DiracDelta(v(i, 1 - s) + v(j, s) + v(k, s));
                                }
                            }
                        }

                        // K
                        double K = 0.0;
                        for (int i = 0; i < 3; i++){
                            for (int j = 0; j < 3; j++){
                                if (i == j)
                                    continue;
                                int k = 3 - i - j;
                                K += v(i, 0) * v(i, 1) * v(j, 0) * v(k, 1) * DiracDelta(v(j, 1) + v(k, 0));
                                K += v(i, 0) * v(i, 1) * v(j, 1) * v(k, 0) * DiracDelta(v(j, 0) + v(k, 1));
                            }
                        }

                        // compute X junction
                        double X = 0.0;
                        for (int i = 0; i < 3; i++){
                            for (int j = 0; j < 3; j++){
                                if (i == j)
                                    continue;
                                int k = 3 - i - j;
                                X += v(i, 0) * v(i, 1) * v(j, 0) * v(j, 1) * DiracDelta(v(k, 0) + v(k, 1));
                            }
                        }

                        // compute T junction
                        double T = 0.0;
                        for (int i = 0; i < 3; i++){
                            for (int j = 0; j < 3; j++){
                                if (i == j)
                                    continue;
                                int k = 3 - i - j;
                                T += v(i, 0) * v(i, 1) * v(j, 0) * DiracDelta(v(j, 1) + v(k, 0) + v(k, 1));
                                T += v(i, 0) * v(i, 1) * v(j, 1) * DiracDelta(v(j, 0) + v(k, 0) + v(k, 1));
                            }
                        }

                        // compute L junction
                        double L = 0.0;
                        for (int i = 0; i < 3; i++){
                            for (int j = 0; j < 3; j++){
                                if (i == j)
                                    continue;
                                int k = 3 - i - j;
                                for (int a = 0; a < 2; a++){
                                    int nota = 1 - a;
                                    for (int b = 0; b < 2; b++){
                                        int notb = 1 - b;
                                        L += v(i, a) * v(j, b) * DiracDelta(v(i, nota) + v(j, notb) + v(k, 0) + v(k, 1));
                                    }
                                }
                            }
                        }

                        //std::cout << " Y-" << Y << " W-" << W << " K-" << K << 
                        //    " X-" << X << " T-" << T << " L-" << L << std::endl; 
                        static const double threshold = 1e-4;
                        if (Y > threshold){
                            junctionWeight += 5.0;
                        }
                        else if (W > threshold){
                            junctionWeight += 5.0;
                        }
                        else if (L > threshold){
                            junctionWeight += 4.0;
                        }
                        else if (K > threshold){
                            junctionWeight += 3.0;
                        }
                        else if (X > threshold){
                            junctionWeight += 5.0;
                        }
                        else if (T > threshold){
                            junctionWeight += 0.1;
                        }
                    }
                    if (lineRelationData.type == LinesNet::LineRelationData::Type::Incidence)
                        junctionWeight += 10.0;

                    auto & line1 = linesNet.lines().data(rd1.asLine.lineHandle).line;
                    auto & line2 = linesNet.lines().data(rd2.asLine.lineHandle).line;

                    auto sampledDir = cam.spatialDirection(lineRelationData.relationCenter);

                    double ratio1 = ComputeDepthRatioOfPointOnSpatialLine(
                        cam.spatialDirection(line1.component.first), 
                        sampledDir, 
                        _globalData.vanishingPoints[line1.claz]);
                    double ratio2 = ComputeDepthRatioOfPointOnSpatialLine(
                        cam.spatialDirection(line2.component.first),
                        sampledDir,
                        _globalData.vanishingPoints[line2.claz]);

                    assert(ratio1 > 0);
                    assert(ratio2 > 0);

                    auto lambda1 = rd1.asLine.etaExpr.expression() * ratio1;
                    auto lambda2 = rd2.asLine.etaExpr.expression() * ratio2;

                    double wLineConnect = 8.0;
                    auto lambdaSquaredDiff = Square(lambda1 - lambda2);
                    con.data.invalidityExpr = lambdaSquaredDiff;
                    con.data.constraintEnergyExpr = (wLineConnect * con.data.invalidityExpr * junctionWeight).cast<double>();
                    
                }
                // line inter view incidence
                else if (con.data.type == ConstraintData::Type::LineInterViewIncidence) {
                    assert(rd1.type == ComponentData::Type::Line && 
                        rd2.type == ComponentData::Type::Line);

                    auto & line1 = lineSpatialAvatars[LineIndex{ rd1.asLine.viewHandle, rd1.asLine.lineHandle }];
                    auto & line2 = lineSpatialAvatars[LineIndex{ rd2.asLine.viewHandle, rd2.asLine.lineHandle }];

                    auto & relationCenter = con.data.asLineInterViewIncidence.relationCenter;

                    double ratio1 = ComputeDepthRatioOfPointOnSpatialLine(
                        line1.component.first, relationCenter, _globalData.vanishingPoints[line1.claz]);
                    double ratio2 = ComputeDepthRatioOfPointOnSpatialLine(
                        line2.component.first, relationCenter, _globalData.vanishingPoints[line2.claz]);

                    assert(ratio1 > 0);
                    assert(ratio2 > 0);

                    auto lambda1 = rd1.asLine.etaExpr.expression() * ratio1;
                    auto lambda2 = rd2.asLine.etaExpr.expression() * ratio2;

                    double wLineConnect = 10.0;
                    auto lambdaSquaredDiff = Square(lambda1 - lambda2);
                    con.data.invalidityExpr = lambdaSquaredDiff;
                    con.data.constraintEnergyExpr = (wLineConnect * con.data.invalidityExpr).cast<double>();

                }
                // region line connectivity
                else if (con.data.type == ConstraintData::Type::RegionLineConnectivity) {

                    assert((rd1.type == ComponentData::Type::Line && 
                        rd2.type == ComponentData::Type::Region ||
                        rd1.type == ComponentData::Type::Region && 
                        rd2.type == ComponentData::Type::Line) &&
                        "invalid component type!");

                    auto prrd1 = &rd1, prrd2 = &rd2;
                    if (rd1.type == ComponentData::Type::Line){
                        std::swap(prrd1, prrd2);
                    }
                    auto & rrd1 = *prrd1; // region
                    auto & rrd2 = *prrd2; // linestruct

                    assert(rrd1.type == ComponentData::Type::Region && 
                        rrd2.type == ComponentData::Type::Line &&
                        "invalid component type!");

                    auto & sampledPoints = con.data.asRegionLineConnectivity.sampledPoints;

                    // get sampled direction expressions
                    std::vector<Expression<Eigen::Vector3d>> sampledDirectionExprs;
                    sampledDirectionExprs.reserve(sampledPoints.size());
                    for (auto & p : sampledPoints) {
                        auto sampledDirectionExpr = composeFunction(graph, [&p]() {
                            return MakeEigenMat(normalize(p));
                        });
                        sampledDirectionExprs.push_back(sampledDirectionExpr);
                    }
                    
                    // get depth expressions of sampled direction
                    std::vector<Expression<double>> lambdasR1(sampledDirectionExprs.size()),
                        lambdasR2(sampledDirectionExprs.size());
                    std::vector<EHandle> lambdaSquaredDiffHandles(sampledDirectionExprs.size());
                    for (int i = 0; i < sampledDirectionExprs.size(); i++) {
                        lambdasR1[i] = dotProd(rrd1.asRegion.thetaExpr.expression(), rrd1.asRegion.thetaExpr.expression()) /
                            dotProd(rrd1.asRegion.thetaExpr.expression(), sampledDirectionExprs[i]);

                        auto & dir = sampledPoints[i]; // we don't need expression here
                        auto & line = _views.data(rrd2.asLine.viewHandle).lineNet->lines().data(rrd2.asLine.lineHandle).line;
                        auto & cam = _views.data(rrd2.asLine.viewHandle).camera;
                        double ratio = ComputeDepthRatioOfPointOnSpatialLine(
                            cam.spatialDirection(line.component.first), dir,
                            _globalData.vanishingPoints[line.claz]);

                        assert(ratio > 0);

                        lambdasR2[i] = (rrd2.asLine.etaExpr.expression() * ratio).cast<double>();
                        lambdaSquaredDiffHandles[i] = Square(lambdasR1[i] - lambdasR2[i]).handle();
                    }
                    EHandle lambdaSquaredDiffSumHandle = HSum<double>(&graph, lambdaSquaredDiffHandles);
                    auto lambdaSquaredDiffSum = graph.as<double>(lambdaSquaredDiffSumHandle);

                    // E_connect
                    double wConnect = 4.0;
                    con.data.invalidityExpr = lambdaSquaredDiffSum / sampledPoints.size();
                    con.data.constraintEnergyExpr = wConnect * con.data.invalidityExpr * sampledPoints.size();                    
                }

            }

            std::cout << "individual energies defined" << std::endl;

            // sum all energy
            std::vector<EHandle> allManhattanEnergyHandles;

            std::map<ConstraintData::Type, std::vector<EHandle>> allConstraintEnergyHandles;
            std::map<ComponentData::Type, std::vector<EHandle>> allReserveScaleEnergyHandles;
            
           
            // add all E_manh and reserve scale energies
            for (auto & rd : _constraints.elements<0>()){
                if (rd.data.type == ComponentData::Type::Region){
                    allManhattanEnergyHandles.push_back(rd.data.asRegion.manhattanEnergyExpr.expression().handle());
                    allReserveScaleEnergyHandles[rd.data.type].push_back(rd.data.reserveScaleEnergyExpr.expression().handle());
                }
                else if (rd.data.type == ComponentData::Type::Line){
                    allReserveScaleEnergyHandles[rd.data.type].push_back(rd.data.reserveScaleEnergyExpr.expression().handle());
                }
            }

            std::cout << "component energies defined" << std::endl;

            // add all EnergyHandles
            for (auto & cd : _constraints.elements<1>()){
                allConstraintEnergyHandles[cd.data.type]
                    .push_back(cd.data.constraintEnergyExpr.expression().handle());
            }

            std::cout << "constraint energies defined" << std::endl;

            // sum all
            DisableableExpression<double> allManhattanEnergySumExpr = 
                graph.as<double>(HSum<double>(&graph, allManhattanEnergyHandles));
            if (allManhattanEnergySumExpr.expression().isInValid()){
                std::cout << "invalid!" << std::endl;
            }
            std::map<ComponentData::Type, DisableableExpression<double>> allReserveScaleEnergySumExpr;
            for (auto & allHandles : allReserveScaleEnergyHandles){
                allReserveScaleEnergySumExpr[allHandles.first] =
                    graph.as<double>(HSum<double>(&graph, allHandles.second));
                if (allReserveScaleEnergySumExpr[allHandles.first].expression().isInValid()){
                    std::cout << "invalid!" << std::endl;
                }
            }
            std::map<ConstraintData::Type, DisableableExpression<double>> allConstraintEnergySumExpr;
            for (auto & allHandles : allConstraintEnergyHandles){
                allConstraintEnergySumExpr[allHandles.first] = 
                    graph.as<double>(HSum<double>(&graph, allHandles.second));
                if (allConstraintEnergySumExpr[allHandles.first].expression().isInValid()){
                    std::cout << "invalid!" << std::endl;
                }
            }

            std::cout << "expressions defined" << std::endl;

            std::cout << "final energy expression defined" << std::endl;

            std::cout << "done building constraints" << std::endl;


            auto startTime = std::chrono::high_resolution_clock::now();

            


            //// optimize lines
            OptimizeConstraintGraphUsingEnergyExpression(
                allConstraintEnergySumExpr[ConstraintData::Type::LineConnectivity].expression() +
                allConstraintEnergySumExpr[ConstraintData::Type::LineInterViewIncidence].expression() +
                allReserveScaleEnergySumExpr[ComponentData::Type::Line].expression(), 
                _constraints, 1e-2, 0.1, 800, ///
                [this, &lineSpatialAvatars, &constraints](int i, double energyVal, double energyValChange) -> bool { 

                if (abs(energyValChange) <= 1e-2){
                    std::vector<Classified<Line3>> rectifiedLines;
                    for (auto & cd : constraints.elements<0>()){
                        if (cd.data.type == ComponentData::Type::Line){
                            auto line = lineSpatialAvatars[LineIndex{ cd.data.asLine.viewHandle, cd.data.asLine.lineHandle }];
                            auto ratio = ComputeDepthRatioOfPointOnSpatialLine(line.component.first,
                                line.component.second, globalData().vanishingPoints[line.claz]);
                            double eta = cd.data.asLine.etaExpr.expression().result();
                            line.component.first = normalize(line.component.first) * eta;
                            line.component.second = normalize(line.component.second) * eta * ratio;
                            rectifiedLines.push_back(line);
                        }
                    }

                    vis::Visualizer3D viz;
                    viz << vis::manip3d::SetWindowName("optimized lines")
                        << vis::manip3d::SetDefaultColor(vis::ColorTag::Yellow)
                        << vis::manip3d::SetColorTableDescriptor(vis::ColorTableDescriptor::RGB)
                        << rectifiedLines
                        << vis::manip3d::AutoSetCamera
                        << vis::manip3d::Show(true);

                    return false;
                }

                return true;


            });

            // 
            std::vector<ComponentHandle> lineCompHandles;
            lineCompHandles.reserve(_constraints.internalElements<0>().size());
            for (auto & cd : _constraints.elements<0>()){
                if (cd.exists && cd.data.type == ComponentData::Type::Line){
                    lineCompHandles.push_back(cd.topo.hd);
                }
            }
            std::vector<ConstraintHandle> lineConsHandles;
            lineConsHandles.reserve(_constraints.internalElements<1>().size());
            for (auto & cd : _constraints.elements<1>()){
                if (cd.exists &&
                    (cd.data.type == ConstraintData::Type::LineConnectivity) ||
                    (cd.data.type == ConstraintData::Type::LineInterViewIncidence)){
                    lineConsHandles.push_back(cd.topo.hd);
                }
            }
            std::vector<ConstraintHandle> mstLineConsHandles;
            auto componentGetter = [this, &constraints](const ConstraintHandle & h) -> std::pair < ComponentHandle, ComponentHandle > {
                return std::make_pair(constraints.topo(h).lowers[0], constraints.topo(h).lowers[1]);
            };
            MinimumSpanningTree(lineCompHandles.begin(), lineCompHandles.end(), 
                lineConsHandles.begin(), lineConsHandles.end(), std::back_inserter(mstLineConsHandles), componentGetter, 
                [this, &constraints](ConstraintHandle a, ConstraintHandle b){
                return constraints.data(a).invalidityExpr.result() < constraints.data(b).invalidityExpr.result();
            });
            // disable all line connection/interview incidences not contained in mst
            for (auto & cd : _constraints.elements<1>()){
                if (cd.exists &&
                    (cd.data.type == ConstraintData::Type::LineConnectivity) ||
                    (cd.data.type == ConstraintData::Type::LineInterViewIncidence)){
                    cd.data.constraintEnergyExpr.disable();
                }
            }
            for (auto & h : mstLineConsHandles){
                _constraints.data(h).constraintEnergyExpr.enable();
            }



            //// optimize lines
            OptimizeConstraintGraphUsingEnergyExpression(
                allConstraintEnergySumExpr[ConstraintData::Type::LineConnectivity].expression() +
                allConstraintEnergySumExpr[ConstraintData::Type::LineInterViewIncidence].expression() +
                allReserveScaleEnergySumExpr[ComponentData::Type::Line].expression(),
                _constraints, 1e-2, 0.01, 800, ///
                [this, &lineSpatialAvatars](int i, double energyVal, double energyValChange) -> bool {

                if (abs(energyValChange) <= 1e-3){
                    return false;
                }
                return true;
            });

            {
                std::vector<Classified<Line3>> rectifiedLines;
                for (auto & cd : _constraints.elements<0>()){
                    if (cd.data.type == ComponentData::Type::Line){
                        auto line = lineSpatialAvatars[LineIndex{ cd.data.asLine.viewHandle, cd.data.asLine.lineHandle }];
                        auto ratio = ComputeDepthRatioOfPointOnSpatialLine(line.component.first,
                            line.component.second, globalData().vanishingPoints[line.claz]);
                        double eta = cd.data.asLine.etaExpr.expression().result();
                        line.component.first = normalize(line.component.first) * eta;
                        line.component.second = normalize(line.component.second) * eta * ratio;
                        rectifiedLines.push_back(line);
                    }
                }

                vis::Visualizer3D viz;
                viz << vis::manip3d::SetWindowName("optimized lines")
                    << vis::manip3d::SetDefaultColor(vis::ColorTag::Yellow)
                    << vis::manip3d::SetColorTableDescriptor(vis::ColorTableDescriptor::RGB)
                    << rectifiedLines
                    << vis::manip3d::AutoSetCamera
                    << vis::manip3d::Show(true);
            }




            //using ComponentTriplet = ConstraintGraph::TripletType<0>;

            //// find connected components of lines based on line connections & line inter-view incidences
            //std::vector<ComponentHandle> lineCompHandles;
            //lineCompHandles.reserve(_constraints.internalElements<0>().size());
            //for (auto & cd : _constraints.elements<0>()){
            //    if (cd.exists && cd.data.type == ComponentData::Type::Line){
            //        lineCompHandles.push_back(cd.topo.hd);
            //    }
            //}
            //auto neighborLinesContainerGetter = [this](const ComponentHandle & ch){
            //    std::vector<ComponentHandle> neighbors;
            //    neighbors.reserve(5);
            //    assert(_constraints.data(ch).type == ComponentData::Type::Line);
            //    auto & topo = _constraints.topo(ch);
            //    for (auto & conh : topo.uppers){
            //        auto & cdata = _constraints.data(conh);
            //        if (cdata.type != ConstraintData::Type::LineConnectivity &&
            //            cdata.type != ConstraintData::Type::LineInterViewIncidence)
            //            continue;
            //        //std::cout << "!!!!!!!!!!!!!!!!!!!" << std::endl;
            //        auto & ctopo = _constraints.topo(conh);
            //        if (ctopo.lowers[0] == ch)
            //            neighbors.push_back(ctopo.lowers[1]);
            //        if (ctopo.lowers[1] == ch)
            //            neighbors.push_back(ctopo.lowers[0]);
            //    }
            //    return neighbors;
            //};
            //std::map<int, std::list<ComponentHandle>> ccid2LineHandles;
            //int nLineCC = ConnectedComponents(lineCompHandles.begin(), lineCompHandles.end(), neighborLinesContainerGetter, 
            //    [this, &ccid2LineHandles](const ComponentHandle & ch, int cid){
            //    _constraints.data(ch).asLine.connectedComponentId = cid;
            //    ccid2LineHandles[cid].push_back(ch);
            //});
            //std::cout << "Connected Components Num of Lines: " << nLineCC << std::endl;
            //for (auto & cc : ccid2LineHandles){
            //    std::cout << "[" << cc.first << "] num: " << cc.second.size();
            //    if (cc.first % 10 == 9)
            //        std::cout << std::endl;
            //}

            //std::vector<EHandle> allCCFirstLineReserveScaleEnergyHandles;
            //// refactor the etaExprs in each line component
            //for (auto & chs : ccid2LineHandles){
            //    assert(!chs.second.empty());
            //    auto firstCh = chs.second.front();
            //    auto & firstLineData = _constraints.data(firstCh).asLine;
            //    allCCFirstLineReserveScaleEnergyHandles.push_back(_constraints.data(firstCh).reserveScaleEnergyExpr.expression().handle());
            //    for (auto ch : chs.second){
            //        auto & currentLineData = _constraints.data(ch).asLine;
            //        if (ch == firstCh){
            //            currentLineData.etaRatio = 1.0;
            //            continue;
            //        }
            //        currentLineData.etaRatio = currentLineData.etaExpr.expression().result() /
            //            firstLineData.etaExpr.expression().result();
            //        auto unitedEtaExpr = firstLineData.etaExpr.expression() * currentLineData.etaRatio;
            //        // replace all eta expressions
            //        currentLineData.etaExpr.expression().replacedWithWhenUsedAsInputs(unitedEtaExpr);
            //        currentLineData.etaExpr.freeze();
            //        //currentLineData.etaExpr = unitedEtaExpr; // replaced!!!
            //    }
            //}
            //auto allCCFirstLineReserveScaleEnergySumExpr =
            //    graph.as<double>(HSum<double>(&graph, allCCFirstLineReserveScaleEnergyHandles));






            //// find connected components of regions based on region overlaps & region connections
            //std::vector<ComponentHandle> regionCompHandles;
            //regionCompHandles.reserve(_constraints.internalElements<0>().size());
            //for (auto & cd : _constraints.elements<0>()){
            //    if (cd.exists && cd.data.type == ComponentData::Type::Region){
            //        regionCompHandles.push_back(cd.topo.hd);
            //    }
            //}
            //auto neighborRegionsContainerGetter = [this](const ComponentHandle & ch){
            //    std::vector<ComponentHandle> neighbors;
            //    neighbors.reserve(5);
            //    assert(_constraints.data(ch).type == ComponentData::Type::Region);
            //    auto & topo = _constraints.topo(ch);
            //    for (auto & conh : topo.uppers){
            //        auto & cdata = _constraints.data(conh);
            //        if (cdata.type != ConstraintData::Type::RegionConnectivity &&
            //            cdata.type != ConstraintData::Type::RegionOverlap)
            //            continue;
            //        if (cdata.type == ConstraintData::Type::RegionOverlap){
            //            if (cdata.asRegionOverlap.overlapRatio < 0.15){
            //                continue;
            //            }
            //        }
            //        if (cdata.type == ConstraintData::Type::RegionConnectivity){
            //            auto & regionCon = cdata.asRegionConnectivity;
            //            auto & boundaryData = _views.data(regionCon.viewHandle).regionNet->regions().data(regionCon.boundaryHandle);
            //            if (!(boundaryData.length > 50 && boundaryData.straightness < 0.01))
            //                continue;
            //        }
            //        ////std::cout << "!!!!!!!!!!!!!!!!!!!" << std::endl;
            //        //auto & ctopo = _constraints.topo(conh);
            //        //if (ctopo.lowers[0] == ch)
            //        //    neighbors.push_back(ctopo.lowers[1]);
            //        //if (ctopo.lowers[1] == ch)
            //        //    neighbors.push_back(ctopo.lowers[0]);
            //    }
            //    return neighbors;
            //};
            //std::map<int, std::list<ComponentHandle>> ccid2RegionHandles;
            //int nRegionCC = ConnectedComponents(regionCompHandles.begin(), regionCompHandles.end(), neighborRegionsContainerGetter,
            //    [this, &ccid2RegionHandles](const ComponentHandle & ch, int cid){
            //    _constraints.data(ch).asRegion.connectedComponentId = cid;
            //    ccid2RegionHandles[cid].push_back(ch);
            //});
            //std::cout << "Connected Components Num of Regions: " << nRegionCC << std::endl;
            //for (auto & cc : ccid2RegionHandles){
            //    std::cout << "[" << cc.first << "] num: " << cc.second.size();
            //    if (cc.first % 10 == 9)
            //        std::cout << std::endl;
            //}


            //std::vector<EHandle> allCCFirstRegionReserveScaleEnergyHandles;
            //// refactor the thetaExprs in each reigon component
            //for (auto & chs : ccid2RegionHandles){
            //    assert(!chs.second.empty());
            //    auto firstCh = chs.second.front();
            //    auto & firstRegionData = _constraints.data(firstCh).asRegion;
            //    allCCFirstRegionReserveScaleEnergyHandles.push_back(_constraints.data(firstCh).reserveScaleEnergyExpr.expression().handle());
            //    for (auto ch : chs.second){
            //        if (ch == firstCh)
            //            continue;
            //        auto & currentRegionData = _constraints.data(ch).asRegion;
            //        auto unitedThetaExpr = firstRegionData.thetaExpr.expression();
            //        currentRegionData.thetaExpr.expression().replacedWithWhenUsedAsInputs(unitedThetaExpr);
            //        currentRegionData.thetaExpr.freeze();
            //    }
            //}
            //auto allCCFirstRegionReserveScaleEnergySumExpr =
            //    graph.as<double>(HSum<double>(&graph, allCCFirstRegionReserveScaleEnergyHandles));
            //
            //// disable region ovelap/connection constraints among the same CC
            //for (auto & con : _constraints.elements<1>()){
            //    auto & ctopo = con.topo;
            //    auto & rd1 = _constraints.data(ctopo.lowers[0]).asRegion;
            //    auto & rd2 = _constraints.data(ctopo.lowers[1]).asRegion;
            //    if (rd1.connectedComponentId == rd2.connectedComponentId)
            //        con.data.constraintEnergyExpr.disable();
            //}
            //std::map<int, std::list<ComponentHandle>> ccid2RegionHandles;
            //int ccid = 0;
            //for (auto & cd : _constraints.elements<0>()){
            //    if (cd.exists && cd.data.type == ComponentData::Type::Region){
            //        ccid2RegionHandles[ccid++] = { cd.topo.hd };
            //    }
            //}




            // reconstruct faces
            OptimizeConstraintGraphUsingEnergyExpression(
                allConstraintEnergySumExpr[ConstraintData::Type::LineConnectivity].expression() +
                allConstraintEnergySumExpr[ConstraintData::Type::LineInterViewIncidence].expression() +
                allReserveScaleEnergySumExpr[ComponentData::Type::Line].expression() +
                //allCCFirstLineReserveScaleEnergySumExpr +
                allConstraintEnergySumExpr[ConstraintData::Type::RegionConnectivity].expression() * 10 +
                allConstraintEnergySumExpr[ConstraintData::Type::RegionLineConnectivity].expression() * 10 +
                allConstraintEnergySumExpr[ConstraintData::Type::RegionOverlap].expression() * 10
                
                , _constraints, 1e-2, 0.0, 5000,
                [this, &lineSpatialAvatars, &startTime, &constraints,
                //nLineCC, &ccid2LineHandles, &ccid2RegionHandles, 
                & allConstraintEnergySumExpr](
                int i, double energyVal, double energyValChange) -> bool {

                std::cout
                    << "region connectivity energy: "
                    << allConstraintEnergySumExpr[ConstraintData::Type::RegionConnectivity].expression().result()
                    << std::endl
                    << "region line connectivity energy: "
                    << allConstraintEnergySumExpr[ConstraintData::Type::RegionLineConnectivity].expression().result()
                    << std::endl
                    << "region overlap energy: "
                    << allConstraintEnergySumExpr[ConstraintData::Type::RegionOverlap].expression().result()
                    << std::endl
                    << "line connectivity energy: "
                    << allConstraintEnergySumExpr[ConstraintData::Type::LineConnectivity].expression().result()
                    << std::endl
                    << "line interview: "
                    << allConstraintEnergySumExpr[ConstraintData::Type::LineInterViewIncidence].expression().result()
                    << std::endl;

                if (i == 0){
                    // count functional eta/thetas
                    int functionalLineEtaNum = 0;
                    int functionalRegionThetaNum = 0;
                    for (auto & cd : constraints.elements<0>()){
                        if (cd.data.type == ComponentData::Type::Line){
                            if (cd.data.asLine.etaExpr.derivativeExpression().isValid())
                                functionalLineEtaNum++;
                        }
                        else if (cd.data.type == ComponentData::Type::Region){
                            if (cd.data.asRegion.thetaExpr.derivativeExpression().isValid())
                                functionalRegionThetaNum++;
                        }
                    }
                    std::cout << "functional line eta num: " << functionalLineEtaNum << std::endl;
                    std::cout << "functional region theta num: " << functionalRegionThetaNum << std::endl;
                    assert(functionalLineEtaNum == nLineCC);
                    //assert(functionalRegionThetaNum == nRegionCC);
                }


                auto curTime = std::chrono::high_resolution_clock::now();
                auto duration = curTime - startTime;
                auto hours = std::chrono::duration_cast<std::chrono::hours>(duration);
                //if (hours.count() < 7) // wait 7 hours
                //    return true;
                if (i % 10 != 0)
                    return true;

                // visualize
                IndexHashMap<RegionIndex, vis::Color> ri2Color;
                for (auto & vd : constraints.elements<0>()){
                    if (vd.data.type == ComponentData::Type::Region){
                        Vec3 theta = deriv::MakeCoreVec(vd.data.asRegion.thetaExpr.expression().result());
                        /*Vec3 theta = deriv::MakeCoreVec(
                            _constraints.data(ccid2RegionHandles[vd.data.asRegion.connectedComponentId].front())
                            .asRegion.thetaExpr.expression().result());*/
                        theta = normalize(theta);
                        Vec3 colorVec = {
                            abs(theta.dot(_globalData.vanishingPoints[0])),
                            abs(theta.dot(_globalData.vanishingPoints[1])),
                            abs(theta.dot(_globalData.vanishingPoints[2]))
                        };
                        colorVec *= 255.0;
                        vis::Color color(colorVec[0], colorVec[1], colorVec[2], 255.0);
                        RegionIndex ri;
                        ri.handle = vd.data.asRegion.regionHandle;
                        ri.viewHandle = vd.data.asRegion.viewHandle;
                        ri2Color[ri] = color;
                    }
                }

                for (auto & vd : _views.elements<0>()){
                    Image im(vd.data.image.size(), CV_8UC4);
                    RegionIndex ri;
                    ri.viewHandle = vd.topo.hd;
                    auto segmentedRegions = vd.data.regionNet->segmentedRegions();
                    for (int y = 0; y < segmentedRegions.rows; y++){
                        for (int x = 0; x < segmentedRegions.cols; x++){
                            ri.handle.id = segmentedRegions.at<int32_t>(y, x);
                            im.at<Vec<uint8_t, 4>>(y, x) = ri2Color[ri];
                        }
                    }
                    vis::Visualizer2D viz;
                    viz.setImage(im);
                    viz.params.winName = "region orientations " + std::to_string(vd.topo.hd.id);
                    viz = viz << vd.data.lineNet->lineSegments();
                    viz << vis::manip2d::Show();
                    //cv::imwrite("region_orientations_" + std::to_string(i) + "_" + 
                    //    std::to_string(vd.topo.hd.id) + ".png", viz.image());
                }

                std::cout << "now visualize spatial lines" << std::endl;

                // show lines
                std::vector<Classified<Line3>> rectifiedLines;
                for (auto & cd : constraints.elements<0>()){
                    if (cd.data.type == ComponentData::Type::Line){
                        auto line = lineSpatialAvatars[LineIndex{ cd.data.asLine.viewHandle, cd.data.asLine.lineHandle }];
                        auto ratio = ComputeDepthRatioOfPointOnSpatialLine(line.component.first,
                            line.component.second, globalData().vanishingPoints[line.claz]);
                        /*double eta = _constraints.data(ccid2LineHandles[cd.data.asLine.connectedComponentId].front())
                            .asLine.etaExpr.expression().result() * cd.data.asLine.etaRatio;*/
                        double eta = cd.data.asLine.etaExpr.expression().result();
                        line.component.first = normalize(line.component.first) * eta;
                        line.component.second = normalize(line.component.second) * eta * ratio;
                        rectifiedLines.push_back(line);
                    }
                }

                vis::Visualizer3D()
                    << vis::manip3d::SetDefaultColor(vis::ColorTag::Yellow)
                    << vis::manip3d::SetColorTableDescriptor(vis::ColorTableDescriptor::RGB)
                    << rectifiedLines
                    << vis::manip3d::AutoSetCamera
                    << vis::manip3d::Show(false);


                std::cout << "now visualize spatial regions" << std::endl;

                // show faces
                std::vector<std::vector<std::pair<Point3, Point2>>> regions;
                regions.reserve(constraints.internalElements<0>().size() * 2);
                for (auto & vd : constraints.elements<0>()){
                    if (vd.data.type == ReconstructionEngine::ComponentData::Type::Region){
                        auto theta = deriv::MakeCoreVec(vd.data.asRegion.thetaExpr.expression().result());
                        /*Vec3 theta = deriv::MakeCoreVec(
                            _constraints.data(ccid2RegionHandles[vd.data.asRegion.connectedComponentId].front())
                            .asRegion.thetaExpr.expression().result());*/
                        auto & cam = views().data(vd.data.asRegion.viewHandle).camera;
                        auto & regionData =
                            views().data(vd.data.asRegion.viewHandle).regionNet->regions().data(vd.data.asRegion.regionHandle);
                        auto & outCountours = regionData.contours.back();
                        std::vector<std::pair<Point3, Point2>> spatialRegion(outCountours.size());
                        for (int i = 0; i < outCountours.size(); i++){
                            Point3 dir = cam.spatialDirection(outCountours[i]);
                            dir /= norm(dir);
                            dir *= (theta.dot(theta) / theta.dot(dir));
                            spatialRegion[i].first = dir;
                            spatialRegion[i].second = Point2(0.0, 0.0);
                        }

                        regions.push_back(spatialRegion);
                        std::reverse(spatialRegion.begin(), spatialRegion.end());
                        regions.push_back(spatialRegion);
                    }
                }

                vis::Visualizer3D()
                    << vis::manip3d::SetDefaultColor(vis::ColorTag::Yellow)
                    << vis::manip3d::SetColorTableDescriptor(vis::ColorTableDescriptor::RGB)
                    << regions
                    << vis::manip3d::AutoSetCamera
                    << vis::manip3d::Show();

                return true;
            });          

        }


        void ReconstructionEngine::reconstructLinesAndFacesII(){

        }

    }
}