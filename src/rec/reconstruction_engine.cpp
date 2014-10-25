
extern "C" {
    #include <gpc.h>
}

#include <glpk.h>
#include <setjmp.h>

#include <Eigen/Core>
#include <Eigen/Dense>
#include <Eigen/Sparse>
#include <Eigen/StdVector>

#include <unsupported/Eigen/NonLinearOptimization>
#include <eiquadprog.hpp>

#include <dlib/matrix.h>
#include <dlib/optimization.h>

#include <GCoptimization.h>

#include "../core/feature.hpp"

#include "../vis/visualize2d.hpp"
#include "../vis/visualize3d.hpp"

#include "reconstruction_engine.hpp"

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
            samplingStepLengthOnRegionBoundaries(40.0),
            samplingStepLengthOnLines(10.0),
            intersectionDistanceThreshold(30), // 30
            incidenceDistanceAlongDirectionThreshold(50), // 50
            incidenceDistanceVerticalDirectionThreshold(4),
            interViewIncidenceAngleAlongDirectionThreshold(M_PI_4 / 4){
        }


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


        void ReconstructionEngine::computeFeatures(ViewHandle h) {
            auto & vd = _views.data(h);
            const Image & im = vd.image;

            // regions
            RegionsNet::Params regionsNetParams;
            regionsNetParams.samplingStepLengthOnBoundary = _params.samplingStepLengthOnRegionBoundaries;
            vd.regionNet = std::make_shared<RegionsNet>(vd.image, regionsNetParams);
            vd.regionNet->buildNetAndComputeGeometricFeatures();
            vd.regionNet->computeImageFeatures();

            // lines
            LinesNet::Params linesNetParams;
            linesNetParams.intersectionDistanceThreshold = _params.intersectionDistanceThreshold;
            linesNetParams.incidenceDistanceVerticalDirectionThreshold = _params.incidenceDistanceVerticalDirectionThreshold;
            linesNetParams.incidenceDistanceAlongDirectionThreshold = _params.incidenceDistanceAlongDirectionThreshold;
            LineSegmentExtractor::Params lsparams;
            lsparams.useLSD = true;
            linesNetParams.lineSegmentExtractor = LineSegmentExtractor(lsparams);
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
                    PixelLoc pixel = PixelLocFromGeoCoord(GeoCoord(p), longitudeDivideNum, latitudeDivideNum);
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
                    Vec3 vecs[] = { vec1, vec1rev, vec2, vec2rev };

                    double score = 0;
                    for (Vec3 & v : vecs){
                        PixelLoc pixel = PixelLocFromGeoCoord(GeoCoord(v), longitudeDivideNum, latitudeDivideNum);
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
                            PixelLoc pixel = PixelLocFromGeoCoord(GeoCoord(v), longitudeDivideNum, latitudeDivideNum);
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
                    Vec3 p3 = vIter->data.camera.spatialDirection(p.value());
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
                    spatialLineSegmentBegin, [&v](const Line2 & line) -> Classified<Line3> {
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

            template <class T, int N>
            inline Line<T, N> NormalizeLine(const Line<T, N> & l) {
                return Line<T, N>(normalize(l.first), normalize(l.second));
            }

            std::vector<int> FillInRectangleWithXs(int extendSize){
                std::vector<int> dx;
                dx.reserve(2 * extendSize + 1);
                for (int a = -extendSize; a <= extendSize; a++) {
                    for (int b = -extendSize; b <= extendSize; b++) {
                        dx.push_back(a);
                    }
                }
                return dx;
            }

            std::vector<int> FillInRectangleWithYs(int extendSize){
                std::vector<int> dy;
                dy.reserve(2 * extendSize + 1);
                for (int a = -extendSize; a <= extendSize; a++) {
                    for (int b = -extendSize; b <= extendSize; b++) {
                        dy.push_back(b);
                    }
                }
                return dy;
            }

        }


        void ReconstructionEngine::recognizeRegionLineRelations() {

            // compute spatial positions of each region
            IndexHashMap<RegionIndex, std::vector<Vec3>>
                regionSpatialContours;
            for (auto & view : _views.elements<0>()) {
                const auto & regions = *view.data.regionNet;
                for (auto & region : regions.regions().elements<0>()) {
                    RegionIndex ri = { view.topo.hd, region.topo.hd };
                    const ReconstructionEngine::ViewData & vd = view.data;
                    const RegionsNet::RegionData & rd = region.data;

                    //assert(!rd.contours.empty() && "Region contour not initialized yet?");
                    std::vector<Vec3> spatialContour;
                    if (!rd.dilatedContours.empty()){
                        for (auto & p : rd.dilatedContours.back()) {
                            auto direction = vd.camera.spatialDirection(p);
                            spatialContour.push_back(direction / norm(direction));
                        }
                    }
                    else{
                        std::cerr << "this region has no dilatedCountour!" << std::endl;
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
            IndexHashMap<std::pair<RegionIndex, RegionIndex>, double> &
                overlappedRegionIndexPairs = _globalData.overlappedRegionIndexPairs;
            overlappedRegionIndexPairs.clear();

            for (auto & rip : regionSpatialContours) {
                auto & ri = rip.first;

                auto & riCountours = regionData(ri).contours;
                if (riCountours.empty()){
                    std::cerr << "this region has no countour!" << std::endl;
                    continue;
                }

                auto & riContour2d = riCountours.front();
                auto & riCamera = _views.data(ri.viewHandle).camera;
                double riArea = regionData(ri).area;
                //double riArea = cv::contourArea(riContour2d);

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
                        //assert(overlapRatio <= 1.5 && "Invalid overlap ratio!");

                        if (overlapRatio > 0.2)
                            overlappedRegionIndexPairs[std::make_pair(relatedRi, ri)] = overlapRatio;
                    }

                    gpc_free_polygon(&relatedRiPoly);
                    gpc_free_polygon(&intersectedPoly);

                    return true;
                });

                gpc_free_polygon(&riPoly);
            }

            //// LINES ////
            // compute spatial normal directions for each line
            IndexHashMap<LineIndex, Classified<Line3>>
                lineSpatialAvatars;
            for (auto & vd : _views.elements<0>()) {
                auto & lines = vd.data.lineNet->lines();
                LineIndex li;
                li.viewHandle = vd.topo.hd;
                auto & cam = vd.data.camera;
                for (auto & ld : lines.elements<0>()) {
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
            for (auto & i : lineSpatialAvatars) {
                linesRTree.insert(i.first);
            }

            // recognize incidence constraints between lines of different views
            IndexHashMap<std::pair<LineIndex, LineIndex>, Vec3> & lineIncidenceRelations 
                = _globalData.lineIncidenceRelationsAcrossViews;
            lineIncidenceRelations.clear();

            for (auto & i : lineSpatialAvatars) {
                auto li = i.first;
                auto & lineData = i.second;
                auto & views = _views;
                linesRTree.search(lookupLineNormal(li), 
                    [this, &li, &lineSpatialAvatars, &views, &lineIncidenceRelations](const LineIndex & relatedLi) -> bool {
                    if (li.viewHandle == relatedLi.viewHandle)
                        return true;
                    if (relatedLi < li) // make sure one relation is stored only once, avoid storing both a-b and b-a
                        return true;

                    auto & line1 = lineSpatialAvatars[li];
                    auto & line2 = lineSpatialAvatars[relatedLi];
                    if (line1.claz != line2.claz) // only incidence relations are recognized here
                        return true;

                    auto normal1 = normalize(line1.component.first.cross(line1.component.second));
                    auto normal2 = normalize(line2.component.first.cross(line2.component.second));
                    auto & vd1 = views.data(li.viewHandle);
                    auto & vd2 = views.data(relatedLi.viewHandle);
                    if (std::min(std::abs(AngleBetweenDirections(normal1, normal2)), std::abs(AngleBetweenDirections(normal1, -normal2))) <
                        vd1.lineNet->params().incidenceDistanceVerticalDirectionThreshold / vd1.camera.focal() +
                        vd2.lineNet->params().incidenceDistanceVerticalDirectionThreshold / vd2.camera.focal()) {

                        auto nearest = DistanceBetweenTwoLines(NormalizeLine(line1.component), NormalizeLine(line2.component));
                        if (AngleBetweenDirections(nearest.second.first.position, nearest.second.second.position) >
                            _params.interViewIncidenceAngleAlongDirectionThreshold) // ignore too far-away relations
                            return true;

                        auto relationCenter = (nearest.second.first.position + nearest.second.second.position) / 2.0;
                        relationCenter /= norm(relationCenter);

                        lineIncidenceRelations[std::make_pair(li, relatedLi)] = relationCenter;
                    }
                    return true;
                });
            }

            // check whether all interview incidences are valid
            IF_DEBUG_USING_VISUALIZERS{
                double maxDist = 0;
                Line3 farthestLine1, farthestLine2;
                for (auto & lir : lineIncidenceRelations) {
                    auto & line1 = lineSpatialAvatars[lir.first.first];
                    auto & line2 = lineSpatialAvatars[lir.first.second];
                    if (line1.claz != line2.claz) {
                        std::cout << "invalid classes!" << std::endl;
                    }
                    auto l1 = NormalizeLine(line1.component);
                    auto l2 = NormalizeLine(line2.component);
                    auto dist = DistanceBetweenTwoLines(l1, l2).first;
                    if (dist > maxDist) {
                        farthestLine1 = l1;
                        farthestLine2 = l2;
                        maxDist = dist;
                    }
                }
                {
                    std::cout << "max dist of interview incidence pair: " << maxDist << std::endl;
                    std::cout << "line1: " << farthestLine1.first << ", " << farthestLine1.second << std::endl;
                    std::cout << "line2: " << farthestLine2.first << ", " << farthestLine2.second << std::endl;
                    auto d = DistanceBetweenTwoLines(farthestLine1, farthestLine2);
                    double angleDist = AngleBetweenDirections(d.second.first.position, d.second.second.position);
                    std::cout << "angle dist: " << angleDist << std::endl;
                }
            }


            // generate sampled points for line-region connections
            IndexHashMap<std::pair<RegionIndex, LineIndex>, std::vector<Vec3>> & regionLineIntersectionSampledPoints 
                = _globalData.regionLineIntersectionSampledPoints;
            regionLineIntersectionSampledPoints.clear();

            static const int extendSize = 7;
            //std::vector<int> dx, dy;
            //dx.reserve(2 * extendSize + 1);
            //dy.reserve(2 * extendSize + 1);
            //for (int a = -extendSize; a <= extendSize; a++) {
            //    for (int b = -extendSize; b <= extendSize; b++) {
            //        dx.push_back(a);
            //        dy.push_back(b);
            //    }
            //}

            static std::vector<int> dx = FillInRectangleWithXs(extendSize);
            static std::vector<int> dy = FillInRectangleWithYs(extendSize);

            for (auto & vd : _views.elements<0>()) {
                RegionIndex ri;
                ri.viewHandle = vd.topo.hd;

                LineIndex li;
                li.viewHandle = vd.topo.hd;

                const Image & segmentedRegions = vd.data.regionNet->segmentedRegions();
                auto & cam = vd.data.camera;

                for (auto & ld : vd.data.lineNet->lines().elements<0>()) {
                    li.handle = ld.topo.hd;

                    auto & line = ld.data.line.component;
                    auto lineDir = normalize(line.direction());
                    double sampleStep = _params.samplingStepLengthOnLines;
                    int sampledNum = static_cast<int>(std::floor(line.length() / sampleStep));

                    for (int i = 0; i < sampledNum; i++) {
                        auto sampledPoint = line.first + lineDir * i * sampleStep;

                        std::set<int32_t> rhids;
                        for (int k = 0; k < dx.size(); k++) {
                            int x = BoundBetween(static_cast<int>(std::round(sampledPoint[0] + dx[k])), 0, segmentedRegions.cols - 1);
                            int y = BoundBetween(static_cast<int>(std::round(sampledPoint[1] + dy[k])), 0, segmentedRegions.rows - 1);
                            PixelLoc p(x, y);
                            rhids.insert(segmentedRegions.at<int32_t>(p));
                        }

                        for (int32_t rhid : rhids) {
                            ri.handle = RegionsNet::RegionHandle(rhid);
                            regionLineIntersectionSampledPoints[std::make_pair(ri, li)]
                                .push_back(normalize(cam.spatialDirection(sampledPoint)));
                        }
                    }
                }
            }



            // compute connected components based on region-region overlaps
            // as merged region indices
            auto overlappedRegionIndicesGetter = [this](const RegionIndex & ri) {
                std::vector<RegionIndex> neighbors;
                for (auto & overlappedRegionPair : _globalData.overlappedRegionIndexPairs){
                    auto overlappingRatio = overlappedRegionPair.second;
                    if (overlappingRatio < 0.2)
                        continue;
                    if (overlappedRegionPair.first.first == ri)
                        neighbors.push_back(overlappedRegionPair.first.second);
                    if (overlappedRegionPair.first.second == ri)
                        neighbors.push_back(overlappedRegionPair.first.first);
                }
                return neighbors;
            };

            // collect all region indices
            std::vector<RegionIndex> regionIndices;
            IndexHashMap<RegionIndex, int> regionIndexToId;
            for (auto & vd : _views.elements<0>()){
                RegionIndex ri;
                ri.viewHandle = vd.topo.hd;
                auto & cam = vd.data.camera;
                for (auto & rd : vd.data.regionNet->regions().elements<0>()){
                    ri.handle = rd.topo.hd;
                    regionIndices.push_back(ri);
                    regionIndexToId[ri] = regionIndices.size() - 1;
                }
            }

            _globalData.regionConnectedComponentIds.clear();
            _globalData.regionConnectedComponentsNum = core::ConnectedComponents(regionIndices.begin(), regionIndices.end(),
                overlappedRegionIndicesGetter, [this](const RegionIndex & ri, int ccid) {
                _globalData.regionConnectedComponentIds[ri] = ccid;
            });

            std::cout << "region ccnum: " << _globalData.regionConnectedComponentsNum << std::endl;




            // compute connected components based on line-line constraints
            auto relatedLineIndicesGetter = [this, &lineIncidenceRelations](const LineIndex & li) {
                std::vector<LineIndex> related;
                // constraints in same view
                auto & lines = _views.data(li.viewHandle).lineNet->lines();
                auto & relationsInSameView = lines.topo(li.handle).uppers;
                for (auto & rh : relationsInSameView) {
                    auto anotherLineHandle = lines.topo(rh).lowers[0];
                    if (lines.data(rh).junctionWeight < 1e-5){
                        continue; // ignore zero weight relations
                    }
                    if (anotherLineHandle == li.handle)
                        anotherLineHandle = lines.topo(rh).lowers[1];
                    related.push_back(LineIndex{ li.viewHandle, anotherLineHandle });
                }
                // incidence constraints across views
                for (auto & interviewIncidence : lineIncidenceRelations) {
                    if (interviewIncidence.first.first == li)
                        related.push_back(interviewIncidence.first.second);
                    else if (interviewIncidence.first.second == li)
                        related.push_back(interviewIncidence.first.first);
                }
                return related;
            };

            // collect all lines
            std::vector<LineIndex> lineIndices;
            IndexHashMap<LineIndex, int> lineIndexToIds;
            for (auto & vd : _views.elements<0>()) {
                LineIndex li;
                li.viewHandle = vd.topo.hd;
                for (auto & ld : vd.data.lineNet->lines().elements<0>()) {
                    li.handle = ld.topo.hd;
                    lineIndices.push_back(li);
                    lineIndexToIds[li] = lineIndices.size() - 1;
                }
            }

            _globalData.lineConnectedComponentIds.clear();
            _globalData.lineConnectedComponentsNum = core::ConnectedComponents(lineIndices.begin(), lineIndices.end(),
                relatedLineIndicesGetter, [this](const LineIndex & li, int ccid) {
                _globalData.lineConnectedComponentIds[li] = ccid;
            });


            std::cout << "line ccnum: " << _globalData.lineConnectedComponentsNum << std::endl;


            IF_DEBUG_USING_VISUALIZERS{
                // visualize connections between regions and lines
                std::unordered_map<ViewHandle, vis::Visualizer2D, HandleHasher<AtLevel<0>>> vizs;
                for (auto & vd : _views.elements<0>()) {
                    //vis::Visualizer2D viz(vd.data.regionNet->image);
                    int height = vd.data.regionNet->image().rows;
                    int width = vd.data.regionNet->image().cols;

                    ImageWithType<Vec3b> coloredOutput(vd.data.regionNet->segmentedRegions().size());
                    vis::ColorTable colors = vis::CreateRandomColorTableWithSize(vd.data.regionNet->regions().internalElements<0>().size());
                    for (int y = 0; y < height; y++) {
                        for (int x = 0; x < width; x++) {
                            coloredOutput(cv::Point(x, y)) =
                                vis::ToVec3b(colors[vd.data.regionNet->segmentedRegions().at<int32_t>(cv::Point(x, y))]);
                        }
                    }
                    vizs[vd.topo.hd].setImage(vd.data.regionNet->image());
                    vizs[vd.topo.hd].params.alphaForNewImage = 0.5;
                    vizs[vd.topo.hd] << coloredOutput;
                    vizs[vd.topo.hd] << vis::manip2d::SetColorTable(vis::ColorTableDescriptor::RGB);
                }

                for (auto & lineIdRi : regionLineIntersectionSampledPoints) {
                    auto & ri = lineIdRi.first.first;
                    auto & li = lineIdRi.first.second;
                    auto & cline2 = _views.data(li.viewHandle).lineNet->lines().data(li.handle).line;
                    auto & cam = _views.data(ri.viewHandle).camera;
                    auto & viz = vizs[ri.viewHandle];

                    viz << vis::manip2d::SetColorTable(vis::ColorTableDescriptor::RGB) << vis::manip2d::SetThickness(3) << cline2;
                    viz << vis::manip2d::SetColor(vis::ColorTag::Black)
                        << vis::manip2d::SetThickness(1);
                    auto & regionCenter = _views.data(ri.viewHandle).regionNet->regions().data(ri.handle).center;
                    for (auto & d : lineIdRi.second) {
                        auto p = cam.screenProjection(d);
                        viz << Line2(regionCenter, p);
                    }
                }

                for (auto & viz : vizs) {
                    viz.second << vis::manip2d::Show();
                }
            }




        }


        void ReconstructionEngine::estimateSpatialLineDepths() {

            using namespace Eigen;
            SparseMatrix<double> A, W;
            VectorXd B;

            // try minimizing ||W(AX-B)||^2

            // collect all lines
            std::vector<LineIndex> lineIndices;
            IndexHashMap<LineIndex, int> lineIndexToIds;
            for (auto & vd : _views.elements<0>()) {
                LineIndex li;
                li.viewHandle = vd.topo.hd;
                for (auto & ld : vd.data.lineNet->lines().elements<0>()) {
                    li.handle = ld.topo.hd;
                    lineIndices.push_back(li);
                    lineIndexToIds[li] = lineIndices.size() - 1;
                }
            }

            // collect all constraints
            std::vector<LineRelationIndex> lineRelationIndices; // constraint indices in same views
            IndexHashMap<LineRelationIndex, int> lineRelationIndexToIds;
            // interview incidence constriants are stored in _globalData

            for (auto & vd : _views.elements<0>()) {
                LineRelationIndex lri;
                lri.viewHandle = vd.topo.hd;
                for (auto & ld : vd.data.lineNet->lines().elements<1>()) {
                    lri.handle = ld.topo.hd;
                    lineRelationIndices.push_back(lri);
                    lineRelationIndexToIds[lri] = lineRelationIndices.size() - 1;
                }
            }

            // pick the first line id in each connected component
            IndexHashSet<LineIndex> firstLineIndexInConnectedComponents;
            std::set<int> ccIdsRecorded;
            for (auto & lineIndexAndItsCCId : _globalData.lineConnectedComponentIds) {
                int ccid = lineIndexAndItsCCId.second;
                if (ccIdsRecorded.find(ccid) == ccIdsRecorded.end()) { // not recorded yet
                    firstLineIndexInConnectedComponents.insert(lineIndexAndItsCCId.first);
                    ccIdsRecorded.insert(ccid);
                }
            }
            static const double constantEtaForFirstLineInEachConnectedComponent = 10.0;
            std::cout << "anchor size: " << firstLineIndexInConnectedComponents.size() << std::endl;
            for (auto & ccId : ccIdsRecorded) {
                std::cout << "ccid: " << ccId << std::endl;
            }


            // setup matrices
            int n = lineIndices.size(); // var num
            int m = lineRelationIndices.size() + _globalData.lineIncidenceRelationsAcrossViews.size();  // cons num

            A.resize(m, n);
            W.resize(m, m);
            B.resize(m);

            // write equations
            int curEquationNum = 0;

            // write intersection/incidence constraint equations in same view
            for (const LineRelationIndex & lri : lineRelationIndices) {
                auto & lrd = lineRelationData(lri);
                auto & relationCenter = lrd.relationCenter;
                //auto & weightDistribution = _views.data(lri.viewHandle).lineNet->lineVotingDistribution();

                auto & topo = _views.data(lri.viewHandle).lineNet->lines().topo(lri.handle);
                auto & camera = _views.data(lri.viewHandle).camera;
                LineIndex li1 = { lri.viewHandle, topo.lowers[0] };
                LineIndex li2 = { lri.viewHandle, topo.lowers[1] };

                int lineId1 = lineIndexToIds[li1];
                int lineId2 = lineIndexToIds[li2];

                auto & line1 = lineData(li1).line;
                auto & line2 = lineData(li2).line;

                auto & vp1 = _globalData.vanishingPoints[line1.claz];
                auto & vp2 = _globalData.vanishingPoints[line2.claz];

                double ratio1 = ComputeDepthRatioOfPointOnSpatialLine(
                    camera.spatialDirection(line1.component.first),
                    camera.spatialDirection(relationCenter), vp1);
                double ratio2 = ComputeDepthRatioOfPointOnSpatialLine(
                    camera.spatialDirection(line2.component.first),
                    camera.spatialDirection(relationCenter), vp2);

                if (!core::Contains(firstLineIndexInConnectedComponents, li1) &&
                    !core::Contains(firstLineIndexInConnectedComponents, li2)) {
                    // eta1 * ratio1 - eta2 * ratio2 = 0
                    A.insert(curEquationNum, lineId1) = ratio1;
                    A.insert(curEquationNum, lineId2) = -ratio2;
                    B(curEquationNum) = 0;
                }
                else if (core::Contains(firstLineIndexInConnectedComponents, li1)) {
                    // const[eta1] * ratio1 - eta2 * ratio2 = 0 -> 
                    // eta2 * ratio2 = const[eta1] * ratio1
                    A.insert(curEquationNum, lineId2) = ratio2;
                    B(curEquationNum) = constantEtaForFirstLineInEachConnectedComponent * ratio1;
                }
                else if (core::Contains(firstLineIndexInConnectedComponents, li2)) {
                    // eta1 * ratio1 - const[eta2] * ratio2 = 0 -> 
                    // eta1 * ratio1 = const[eta2] * ratio2
                    A.insert(curEquationNum, lineId1) = ratio1;
                    B(curEquationNum) = constantEtaForFirstLineInEachConnectedComponent * ratio2;
                }

                // set junction weights
                W.insert(curEquationNum, curEquationNum) = lrd.junctionWeight;

                curEquationNum++;
            }

            // write inter-view incidence constraints
            for (auto & lineIncidenceAcrossView : _globalData.lineIncidenceRelationsAcrossViews) {
                auto & li1 = lineIncidenceAcrossView.first.first;
                auto & li2 = lineIncidenceAcrossView.first.second;
                auto & relationCenter = lineIncidenceAcrossView.second;

                auto & camera1 = _views.data(li1.viewHandle).camera;
                auto & camera2 = _views.data(li2.viewHandle).camera;

                int lineId1 = lineIndexToIds[li1];
                int lineId2 = lineIndexToIds[li2];

                auto & line1 = lineData(li1).line;
                auto & line2 = lineData(li2).line;

                auto & vp1 = _globalData.vanishingPoints[line1.claz];
                auto & vp2 = _globalData.vanishingPoints[line2.claz];

                double ratio1 = ComputeDepthRatioOfPointOnSpatialLine(
                    normalize(camera1.spatialDirection(line1.component.first)),
                    normalize(relationCenter), vp1);
                double ratio2 = ComputeDepthRatioOfPointOnSpatialLine(
                    normalize(camera2.spatialDirection(line2.component.first)),
                    normalize(relationCenter), vp2);

                if (ratio1 == 0.0 || ratio2 == 0.0) {
                    std::cout << "!!!!!!!ratio is zero!!!!!!!!" << std::endl;
                }

                if (!core::Contains(firstLineIndexInConnectedComponents, li1) &&
                    !core::Contains(firstLineIndexInConnectedComponents, li2)) {
                    // eta1 * ratio1 - eta2 * ratio2 = 0
                    A.insert(curEquationNum, lineId1) = ratio1;
                    A.insert(curEquationNum, lineId2) = -ratio2;
                    B(curEquationNum) = 0;
                } else if (core::Contains(firstLineIndexInConnectedComponents, li1)) {
                    // const[eta1] * ratio1 - eta2 * ratio2 = 0 -> 
                    // eta2 * ratio2 = const[eta1] * ratio1
                    A.insert(curEquationNum, lineId2) = ratio2;
                    B(curEquationNum) = constantEtaForFirstLineInEachConnectedComponent * ratio1;
                } else if (core::Contains(firstLineIndexInConnectedComponents, li2)) {
                    // eta1 * ratio1 - const[eta2] * ratio2 = 0 -> 
                    // eta1 * ratio1 = const[eta2] * ratio2
                    A.insert(curEquationNum, lineId1) = ratio1;
                    B(curEquationNum) = constantEtaForFirstLineInEachConnectedComponent * ratio2;
                }

                double junctionWeight = 5.0;
                W.insert(curEquationNum, curEquationNum) = junctionWeight;

                curEquationNum++;
            }

            // solve the equation system
            static const bool useWeights = true;

            VectorXd X;
            SparseQR<SparseMatrix<double>, COLAMDOrdering<int>> solver;
            static_assert(!(SparseMatrix<double>::IsRowMajor), "COLAMDOrdering only supports column major");
            SparseMatrix<double> WA = W * A;
            A.makeCompressed();
            WA.makeCompressed();
            solver.compute(useWeights ? WA : A);
            if (solver.info() != Success) {
                assert(0);
                std::cout << "computation error" << std::endl;
                return;
            }
            VectorXd WB = W * B;
            X = solver.solve(useWeights ? WB : B);
            if (solver.info() != Success) {
                assert(0);
                std::cout << "solving error" << std::endl;
                return;
            }

            // fill back all etas
            int k = 0;
            for (int i = 0; i < lineIndices.size(); i++) {
                auto & li = lineIndices[i];
                double eta = X(i);
                if (firstLineIndexInConnectedComponents.find(li) != firstLineIndexInConnectedComponents.end()) { // is first of a cc
                    eta = constantEtaForFirstLineInEachConnectedComponent;
                    std::cout << "is the " << (++k) << "-th anchor!" << std::endl;
                }
                auto & line2 = _views.data(li.viewHandle).lineNet->lines().data(li.handle).line;                
                auto & camera = _views.data(li.viewHandle).camera;
                Line3 line3 = { 
                    normalize(camera.spatialDirection(line2.component.first)),
                    normalize(camera.spatialDirection(line2.component.second))
                };

                std::cout << "eta: " << eta << " --- " << "ccid: " << _globalData.lineConnectedComponentIds[li] << std::endl;
                
                double resizeScale = eta / norm(line3.first);
                line3.first *= resizeScale;
                line3.second *= (resizeScale * 
                    ComputeDepthRatioOfPointOnSpatialLine(line3.first, line3.second, _globalData.vanishingPoints[line2.claz]));

                _globalData.reconstructedLines[li] = line3;
            }



            // visualize ccids
            // display reconstructed lines
            IF_DEBUG_USING_VISUALIZERS{
                vis::Visualizer3D viz;
                viz << vis::manip3d::SetBackgroundColor(vis::ColorTag::White)
                    << vis::manip3d::SetDefaultColorTable(vis::CreateRandomColorTableWithSize(_globalData.lineConnectedComponentsNum))
                    << vis::manip3d::SetDefaultLineWidth(2.0);
                for (auto & l : _globalData.reconstructedLines) {
                    viz << core::ClassifyAs(NormalizeLine(l.second), _globalData.lineConnectedComponentIds[l.first]);
                }
                viz << vis::manip3d::SetDefaultLineWidth(4.0);
                for (auto & c : _globalData.lineIncidenceRelationsAcrossViews) {
                    auto & line1 = _globalData.reconstructedLines[c.first.first];
                    auto & line2 = _globalData.reconstructedLines[c.first.second];
                    auto nearest = DistanceBetweenTwoLines(NormalizeLine(line1), NormalizeLine(line2));
                    viz << vis::manip3d::SetDefaultForegroundColor(vis::ColorTag::Black)
                        << Line3(nearest.second.first.position, nearest.second.second.position);
                }
                viz << vis::manip3d::SetWindowName("not-yet-reconstructed lines with ccids");
                viz << vis::manip3d::Show(false, true);
            }

            IF_DEBUG_USING_VISUALIZERS{
                vis::Visualizer3D viz;
                viz << vis::manip3d::SetBackgroundColor(vis::ColorTag::White)
                    << vis::manip3d::SetDefaultColorTable(vis::CreateRandomColorTableWithSize(_globalData.lineConnectedComponentsNum))
                    << vis::manip3d::SetDefaultLineWidth(4.0);
                for (auto & l : _globalData.reconstructedLines) {
                    viz << core::ClassifyAs(l.second, _globalData.lineConnectedComponentIds[l.first]);
                }
                viz << vis::manip3d::SetWindowName("reconstructed lines with ccids");
                viz << vis::manip3d::Show(false, true);
            }

            IF_DEBUG_USING_VISUALIZERS{ // show interview constraints
                vis::Visualizer3D viz;
                viz << vis::manip3d::SetBackgroundColor(vis::ColorTag::White)
                    << vis::manip3d::SetDefaultColorTable(vis::CreateRandomColorTableWithSize(_globalData.lineConnectedComponentsNum))
                    << vis::manip3d::SetDefaultLineWidth(2.0);
                for (auto & l : _globalData.reconstructedLines) {
                    viz << core::ClassifyAs(l.second, _globalData.lineConnectedComponentIds[l.first]);
                }
                viz << vis::manip3d::SetDefaultLineWidth(4.0);
                for (auto & c : _globalData.lineIncidenceRelationsAcrossViews) {
                    auto & line1 = _globalData.reconstructedLines[c.first.first];
                    auto & line2 = _globalData.reconstructedLines[c.first.second];
                    auto nearest = DistanceBetweenTwoLines(line1, line2);
                    viz << vis::manip3d::SetDefaultForegroundColor(vis::ColorTag::Black)
                        << Line3(nearest.second.first.position, nearest.second.second.position);
                }
                viz << vis::manip3d::SetWindowName("reconstructed lines with interview constraints");
                viz << vis::manip3d::Show(true, true);
            }



        }

        namespace {

            template <class FunctorT>
            struct DataCostFunctorWrapper : GCoptimization::DataCostFunctor{
                inline DataCostFunctorWrapper(FunctorT && f) : fun(std::forward<FunctorT>(f)) {}
                virtual GCoptimization::EnergyTermType compute(GCoptimization::SiteID s, GCoptimization::LabelID l) override {
                    return fun(s, l);
                }
                FunctorT fun;
            };
            template <class FunctorT>
            inline DataCostFunctorWrapper<FunctorT> * AllocDataCostFunctor(FunctorT && f) {
                return new DataCostFunctorWrapper<FunctorT>(std::forward<FunctorT>(f));
            }

            template <class FunctorT>
            struct SmoothCostFunctorWrapper : GCoptimization::SmoothCostFunctor {
                inline SmoothCostFunctorWrapper(FunctorT && f) : fun(std::forward<FunctorT>(f)){}
                virtual GCoptimization::EnergyTermType compute(
                    GCoptimization::SiteID s1, GCoptimization::SiteID s2,
                    GCoptimization::LabelID l1, GCoptimization::LabelID l2) override {
                    return fun(s1, s2, l1, l2);
                }
                FunctorT fun;
            };
            template <class FunctorT>
            inline SmoothCostFunctorWrapper<FunctorT> * AllocSmoothCostFunctor(FunctorT && f) {
                return new SmoothCostFunctorWrapper<FunctorT>(std::forward<FunctorT>(f));
            }

            inline Point2 ToPoint2(const PixelLoc & p) {
                return Point2(p.x, p.y);
            }

        }

        void ReconstructionEngine::initializeRegionOrientations() {

            // initialize region data
            std::vector<RegionIndex> regionIndices;
            std::map<RegionIndex, int> regionIndexToGraphSiteId;

            for (auto & vd : _views.elements<0>()){
                RegionIndex ri;
                ri.viewHandle = vd.topo.hd;
                auto & cam = vd.data.camera;
                for (auto & rd : vd.data.regionNet->regions().elements<0>()){
                    ri.handle = rd.topo.hd;
                    regionIndices.push_back(ri);
                    regionIndexToGraphSiteId[ri] = regionIndices.size() - 1;
                }
            }


            enum RegionOrientationType {
                VP0 = 0,
                VP1 = 1,
                VP2 = 2,
                //Planar = 3,
                //Other = 4,
                RegionOrientationTypeNum
            };

            GCoptimizationGeneralGraph graph(regionIndices.size(), RegionOrientationTypeNum);

            // set neighbors
            // region boundaries
            for (auto & vd : _views.elements<0>()){
                for (auto & bd : vd.data.regionNet->regions().elements<1>()){
                    RegionIndex ri1 = { vd.topo.hd, bd.topo.lowers[0] };
                    RegionIndex ri2 = { vd.topo.hd, bd.topo.lowers[1] };
                    auto siteId1 = regionIndexToGraphSiteId[ri1];
                    auto siteId2 = regionIndexToGraphSiteId[ri2];
                    graph.setNeighbors(siteId1, siteId2);
                }
            }
            // region overlaps
            for (auto & rp : _globalData.overlappedRegionIndexPairs){
                auto overlappingRatio = rp.second;
                if (overlappingRatio < 0.05)
                    continue;
                auto siteId1 = regionIndexToGraphSiteId[rp.first.first];
                auto siteId2 = regionIndexToGraphSiteId[rp.first.second];
                graph.setNeighbors(siteId1, siteId2);
            }

            
            // data costs for regions
            // store region costs corresponding to orientations
            IndexHashMap<RegionIndex, std::array<double, 3>> regionOrientationCosts;
            for (auto & ri : regionIndices){
                regionOrientationCosts[ri] = { { 0.0, 0.0, 0.0 } };
            }            
            for (auto & regionLine : _globalData.regionLineIntersectionSampledPoints){
                auto & ri = regionLine.first.first;
                auto & nearbyLi = regionLine.first.second;
                int nearbyLineClass = lineData(nearbyLi).line.claz;
                if (nearbyLineClass == -1)
                    continue;
                size_t sampledPointsNum = regionLine.second.size();
                regionOrientationCosts[ri][nearbyLineClass] += sampledPointsNum;
            }
            for (auto & c : regionOrientationCosts) {
                for (int i = 0; i < 3; i++) {
                    for (int j = i + 1; j < 3; j++) {
                        int k = 0 + 1 + 2 - i - j;
                        c.second[i] = 1.0 - Gaussian(c.second[i], 7.0);
                    }
                }
            }


            static const int scaleFactor = 1000;

            graph.setDataCostFunctor(AllocDataCostFunctor([this, &regionIndices, &regionIndexToGraphSiteId, &regionOrientationCosts](
                GCoptimization::SiteID s, GCoptimization::LabelID l){
                // nearby lines
                auto & ri = regionIndices[s];
                if (l < 3){
                    return regionOrientationCosts[ri][l] * scaleFactor;
                }
                else{
                    NOT_IMPLEMENTED_YET();
                }
            }));


            // smooth costs for region relations

            // store region folding costs
            IndexHashMap<std::pair<RegionIndex, RegionIndex>, std::array<double, 3>> 
                regionFoldingCosts;
            for (auto & vd : _views.elements<0>()){
                HPoint2 vps[] = {
                    vd.data.camera.screenProjectionInHPoint(_globalData.vanishingPoints[0]),
                    vd.data.camera.screenProjectionInHPoint(_globalData.vanishingPoints[1]),
                    vd.data.camera.screenProjectionInHPoint(_globalData.vanishingPoints[2])
                };
                for (int i = 0; i < 3; i++) {
                    auto vpp = vps[i].value();
                    std::cout << "vp[" << i << "] = (" << vpp[0] << "," << vpp[1] << ")" << std::endl;
                }
                auto & regions = vd.data.regionNet->regions();
                for (auto & bd : regions.elements<1>()){
                    if (bd.data.sampledPoints.empty())
                        continue;

                    RegionIndex ri1 = { vd.topo.hd, bd.topo.lowers[0] };
                    RegionIndex ri2 = { vd.topo.hd, bd.topo.lowers[1] };
                    Point2 sampledPointsCenter(0, 0);
                    int num = 0;
                    for (auto & ps : bd.data.sampledPoints) {
                        for (auto & p : ps) {
                            sampledPointsCenter = sampledPointsCenter + p;
                            num++;
                        }
                    }
                    sampledPointsCenter /= num;
                    
                    for (int i = 0; i < 3; i++) {
                        Vec2 midToVP = (vps[i] - HPoint2(sampledPointsCenter)).numerator;
                        Vec2 edgeDir = bd.data.fittedLine.direction;
                        double angle = AngleBetweenUndirectedVectors(midToVP, edgeDir);
                        double cost = 1.0 - Gaussian(angle, M_PI / 16.0) * BoundBetween(bd.data.straightness, 0.0, 1.0);
                        cost = cost < 0.1 ? 0 : cost;
                        regionFoldingCosts[std::make_pair(ri1, ri2)][i] = cost;
                        regionFoldingCosts[std::make_pair(ri2, ri1)][i] = cost;
                    }

                   //std::cout << "region folding cost: " 
                   //     << regionFoldingCosts[std::make_pair(ri1, ri2)][0] << ", "
                   //     << regionFoldingCosts[std::make_pair(ri1, ri2)][1] << ", "
                   //     << regionFoldingCosts[std::make_pair(ri1, ri2)][2] << std::endl;
                }
            }



            graph.setSmoothCostFunctor(AllocSmoothCostFunctor([this, &regionIndices, &regionIndexToGraphSiteId, &regionFoldingCosts](
                GCoptimization::SiteID s1, GCoptimization::SiteID s2,
                GCoptimization::LabelID l1, GCoptimization::LabelID l2){
                auto & ri1 = regionIndices[s1];
                auto & ri2 = regionIndices[s2];
                if (ri1.viewHandle == ri2.viewHandle) { // region boundary
                    if (l1 < 3 && l2 < 3){ // folding
                        if (l1 == l2) { // same orientation cost
                            return 0; // TODO 
                        }
                        int foldOrientation = 0 + 1 + 2 - l1 - l2;
                        assert(regionFoldingCosts.find(std::make_pair(ri1, ri2)) != regionFoldingCosts.end());
                        return (int)regionFoldingCosts[std::make_pair(ri1, ri2)][foldOrientation] * 1 * scaleFactor;
                    }
                    else{
                        NOT_IMPLEMENTED_YET();
                    }
                }
                else { // region overlap
                    return l1 == l2 ? 0 : 1 * scaleFactor;
                }
            }));


            IF_DEBUG_USING_VISUALIZERS{

                // visualize data costs and folding costs
                for (auto & vd : _views.elements<0>()){
                    ImageWithType<Vec<uint8_t, 3>> orientationImage = ImageWithType<Vec<uint8_t, 3>>::zeros(vd.data.image.size());
                    for (int y = 0; y < vd.data.image.rows; y++){
                        for (int x = 0; x < vd.data.image.cols; x++){
                            auto regionId = vd.data.regionNet->segmentedRegions().at<int32_t>(PixelLoc(x, y));
                            RegionIndex ri = { vd.topo.hd, RegionsNet::RegionHandle(regionId) };
                            auto costs = regionOrientationCosts[ri];
                            Vec<uint8_t, 3> color(255 - costs[2] * 255, 255 - costs[1] * 255, 255 - costs[0] * 255);
                            orientationImage(PixelLoc(x, y)) = color;
                        }
                    }

                    vis::Visualizer2D viz(orientationImage);

                    for (auto & bd : vd.data.regionNet->regions().elements<1>()){
                        RegionIndex ri1 = { vd.topo.hd, bd.topo.lowers[0] };
                        RegionIndex ri2 = { vd.topo.hd, bd.topo.lowers[1] };
                        auto costs = regionFoldingCosts[std::make_pair(ri1, ri2)];
                        vis::Color color(255 - costs[2] * 255, 255 - costs[1] * 255, 255 - costs[0] * 255);
                        auto & edges = bd.data.edges;
                        for (auto & e : edges){
                            for (int i = 0; i < e.size() - 1; i++){
                                viz.params.color = color;
                                viz.params.thickness = 2;
                                viz = viz << Line<int, 2>(e[i], e[i + 1]);
                            }
                        }
                    }

                    viz << vis::manip2d::Show();
                }

            }


            std::cout << "energy before graph-cut: " << graph.compute_energy() << std::endl;
            graph.expansion();
            graph.swap();
            std::cout << "energy after graph-cut: " << graph.compute_energy() << std::endl;

            for (int i = 0; i < regionIndices.size(); i++){
                _globalData.regionOrientations[regionIndices[i]] = graph.whatLabel(i);
            }

            IF_DEBUG_USING_VISUALIZERS{

                static const vis::Color colors[] = {
                    vis::ColorFromTag(vis::ColorTag::Red),
                    vis::ColorFromTag(vis::ColorTag::Green),
                    vis::ColorFromTag(vis::ColorTag::Blue),
                    vis::ColorFromTag(vis::ColorTag::Yellow),
                    vis::ColorFromTag(vis::ColorTag::White)
                };

                // visualize result region labels
                for (auto & vd : _views.elements<0>()) {
                    RegionIndex ri;
                    ri.viewHandle = vd.topo.hd;

                    auto & regions = *vd.data.regionNet;
                    int width = regions.segmentedRegions().cols;
                    int height = regions.segmentedRegions().rows;
                    ImageWithType<cv::Vec<uint8_t, 3>> coloredOutput = ImageWithType<cv::Vec<uint8_t, 3>>::zeros(height, width);
                    for (int y = 0; y < height; y++) {
                        for (int x = 0; x < width; x++) {
                            auto regionId = regions.segmentedRegions().at<int32_t>(cv::Point(x, y));
                            ri.handle.id = regionId;
                            auto & color = colors[_globalData.regionOrientations[ri]];
                            coloredOutput(cv::Point(x, y)) =
                                Vec<uint8_t, 3>((uint8_t)color[0], (uint8_t)color[1], (uint8_t)color[2]);
                        }
                    }

                    vis::Visualizer2D(coloredOutput) << vis::manip2d::Show();
                }




            }

        }



        void ReconstructionEngine::estimateRegionPlanes() {

            // collect all region indices
            std::vector<RegionIndex> regionIndices;
            IndexHashMap<RegionIndex, int> regionIndexToId;
            for (auto & vd : _views.elements<0>()){
                RegionIndex ri;
                ri.viewHandle = vd.topo.hd;
                auto & cam = vd.data.camera;
                for (auto & rd : vd.data.regionNet->regions().elements<0>()){
                    ri.handle = rd.topo.hd;
                    regionIndices.push_back(ri);
                    regionIndexToId[ri] = regionIndices.size() - 1;
                }
            }

            // initialize region planes
            _globalData.regionConnectedComponentPlanes = std::vector<Plane3>(_globalData.regionConnectedComponentsNum);
            for (auto & ri : regionIndices){
                auto & rd = regionData(ri);
                Vec3 centerDir = _views.data(ri.viewHandle).camera.spatialDirection(rd.center);
                int regionCCId = _globalData.regionConnectedComponentIds[ri];
                _globalData.regionConnectedComponentPlanes[regionCCId].anchor = normalize(centerDir) * 10;
                _globalData.regionConnectedComponentPlanes[regionCCId].normal = normalize(centerDir);
            }
            // initialize line depth factors
            _globalData.lineConnectedComponentDepthFactors = std::vector<double>(_globalData.lineConnectedComponentsNum, 1.0);


            auto displayAll = [&](int highlightedLineCCId, int highlightedRegionCCId,
                const core::MaxHeap<int> & lineCCIdsNotCheckedYet, 
                const core::MaxHeap<int> & regionCCIdsNotCheckedYet){
                
                IF_DEBUG_USING_VISUALIZERS{
                    std::vector<Line3> linesRepresentingSampledPoints;

                    // line-region connections
                    for (auto & pp : _globalData.regionLineIntersectionSampledPoints){
                        RegionIndex ri = pp.first.first;
                        LineIndex li = pp.first.second;
                        int regionCCId = _globalData.regionConnectedComponentIds[ri];
                        int lineCCId = _globalData.lineConnectedComponentIds[li];
                        if (regionCCIdsNotCheckedYet.contains(regionCCId) || lineCCIdsNotCheckedYet.contains(lineCCId))
                            continue;

                        auto line = _globalData.reconstructedLines[li];
                        double depthFactor = _globalData.lineConnectedComponentDepthFactors[lineCCId];
                        line.first *= depthFactor;
                        line.second *= depthFactor;

                        const std::vector<Vec3> & selectedSampledPoints = pp.second;

                        for (const Vec3 & sampleRay : selectedSampledPoints){
                            Point3 pointOnLine = DistanceBetweenTwoLines(InfiniteLine3(Point3(0, 0, 0), sampleRay), line.infinieLine()).second.second;
                            Point3 pointOnRegion = IntersectionOfLineAndPlane(InfiniteLine3(Point3(0, 0, 0), sampleRay),
                                _globalData.regionConnectedComponentPlanes[regionCCId]).position;
                            linesRepresentingSampledPoints.emplace_back(pointOnLine, pointOnRegion);
                        }
                    }

                    // paint regions
                    std::vector<vis::SpatialProjectedPolygon> spps, highlightedSpps;
                    spps.reserve(regionIndices.size());
                    static const int stepSize = 10;

                    for (auto & ri : regionIndices){
                        vis::SpatialProjectedPolygon spp;
                        int regionCCId = _globalData.regionConnectedComponentIds[ri];
                        if (regionCCIdsNotCheckedYet.contains(regionCCId)) // igore not reconstructed regions 
                            continue;

                        spp.plane = _globalData.regionConnectedComponentPlanes[regionCCId];
                        auto & rd = regionData(ri);
                        if (rd.contours.back().size() < 3)
                            continue;

                        spp.corners.reserve(rd.contours.back().size() / double(stepSize));
                        auto & cam = _views.data(ri.viewHandle).camera;

                        PixelLoc lastPixel;
                        for (int i = 0; i < rd.contours.back().size(); i++){
                            if (spp.corners.empty()){
                                spp.corners.push_back(cam.spatialDirection(ToPoint2(rd.contours.back()[i])));
                                lastPixel = rd.contours.back()[i];
                            }
                            else {
                                if (Distance(lastPixel, rd.contours.back()[i]) >= stepSize){
                                    spp.corners.push_back(cam.spatialDirection(ToPoint2(rd.contours.back()[i])));
                                    lastPixel = rd.contours.back()[i];
                                }
                            }
                        }

                        spp.projectionCenter = cam.eye();
                        if (spp.corners.size() > 3){
                            spps.push_back(spp);
                            if (_globalData.regionConnectedComponentIds[ri] == highlightedRegionCCId){
                                highlightedSpps.push_back(spp);
                            }
                        }
                    }

                    vis::Visualizer3D viz;
                    viz << vis::manip3d::SetBackgroundColor(vis::ColorTag::White)
                        << vis::manip3d::SetDefaultLineWidth(1.0)
                        << vis::manip3d::SetDefaultForegroundColor(vis::ColorTag::DimGray)
                        << linesRepresentingSampledPoints
                        << vis::manip3d::SetDefaultLineWidth(5.0);

                    viz << vis::manip3d::SetDefaultColorTable(vis::CreateRandomColorTableWithSize(_globalData.lineConnectedComponentsNum));

                    // paint lines
                    std::vector<Line3> highlightedLines;
                    for (auto & l : _globalData.reconstructedLines) {
                        int lineCCId = _globalData.lineConnectedComponentIds[l.first];
                        if (lineCCIdsNotCheckedYet.contains(lineCCId))
                            continue;
                        auto line = l.second;
                        double depthFactor = _globalData.lineConnectedComponentDepthFactors[lineCCId];
                        line.first *= depthFactor;
                        line.second *= depthFactor;
                        if (lineCCId == highlightedLineCCId)
                            highlightedLines.push_back(line);
                        viz << core::ClassifyAs(line, lineCCId);
                    }

                    viz << vis::manip3d::Begin(spps)
                        << vis::manip3d::SetTexture(_globalData.panorama)
                        << vis::manip3d::End
                        << vis::manip3d::SetDefaultLineWidth(6.0)
                        << vis::manip3d::SetDefaultForegroundColor(vis::ColorTag::Black)
                        << BoundingBoxOfContainer(highlightedSpps)
                        << BoundingBoxOfContainer(highlightedLines)
                        << vis::manip3d::SetWindowName("initial region planes and reconstructed lines")
                        << vis::manip3d::Show();
                }
            };

            // collect line indices in same ccs
            std::vector<IndexHashSet<LineIndex>> lineIndicesInCCs(_globalData.lineConnectedComponentsNum);
            for (auto & ccid : _globalData.lineConnectedComponentIds){
                lineIndicesInCCs[ccid.second].insert(ccid.first);
            }

            // collect region indices in same merge regions
            std::vector<IndexHashSet<RegionIndex>> regionIndicesInCCs(_globalData.regionConnectedComponentsNum);
            for (auto & ccid : _globalData.regionConnectedComponentIds){
                regionIndicesInCCs[ccid.second].insert(ccid.first);
            }

            
            struct RegionCCRecInfo { // information for reconstruction of regions
                double areaRatio;
                Rational samplePointsNumWithOtherRegionsAnchoredRatio;
                Rational samplePointsNumWithLinesAnchoredRatio;
                inline RegionCCRecInfo() 
                    : samplePointsNumWithOtherRegionsAnchoredRatio(0, 0),
                    samplePointsNumWithLinesAnchoredRatio(0, 0),
                    areaRatio(0) {
                }
                inline bool reconstructible() const {
                    return samplePointsNumWithLinesAnchoredRatio.denominator > 0 || 
                        samplePointsNumWithOtherRegionsAnchoredRatio.denominator > 0;
                }
                inline double priority() const { 
                    double samplePointsCompleteRatioWithOtherRegions =
                        (samplePointsNumWithOtherRegionsAnchoredRatio.denominator == 0 ? 0 :
                        samplePointsNumWithOtherRegionsAnchoredRatio.value());
                    double samplePointsCompleteRatioWithLines =
                        (samplePointsNumWithLinesAnchoredRatio.denominator == 0 ? 0 :
                        samplePointsNumWithLinesAnchoredRatio.value());
                    return samplePointsCompleteRatioWithOtherRegions * 0.7
                        + samplePointsCompleteRatioWithLines * 0.29
                        + areaRatio * 0.01;
                }
            };
            struct LineCCRecInfo { // information for reconstruction of line ccs
                double linesNumRatio;
                Rational samplePointsNumWithRegionsAnchoredRatio;
                inline LineCCRecInfo() : samplePointsNumWithRegionsAnchoredRatio(0, 0), linesNumRatio(0) {}
                inline bool reconstructible() const { return true; }
                inline double priority() const {
                    double samplePointsCompleteRatio = 
                        samplePointsNumWithRegionsAnchoredRatio.denominator == 0 ?
                        0 : samplePointsNumWithRegionsAnchoredRatio.value();
                    return samplePointsCompleteRatio * 0.9 + linesNumRatio * 0.1; // the advantage over regions
                }
            };



            /// count sample point for each region
            std::vector<RegionCCRecInfo> regionCCRecInfos(_globalData.regionConnectedComponentsNum); // [boundary sample points, line sample points]
            // collect region boundary sample point num
            for (auto & v : _views.elements<0>()){
                for (auto & b : v.data.regionNet->regions().elements<1>()){
                    int n = 0;
                    for (auto & pts : b.data.sampledPoints){
                        n += pts.size();
                    }
                    auto ri1 = RegionIndex{ v.topo.hd, b.topo.lowers[0] };
                    auto ri2 = RegionIndex{ v.topo.hd, b.topo.lowers[1] };                    
                    regionCCRecInfos[_globalData.regionConnectedComponentIds[ri1]]
                        .samplePointsNumWithOtherRegionsAnchoredRatio.denominator += n;
                    regionCCRecInfos[_globalData.regionConnectedComponentIds[ri2]]
                        .samplePointsNumWithOtherRegionsAnchoredRatio.denominator += n;
                }
                for (auto & r : v.data.regionNet->regions().elements<0>()){
                    auto ri = RegionIndex{ v.topo.hd, r.topo.hd };
                    regionCCRecInfos[_globalData.regionConnectedComponentIds[ri]].areaRatio += 
                        r.data.area / double(v.data.image.cols * v.data.image.rows);
                }
            }
            // collect region line sample point num
            for (auto & pp : _globalData.regionLineIntersectionSampledPoints){
                regionCCRecInfos[_globalData.regionConnectedComponentIds[pp.first.first]]
                    .samplePointsNumWithLinesAnchoredRatio.denominator += pp.second.size();
            }


            /// count sample point for each line CC
            std::vector<LineCCRecInfo> lineCCRecInfos(_globalData.lineConnectedComponentsNum);
            for (auto & pp : _globalData.regionLineIntersectionSampledPoints){
                auto & lineCCRecInfo = lineCCRecInfos[_globalData.lineConnectedComponentIds[pp.first.second]];
                lineCCRecInfo.samplePointsNumWithRegionsAnchoredRatio.denominator += pp.second.size();
            }
            for (int i = 0; i < _globalData.lineConnectedComponentsNum; i++){
                lineCCRecInfos[i].linesNumRatio = lineIndicesInCCs[i].size() / double(_globalData.reconstructedLines.size());
            }


            // maximum heap for checking            
            std::vector<int> lineCCIds(_globalData.lineConnectedComponentsNum);
            std::iota(lineCCIds.begin(), lineCCIds.end(), 0);
            core::MaxHeap<int> lineCCIdsForChecking(lineCCIds.begin(), lineCCIds.end(), 
                [&lineCCRecInfos](int ccId) -> double {
                return lineCCRecInfos[ccId].priority();
            });
            std::vector<int> regionCCIds(_globalData.regionConnectedComponentsNum);
            std::iota(regionCCIds.begin(), regionCCIds.end(), 0);
            core::MaxHeap<int> regionCCIdsForChecking(regionCCIds.begin(), regionCCIds.end(),
                [&regionCCRecInfos](int ccId) -> double{
                return regionCCRecInfos[ccId].priority();
            });


            // data for reserving consistency in reconstruction
            // TODO
            std::vector<double> pointAnchors[3];

            int iterCount = 0;

            //// begin iteration            
            while (!(regionCCIdsForChecking.empty() && lineCCIdsForChecking.empty())) {
                iterCount++;

                double regionTopScore = regionCCIdsForChecking.empty() ? 0.0 : regionCCIdsForChecking.topScore();
                double lineCCTopScore = lineCCIdsForChecking.empty() ? 0.0 : lineCCIdsForChecking.topScore();
                if (regionTopScore < lineCCTopScore){
                    // line ccid to check
                    int lineCCId = lineCCIdsForChecking.top();
                    lineCCIdsForChecking.pop();
                    if (!lineCCRecInfos[lineCCId].reconstructible()){
                        continue;
                    }

                    // compute depth of this line cc
                    std::cout << "processing line cc: " << lineCCId << " - " << lineCCTopScore << std::endl;
                    // collect sample points with related checked regions
                    std::vector<double> candidateDepthFactors;
                    for (auto & pp : _globalData.regionLineIntersectionSampledPoints){
                        RegionIndex ri = pp.first.first;
                        LineIndex li = pp.first.second;
                        int thisRegionCCId = _globalData.regionConnectedComponentIds[ri];
                        int thisLineCCId = _globalData.lineConnectedComponentIds[li];
                        if (thisLineCCId != lineCCId){
                            continue;
                        }
                        if (!core::Contains(regionCCIdsForChecking, thisRegionCCId)){ // checked already, use as anchors
                            auto & samplePoints = pp.second;
                            auto & line = _globalData.reconstructedLines[li];
                            for (auto & p : samplePoints){
                                double depthVarOnLine = norm(DistanceBetweenTwoLines(line.infinieLine(), InfiniteLine3(Point3(0, 0, 0), p))
                                    .second.second);
                                double depthValueOnRegion = norm(IntersectionOfLineAndPlane(InfiniteLine3(Point3(0, 0, 0), p),
                                    _globalData.regionConnectedComponentPlanes[thisRegionCCId]).position);
                                if (IsInfOrNaN(depthVarOnLine) || IsInfOrNaN(depthValueOnRegion))
                                    continue;
                                candidateDepthFactors.push_back((depthValueOnRegion / depthVarOnLine));
                            }
                        }
                    }
                    if (candidateDepthFactors.empty()){ // isolated or first
                        // TODO
                        _globalData.lineConnectedComponentDepthFactors[lineCCId] = 1;
                    }
                    else {
                        double bestVotes = -1;
                        double bestDepthFactor = -1;
                        for (const double & df : candidateDepthFactors){
                            double vote = 0;
                            for (const double & df2 : candidateDepthFactors){
                                vote += Gaussian(df - df2, 0.01);
                            }
                            if (!IsInfOrNaN(vote) && vote > bestVotes){
                                bestDepthFactor = df;
                                bestVotes = vote;
                            }
                        }
                        assert(bestDepthFactor > 0);
                        _globalData.lineConnectedComponentDepthFactors[lineCCId] = bestDepthFactor;
                    }


                    // update related region scores
                    for (auto & pp : _globalData.regionLineIntersectionSampledPoints){
                        LineIndex li = pp.first.second;
                        RegionIndex ri = pp.first.first;
                        int thisLineCCId = _globalData.lineConnectedComponentIds[li];
                        int thisRegionCCId = _globalData.regionConnectedComponentIds[ri];
                        if (!core::Contains(regionCCIdsForChecking, thisRegionCCId)){ // checked already
                            continue;
                        }
                        if (thisLineCCId == lineCCId){ // related 
                            regionCCRecInfos[thisRegionCCId].samplePointsNumWithLinesAnchoredRatio.numerator += pp.second.size();
                            regionCCIdsForChecking.increaseScoreTo(thisRegionCCId, regionCCRecInfos[thisRegionCCId].priority());
                        }
                    }


                    displayAll(lineCCId, -1, lineCCIdsForChecking, regionCCIdsForChecking);

                }
                else if(lineCCTopScore < regionTopScore){
                    // region index to check
                    int regionCCId = regionCCIdsForChecking.top();
                    regionCCIdsForChecking.pop();
                    if (!regionCCRecInfos[regionCCId].reconstructible()){
                        continue;
                    }

                    // compute equation of this region
                    std::cout << "processing region cc: " << regionCCId  << " - " << regionTopScore << std::endl;
                    IndexHashMap<RegionIndex, std::vector<Point3>> anchorsOnRelatedRegions;
                    IndexHashMap<LineIndex, std::vector<Point3>> anchorsOnRelatedLines;

                    // collect roots of candidate planes
                    std::vector<Point3> candidatePlaneRoots;
                    static const double planeRootDistThres = 0.01;
                    DefaultInfluenceBoxFunctor<Point3> influenceBoxFun(planeRootDistThres);
                    RTreeWrapper<Point3, decltype(influenceBoxFun)> candidatePlaneRootsRTree(influenceBoxFun); // rtree for detecting duplicates

                    // collect anchors on related regions
                    for (auto & ri : regionIndicesInCCs[regionCCId]){
                        auto & cam = _views.data(ri.viewHandle).camera;
                        auto & regions = _views.data(ri.viewHandle).regionNet->regions();
                        auto & bdis = regions.topo(ri.handle).uppers;
                        for (auto & bdi : bdis){
                            auto r1 = regions.topo(bdi).lowers.front();
                            auto r2 = regions.topo(bdi).lowers.back();
                            assert(r1 != r2);
                            RegionIndex relatedRi;
                            relatedRi.viewHandle = ri.viewHandle;
                            
                            if (r1 == ri.handle)
                                relatedRi.handle = r2;
                            else if (r2 == ri.handle)
                                relatedRi.handle = r1;
                            else{
                                assert(0 && "impossible!");
                            }
                            int relatedRegionCCId = _globalData.regionConnectedComponentIds[relatedRi];
                            if (relatedRegionCCId == regionCCId)
                                continue;
                            if (core::Contains(regionCCIdsForChecking, relatedRegionCCId)){
                                continue;
                            }

                            const Plane3 & relatedRegionCCPlane = _globalData.regionConnectedComponentPlanes[relatedRegionCCId];
                            if (regions.data(bdi).straightness < 0.8){ // collect the plane of neighbor regions with rough boundaries
                                auto root = relatedRegionCCPlane.root();
                                if (!HasValue(root, IsInfOrNaN<double>) &&
                                    !candidatePlaneRootsRTree.contains(root, [](const Point3 & a, const Point3 & b){return Distance(a, b) < planeRootDistThres; })){
                                    candidatePlaneRoots.push_back(root);
                                    candidatePlaneRootsRTree.insert(root);
                                }
                            }

                            for (auto & pts : regions.data(bdi).sampledPoints){
                                for (auto & p : pts){
                                    Vec3 pdir = cam.spatialDirection(p);
                                    auto pOnPlane = IntersectionOfLineAndPlane(InfiniteLine3(Point3(0, 0, 0), pdir), relatedRegionCCPlane).position;
                                    anchorsOnRelatedRegions[relatedRi].push_back(pOnPlane);
                                }
                            }
                        }
                    }

                    // collect anchors on related lines
                    for (auto & pp : _globalData.regionLineIntersectionSampledPoints){
                        LineIndex li = pp.first.second;
                        RegionIndex ri = pp.first.first;
                        int thisLineCCId = _globalData.lineConnectedComponentIds[li];
                        int thisRegionCCId = _globalData.regionConnectedComponentIds[ri];
                        if (core::Contains(lineCCIdsForChecking, thisLineCCId)){ // ignore unchecked lines
                            continue;
                        }
                        if (thisRegionCCId == regionCCId){ // related
                            auto & samplePoints = pp.second;
                            auto & line = _globalData.reconstructedLines[li];
                            for (auto & p : samplePoints){
                                auto pOnLine = DistanceBetweenTwoLines(line.infinieLine(), InfiniteLine3(Point3(0, 0, 0), p))
                                    .second.second;
                                anchorsOnRelatedLines[li].push_back(pOnLine * _globalData.lineConnectedComponentDepthFactors[thisLineCCId]);
                            }
                        }
                    }

                    int anchorsOnRelatedRegionsNum = 0;
                    for (auto & as : anchorsOnRelatedRegions)
                        anchorsOnRelatedRegionsNum += as.second.size();
                    int anchorsOnRelatedLinesNum = 0;
                    for (auto & as : anchorsOnRelatedLines)
                        anchorsOnRelatedLinesNum += as.second.size();
                    assert(anchorsOnRelatedRegionsNum == regionCCRecInfos[regionCCId].samplePointsNumWithOtherRegionsAnchoredRatio.numerator);
                    assert(anchorsOnRelatedLinesNum == regionCCRecInfos[regionCCId].samplePointsNumWithLinesAnchoredRatio.numerator);
                    std::cout << "anchors num for this region: " << anchorsOnRelatedRegionsNum << "  -  " << anchorsOnRelatedLinesNum << std::endl;
                    
                    assert(!anchorsOnRelatedRegions.empty() || !anchorsOnRelatedLines.empty());

                    // collect roots of planes formed by multiple related instances
                    std::vector<Point3> allAnchors;
                    allAnchors.reserve(anchorsOnRelatedRegions.size() * 2 + anchorsOnRelatedLines.size() * 2);
                    for (auto & as : anchorsOnRelatedRegions){
                        assert(!as.second.empty());
                        if (as.second.size() == 1)
                            allAnchors.push_back(as.second.front());
                        else{
                            allAnchors.push_back(as.second.front());
                            allAnchors.push_back(as.second.back());
                        }
                    }
                    for (auto & as : anchorsOnRelatedLines){
                        assert(!as.second.empty());
                        if (as.second.size() == 1)
                            allAnchors.push_back(as.second.front());
                        else{
                            allAnchors.push_back(as.second.front());
                            allAnchors.push_back(as.second.back());
                        }
                    }

                    for (int i = 0; i < allAnchors.size(); i++){
                        for (int j = i + 1; j < allAnchors.size(); j++){
                            for (int k = j + 1; k < allAnchors.size(); k++){
                                auto plane = Plane3From3Points(allAnchors[i], allAnchors[j], allAnchors[k]);
                                auto root = plane.root();
                                if (norm(root) < 3)
                                    continue;
                                //double vpFactors[] = { 
                                //    root.dot(normalize(_globalData.vanishingPoints[0])), 
                                //    root.dot(normalize(_globalData.vanishingPoints[1])), 
                                //    root.dot(normalize(_globalData.vanishingPoints[2])) 
                                //};
                                //double vpFactorsMax = std::max({ abs(vpFactors[0]), abs(vpFactors[1]), abs(vpFactors[2]) });
                                //for (auto & vpFactor : vpFactors){
                                //    if (abs(vpFactor) < 0.2 * vpFactorsMax)
                                //        vpFactor = 0.0;
                                //}
                                ////std::cout << "raw root: " << root;
                                //root = vpFactors[0] * normalize(_globalData.vanishingPoints[0])
                                //    + vpFactors[1] * normalize(_globalData.vanishingPoints[1])
                                //    + vpFactors[2] * normalize(_globalData.vanishingPoints[2]);
                                //std::cout << "refined root: " << root << std::endl;
                                if (!HasValue(root, IsInfOrNaN<double>) &&
                                    !candidatePlaneRootsRTree.contains(root, [](const Point3 & a, const Point3 & b){return Distance(a, b) < planeRootDistThres; })){
                                    candidatePlaneRoots.push_back(root);
                                    candidatePlaneRootsRTree.insert(root);
                                }
                            }
                        }
                    }

                    if (candidatePlaneRoots.empty()){
                        std::cerr << "can not determine the plane of this region" << std::endl;
                        continue;
                    }


                    // vote for planes
                    static const double distFromPointToPlaneThres = 0.5;
                    int bestPlaneId = -1;
                    double bestVote = -1;
                    for (int i = 0; i < candidatePlaneRoots.size(); i++){
                        auto & root = candidatePlaneRoots[i];
                        double vote = 0;
                        for (auto & as : anchorsOnRelatedRegions){
                            for (auto & a : as.second){
                                double distanceToPlane = abs((a - root).dot(normalize(root)));
                                if (distanceToPlane > distFromPointToPlaneThres)
                                    continue;
                                vote += Gaussian(distanceToPlane, distFromPointToPlaneThres);
                            }
                        }
                        for (auto & as : anchorsOnRelatedLines){
                            for (auto & a : as.second){
                                double distanceToPlane = abs((a - root).dot(normalize(root)));
                                if (distanceToPlane > distFromPointToPlaneThres)
                                    continue;
                                vote += Gaussian(distanceToPlane, distFromPointToPlaneThres);
                            }
                        }
                        if (vote > bestVote){
                            bestVote = vote;
                            bestPlaneId = i;
                        }
                    }
                    assert(bestPlaneId != -1);
                    // get best plane
                    _globalData.regionConnectedComponentPlanes[regionCCId] = 
                        Plane3(candidatePlaneRoots[bestPlaneId], candidatePlaneRoots[bestPlaneId]);


                    //displayAll(regionCCId, lineCCIdsForChecking, regionCCIdsForChecking);



                    // update scores of related line ccs
                    for (auto & pp : _globalData.regionLineIntersectionSampledPoints){
                        LineIndex li = pp.first.second;
                        RegionIndex ri = pp.first.first;
                        int thisLineCCId = _globalData.lineConnectedComponentIds[li];
                        int thisRegionCCId = _globalData.regionConnectedComponentIds[ri];
                        if (!core::Contains(lineCCIdsForChecking, thisLineCCId)){ // checked already
                            continue;
                        }
                        if (thisRegionCCId == regionCCId){ // related
                            lineCCRecInfos[thisLineCCId].samplePointsNumWithRegionsAnchoredRatio.numerator += pp.second.size();
                            lineCCIdsForChecking.increaseScoreTo(thisLineCCId, lineCCRecInfos[thisLineCCId].priority());
                        }
                    }
                    // update scores of related region ccs 
                    // (via boundray! since overlapping constraints are already used to construct CCs)
                    for (auto & v : _views.elements<0>()){
                        for (auto & b : v.data.regionNet->regions().elements<1>()){
                            auto ri1 = RegionIndex{ v.topo.hd, b.topo.lowers[0] };
                            auto ri2 = RegionIndex{ v.topo.hd, b.topo.lowers[1] };
                            int thisRegionCCId1 = _globalData.regionConnectedComponentIds[ri1];
                            int thisRegionCCId2 = _globalData.regionConnectedComponentIds[ri2];
                            if (thisRegionCCId1 == regionCCId && core::Contains(regionCCIdsForChecking, thisRegionCCId2)){ // related
                                for (auto & pts : b.data.sampledPoints){
                                    regionCCRecInfos[thisRegionCCId2].samplePointsNumWithOtherRegionsAnchoredRatio.numerator += pts.size();
                                }
                                regionCCIdsForChecking.increaseScoreTo(thisRegionCCId2, regionCCRecInfos[thisRegionCCId2].priority());
                            }
                            else if (thisRegionCCId2 == regionCCId && core::Contains(regionCCIdsForChecking, thisRegionCCId1)){ // related
                                for (auto & pts : b.data.sampledPoints){
                                    regionCCRecInfos[thisRegionCCId1].samplePointsNumWithOtherRegionsAnchoredRatio.numerator += pts.size();
                                }
                                regionCCIdsForChecking.increaseScoreTo(thisRegionCCId1, regionCCRecInfos[thisRegionCCId1].priority());
                            }
                        }
                    }

                    //if (iterCount % 50 == 0)
                    //    displayAll(-1, regionCCId, lineCCIdsForChecking, regionCCIdsForChecking);
                }
                

                
            }


            //displayAll(-1, lineCCIdsForChecking, regionCCIdsForChecking);
           


        }












    }
}