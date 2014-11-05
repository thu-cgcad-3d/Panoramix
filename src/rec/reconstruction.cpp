extern "C" {
    #include <gpc.h>
}

#include <random>
#include <thread>

#include <Eigen/Core>
#include <Eigen/Dense>
#include <Eigen/Sparse>
#include <Eigen/StdVector>

#include <unsupported/Eigen/NonLinearOptimization>

#include "../core/debug.hpp"
#include "../core/utilities.hpp"
#include "../core/containers.hpp"
#include "../core/algorithms.hpp"

#include "../vis/visualize2d.hpp"
#include "../vis/visualize3d.hpp"

#include "reconstruction.hpp"

namespace panoramix {
    namespace rec {

        View<PanoramicCamera> CreatePanoramicView(const Image & panorama) {
            return View<PanoramicCamera>{panorama, PanoramicCamera(panorama.cols / M_PI / 2.0)};
        }

        std::vector<View<PerspectiveCamera>> PerspectiveSampling(const View<PanoramicCamera> & panoView,
            const std::vector<PerspectiveCamera> & cameras) {

            std::vector<View<PerspectiveCamera>> views(cameras.size());
            for (int i = 0; i < cameras.size(); i++){
                views[i].camera = cameras[i];
                views[i].image = MakeCameraSampler(views[i].camera, panoView.camera)(panoView.image);
            }

            return views;
        }


        std::pair<RegionsNet, LinesNet> InitializeFeatureNets(const View<PerspectiveCamera> & view,
            double samplingStepLengthOnRegionBoundaries,
            double intersectionDistanceThreshold,
            double incidenceDistanceVerticalDirectionThreshold,
            double incidenceDistanceAlongDirectionThreshold) {

            // regions
            RegionsNet::Params regionsNetParams;
            regionsNetParams.samplingStepLengthOnBoundary = samplingStepLengthOnRegionBoundaries;
            RegionsNet regionsNet(view.image, regionsNetParams);
            regionsNet.buildNetAndComputeGeometricFeatures();
            regionsNet.computeImageFeatures();

            // lines
            LinesNet::Params linesNetParams;
            linesNetParams.intersectionDistanceThreshold = intersectionDistanceThreshold;
            linesNetParams.incidenceDistanceVerticalDirectionThreshold = incidenceDistanceVerticalDirectionThreshold;
            linesNetParams.incidenceDistanceAlongDirectionThreshold = incidenceDistanceAlongDirectionThreshold;
            LineSegmentExtractor::Params lsparams;
            lsparams.useLSD = true;
            linesNetParams.lineSegmentExtractor = LineSegmentExtractor(lsparams);
            LinesNet linesNet(view.image, linesNetParams);

            // move as pair
            return std::make_pair(std::move(regionsNet), std::move(linesNet));
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
                double angleThreshold = M_PI / 3, double sigma = 0.1) {

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

            // region/line data getter
            inline const RegionsNet::RegionData & GetData(const RegionIndex & i, const std::vector<RegionsNet> & nets){
                return nets[i.viewId].regions().data(i.handle);
            }
            inline const RegionsNet::BoundaryData & GetData(const RegionBoundaryIndex & i, const std::vector<RegionsNet> & nets){
                return nets[i.viewId].regions().data(i.handle);
            }
            inline const LinesNet::LineData & GetData(const LineIndex & i, const std::vector<LinesNet> & nets){
                return nets[i.viewId].lines().data(i.handle);
            }
            inline const LinesNet::LineRelationData & GetData(const LineRelationIndex & i, const std::vector<LinesNet> & nets){
                return nets[i.viewId].lines().data(i.handle);
            }
        }

        std::array<Vec3, 3> EstimateVanishingPointsAndClassifyLines(const std::vector<View<PerspectiveCamera>> & views,
            std::vector<LinesNet> & linesNets) {

            assert(views.size() == linesNets.size() && "num of views and linesNets mismatched!");

            // collect line intersections
            size_t lineIntersectionsNum = 0;
            for (const auto & linesNet : linesNets) // count intersecton num
                lineIntersectionsNum += linesNet.lineSegmentIntersections().size();
            std::vector<Vec3> intersections;
            intersections.reserve(lineIntersectionsNum);
            for (int i = 0; i < views.size(); i++){ // projection 2d intersections to global GeoCoord
                for(const auto & p : linesNets[i].lineSegmentIntersections()){
                    Vec3 p3 = views[i].camera.spatialDirection(p.value());
                    intersections.push_back(p3 / norm(p3)); // normalized
                }
            }

            // find vanishing points;
            auto vanishingPoints = FindVanishingPoints(intersections);

            // add spatial line segments from line segments of all views
            size_t spatialLineSegmentsNum = 0;
            for (auto & linesNet : linesNets)
                spatialLineSegmentsNum += linesNet.lineSegments().size();
            std::vector<Classified<Line3>> spatialLineSegments;
            spatialLineSegments.reserve(spatialLineSegmentsNum);
            for (int i = 0; i < views.size(); i++){
                for (const auto & line : linesNets[i].lineSegments()) {
                    auto & p1 = line.first;
                    auto & p2 = line.second;
                    auto pp1 = views[i].camera.spatialDirection(p1);
                    auto pp2 = views[i].camera.spatialDirection(p2);
                    Classified<Line3> cline3;
                    cline3.claz = -1;
                    cline3.component = Line3{ pp1, pp2 };
                    spatialLineSegments.push_back(cline3);
                }
            }

            // classify lines
            ClassifyLines(vanishingPoints, spatialLineSegments);

            // build lines net and compute features
            auto spatialLineSegmentBegin = spatialLineSegments.begin();
            for (int i = 0; i < views.size(); i++){
                std::array<HPoint2, 3> projectedVPs;
                for (int j = 0; j < 3; j++){
                    projectedVPs[j] = views[i].camera.screenProjectionInHPoint(vanishingPoints[j]);
                }
                std::vector<int> lineClasses(linesNets[i].lineSegments().size());
                for (auto & lineClass : lineClasses){
                    lineClass = spatialLineSegmentBegin->claz;
                    ++spatialLineSegmentBegin;
                }
                linesNets[i].buildNetAndComputeFeaturesUsingVanishingPoints(projectedVPs, lineClasses);
            }

            return vanishingPoints;

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

        void RecognizeRegionLineConstraints(const std::vector<View<PerspectiveCamera>> & views,
            const std::vector<RegionsNet> & regionsNets, const std::vector<LinesNet> & linesNets,
            ComponentIndexHashMap<std::pair<RegionIndex, RegionIndex>, double> & regionOverlappings,
            ComponentIndexHashMap<std::pair<RegionIndex, LineIndex>, std::vector<Vec3>> & regionLineConnections,
            ComponentIndexHashMap<std::pair<LineIndex, LineIndex>, Vec3> & interViewLineIncidences,
            double interViewIncidenceAngleAlongDirectionThreshold,
            double samplingStepLengthOnLines){

            assert(views.size() == regionsNets.size());
            assert(views.size() == linesNets.size());

            // compute spatial positions of each region
            ComponentIndexHashMap<RegionIndex, std::vector<Vec3>>
                regionSpatialContours;
            for (int i = 0; i < views.size(); i++) {
                const auto & regions = regionsNets[i];
                for (auto & region : regions.regions().elements<0>()) {
                    RegionIndex ri = { i, region.topo.hd };
                    const RegionsNet::RegionData & rd = region.data;
                    std::vector<Vec3> spatialContour;
                    if (!rd.dilatedContours.empty()){
                        for (auto & p : rd.dilatedContours.back()) {
                            auto direction = views[i].camera.spatialDirection(p);
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
            regionOverlappings.clear();

            for (auto & rip : regionSpatialContours) {
                auto & ri = rip.first;

                auto & riCountours = GetData(ri, regionsNets).contours;
                if (riCountours.empty()){
                    std::cerr << "this region has no countour!" << std::endl;
                    continue;
                }

                auto & riContour2d = riCountours.front();
                auto & riCamera = views[ri.viewId].camera;
                double riArea = GetData(ri, regionsNets).area;
                //double riArea = cv::contourArea(riContour2d);

                gpc_polygon riPoly;
                ConvertToGPCPolygon(riContour2d, riPoly);

                regionsRTree.search(lookupRegionBB(ri),
                    [&ri, &riContour2d, &riPoly, &riCamera, riArea, &regionOverlappings, &regionSpatialContours](
                    const RegionIndex & relatedRi) {

                    if (ri.viewId == relatedRi.viewId) {
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
                            regionOverlappings[std::make_pair(relatedRi, ri)] = overlapRatio;
                    }

                    gpc_free_polygon(&relatedRiPoly);
                    gpc_free_polygon(&intersectedPoly);

                    return true;
                });

                gpc_free_polygon(&riPoly);
            }

            //// LINES ////
            // compute spatial normal directions for each line
            ComponentIndexHashMap<LineIndex, Classified<Line3>>
                lineSpatialAvatars;
            for (int i = 0; i < views.size(); i++) {
                auto & lines = linesNets[i].lines();
                LineIndex li;
                li.viewId = i;
                auto & cam = views[i].camera;
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
            interViewLineIncidences.clear();

            for (auto & i : lineSpatialAvatars) {
                auto li = i.first;
                auto & lineData = i.second;
                linesRTree.search(lookupLineNormal(li),
                    [&interViewIncidenceAngleAlongDirectionThreshold, &li, &lineSpatialAvatars, &views, &linesNets, &interViewLineIncidences](const LineIndex & relatedLi) -> bool {
                    if (li.viewId == relatedLi.viewId)
                        return true;
                    if (relatedLi < li) // make sure one relation is stored only once, avoid storing both a-b and b-a
                        return true;

                    auto & line1 = lineSpatialAvatars[li];
                    auto & line2 = lineSpatialAvatars[relatedLi];
                    if (line1.claz != line2.claz) // only incidence relations are recognized here
                        return true;

                    auto normal1 = normalize(line1.component.first.cross(line1.component.second));
                    auto normal2 = normalize(line2.component.first.cross(line2.component.second));

                    if (std::min(std::abs(AngleBetweenDirections(normal1, normal2)), std::abs(AngleBetweenDirections(normal1, -normal2))) <
                        linesNets[li.viewId].params().incidenceDistanceVerticalDirectionThreshold / views[li.viewId].camera.focal() +
                        linesNets[relatedLi.viewId].params().incidenceDistanceVerticalDirectionThreshold / views[relatedLi.viewId].camera.focal()) {

                        auto nearest = DistanceBetweenTwoLines(NormalizeLine(line1.component), NormalizeLine(line2.component));
                        if (AngleBetweenDirections(nearest.second.first.position, nearest.second.second.position) >
                            interViewIncidenceAngleAlongDirectionThreshold) // ignore too far-away relations
                            return true;

                        auto relationCenter = (nearest.second.first.position + nearest.second.second.position) / 2.0;
                        relationCenter /= norm(relationCenter);

                        interViewLineIncidences[std::make_pair(li, relatedLi)] = relationCenter;
                    }
                    return true;
                });
            }

            // check whether all interview incidences are valid
            IF_DEBUG_USING_VISUALIZERS{
                double maxDist = 0;
                Line3 farthestLine1, farthestLine2;
                for (auto & lir : interViewLineIncidences) {
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
            regionLineConnections.clear();

            static const int OPT_ExtendSize = 2;

            static std::vector<int> dx = FillInRectangleWithXs(OPT_ExtendSize);
            static std::vector<int> dy = FillInRectangleWithYs(OPT_ExtendSize);

            for (int i = 0; i < views.size(); i++) {
                RegionIndex ri;
                ri.viewId = i;

                LineIndex li;
                li.viewId = i;

                const Image & segmentedRegions = regionsNets[i].segmentedRegions();
                auto & cam = views[i].camera;

                for (auto & ld : linesNets[i].lines().elements<0>()) {
                    li.handle = ld.topo.hd;

                    auto & line = ld.data.line.component;
                    auto lineDir = normalize(line.direction());
                    double sampleStep = samplingStepLengthOnLines;
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
                            regionLineConnections[std::make_pair(ri, li)]
                                .push_back(normalize(cam.spatialDirection(sampledPoint)));
                        }
                    }
                }
            }

        }



        namespace {

            void CollectIndices(const std::vector<View<PerspectiveCamera>> & views,
                const std::vector<RegionsNet> & regionsNets,
                std::vector<RegionIndex> & regionIndices,
                ComponentIndexHashMap<RegionIndex, int> & regionIndexToId){
                regionIndices.clear();
                regionIndexToId.clear();
                for (int i = 0; i < views.size(); i++){
                    RegionIndex ri;
                    ri.viewId = i;
                    for (auto & rd : regionsNets[i].regions().elements<0>()){
                        ri.handle = rd.topo.hd;
                        regionIndices.push_back(ri);
                        regionIndexToId[ri] = regionIndices.size() - 1;
                    }
                }
            }

            void CollectIndices(const std::vector<View<PerspectiveCamera>> & views,
                const std::vector<LinesNet> & linesNets,
                std::vector<LineIndex> & lineIndices,
                ComponentIndexHashMap<LineIndex, int> & lineIndexToIds){
                lineIndices.clear();
                lineIndexToIds.clear();
                for (int i = 0; i < views.size(); i++) {
                    LineIndex li;
                    li.viewId = i;
                    for (auto & ld : linesNets[i].lines().elements<0>()) {
                        li.handle = ld.topo.hd;
                        lineIndices.push_back(li);
                        lineIndexToIds[li] = lineIndices.size() - 1;
                    }
                }
            }

        }

        static const double MinimumJunctionWeght = 1e-5;

        void ComputeConnectedComponentsUsingRegionLineConstraints(const std::vector<View<PerspectiveCamera>> & views,
            const std::vector<RegionsNet> & regionsNets, const std::vector<LinesNet> & linesNets,
            const ComponentIndexHashMap<std::pair<RegionIndex, RegionIndex>, double> & regionOverlappings,
            const ComponentIndexHashMap<std::pair<RegionIndex, LineIndex>, std::vector<Vec3>> & regionLineConnections,
            const ComponentIndexHashMap<std::pair<LineIndex, LineIndex>, Vec3> & interViewLineIncidences,
            int & regionConnectedComponentsNum, ComponentIndexHashMap<RegionIndex, int> & regionConnectedComponentIds,
            int & lineConnectedComponentsNum, ComponentIndexHashMap<LineIndex, int> & lineConnectedComponentIds){

            assert(views.size() == regionsNets.size());
            assert(views.size() == linesNets.size());

            int n = views.size();

            // compute connected components based on region-region overlaps
            // as merged region indices
            auto overlappedRegionIndicesGetter = [&](const RegionIndex & ri) {
                std::vector<RegionIndex> neighbors;
                for (auto & overlappedRegionPair : regionOverlappings){
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
            ComponentIndexHashMap<RegionIndex, int> regionIndexToId;
            CollectIndices(views, regionsNets, regionIndices, regionIndexToId);

            regionConnectedComponentIds.clear();
            regionConnectedComponentsNum = core::ConnectedComponents(regionIndices.begin(), regionIndices.end(),
                overlappedRegionIndicesGetter, [&regionConnectedComponentIds](const RegionIndex & ri, int ccid) {
                regionConnectedComponentIds[ri] = ccid;
            });

            std::cout << "region ccnum: " << regionConnectedComponentsNum << std::endl;




            // compute connected components based on line-line constraints
            auto relatedLineIndicesGetter = [&](const LineIndex & li) {
                std::vector<LineIndex> related;
                // constraints in same view
                auto & lines = linesNets[li.viewId].lines();
                auto & relationsInSameView = lines.topo(li.handle).uppers;
                for (auto & rh : relationsInSameView) {
                    auto anotherLineHandle = lines.topo(rh).lowers[0];
                    if (lines.data(rh).junctionWeight < MinimumJunctionWeght){
                        continue; // ignore zero weight relations
                    }
                    if (anotherLineHandle == li.handle)
                        anotherLineHandle = lines.topo(rh).lowers[1];
                    related.push_back(LineIndex{ li.viewId, anotherLineHandle });
                }
                // incidence constraints across views
                for (auto & interviewIncidence : interViewLineIncidences) {
                    if (interviewIncidence.first.first == li)
                        related.push_back(interviewIncidence.first.second);
                    else if (interviewIncidence.first.second == li)
                        related.push_back(interviewIncidence.first.first);
                }
                return related;
            };

            // collect all lines
            std::vector<LineIndex> lineIndices;
            ComponentIndexHashMap<LineIndex, int> lineIndexToIds;
            CollectIndices(views, linesNets, lineIndices, lineIndexToIds);

            lineConnectedComponentIds.clear();
            lineConnectedComponentsNum = core::ConnectedComponents(lineIndices.begin(), lineIndices.end(),
                relatedLineIndicesGetter, [&lineConnectedComponentIds](const LineIndex & li, int ccid) {
                lineConnectedComponentIds[li] = ccid;
            });


            std::cout << "line ccnum: " << lineConnectedComponentsNum << std::endl;


            IF_DEBUG_USING_VISUALIZERS{
                // visualize connections between regions and lines
                std::unordered_map<int, vis::Visualizer2D, HandleHasher<AtLevel<0>>> vizs;
                for (int i = 0; i < n; i++) {
                    //vis::Visualizer2D viz(vd.data.regionNet->image);
                    int height = views[i].image.rows;
                    int width = views[i].image.cols;

                    ImageWithType<Vec3b> coloredOutput(regionsNets[i].segmentedRegions().size());
                    vis::ColorTable colors = vis::CreateRandomColorTableWithSize(regionsNets[i].regions().internalElements<0>().size());
                    for (int y = 0; y < height; y++) {
                        for (int x = 0; x < width; x++) {
                            coloredOutput(cv::Point(x, y)) =
                                vis::ToVec3b(colors[regionsNets[i].segmentedRegions().at<int32_t>(cv::Point(x, y))]);
                        }
                    }
                    vizs[i].setImage(views[i].image);
                    vizs[i].params.alphaForNewImage = 0.5;
                    vizs[i] << coloredOutput;
                    vizs[i] << vis::manip2d::SetColorTable(vis::ColorTableDescriptor::RGB);
                }

                for (auto & lineIdRi : regionLineConnections) {
                    auto & ri = lineIdRi.first.first;
                    auto & li = lineIdRi.first.second;
                    auto & cline2 = linesNets[li.viewId].lines().data(li.handle).line;
                    auto & cam = views[ri.viewId].camera;
                    auto & viz = vizs[ri.viewId];

                    viz << vis::manip2d::SetColorTable(vis::ColorTableDescriptor::RGB) << vis::manip2d::SetThickness(3) << cline2;
                    viz << vis::manip2d::SetColor(vis::ColorTag::Black)
                        << vis::manip2d::SetThickness(1);
                    auto & regionCenter = regionsNets[ri.viewId].regions().data(ri.handle).center;
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



        namespace {

            void EstimateSpatialLineDepthsOnce(const std::vector<View<PerspectiveCamera>> & views,
                const std::vector<LinesNet> & linesNets,
                const std::array<Vec3, 3> & vanishingPoints,
                const std::vector<LineIndex> & lineIndices,
                const std::vector<LineRelationIndex> & lineRelationIndices,
                const ComponentIndexHashMap<std::pair<LineIndex, LineIndex>, Vec3> & interViewLineIncidences,
                int lineConnectedComponentsNum, const ComponentIndexHashMap<LineIndex, int> & lineConnectedComponentIds,
                ComponentIndexHashMap<LineIndex, Line3> & reconstructedLines,
                double constantEtaForFirstLineInEachConnectedComponent, 
                bool useWeights){

                ComponentIndexHashMap<LineIndex, int> lineIndexToIds;
                for (int i = 0; i < lineIndices.size(); i++)
                    lineIndexToIds[lineIndices[i]] = i;

                using namespace Eigen;
                SparseMatrix<double> A, W;
                VectorXd B;

                // try minimizing ||W(AX-B)||^2

                // pick the first line id in each connected component
                ComponentIndexHashSet<LineIndex> firstLineIndexInConnectedComponents;
                std::set<int> ccIdsRecorded;
                for (auto & lineIndexAndItsCCId : lineConnectedComponentIds) {
                    int ccid = lineIndexAndItsCCId.second;
                    if (ccIdsRecorded.find(ccid) == ccIdsRecorded.end()) { // not recorded yet
                        firstLineIndexInConnectedComponents.insert(lineIndexAndItsCCId.first);
                        ccIdsRecorded.insert(ccid);
                    }
                }

                std::cout << "anchor size: " << firstLineIndexInConnectedComponents.size() << std::endl;
                for (auto & ccId : ccIdsRecorded) {
                    std::cout << "ccid: " << ccId << std::endl;
                }


                // setup matrices
                int n = lineIndices.size(); // var num
                int m = lineRelationIndices.size() + interViewLineIncidences.size();  // cons num

                A.resize(m, n);
                W.resize(m, m);
                B.resize(m);

                // write equations
                int curEquationNum = 0;

                // write intersection/incidence constraint equations in same view
                for (const LineRelationIndex & lri : lineRelationIndices) {
                    auto & lrd = GetData(lri, linesNets);
                    auto & relationCenter = lrd.relationCenter;
                    //auto & weightDistribution = _views.data(lri.viewHandle).lineNet->lineVotingDistribution();

                    auto & topo = linesNets[lri.viewId].lines().topo(lri.handle);
                    auto & camera = views[lri.viewId].camera;
                    LineIndex li1 = { lri.viewId, topo.lowers[0] };
                    LineIndex li2 = { lri.viewId, topo.lowers[1] };

                    int lineId1 = lineIndexToIds[li1];
                    int lineId2 = lineIndexToIds[li2];

                    auto & line1 = GetData(li1, linesNets).line;
                    auto & line2 = GetData(li2, linesNets).line;

                    auto & vp1 = vanishingPoints[line1.claz];
                    auto & vp2 = vanishingPoints[line2.claz];

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
                    W.insert(curEquationNum, curEquationNum) = lrd.junctionWeight < MinimumJunctionWeght ? 0.0 : lrd.junctionWeight;

                    curEquationNum++;
                }

                // write inter-view incidence constraints
                for (auto & lineIncidenceAcrossView : interViewLineIncidences) {
                    auto & li1 = lineIncidenceAcrossView.first.first;
                    auto & li2 = lineIncidenceAcrossView.first.second;
                    auto & relationCenter = lineIncidenceAcrossView.second;

                    auto & camera1 = views[li1.viewId].camera;
                    auto & camera2 = views[li2.viewId].camera;

                    int lineId1 = lineIndexToIds[li1];
                    int lineId2 = lineIndexToIds[li2];

                    auto & line1 = GetData(li1, linesNets).line;
                    auto & line2 = GetData(li2, linesNets).line;

                    auto & vp1 = vanishingPoints[line1.claz];
                    auto & vp2 = vanishingPoints[line2.claz];

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

                    double junctionWeight = 5.0;
                    W.insert(curEquationNum, curEquationNum) = junctionWeight;

                    curEquationNum++;
                }

                // solve the equation system
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
                    auto & line2 = linesNets[li.viewId].lines().data(li.handle).line;
                    auto & camera = views[li.viewId].camera;
                    Line3 line3 = {
                        normalize(camera.spatialDirection(line2.component.first)),
                        normalize(camera.spatialDirection(line2.component.second))
                    };

                    //std::cout << "eta: " << eta << " --- " << "ccid: " << lineConnectedComponentIds.at(li) << std::endl;

                    double resizeScale = eta / norm(line3.first);
                    line3.first *= resizeScale;
                    line3.second *= (resizeScale *
                        ComputeDepthRatioOfPointOnSpatialLine(line3.first, line3.second, vanishingPoints[line2.claz]));

                    reconstructedLines[li] = line3;
                }


            }

        }

        void EstimateSpatialLineDepths(const std::vector<View<PerspectiveCamera>> & views,
            const std::vector<LinesNet> & linesNets,
            const std::array<Vec3, 3> & vanishingPoints,
            const ComponentIndexHashMap<std::pair<LineIndex, LineIndex>, Vec3> & interViewLineIncidences,
            int lineConnectedComponentsNum, const ComponentIndexHashMap<LineIndex, int> & lineConnectedComponentIds,
            ComponentIndexHashMap<LineIndex, Line3> & reconstructedLines,
            double constantEtaForFirstLineInEachConnectedComponent,
            bool twiceEstimation){

            assert(views.size() == linesNets.size());

            // collect all lines
            std::vector<LineIndex> lineIndices;
            ComponentIndexHashMap<LineIndex, int> lineIndexToIds;
            CollectIndices(views, linesNets, lineIndices, lineIndexToIds);

            // collect all same view constraints
            std::vector<LineRelationIndex> lineRelationIndices; // constraint indices in same views
            for (int i = 0; i < views.size(); i++) {
                LineRelationIndex lri;
                lri.viewId = i;
                for (auto & ld : linesNets[i].lines().elements<1>()) {
                    lri.handle = ld.topo.hd;
                    lineRelationIndices.push_back(lri);
                }
            }

            // reconstruct
            ComponentIndexHashMap<LineIndex, Line3> reconstructedLinesOriginal;
            EstimateSpatialLineDepthsOnce(views, linesNets, vanishingPoints, lineIndices, lineRelationIndices, 
                interViewLineIncidences, lineConnectedComponentsNum, lineConnectedComponentIds, 
                reconstructedLinesOriginal, constantEtaForFirstLineInEachConnectedComponent, true);

            if (!twiceEstimation){
                reconstructedLines = std::move(reconstructedLinesOriginal);
                return;
            }

            // store all line constraints homogeneously
            struct ConstraintBetweenLines{
                enum { InnerView, InterView } type;
                LineRelationIndex lineRelationIndex;
                std::pair<LineIndex, LineIndex> linePairIndex;
                double distance;
            };
            std::vector<ConstraintBetweenLines> homogeneousConstaints;
            homogeneousConstaints.reserve(lineRelationIndices.size() + interViewLineIncidences.size());
            for (auto & lri : lineRelationIndices){
                int viewId = lri.viewId;
                // ignore too light constraints, same in ComputeConnectedComponentsUsingRegionLineConstraints
                if (GetData(lri, linesNets).junctionWeight < MinimumJunctionWeght) 
                    continue;
                auto lineHandles = linesNets.at(viewId).lines().topo(lri.handle).lowers;
                auto & line1 = reconstructedLinesOriginal[LineIndex{ viewId, lineHandles[0] }];
                auto & line2 = reconstructedLinesOriginal[LineIndex{ viewId, lineHandles[1] }];
                auto nearestPoints = DistanceBetweenTwoLines(line1.infinieLine(), line2.infinieLine()).second;
                auto c = (nearestPoints.first + nearestPoints.second) / 2.0;
                double distance = abs((nearestPoints.first - nearestPoints.second).dot(normalize(c))) 
                    / constantEtaForFirstLineInEachConnectedComponent;
                homogeneousConstaints.push_back(ConstraintBetweenLines{ 
                    ConstraintBetweenLines::InnerView, lri, std::pair<LineIndex, LineIndex>(), distance 
                });
            }
            for (auto & ivl : interViewLineIncidences){
                auto & line1 = reconstructedLinesOriginal[ivl.first.first];
                auto & line2 = reconstructedLinesOriginal[ivl.first.second];
                auto nearestPoints = DistanceBetweenTwoLines(line1.infinieLine(), line2.infinieLine()).second;
                auto c = (nearestPoints.first + nearestPoints.second) / 2.0;
                double distance = abs((nearestPoints.first - nearestPoints.second).dot(normalize(c)))
                    / constantEtaForFirstLineInEachConnectedComponent;
                homogeneousConstaints.push_back(ConstraintBetweenLines{ 
                    ConstraintBetweenLines::InterView, LineRelationIndex(), ivl.first, distance 
                });
            }

            std::cout << "original line constraints num = " << homogeneousConstaints.size() << std::endl;
            std::vector<int> constraintIds(homogeneousConstaints.size());
            std::iota(constraintIds.begin(), constraintIds.end(), 0);

            // minimum spanning tree
            auto edgeVertsGetter = [&homogeneousConstaints, &linesNets](int cid) -> std::pair<LineIndex, LineIndex> {
                std::pair<LineIndex, LineIndex> verts;
                auto & c = homogeneousConstaints[cid];
                if (c.type == ConstraintBetweenLines::InnerView){
                    int viewId = c.lineRelationIndex.viewId;
                    auto lineHandles = linesNets.at(viewId).lines()
                        .topo(c.lineRelationIndex.handle).lowers;
                    verts = std::make_pair(LineIndex{ viewId, lineHandles[0] }, LineIndex{ viewId, lineHandles[1] });
                }
                else if (c.type == ConstraintBetweenLines::InterView){
                    verts = c.linePairIndex;
                }
                return verts;
            };

            std::vector<int> reservedHomogeneousConstaintsIds;
            reservedHomogeneousConstaintsIds.reserve(homogeneousConstaints.size() / 2);
            core::MinimumSpanningTree(lineIndices.begin(), lineIndices.end(), 
                constraintIds.begin(), constraintIds.end(),
                std::back_inserter(reservedHomogeneousConstaintsIds), edgeVertsGetter,
                [&homogeneousConstaints](int cid1, int cid2)->bool {
                return homogeneousConstaints[cid1].distance < homogeneousConstaints[cid2].distance;
            });

            std::cout << "line constraints num after MST = " << reservedHomogeneousConstaintsIds.size() << std::endl;


            // build trimmed line relation indices and inter-view-incidences
            std::vector<LineRelationIndex> trimmedLineRelationIndices;
            trimmedLineRelationIndices.reserve(reservedHomogeneousConstaintsIds.size() / 2);
            ComponentIndexHashMap<std::pair<LineIndex, LineIndex>, Vec3> trimmedInterViewLineIncidences;
            for (int i : reservedHomogeneousConstaintsIds){
                auto & c = homogeneousConstaints[i];
                if (c.type == ConstraintBetweenLines::InnerView){
                    trimmedLineRelationIndices.push_back(c.lineRelationIndex);
                }
                else if (c.type == ConstraintBetweenLines::InterView){
                    trimmedInterViewLineIncidences.emplace(c.linePairIndex, interViewLineIncidences.at(c.linePairIndex));
                }
            }

            // reconstruct again
            EstimateSpatialLineDepthsOnce(views, linesNets, vanishingPoints, lineIndices, trimmedLineRelationIndices,
                trimmedInterViewLineIncidences, lineConnectedComponentsNum, lineConnectedComponentIds,
                reconstructedLines, constantEtaForFirstLineInEachConnectedComponent, false);


            // visualize ccids
            // display reconstructed lines
            IF_DEBUG_USING_VISUALIZERS{
                vis::Visualizer3D viz;
                viz << vis::manip3d::SetBackgroundColor(vis::ColorTag::White)
                    << vis::manip3d::SetDefaultColorTable(vis::CreateRandomColorTableWithSize(lineConnectedComponentsNum))
                    << vis::manip3d::SetDefaultLineWidth(2.0);
                for (auto & l : reconstructedLines) {
                    viz << core::ClassifyAs(NormalizeLine(l.second), lineConnectedComponentIds.at(l.first));
                }
                viz << vis::manip3d::SetDefaultLineWidth(4.0);
                for (auto & c : interViewLineIncidences) {
                    auto & line1 = reconstructedLines[c.first.first];
                    auto & line2 = reconstructedLines[c.first.second];
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
                    << vis::manip3d::SetDefaultColorTable(vis::CreateRandomColorTableWithSize(lineConnectedComponentsNum))
                    << vis::manip3d::SetDefaultLineWidth(4.0);
                for (auto & l : reconstructedLinesOriginal) {
                    viz << core::ClassifyAs(l.second, lineConnectedComponentIds.at(l.first));
                }
                viz << vis::manip3d::SetWindowName("reconstructed lines with ccids, 1st time");
                viz << vis::manip3d::Show(false, true);
            }
            IF_DEBUG_USING_VISUALIZERS{
                vis::Visualizer3D viz;
                viz << vis::manip3d::SetBackgroundColor(vis::ColorTag::White)
                    << vis::manip3d::SetDefaultColorTable(vis::CreateRandomColorTableWithSize(lineConnectedComponentsNum))
                    << vis::manip3d::SetDefaultLineWidth(4.0);
                for (auto & l : reconstructedLines) {
                    viz << core::ClassifyAs(l.second, lineConnectedComponentIds.at(l.first));
                }
                viz << vis::manip3d::SetWindowName("reconstructed lines with ccids, 2nd time");
                viz << vis::manip3d::Show(false, true);
            }

            IF_DEBUG_USING_VISUALIZERS{ // show interview constraints
                vis::Visualizer3D viz;
                viz << vis::manip3d::SetBackgroundColor(vis::ColorTag::White)
                    << vis::manip3d::SetDefaultColorTable(vis::CreateRandomColorTableWithSize(lineConnectedComponentsNum))
                    << vis::manip3d::SetDefaultLineWidth(2.0);
                for (auto & l : reconstructedLines) {
                    viz << core::ClassifyAs(l.second, lineConnectedComponentIds.at(l.first));
                }
                viz << vis::manip3d::SetDefaultLineWidth(4.0);
                for (auto & c : interViewLineIncidences) {
                    auto & line1 = reconstructedLines[c.first.first];
                    auto & line2 = reconstructedLines[c.first.second];
                    auto nearest = DistanceBetweenTwoLines(line1, line2);
                    viz << vis::manip3d::SetDefaultForegroundColor(vis::ColorTag::Black)
                        << Line3(nearest.second.first.position, nearest.second.second.position);
                }
                viz << vis::manip3d::SetWindowName("reconstructed lines with interview constraints");
                viz << vis::manip3d::Show(true, true);
            }

        }



        // display options
        static const bool OPT_DisplayMessages = false;
        static const bool OPT_DisplayOnEachTrial = false;
        static const bool OPT_DisplayOnEachLineCCRegonstruction = false;
        static const bool OPT_DisplayOnEachRegionRegioncstruction = false;
        static const bool OPT_DisplayOnEachIteration = false;
        static const int OPT_DisplayOnEachIterationInterval = 500;
        static const bool OPT_DisplayAtLast = true;


        // algorithm options
        static const bool OPT_OnlyConsiderManhattanPlanes = true;
        static const bool OPT_IgnoreTooSkewedPlanes = true;
        static const bool OPT_IgnoreTooFarAwayPlanes = true;
        static const int OPT_RandomRepeatNum = 100;
        static const int OPT_MaxSolutionNumForEachCC = 1;
        static const int OPT_NumThreads = 3;


        namespace {

            inline Point2 ToPoint2(const PixelLoc & p) {
                return Point2(p.x, p.y);
            }

            double ComputeVisualAreaOfDirections(const Plane3 & tplane, const Vec3 & x, const Vec3 & y, 
                const std::vector<Vec3> & dirs, bool convexify){
                if (dirs.size() <= 2)
                    return 0.0;
                std::vector<Point2f> pointsOnPlane(dirs.size());
                static const Point3 zeroPoint(0, 0, 0);
                for (int i = 0; i < dirs.size(); i++){
                    auto pOnPlane = IntersectionOfLineAndPlane(InfiniteLine3(zeroPoint, dirs[i]), tplane).position;
                    auto pOnPlaneOffsetted = pOnPlane - tplane.anchor;
                    pointsOnPlane[i] = Point2f(pOnPlane.dot(x), pOnPlane.dot(y));
                }
                // compute convex hull and contour area
                if (convexify){
                    cv::convexHull(pointsOnPlane, pointsOnPlane, false, true);
                }
                return cv::contourArea(pointsOnPlane);
            }

            // reconstruction context
            struct RecContext {
                const std::vector<View<PerspectiveCamera>> & views;
                const std::vector<RegionsNet> & regionsNets;
                const std::vector<LinesNet> & linesNets;
                const std::array<Vec3, 3> & vanishingPoints;
                const ComponentIndexHashMap<std::pair<RegionIndex, RegionIndex>, double> & regionOverlappings;
                const ComponentIndexHashMap<std::pair<RegionIndex, LineIndex>, std::vector<Vec3>> & regionLineConnections;
                const ComponentIndexHashMap<std::pair<LineIndex, LineIndex>, Vec3> & interViewLineIncidences;
                int regionConnectedComponentsNum; 
                const ComponentIndexHashMap<RegionIndex, int> & regionConnectedComponentIds;
                int lineConnectedComponentsNum;
                const ComponentIndexHashMap<LineIndex, int> & lineConnectedComponentIds;
                const ComponentIndexHashMap<LineIndex, Line3> & reconstructedLines;
                const ComponentIndexHashMap<RegionIndex, Plane3> & reconstructedPlanes;
                const Image & globalTexture;
                const Box3 & initialBoundingBox;
            };


            // choice
            struct Choice {
                bool isRegionCC;
                inline bool isLineCC() const { return !isRegionCC; }
                int ccId; // regionCCId or lineCCId
                int choiceId;
            };





            // region cc reconstruction information
            struct RegionCCRecInfo {
                int regionCCId;
                ComponentIndexHashSet<RegionIndex> regionIndices;
                Plane3 tangentialPlane;
                Vec3 xOnTangentialPlane, yOnTangentialPlane;
                double regionVisualArea;
                double regionConvexContourVisualArea;
                Rational anchoredConnectionsWithOtherRegions;
                Rational anchoredConnectionsWithLines;
                std::vector<Point3> anchoredConnections;

                struct PlaneConfidenceData {
                    Plane3 plane;
                    std::vector<int> inlierAnchors;
                    double regionInlierAnchorsConvexContourVisualArea;
                    double regionInlierAnchorsDistanceVotesSum;
                };
                using PlaneConfidenceMap = VecMap<double, 3, PlaneConfidenceData>;
                PlaneConfidenceMap candidatePlanesByRoot;
            };

            // initialize region cc reconstruction information
            RegionCCRecInfo CreateRegionCCRecInfo(int regionCCId, const RecContext & context){
                assert(context.views.size() == context.regionsNets.size());

                RegionCCRecInfo rci;
                rci.regionCCId = regionCCId;
                rci.candidatePlanesByRoot = RegionCCRecInfo::PlaneConfidenceMap(0.05);

                // collect region indices
                for (auto & rcc : context.regionConnectedComponentIds){
                    if (rcc.second == regionCCId)
                        rci.regionIndices.insert(rcc.first);
                }

                // locate tangential coordinates
                std::vector<Vec3> outerContourDirections;
                Vec3 regionsCenterDirection(0, 0, 0);
                for (auto & ri : rci.regionIndices){
                    auto & cam = context.views[ri.viewId].camera;
                    regionsCenterDirection += normalize(cam.spatialDirection(GetData(ri, context.regionsNets).center));
                    auto & regionOuterContourPixels = GetData(ri, context.regionsNets).contours.back();
                    for (auto & pixel : regionOuterContourPixels){
                        outerContourDirections.push_back(cam.spatialDirection(pixel));
                    }
                }
                regionsCenterDirection /= norm(regionsCenterDirection);
                rci.tangentialPlane = Plane3(regionsCenterDirection, regionsCenterDirection);
                std::tie(rci.xOnTangentialPlane, rci.yOnTangentialPlane) = 
                    ProposeXYDirectionsFromZDirection(rci.tangentialPlane.normal);

                // compute visual areas
                rci.regionVisualArea = ComputeVisualAreaOfDirections(rci.tangentialPlane, 
                    rci.xOnTangentialPlane, rci.yOnTangentialPlane,
                    outerContourDirections, false);
                rci.regionConvexContourVisualArea = ComputeVisualAreaOfDirections(rci.tangentialPlane,
                    rci.xOnTangentialPlane, rci.yOnTangentialPlane,
                    outerContourDirections, true);
                
                // calc sample points num with other regions
                rci.anchoredConnectionsWithOtherRegions = Rational(0.0, 0.0);
                for (auto & ri : rci.regionIndices){
                    auto & regions = context.regionsNets[ri.viewId].regions();
                    auto & bdHandles = regions.topo(ri.handle).uppers;
                    for (auto bh : bdHandles){
                        for (auto & pts : regions.data(bh).sampledPoints){
                            rci.anchoredConnectionsWithOtherRegions.denominator += pts.size();
                        }
                    }
                }

                // calc sample points num with other lines
                rci.anchoredConnectionsWithLines = Rational(0.0, 0.0);
                for (auto & p : context.regionLineConnections){
                    auto & ri = p.first.first;
                    if (context.regionConnectedComponentIds.at(ri) == rci.regionCCId){
                        rci.anchoredConnectionsWithLines.denominator += p.second.size();
                    }
                }

                return rci;
            }
            
            void UpdateRegionCCUsingLastInsertedAnchor(RegionCCRecInfo & rci, const RecContext & context){
                // add new candidate plane root
                double scale = context.initialBoundingBox.outerSphere().radius;
                for (auto & vp : context.vanishingPoints){
                    Plane3 plane(rci.anchoredConnections.back(), vp);
                    if (OPT_IgnoreTooSkewedPlanes){
                        if (norm(plane.root()) <= scale / 5.0)
                            continue;
                    }
                    if (OPT_IgnoreTooFarAwayPlanes){
                        bool valid = true;
                        for (auto & ri : rci.regionIndices){
                            if (!valid)
                                break;
                            auto & rd = GetData(ri, context.regionsNets);
                            if (rd.contours.back().size() < 3)
                                continue;
                            auto & cam = context.views[ri.viewId].camera;
                            for (int i = 0; i < rd.contours.back().size(); i++){
                                auto dir = cam.spatialDirection(ToPoint2(rd.contours.back()[i]));
                                auto intersectionOnPlane = IntersectionOfLineAndPlane(InfiniteLine3(Point3(0, 0, 0), dir), plane).position;
                                if (norm(intersectionOnPlane) > scale * 5.0){
                                    valid = false;
                                    break;
                                }
                            }
                        }
                        if (!valid)
                            continue;
                    }

                    static const double distFromPointToPlaneThres = scale / 12.0;

                    // insert new root data
                    auto & pcd = rci.candidatePlanesByRoot[plane.root()];
                    pcd.plane = plane;

                    // collect distance votes
                    double distVotes = 0.0;
                    std::vector<Vec3> nearbyAnchors;
                    for (int i = 0; i < rci.anchoredConnections.size(); i++){
                        double distanceToPlane = plane.distanceTo(rci.anchoredConnections[i]);
                        if (distanceToPlane > distFromPointToPlaneThres)
                            continue;
                        distVotes += Gaussian(distanceToPlane, distFromPointToPlaneThres);
                        pcd.inlierAnchors.push_back(i);
                        nearbyAnchors.push_back(rci.anchoredConnections[i]);
                    }
                    pcd.regionInlierAnchorsDistanceVotesSum = distVotes;
                    pcd.regionInlierAnchorsConvexContourVisualArea = 
                        ComputeVisualAreaOfDirections(rci.tangentialPlane, 
                        rci.xOnTangentialPlane, rci.yOnTangentialPlane, nearbyAnchors, true);
                }
            }

            inline void InsertLineAnchorToRegionCC(RegionCCRecInfo & rci, 
                const LineIndex & li, const Point3 & anchor,
                const RecContext & context){                
                rci.anchoredConnections.push_back(anchor);
                rci.anchoredConnectionsWithLines.numerator++;
                UpdateRegionCCUsingLastInsertedAnchor(rci, context);
            }

            inline void InsertRegionAnchorToRegionCC(RegionCCRecInfo & rci, 
                const RegionIndex & ri, const Point3 & anchor,
                const RecContext & context){
                rci.anchoredConnections.push_back(anchor);
                rci.anchoredConnectionsWithOtherRegions.numerator++;
                UpdateRegionCCUsingLastInsertedAnchor(rci, context);
            }

            void RegisterCurrentChoicesAndProbabilities(const RegionCCRecInfo & rci, 
                std::vector<Choice> & choices, std::vector<double> & probabilities, 
                const RecContext & context, double baseProb, int maxChoiceNum){

                std::vector<Scored<Choice>> newChoices;

                int planeId = 0;
                double maxMeanVote = 0.0;
                for (auto & c : rci.candidatePlanesByRoot){
                    maxMeanVote = std::max(maxMeanVote, c.second.regionInlierAnchorsDistanceVotesSum / c.second.inlierAnchors.size());
                }
                for (auto & c : rci.candidatePlanesByRoot){
                    auto & candidatePlaneData = c.second;
                    Choice choice = { true, rci.regionCCId, planeId++ };
                    auto completenessWithOtherRegions = rci.anchoredConnectionsWithOtherRegions.value(0.0);
                    auto completenessWithOtherLines = rci.anchoredConnectionsWithLines.value(0.0);
                    auto inlierOccupationRatio =
                        candidatePlaneData.regionInlierAnchorsConvexContourVisualArea / rci.regionConvexContourVisualArea;
                    //assert(candidatePlaneData.regionInlierAnchorsConvexContourVisualArea <= rci.regionConvexContourVisualArea);
                    double probability =
                        ((completenessWithOtherRegions * 0.7 + completenessWithOtherLines * 0.3) *
                        (inlierOccupationRatio > 0.3 ? 1.0 : 0.01)) *
                        c.second.regionInlierAnchorsDistanceVotesSum / c.second.inlierAnchors.size();

                    assert(c.second.regionInlierAnchorsDistanceVotesSum / c.second.inlierAnchors.size() <= 1.0);
                    newChoices.push_back(ScoreAs(choice, probability + baseProb));
                }

                std::sort(newChoices.begin(), newChoices.end(), std::greater<void>());
                for (int i = 0; i < std::min<int>(maxChoiceNum, newChoices.size()); i++){
                    choices.push_back(newChoices[i].component);
                    probabilities.push_back(newChoices[i].score);
                }
            }

            Plane3 ChoosePrediction(const std::vector<RegionCCRecInfo> & regionCCRecInfos, const Choice & choice){
                assert(choice.isRegionCC);
                return (regionCCRecInfos[choice.ccId].candidatePlanesByRoot.begin() + choice.choiceId)->second.plane;
            }

            



            // line cc reconstruction information
            struct LineCCRecInfo {
                int lineCCId;
                ComponentIndexHashSet<LineIndex> lineIndices;
                Rational anchoredConnectionsWithOtherRegions;
                std::vector<Point3> anchoredConnections;
                using DepthConfidenceMap = VecMap<double, 1, double>;
                DepthConfidenceMap candidateDepthFactors;
            };

            // initialize line cc reconstruction information
            LineCCRecInfo CreateLineCCRecInfo(int lineCCid, const RecContext & context){
                assert(context.views.size() == context.linesNets.size());
                LineCCRecInfo lci;
                lci.lineCCId = lineCCid;
                lci.candidateDepthFactors = LineCCRecInfo::DepthConfidenceMap(1e-4);
                lci.candidateDepthFactors[1.0] = 0.1; // initialize with a 1.0 depth for start
                // collect lis
                for (auto & p : context.lineConnectedComponentIds){
                    if (p.second == lineCCid){
                        lci.lineIndices.insert(p.first);
                    }
                }
                // calc sample points num with regions
                lci.anchoredConnectionsWithOtherRegions = Rational(0, 0);
                for (auto & p : context.regionLineConnections){
                    auto & li = p.first.second;
                    if (context.lineConnectedComponentIds.at(li) == lineCCid){
                        lci.anchoredConnectionsWithOtherRegions.denominator += p.second.size();
                    }
                }
                return lci;
            }

            void UpdateLineCCUsingLastInsertedAnchor(LineCCRecInfo & lci, 
                const RegionIndex & ri, const LineIndex & li, 
                const RecContext & context){
                auto & line = context.reconstructedLines.at(li);
                auto & anchor = lci.anchoredConnections.back();
                double depthVarOnLine = norm(DistanceBetweenTwoLines(line.infinieLine(), InfiniteLine3(Point3(0, 0, 0), anchor))
                    .second.second);
                double depthValueOnRegion = norm(anchor);
                if (!IsInfOrNaN(depthVarOnLine) && !IsInfOrNaN(depthValueOnRegion))
                    lci.candidateDepthFactors[depthValueOnRegion / depthVarOnLine] += 1.0;
            }

            inline void InsertRegionAnchorToLineCC(LineCCRecInfo & lci,
                const RegionIndex & ri, const LineIndex & li, const Point3 & anchor,
                const RecContext & context){
                lci.anchoredConnections.push_back(anchor);
                lci.anchoredConnectionsWithOtherRegions.numerator++; 
                UpdateLineCCUsingLastInsertedAnchor(lci, ri, li, context);
            }

            void RegisterCurrentChoicesAndProbabilities(const LineCCRecInfo & lci,
                std::vector<Choice> & choices, std::vector<double> & probabilities,
                const RecContext & context, double baseProb, int maxChoiceNum){
                std::vector<Scored<Choice>> newChoices;
                int depthId = 0;
                double maxVote = 0.0;
                for (auto & c : lci.candidateDepthFactors)
                    maxVote = maxVote < c.second ? c.second : maxVote;
                for (auto & c : lci.candidateDepthFactors){
                    auto & candidateDepthVote = c.second;
                    Choice choice = { false, lci.lineCCId, depthId++ };
                    auto completenessWithOtherRegions = lci.anchoredConnectionsWithOtherRegions.value(0.0);
                    double probability = (completenessWithOtherRegions * 0.9
                        + double(lci.lineIndices.size()) / context.reconstructedLines.size() * 0.1)
                        * candidateDepthVote / maxVote;
                    newChoices.push_back(ScoreAs(choice, probability + baseProb));
                }

                std::sort(newChoices.begin(), newChoices.end(), std::greater<void>());
                for (int i = 0; i < std::min<int>(maxChoiceNum, newChoices.size()); i++){
                    choices.push_back(newChoices[i].component);
                    probabilities.push_back(newChoices[i].score);
                }
            }

            double ChoosePrediction(const std::vector<LineCCRecInfo> & lineCCRecInfos, const Choice & choice){
                assert(!choice.isRegionCC);
                return (lineCCRecInfos[choice.ccId].candidateDepthFactors.begin() + choice.choiceId)->first[0];
            }


            void DisplayReconstruction(int highlightedRegionCCId, int highlightedLineCCId,
                const std::set<int> & regionCCIdsNotDeterminedYet, const std::set<int> & lineCCIdsNotDeterminedYet,
                const std::vector<Plane3> & regionConnectedComponentPlanes,
                const std::vector<double> & lineConnectedComponentDepthFactors,
                const RecContext & context){

                std::vector<Line3> linesRepresentingSampledPoints;

                // line-region connections
                for (auto & pp : context.regionLineConnections){
                    RegionIndex ri = pp.first.first;
                    LineIndex li = pp.first.second;
                    int regionCCId = context.regionConnectedComponentIds.at(ri);
                    int lineCCId = context.lineConnectedComponentIds.at(li);
                    if (core::Contains(regionCCIdsNotDeterminedYet, regionCCId) || 
                        core::Contains(lineCCIdsNotDeterminedYet, lineCCId))
                        continue;

                    auto line = context.reconstructedLines.at(li);
                    double depthFactor = lineConnectedComponentDepthFactors[lineCCId];
                    line.first *= depthFactor;
                    line.second *= depthFactor;

                    const std::vector<Vec3> & selectedSampledPoints = pp.second;

                    for (const Vec3 & sampleRay : selectedSampledPoints){
                        Point3 pointOnLine = DistanceBetweenTwoLines(InfiniteLine3(Point3(0, 0, 0), sampleRay), line.infinieLine()).second.second;
                        Point3 pointOnRegion = IntersectionOfLineAndPlane(InfiniteLine3(Point3(0, 0, 0), sampleRay),
                            regionConnectedComponentPlanes[regionCCId]).position;
                        linesRepresentingSampledPoints.emplace_back(pointOnLine, pointOnRegion);
                    }
                }


                // paint regions
                std::vector<vis::SpatialProjectedPolygon> spps, highlightedSpps;
                spps.reserve(context.regionConnectedComponentIds.size());
                static const int stepSize = 10;

                for (auto & r : context.regionConnectedComponentIds){
                    auto & ri = r.first;
                    vis::SpatialProjectedPolygon spp;
                    int regionCCId = context.regionConnectedComponentIds.at(ri);
                    if (core::Contains(regionCCIdsNotDeterminedYet, regionCCId)) // igore not reconstructed regions 
                        continue;

                    spp.plane = regionConnectedComponentPlanes[regionCCId];
                    auto & rd = GetData(ri, context.regionsNets);
                    if (rd.contours.back().size() < 3)
                        continue;

                    spp.corners.reserve(rd.contours.back().size() / double(stepSize));
                    auto & cam = context.views[ri.viewId].camera;

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
                        if (context.regionConnectedComponentIds.at(ri) == highlightedRegionCCId){
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

                viz << vis::manip3d::SetDefaultColorTable(vis::CreateRandomColorTableWithSize(context.lineConnectedComponentsNum));

                // paint lines
                std::vector<Line3> highlightedLines;
                for (auto & l : context.reconstructedLines) {
                    int lineCCId = context.lineConnectedComponentIds.at(l.first);
                    if (core::Contains(lineCCIdsNotDeterminedYet, lineCCId))
                        continue;
                    auto line = l.second;
                    double depthFactor = lineConnectedComponentDepthFactors[lineCCId];
                    line.first *= depthFactor;
                    line.second *= depthFactor;
                    if (lineCCId == highlightedLineCCId)
                        highlightedLines.push_back(line);
                    viz << core::ClassifyAs(line, lineCCId);
                }

                viz << vis::manip3d::SetBackgroundColor(vis::ColorTag::White)
                    << vis::manip3d::Begin(spps)
                    << vis::manip3d::SetTexture(context.globalTexture)
                    << vis::manip3d::End
                    << vis::manip3d::SetDefaultLineWidth(6.0)
                    << vis::manip3d::SetDefaultForegroundColor(vis::ColorTag::Black)
                    << BoundingBoxOfContainer(highlightedSpps)
                    << BoundingBoxOfContainer(highlightedLines)
                    << vis::manip3d::SetWindowName("initial region planes and reconstructed lines")
                    << vis::manip3d::Show();

            }

        }


        void EstimateSpatialRegionPlanes(const std::vector<View<PerspectiveCamera>> & views,
            const std::vector<RegionsNet> & regionsNets, const std::vector<LinesNet> & linesNets,
            const std::array<Vec3, 3> & vanishingPoints,
            const ComponentIndexHashMap<std::pair<RegionIndex, RegionIndex>, double> & regionOverlappings,
            const ComponentIndexHashMap<std::pair<RegionIndex, LineIndex>, std::vector<Vec3>> & regionLineConnections,
            const ComponentIndexHashMap<std::pair<LineIndex, LineIndex>, Vec3> & interViewLineIncidences,
            int regionConnectedComponentsNum, const ComponentIndexHashMap<RegionIndex, int> & regionConnectedComponentIds,
            int lineConnectedComponentsNum, const ComponentIndexHashMap<LineIndex, int> & lineConnectedComponentIds,
            ComponentIndexHashMap<LineIndex, Line3> & reconstructedLines,
            ComponentIndexHashMap<RegionIndex, Plane3> & reconstructedPlanes,
            const Image & globalTexture){

            std::cout << "invoking " << __FUNCTION__ << std::endl;

            Box3 bbox;
            for (auto & l : reconstructedLines)
                bbox |= BoundingBox(l.second);
            double scale = bbox.outerSphere().radius;

            const RecContext context = { views, regionsNets, linesNets, vanishingPoints,
                regionOverlappings, regionLineConnections, interViewLineIncidences,
                regionConnectedComponentsNum, regionConnectedComponentIds,
                lineConnectedComponentsNum, lineConnectedComponentIds,
                reconstructedLines, reconstructedPlanes, globalTexture,
                bbox
            };


            // initial cc ids status
            std::vector<int> regionCCIds(regionConnectedComponentsNum);
            std::iota(regionCCIds.begin(), regionCCIds.end(), 0);
            std::set<int> initialRegionCCIdsNotDeterminedYet(regionCCIds.begin(), regionCCIds.end());

            std::vector<int> lineCCIds(lineConnectedComponentsNum);
            std::iota(lineCCIds.begin(), lineCCIds.end(), 0);
            std::set<int> initialLineCCIdsNotDeterminedYet(lineCCIds.begin(), lineCCIds.end());


            // initialize region planes
            std::vector<Plane3> initialRegionConnectedComponentPlanes(regionConnectedComponentsNum);
            for (auto & r : regionConnectedComponentIds){
                auto & ri = r.first;
                auto & rd = GetData(ri, regionsNets);
                Vec3 centerDir = views[ri.viewId].camera.spatialDirection(rd.center);
                int regionCCId = regionConnectedComponentIds.at(ri);
                initialRegionConnectedComponentPlanes[regionCCId].anchor = normalize(centerDir) * scale;
                initialRegionConnectedComponentPlanes[regionCCId].normal = normalize(centerDir);
            }
            // initialize line depth factors
            std::vector<double> initialLineConnectedComponentDepthFactors(lineConnectedComponentsNum, 1.0);


            // try randomly
            int trialNum = OPT_RandomRepeatNum;

            // start MC reasoning
            std::vector<Scored<std::pair<std::vector<Plane3>, std::vector<double>>>>
                candidates(trialNum,
                ScoreAs(std::make_pair(initialRegionConnectedComponentPlanes, initialLineConnectedComponentDepthFactors),
                0.0));

            auto task = [&candidates, &context, &initialRegionCCIdsNotDeterminedYet, &initialLineCCIdsNotDeterminedYet](int t){
                std::cout << "task: " << t << std::endl;
                auto & candidate = candidates[t];

                // random engine initialized
                std::random_device rd;
                std::default_random_engine gen(rd());

                // build initial reconstruction info data table
                std::vector<RegionCCRecInfo> regionCCRecInfos;
                regionCCRecInfos.reserve(context.regionConnectedComponentsNum);
                for (int i = 0; i < context.regionConnectedComponentsNum; i++)
                    regionCCRecInfos.push_back(CreateRegionCCRecInfo(i, context));

                std::vector<LineCCRecInfo> lineCCRecInfos;
                lineCCRecInfos.reserve(context.lineConnectedComponentsNum);
                for (int i = 0; i < context.lineConnectedComponentsNum; i++)
                    lineCCRecInfos.push_back(CreateLineCCRecInfo(i, context));


                // undetermined checkers
                std::set<int> regionCCIdsNotDeterminedYet = initialRegionCCIdsNotDeterminedYet;
                std::set<int> lineCCIdsNotDeterminedYet = initialLineCCIdsNotDeterminedYet;

                // choices and probabilities
                std::vector<Choice> choices;
                std::vector<double> choiceProbabilities;
                choices.reserve(regionCCIdsNotDeterminedYet.size() + lineCCIdsNotDeterminedYet.size());
                choiceProbabilities.reserve(regionCCIdsNotDeterminedYet.size() + lineCCIdsNotDeterminedYet.size());

                // start expansion
                std::cout << "start expansion" << std::endl;
                while ((regionCCIdsNotDeterminedYet.size() + lineCCIdsNotDeterminedYet.size()) > 0){

                    // collect choices and probabilities
                    choices.clear();
                    choiceProbabilities.clear();
                    for (int rid : regionCCIdsNotDeterminedYet){
                        RegisterCurrentChoicesAndProbabilities(regionCCRecInfos[rid], choices, choiceProbabilities, context, 
                            1e-5, 1);
                    }
                    for (int lid : lineCCIdsNotDeterminedYet){
                        RegisterCurrentChoicesAndProbabilities(lineCCRecInfos[lid], choices, choiceProbabilities, context, 
                            1e-5, 1);
                    }

                    assert(choices.size() == choiceProbabilities.size());

                    if (std::accumulate(choiceProbabilities.begin(), choiceProbabilities.end(), 0.0) == 0.0){
                        std::cerr << "all zero probabilities!" << std::endl;
                        break;
                    }

                    // a workaround constructor since VS2013 lacks the range iterator constructor for std::discrete_distribution
                    int ord = 0;
                    std::discrete_distribution<int> distribution(choiceProbabilities.size(),
                        0.0, 1000.0,
                        [&choiceProbabilities, &ord](double){
                        return choiceProbabilities[ord++];
                    });

                    // made choice
                    const Choice & choice = choices[distribution(gen)];
                    if (choice.isRegionCC){ // a region cc is chosen
                        if (OPT_DisplayMessages)
                            std::cout << "chosen unit - region cc: " << choice.ccId << "  plane chosen: " << choice.choiceId << std::endl;
                        candidate.component.first[choice.ccId] = ChoosePrediction(regionCCRecInfos, choice);
                        regionCCIdsNotDeterminedYet.erase(choice.ccId);

                        // update related line ccs
                        for (auto & pp : context.regionLineConnections){
                            LineIndex li = pp.first.second;
                            RegionIndex ri = pp.first.first;
                            int thisLineCCId = context.lineConnectedComponentIds.at(li);
                            int thisRegionCCId = context.regionConnectedComponentIds.at(ri);
                            if (!core::Contains(lineCCIdsNotDeterminedYet, thisLineCCId)){ // checked already
                                continue;
                            }
                            if (thisRegionCCId == choice.ccId){ // related
                                auto & plane = candidate.component.first[thisRegionCCId];
                                auto & samplePoints = pp.second;
                                for (auto & p : samplePoints){
                                    auto anchor = IntersectionOfLineAndPlane(InfiniteLine3(Point3(0, 0, 0), p), plane).position;
                                    InsertRegionAnchorToLineCC(lineCCRecInfos[thisLineCCId], ri, li, anchor, context);
                                }
                            }
                        }

                        // update related region ccs 
                        for (int i = 0; i < context.views.size(); i++){
                            for (auto & b : context.regionsNets[i].regions().elements<1>()){
                                auto ri1 = RegionIndex{ i, b.topo.lowers[0] };
                                auto ri2 = RegionIndex{ i, b.topo.lowers[1] };
                                int thisRegionCCId1 = context.regionConnectedComponentIds.at(ri1);
                                int thisRegionCCId2 = context.regionConnectedComponentIds.at(ri2);
                                int determinedRegionCCId = choice.ccId, notDeterminedRegionCCId;
                                RegionIndex determinedRi;
                                if (thisRegionCCId1 == choice.ccId && core::Contains(regionCCIdsNotDeterminedYet, thisRegionCCId2)){
                                    determinedRi = ri1;
                                    notDeterminedRegionCCId = thisRegionCCId2;
                                }
                                else if (thisRegionCCId2 == choice.ccId && core::Contains(regionCCIdsNotDeterminedYet, thisRegionCCId1)){
                                    determinedRi = ri2;
                                    notDeterminedRegionCCId = thisRegionCCId1;
                                }
                                else{
                                    continue;
                                }
                                auto & plane = candidate.component.first[determinedRegionCCId];
                                auto & cam = context.views[i].camera;
                                for (auto & pts : b.data.sampledPoints){
                                    for (auto & p : pts){
                                        auto dir = cam.spatialDirection(p);
                                        auto anchor = IntersectionOfLineAndPlane(InfiniteLine3(Point3(0, 0, 0), dir), plane).position;
                                        InsertRegionAnchorToRegionCC(regionCCRecInfos[notDeterminedRegionCCId], determinedRi, anchor, context);
                                    }
                                }
                            }
                        }

                    }
                    else{ // a line cc is chosen
                        if (OPT_DisplayMessages)
                            std::cout << "chosen unit - line cc: " << choice.ccId << "  depthfactor chosen: " << choice.choiceId << std::endl;
                        candidate.component.second[choice.ccId] = ChoosePrediction(lineCCRecInfos, choice);
                        lineCCIdsNotDeterminedYet.erase(choice.ccId);

                        // update related unchecked region anchors
                        for (auto & pp : context.regionLineConnections){
                            LineIndex li = pp.first.second;
                            RegionIndex ri = pp.first.first;
                            int thisLineCCId = context.lineConnectedComponentIds.at(li);
                            int thisRegionCCId = context.regionConnectedComponentIds.at(ri);
                            if (!core::Contains(regionCCIdsNotDeterminedYet, thisRegionCCId)){ // checked already, ignore
                                continue;
                            }
                            if (thisLineCCId == choice.ccId){ // related 
                                auto & samplePoints = pp.second;
                                const auto & line = context.reconstructedLines.at(li);
                                for (auto & p : samplePoints){ // insert anchors into the rec info of the related region
                                    auto pOnLine = DistanceBetweenTwoLines(line.infinieLine(), InfiniteLine3(Point3(0, 0, 0), p))
                                        .second.second;
                                    InsertLineAnchorToRegionCC(regionCCRecInfos[thisRegionCCId], li,
                                        pOnLine * candidate.component.second[thisLineCCId], context);
                                }
                            }
                        }
                    }

                } // while
                std::cout << "expansion done" << std::endl;


                // score this candidate
                double distanceSumOfRegionRegionConnections = 0.0;
                double distanceSumOfRegionLineConnections = 0.0;

                auto & regionConnectedComponentPlanes = candidate.component.first;
                auto & lineConnectedComponentDepthFactors = candidate.component.second;

                // region region
                for (int i = 0; i < context.views.size(); i++){
                    auto & cam = context.views[i].camera;
                    for (auto & b : context.regionsNets[i].regions().elements<1>()){
                        auto ri1 = RegionIndex{ i, b.topo.lowers[0] };
                        auto ri2 = RegionIndex{ i, b.topo.lowers[1] };
                        int regionCCId1 = context.regionConnectedComponentIds.at(ri1);
                        int regionCCId2 = context.regionConnectedComponentIds.at(ri2);
                        for (auto & pts : b.data.sampledPoints){
                            for (auto & p : pts){
                                auto dir = cam.spatialDirection(p);
                                auto anchor1 = IntersectionOfLineAndPlane(InfiniteLine3(Point3(0, 0, 0), dir),
                                    regionConnectedComponentPlanes[regionCCId1]).position;
                                auto anchor2 = IntersectionOfLineAndPlane(InfiniteLine3(Point3(0, 0, 0), dir),
                                    regionConnectedComponentPlanes[regionCCId2]).position;
                                double distance = Distance(anchor1, anchor2);
                                distanceSumOfRegionRegionConnections += distance;
                            }
                        }
                    }
                }
                // region line
                for (auto & pp : context.regionLineConnections){
                    LineIndex li = pp.first.second;
                    RegionIndex ri = pp.first.first;
                    int lineCCId = context.lineConnectedComponentIds.at(li);
                    int regionCCId = context.regionConnectedComponentIds.at(ri);
                    auto & samplePoints = pp.second;
                    auto line = context.reconstructedLines.at(li);
                    line.first *= lineConnectedComponentDepthFactors[lineCCId];
                    line.second *= lineConnectedComponentDepthFactors[lineCCId];
                    for (auto & p : samplePoints){ // insert anchors into the rec info of the related region
                        auto pOnLine = DistanceBetweenTwoLines(line.infinieLine(), InfiniteLine3(Point3(0, 0, 0), p))
                            .second.second;
                        auto pOnRegion = IntersectionOfLineAndPlane(InfiniteLine3(Point3(0, 0, 0), p),
                            regionConnectedComponentPlanes[regionCCId]).position;
                        double distance = Distance(pOnLine, pOnRegion);
                        distanceSumOfRegionLineConnections += distance;
                    }
                }

                std::cout << "distance sum of region-region connections: " << distanceSumOfRegionRegionConnections << std::endl;
                std::cout << "distance sum of region-line connections: " << distanceSumOfRegionLineConnections << std::endl;

                candidate.score = -(distanceSumOfRegionRegionConnections + distanceSumOfRegionLineConnections);

                if (OPT_DisplayOnEachTrial){
                    IF_DEBUG_USING_VISUALIZERS{
                        DisplayReconstruction(-1, -1, {}, {},
                        candidate.component.first, candidate.component.second, context);
                    }
                }

            }; // task

            // run tasks
            int parallelNum = OPT_NumThreads;
            for (int i = 0; i < trialNum; i += parallelNum){
                std::vector<std::thread> parallelThreads(parallelNum);
                for (int t = i; t < std::min(trialNum, i + parallelNum); t++){
                    parallelThreads[t-i] = std::thread(task, t);
                }
                for (int t = i; t < std::min(trialNum, i + parallelNum); t++){
                    parallelThreads[t-i].join();
                }
            }


            // select best candidate

            auto result = std::max_element(candidates.begin(), candidates.end())->component;

            // visualize result of this task
            if (OPT_DisplayAtLast){
                IF_DEBUG_USING_VISUALIZERS{
                    DisplayReconstruction(-1, -1, {}, {},
                    result.first, result.second, context);
                }
            }



        }

    }
}