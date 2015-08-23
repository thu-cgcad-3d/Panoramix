#include "../core/containers.hpp"

#include "rl_graph_modeling.hpp"

namespace pano {

    namespace experimental {

        std::vector<SectionalPiece> MakeSectionalPieces(const HandledTable<RegionHandle, std::vector<Polygon3>> & polygons,
            const Plane3 & cutplane, bool removeIllPosedPolygons) {

            if (removeIllPosedPolygons) {
                /// todo
            }

            Box3 bboxAll;
            for (auto & p : polygons) {
                bboxAll |= BoundingBoxOfContainer(p);
            }
            double gsize = std::cbrt(bboxAll.size(0) * bboxAll.size(1) * bboxAll.size(2));

            std::vector<SectionalPiece> segments;
            std::vector<double> startAngles, endAngles;

            Vec3 x, y;
            std::tie(x, y) = ProposeXYDirectionsFromZDirection(cutplane.normal);
            Point3 original = cutplane.root();

            RTree<Box2, int> segmentIdRTree;
            Box2 bboxOnCutplane;

            for (auto it = polygons.begin(); it != polygons.end(); ++it) {
                RegionHandle rh = it.hd();

                for (const Polygon3 & polygon : *it) {
                    if (polygon.corners.size() <= 2)
                        continue;

                    auto plane = polygon.plane();

                    if (core::IsFuzzyParallel(plane.normal, cutplane.normal, 0.01))
                        continue;

                    Vec3 along = normalize(cutplane.normal.cross(plane.normal)); // intersection line direction
                    std::vector<Scored<std::pair<Point3, Point2>>> cutpoints;

                    for (int i = 1; i <= polygon.corners.size(); i++) {
                        auto & lastP = polygon.corners[i - 1];
                        auto & p = polygon.corners[i % polygon.corners.size()];
                        double lastDist = cutplane.signedDistanceTo(lastP);
                        double dist = cutplane.signedDistanceTo(p);

                        if (lastDist < 0 && dist > 0 || lastDist > 0 && dist < 0) { // intersection!
                            Point3 intersection = (lastP * abs(dist) + p * abs(lastDist)) / (abs(dist) + abs(lastDist));
                            double order = intersection.dot(along);
                            Point2 projOnCutPlane((intersection - original).dot(x), (intersection - original).dot(y));
                            cutpoints.push_back(ScoreAs(std::make_pair(intersection, projOnCutPlane), order));
                        }
                    }

                    if (cutpoints.empty()) {
                        continue;
                    }
                    assert(cutpoints.size() % 2 == 0);
                    if (cutpoints.size() % 2 != 0) {
                        std::cout << "ODD cutpoints num!!!" << std::endl;
                    }

                    // chain the segments
                    std::sort(cutpoints.begin(), cutpoints.end());
                    for (int i = 0; i < cutpoints.size(); i += 2) {
                        auto box2 = BoundingBox(cutpoints[i].component.second) | BoundingBox(cutpoints[i + 1].component.second);
                        bboxOnCutplane |= box2;

                        SectionalPiece segment;
                        segment.rh = rh;
                        segment.range.first = cutpoints[i].component.first;
                        segment.range.second = cutpoints[i + 1].component.first;
                        segments.push_back(segment);
                        segmentIdRTree.insert(box2, segments.size() - 1);

                        startAngles.push_back(SignedAngleBetweenDirections(Vec2(1, 0), cutpoints[i].component.second));
                        endAngles.push_back(SignedAngleBetweenDirections(Vec2(1, 0), cutpoints[i + 1].component.second));
                    }
                }
            }

            if (segments.empty())
                return segments;

            if (sqrt(bboxOnCutplane.size(0) * bboxOnCutplane.size(1)) < gsize * 0.1)
                return std::vector<SectionalPiece>();

            std::vector<int> ids(segments.size());
            std::iota(ids.begin(), ids.end(), 0);
            std::sort(ids.begin(), ids.end(), [&startAngles](int id1, int id2) {return startAngles[id1] < startAngles[id2]; });

            std::vector<int> ids2(1, ids.front());
            ids2.reserve(ids.size());
            for (int i = 1; i < ids.size(); i++) {
                int id = ids[i];
                Radian lastEndAngle = endAngles[ids2.back()];
                Radian thisStartAngle = startAngles[id];
                if (thisStartAngle <= lastEndAngle) {
                    continue;
                }
                ids2.push_back(id);
            }

            std::vector<SectionalPiece> segs;
            segs.reserve(segments.size());
            for (int id : ids2) {
                segs.push_back(std::move(segments[id]));
            }

            return segs;
        }



        Chain3 MakeChain(const std::vector<SectionalPiece> & pieces, bool closed) {
            Chain3 chain;
            chain.closed = closed;
            chain.points.reserve(pieces.size() * 2);
            for (auto & seg : pieces) {
                chain.points.push_back(seg.range.first);
                chain.points.push_back(seg.range.second);
            }

            if (!closed) {

            }

            return chain;
        }





        LayeredShape3 OptimizedModel(const std::vector<SectionalPiece> & pieces, 
            const RLGraph & mg, const RLGraphControls & controls) {

            NOT_IMPLEMENTED_YET();




        }












        std::pair<double, double> EstimateEffectiveRangeAlongDirection(
            const HandledTable<RegionHandle, std::vector<Polygon3>> & polygons,
            const Vec3 & direction, double stepLen, double minEffectiveAreaRatio,
            double gamma1, double gamma2) {

            double minv = std::numeric_limits<double>::max();
            double maxv = std::numeric_limits<double>::lowest();

            auto ndir = normalize(direction);
            for (auto it = polygons.begin(); it != polygons.end(); ++it) {
                for (auto & ply : *it) {
                    for (auto & p : ply.corners) {
                        double pos = p.dot(ndir);
                        if (minv > pos) minv = pos;
                        if (maxv < pos) maxv = pos;
                    }
                }
            }

            assert(maxv - minv > 0);

            std::vector<double> areas;
            std::vector<double> dareas;
            for (double x = minv; x <= maxv; x += stepLen) {
                Plane3 cutplane(ndir * x, ndir);
                auto chain = MakeChain(MakeSectionalPieces(polygons, cutplane));

                // optimize chain


                double area = chain.size() < 3 ? 0.0 : Area(chain);
                assert(area >= 0.0);
                dareas.push_back(area - (areas.empty() ? 0.0 : areas.back()));
                areas.push_back(area);
            }

            std::pair<double, double> range(0, 0);
            if (areas.size() < 2) {
                return range;
            }
            dareas[0] = dareas[1];

            size_t maxCutAreaId = std::max_element(areas.begin(), areas.end()) - areas.begin();
            double maxCutArea = areas[maxCutAreaId];

            int start = 0;
            while (start < areas.size() && areas[start] < minEffectiveAreaRatio * maxCutArea) ++start;
            int end = areas.size() - 1;
            while (end >= 0 && areas[end] < minEffectiveAreaRatio * maxCutArea) --end;

            std::cout << "all: " << areas.size() << std::endl;
            std::cout << "start: " << start << std::endl;
            std::cout << "end: " << end << std::endl;

            if (start >= end)
                return range;
            if (start > maxCutAreaId || end < maxCutAreaId)
                return range;

            range.first = start * stepLen + minv;
            range.second = end * stepLen + minv;
            //for (int i = start; i < maxCutAreaId; i++){
            //    if (i != 0 && (dareas[i] - dareas[i - 1]) / maxCutArea * ((maxv - minv) / stepLen) < -gamma1){
            //        range.first = i * stepLen + minv;
            //        break;
            //    }
            //}

            //for (int i = end; i > maxCutAreaId; i--){
            //    if (i != areas.size() - 1 && (dareas[i] - dareas[i + 1]) / maxCutArea * ((maxv - minv) / stepLen) < -gamma2){
            //        range.second = i * stepLen + minv;
            //        break;
            //    }
            //}

            return range;

        }




    }

}