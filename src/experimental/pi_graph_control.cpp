#include "../core/clock.hpp"
#include "../core/utility.hpp"

#include "pi_graph_control.hpp"

namespace pano {
    namespace experimental {   


        std::vector<int> CollectSegsIntersectingDirection(const Vec3 & direction,
            bool alsoConsiderBackward, const PIGraph & mg, double rangeAngle) {
            // find peaky regions
            std::vector<int> peakySegs;
            for (int i = 0; i < mg.nsegs; i++) {
                double radiusAngle = 0.0;
                for (int bnd : mg.seg2bnds[i]) {
                    for (int bp : mg.bnd2bndPieces[bnd]) {
                        for (auto & c : mg.bndPiece2dirs[bp]) {
                            double angle = AngleBetweenDirections(mg.seg2center[i], c);
                            if (angle > radiusAngle) {
                                radiusAngle = angle;
                            }
                        }
                    }
                }

                if ((alsoConsiderBackward ? AngleBetweenUndirectedVectors(direction, mg.seg2center[i]) :
                    AngleBetweenDirections(direction, mg.seg2center[i])) > radiusAngle) {
                    continue;
                }

                float ppcFocal = 100.0f;
                int ppcSize = 2 * radiusAngle * ppcFocal + 10;
                Vec3 x;
                std::tie(x, std::ignore) = ProposeXYDirectionsFromZDirection(mg.seg2center[i]);
                PartialPanoramicCamera ppc(ppcSize, ppcSize, ppcFocal, Point3(0, 0, 0), mg.seg2center[i], x);
                Imageub mask = Imageub::zeros(ppc.screenSize());

                // project contours to ppc
                auto & contours = mg.seg2contours[i];
                std::vector<std::vector<Point2i>> contourProjs(contours.size());
                for (int k = 0; k < contours.size(); k++) {
                    auto & contourProj = contourProjs[k];
                    contourProj.reserve(contours[k].size());
                    for (auto & d : contours[k]) {
                        contourProj.push_back(ecast<int>(ppc.toScreen(d)));
                    }
                }
                cv::fillPoly(mask, contourProjs, (uint8_t)1);

                // intersection test
                auto p = ToPixel(ppc.toScreen(direction));
                auto p2 = ToPixel(ppc.toScreen(-direction));

                int dilateSize = ppcFocal * rangeAngle;
                bool intersected = false;
                for (int x = -dilateSize; x <= dilateSize; x++) {
                    if (intersected)
                        break;
                    for (int y = -dilateSize; y <= dilateSize; y++) {
                        if (intersected)
                            break;
                        auto pp = Pixel(p.x + x, p.y + y);
                        if (Contains(mask, pp) && mask(pp)) {
                            peakySegs.push_back(i);
                            intersected = true;
                        }
                        if (alsoConsiderBackward) {
                            auto pp2 = Pixel(p2.x + x, p2.y + y);
                            if (Contains(mask, pp2) && mask(pp2)) {
                                peakySegs.push_back(i);
                                intersected = true;
                            }
                        }
                    }
                }

            }

            return peakySegs;

        }



        namespace {

            // returns false if confliction occurs
            bool MakeRegionPlaneUsable(int seg, bool usable, PIGraph & mg) {
                auto & p = mg.seg2control[seg];
                if (p.used == usable)
                    return true;
                if (p.used && !usable) {
                    p.used = false;
                    p.orientationClaz = p.orientationNotClaz = -1;
                    return true;
                }
                p.orientationClaz = p.orientationNotClaz = -1;
                return true;
            }


            // returns false if confliction occurs
            bool MakeRegionPlaneToward(int seg, int normalVPId, PIGraph & mg) {
                auto & p = mg.seg2control[seg];
                if (!p.used)
                    return true;
                assert(normalVPId != -1);
                if (p.orientationClaz != -1) {
                    if (p.orientationClaz != normalVPId)
                        return false;
                    return true;
                }
                if (p.orientationNotClaz == -1) {
                    p.orientationClaz = normalVPId;
                    return true;
                }
                auto & dir = mg.vps[p.orientationNotClaz];
                if (IsFuzzyPerpendicular(mg.vps[normalVPId], dir)) {
                    p.orientationClaz = normalVPId;
                    p.orientationNotClaz = -1;
                    return true;
                }
                return false;
            }

            // returns false if confliction occurs
            bool MakeRegionPlaneAlsoAlong(int seg, int alongVPId, PIGraph & mg) {
                auto & p = mg.seg2control[seg];
                if (!p.used)
                    return true;
                assert(alongVPId != -1);
                auto & dir = mg.vps[alongVPId];
                if (p.orientationClaz != -1) {
                    auto & normal = mg.vps[p.orientationClaz];
                    return IsFuzzyPerpendicular(normal, dir);
                }
                if (p.orientationNotClaz == -1) {
                    p.orientationNotClaz = alongVPId;
                    return true;
                }
                if (p.orientationNotClaz == alongVPId)
                    return true;

                auto newNormal = dir.cross(mg.vps[p.orientationNotClaz]);
                double minAngle = M_PI;
                for (int i = 0; i < mg.vps.size(); i++) {
                    double angle = AngleBetweenUndirectedVectors(mg.vps[i], newNormal);
                    if (angle < minAngle) {
                        p.orientationClaz = i;
                        minAngle = angle;
                    }
                }
                if (p.orientationClaz != -1) {
                    p.orientationNotClaz = -1;
                    return true;
                }
                return false;
            }

        }


        void AttachPrincipleDirectionConstraints(PIGraph & mg, double angle) {
            SetClock();
            std::vector<std::vector<int>> peakySegs(mg.vps.size());
            for (int i = 0; i < mg.vps.size(); i++) {
                peakySegs[i] = CollectSegsIntersectingDirection(mg.vps[i], true, mg, angle);
            }
            for (int i = 0; i < mg.vps.size(); i++) {
                auto & vp = mg.vps[i];
                auto & segs = peakySegs[i];
                for (auto seg : segs) {
                    if (!mg.seg2control[seg].used) {
                        continue;
                    }
                    MakeRegionPlaneToward(seg, i, mg);
                    for (int lp : mg.seg2linePieces[seg]) {
                        if (mg.linePiece2samples[lp].size() < 2)
                            continue;
                        int & lc = mg.lines[mg.linePiece2line[lp]].claz;
                        if (lc == i) {
                            lc = -1;
                        }
                    }
                }
            }
        }

        void AttachWallConstraints(PIGraph & mg, double rangeAngle, int vertVPId) {
            SetClock();
            auto & vertical = mg.vps[vertVPId];

            std::vector<int> horizontalSegs;
            for (int seg = 0; seg < mg.nsegs; seg ++) {
                if (!mg.seg2control[seg].used)
                    continue;
                auto & contours = mg.seg2contours[seg];
                bool intersected = false;
                for (auto & cs : contours) {
                    if (intersected)
                        break;
                    for (auto & c : cs) {
                        double angle = M_PI_2 - AngleBetweenUndirectedVectors(c, vertical);
                        if (angle <= rangeAngle) {
                            intersected = true;
                            break;
                        }
                    }
                }

                if (intersected) {
                    horizontalSegs.push_back(seg);
                }
            }
            for (auto h : horizontalSegs) {
                MakeRegionPlaneAlsoAlong(h, vertVPId, mg);
            }
        }

        void DisableTopSeg(PIGraph & mg) {

        }

        void DisableBottomSeg(PIGraph & mg) {

        }

        void DisableInvalidConstraints(PIGraph & mg) {

        }

        void AttachGCConstraints(PIGraph & mg, const Image5d & gc) {

        }


        void DetectAndApplyOcclusions(PIGraph & mg) {

        }



    }
}