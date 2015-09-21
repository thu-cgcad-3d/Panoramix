#include "../core/image.hpp"
#include "../core/utility.hpp"
#include "../core/containers.hpp"

#include "rl_graph.hpp"
#include "pi_graph.hpp"

namespace pano {

    namespace experimental {

     

        namespace {


            template <class T>
            inline size_t ElementsNum(const std::vector<T> & v) {
                return v.size();
            }

            template <class T>
            inline size_t ElementsNum(const std::vector<std::vector<T>> & v) {
                size_t n = 0;
                for (auto & vv : v) {
                    n += vv.size();
                }
                return n;
            }

            template <class T>
            inline double ChainLength(const std::vector<T> & v) {
                double len = 0;
                for (int i = 1; i < v.size(); i++) {
                    len += Distance(v[i - 1], v[i]);
                }
                return len;
            }

            template <class T>
            inline double ChainLength(const std::vector<std::vector<T>> & v) {
                double len = 0;
                for (const auto & vv : v) {
                    for (int i = 1; i < vv.size(); i++) {
                        len += Distance(vv[i - 1], vv[i]);
                    }
                }
                return len;
            }
            
            template <class T>
            inline std::pair<T, T> MakeOrderedPair(const T & a, const T & b) {
                return a < b ? std::make_pair(a, b) : std::make_pair(b, a);
            }

            bool AllAlong(const std::vector<std::vector<Vec3>> & pts, const Vec3 & from, const Vec3 & to, double angleThres) {
                auto n = normalize(from.cross(to));
                return std::all_of(pts.begin(), pts.end(), [&n, &angleThres](const std::vector<Vec3> & ps) {
                    return std::all_of(ps.begin(), ps.end(), [&n, &angleThres](const Vec3 & p) {
                        return abs(M_PI_2 - AngleBetweenDirections(n, p)) < angleThres;
                    });
                });
            }

            bool AllAlong(const std::vector<Vec3> & pts, const Vec3 & from, const Vec3 & to, double angleThres) {
                auto n = normalize(from.cross(to));
                return std::all_of(pts.begin(), pts.end(), [&n, angleThres](const Vec3 & p) {
                    return abs(M_PI_2 - AngleBetweenDirections(n, p)) < angleThres;
                });
            }

            struct ComparePixel {
                inline bool operator ()(const Pixel & a, const Pixel & b) const {
                    if (a.x != b.x)
                        return a.x < b.x;
                    return a.y < b.y;
                }
            };

        }

        PIGraph BuildPIGraph(const PanoramicView & view, const std::vector<Vec3> & vps, int verticalVPId,
            const Imagei & segs, const std::vector<Classified<Line3>> & lines, 
            double bndPieceSplitAngleThres, 
            double bndPieceClassifyAngleThres,
            double bndPieceBoundToLineAngleThres,
            double intersectionAngleThreshold,
            double incidenceAngleAlongDirectionThreshold,
            double incidenceAngleVerticalDirectionThreshold) {

            PIGraph mg;
            mg.view = view;
            assert(vps.size() >= 3);
            mg.vps = std::vector<Vec3>(vps.begin(), vps.begin() + 3);
            mg.verticalVPId = verticalVPId;

            mg.segs = segs;
            assert(IsDenseSegmentation(segs));
            int width = segs.cols;
            int height = segs.rows;
            int nsegs = MinMaxValOfImage(segs).second + 1;

            // init segs
            mg.nsegs = nsegs;
            mg.seg2area.resize(nsegs);
            mg.seg2bnds.resize(nsegs);
            mg.seg2center.resize(nsegs);
            mg.seg2control.resize(nsegs);
            mg.seg2linePieces.resize(nsegs);
            mg.seg2plane.resize(nsegs);
            mg.seg2contours.resize(nsegs);

            double fullArea = 0.0;
            for (auto it = segs.begin(); it != segs.end(); ++it) {
                double weight = cos((it.pos().y - height / 2.0) / height * M_PI);
                mg.seg2area[*it] += weight;
                fullArea += weight;
            }
            for (int i = 0; i < nsegs; i++) {
                mg.seg2area[i] /= fullArea;
                auto & control = mg.seg2control[i];
                control.orientationClaz = control.orientationNotClaz = -1;
                control.used = true;
                auto & center = mg.seg2center[i];
                mg.seg2plane[i].normal = - center;
                mg.seg2plane[i].anchor = center;
            }
            for (int i = 0; i < nsegs; i++) {
                Image regionMask = (segs == i);

                // find contour of the region
                std::vector<std::vector<Pixel>> contours;
                cv::findContours(regionMask, contours, CV_RETR_EXTERNAL, CV_CHAIN_APPROX_SIMPLE); // CV_RETR_EXTERNAL: get only the outer contours
                if (contours.empty()) {
                    continue;
                }

                Vec3 centerDirection(0, 0, 0);
                std::vector<Vec3> directions;
                directions.reserve(ElementsNum(contours));
                for (auto & cs : contours) {
                    for (auto & c : cs) {
                        directions.push_back(normalize(view.camera.toSpace(c)));
                        centerDirection += directions.back();
                    }
                }
                centerDirection /= norm(centerDirection);
                // get max angle distance from center direction
                double radiusAngle = 0.0;
                for (auto & d : directions) {
                    double a = AngleBetweenDirections(centerDirection, d);
                    if (radiusAngle < a) {
                        radiusAngle = a;
                    }
                }

                // perform a more precise sample !
                int newSampleSize = view.camera.focal() * radiusAngle * 2 + 2;
                PartialPanoramicCamera sCam(newSampleSize, newSampleSize, view.camera.focal(), view.camera.eye(), centerDirection,
                    ProposeXYDirectionsFromZDirection(centerDirection).second);
                Imagei sampledSegmentedRegions = MakeCameraSampler(sCam, view.camera)(segs);

                // collect better contours
                contours.clear();
                regionMask = (sampledSegmentedRegions == i);
                cv::findContours(regionMask, contours, CV_RETR_EXTERNAL, CV_CHAIN_APPROX_SIMPLE); // CV_RETR_EXTERNAL: get only the outer contours
                std::sort(contours.begin(), contours.end(),
                    [](const std::vector<Pixel> & ca, const std::vector<Pixel> & cb) {return ca.size() > cb.size(); });

                auto iter = std::find_if(contours.begin(), contours.end(), [](const std::vector<Pixel> & c) {return c.size() <= 2; });
                contours.erase(iter, contours.end());

                mg.seg2contours[i].resize(contours.size());
                mg.seg2center[i] = Origin();
                for (int j = 0; j < contours.size(); j++) {
                    auto & cs = mg.seg2contours[i][j];
                    cs.reserve(contours[j].size());
                    for (auto & p : contours[j]) {
                        cs.push_back(normalize(sCam.toSpace(p)));
                        mg.seg2center[i] += cs.back();
                    }
                }
                mg.seg2center[i] /= norm(mg.seg2center[i]);
            }


            // init lines
            mg.lines = lines;
            int nlines = lines.size();
            mg.line2linePieces.resize(nlines);
            mg.line2lineRelations.resize(nlines);
            mg.line2reconstructed.resize(nlines);
            for (int i = 0; i < nlines; i++) {
                mg.line2reconstructed[i] = mg.lines[i].component;
            }

            std::map<std::set<int>, std::vector<int>> segs2juncs;
            std::vector<std::vector<int>> seg2juncs(nsegs);
            std::map<std::pair<int, int>, std::set<Pixel, ComparePixel>> segpair2pixels;

            std::map<Pixel, std::set<int>, ComparePixel> pixel2segs;
            std::map<Pixel, int, ComparePixel> pixel2junc;

            std::vector<Pixel> juncPositions;
            std::vector<std::vector<int>> junc2segs;

            std::cout << "recording pixels" << std::endl;
            for (auto it = segs.begin(); it != segs.end(); ++it) {
                auto p = it.pos();

                // seg ids related
                std::set<int> idset = {
                    segs(p),
                    segs(Pixel((p.x + 1) % segs.cols, p.y))
                };
                if (p.y < height - 1) { // note that the top/bottom borders cannot be crossed!
                    idset.insert(segs(Pixel(p.x, p.y + 1)));
                    idset.insert(segs(Pixel((p.x + 1) % segs.cols, p.y + 1)));
                }

                if (idset.size() <= 1) {
                    continue;
                }

                // meet a bnd or a junc (and bnd) pixel
                pixel2segs[p] = std::move(idset);
                auto & relatedsegs = pixel2segs[p];

                // register this pixel as a bnd candidate for bnd of related segids
                for (auto ii = relatedsegs.begin(); ii != relatedsegs.end(); ++ii) {
                    for (auto jj = std::next(ii); jj != relatedsegs.end(); ++jj) {
                        segpair2pixels[MakeOrderedPair(*ii, *jj)].insert(p);
                    }
                }

                // meet a junc
                if (relatedsegs.size() >= 3) {
                    // create a new junc
                    juncPositions.push_back(p);
                    int newjuncid = juncPositions.size() - 1;
                    pixel2junc[p] = newjuncid;
                    segs2juncs[relatedsegs].push_back(newjuncid);
                    for (int segid : relatedsegs) {
                        seg2juncs[segid].push_back(newjuncid);
                    }
                    junc2segs.emplace_back(relatedsegs.begin(), relatedsegs.end());
                }
            }

            // now we have juncs
            mg.junc2positions.resize(juncPositions.size());
            for (int i = 0; i < juncPositions.size(); i++) {
                mg.junc2positions[i] = normalize(view.camera.toSpace(juncPositions[i]));
            }
            mg.junc2bnds.resize(juncPositions.size());


            std::cout << "connecting boundaries" << std::endl;


            // connect different junctions using allbndpixels
            // and thus generate seperated bnds
            std::vector<std::vector<Pixel>> bndPixels;
            for (int i = 0; i < juncPositions.size(); i++) {
                auto & juncpos = juncPositions[i];
                auto & relatedSegIds = junc2segs[i];

                for (int ii = 0; ii < relatedSegIds.size(); ii++) {
                    for (int jj = ii + 1; jj < relatedSegIds.size(); jj++) {
                        int segi = relatedSegIds[ii];
                        int segj = relatedSegIds[jj];

                        const auto & pixelsForThisSegPair = segpair2pixels.at(MakeOrderedPair(segi, segj));
                        std::vector<Pixel> pixelsForThisBnd;
                        int saySegIIsOnLeft = 0, saySegIIsOnRight = 0;

                        std::set<Pixel, ComparePixel> visitedPixels;

                        // use BFS
                        std::queue<Pixel> Q;
                        Q.push(juncpos);
                        visitedPixels.insert(juncpos);

                        while (!Q.empty()) {
                            auto curp = Q.front();
                            pixelsForThisBnd.push_back(curp);

                            // find another junc!
                            if (Contains(pixel2junc, curp) && i < pixel2junc.at(curp) && pixelsForThisBnd.size() > 1) {
                                int tojuncid = pixel2junc.at(curp);

                                // make a new bnd!
                                bndPixels.push_back(std::move(pixelsForThisBnd));
                                int newbndid = bndPixels.size() - 1;
                                mg.bnd2juncs.emplace_back(i, tojuncid);

                                if (saySegIIsOnLeft < saySegIIsOnRight) {
                                    mg.bnd2segs.emplace_back(segj, segi); // left, right
                                } else {
                                    mg.bnd2segs.emplace_back(segi, segj);
                                }
                                saySegIIsOnLeft = saySegIIsOnRight = 0;

                                mg.junc2bnds[i].push_back(newbndid);
                                mg.junc2bnds[tojuncid].push_back(newbndid);

                                mg.seg2bnds[segi].push_back(newbndid);
                                mg.seg2bnds[segj].push_back(newbndid);

                                break;
                            }

                            Q.pop();

                            // * - * - *
                            // | d | a |
                            // * -[*]- *
                            // | c | b |
                            // * - * - *
                            static const int dxs[] = { 1, 0, -1, 0 };
                            static const int dys[] = { 0, 1, 0, -1 };

                            static const int leftdxs[] = { 1, 1, 0, 0 };
                            static const int leftdys[] = { 0, 1, 1, 0 };

                            static const int rightdxs[] = { 1, 0, 0, 1 };
                            static const int rightdys[] = { 1, 1, 0, 0 };

                            for (int k = 0; k < 4; k++) {
                                auto nextp = curp + Pixel(dxs[k], dys[k]);
                                nextp.x = (nextp.x + width) % width;

                                if (nextp.y >= height || nextp.y < 0) { // note that the top/bottom borders cannot be crossed!
                                    continue;
                                }
                                if (!Contains(pixelsForThisSegPair, nextp)) {
                                    continue;
                                }
                                if (Contains(visitedPixels, nextp)) {
                                    continue;
                                }

                                auto rightp = curp + Pixel(rightdxs[k], rightdys[k]);
                                rightp.x = (rightp.x + width) % width;
                                auto leftp = curp + Pixel(leftdxs[k], leftdys[k]);
                                leftp.x = (leftp.x + width) % width;
                                if (Contains(segs, rightp) && Contains(segs, leftp)) {
                                    if (segs(rightp) == segi && segs(leftp) == segj) {
                                        saySegIIsOnRight++;
                                    } else if (segs(rightp) == segj && segs(leftp) == segi) {
                                        saySegIIsOnLeft++;
                                    } else {
                                        continue;
                                    }
                                }
                                Q.push(nextp);
                                visitedPixels.insert(nextp);
                            }
                        }
                    }
                }            
            }

            // split bnds into pieces
            std::cout << "splitting boundaries" << std::endl;
            std::vector<std::vector<Vec3>> bnd2SmoothedPostions(bndPixels.size());
            for (int i = 0; i < bndPixels.size(); i++) {
                assert(bndPixels[i].size() >= 2);
                for (int j = 0; j < bndPixels[i].size(); j+=2) {
                    bnd2SmoothedPostions[i].push_back(normalize(view.camera.toSpace(bndPixels[i][j])));
                }
                if (bnd2SmoothedPostions[i].size() == 1) {
                    bnd2SmoothedPostions[i].push_back(normalize(view.camera.toSpace(bndPixels[i].back())));
                }
            }
            mg.bnd2bndPieces.resize(mg.bnd2segs.size());
            for (int i = 0; i < bnd2SmoothedPostions.size(); i++) {
                const auto & dirs = bnd2SmoothedPostions[i];
                assert(dirs.size() > 0);

                std::vector<Vec3> curPiece = { dirs.front() };
                for (int j = 1; j <= dirs.size(); j++) {
                    if (j < dirs.size() && AllAlong(curPiece, curPiece.front(), dirs[j], bndPieceSplitAngleThres)) {
                        curPiece.push_back(dirs[j]);
                    } else {
                        if (curPiece.size() >= 2) {
                            mg.bndPiece2dirs.push_back(std::move(curPiece));
                            mg.bndPiece2bnd.push_back(i);
                            int bndPieceId = mg.bndPiece2dirs.size() - 1;
                            mg.bnd2bndPieces[i].push_back(bndPieceId);
                            curPiece.clear();
                        }                        
                        if (j < dirs.size()) {
                            curPiece.push_back(dirs[j]);
                        }
                    }
                }
            }

            // bndPieces properties
            // classes
            std::cout << "classifying boundary pieces" << std::endl;
            mg.bndPiece2classes.resize(mg.bndPiece2dirs.size(), -1);
            for (int i = 0; i < mg.bndPiece2dirs.size(); i++) {
                auto & piece = mg.bndPiece2dirs[i];
                assert(piece.size() > 1);
                Vec3 center = piece[piece.size() / 2];
                mg.bndPiece2classes[i] = -1;
                assert(vps.size() >= 3);
                for (int j = 0; j < 3; j++) {
                    if (AllAlong(piece, center, vps[j], bndPieceClassifyAngleThres)) {
                        mg.bndPiece2classes[i] = j;
                        break;
                    }
                }
            }
            // lengths
            mg.bndPiece2length.resize(mg.bndPiece2dirs.size(), 0);
            for (int i = 0; i < mg.bndPiece2dirs.size(); i++) {
                auto & piece = mg.bndPiece2dirs[i];
                for (int j = 1; j < piece.size(); j++) {
                    mg.bndPiece2length[i] += 
                        AngleBetweenDirections(piece[j-1], piece[j]);
                }
            }
            mg.bndPiece2linePieces.resize(mg.bndPiece2dirs.size());
            mg.bndPiece2occlusion.resize(mg.bndPiece2dirs.size(), OcclusionRelation::Unknown);


            // register bndPiece dirs in RTree
            RTreeMap<Vec3, int> bndPieceRTree;
            for (int i = 0; i < mg.bndPiece2dirs.size(); i++) {
                for (auto & d : mg.bndPiece2dirs[i]) {
                    bndPieceRTree.emplace(normalize(d), i);
                }
            }
            RTreeMap<Vec3, int> lineRTree;
            static const double lineSampleAngle = DegreesToRadians(0.5);
            std::vector<std::vector<Vec3>> lineSamples(mg.lines.size());
            for (int i = 0; i < mg.lines.size(); i++) {
                auto & line = mg.lines[i];
                double angle = AngleBetweenDirections(line.component.first, line.component.second);
                for (double a = 0; a <= angle; a += lineSampleAngle) {
                    Vec3 sample = normalize(RotateDirection(line.component.first, line.component.second, a));
                    lineSamples[i].push_back(sample);
                    lineRTree.emplace(sample, i);
                }
            }           

            // split lines to linePieces
            for (int i = 0; i < mg.lines.size(); i++) {
                Vec3 lineRotNormal = normalize(mg.lines[i].component.first.cross(mg.lines[i].component.second));
                auto & samples = lineSamples[i];
                int lastDetectedBndPiece = -1;
                int lastDetectedSeg = -1;
                std::vector<Vec3> collectedSamples;
                for (int j = 0; j <= samples.size(); j++) {
                    double minAngleDist = bndPieceBoundToLineAngleThres;
                    int nearestBndPiece = -1;
                    int nearestSeg = -1;
                    if (j < samples.size()) {
                        auto & d = samples[j];
                        bndPieceRTree.search(BoundingBox(d).expand(bndPieceBoundToLineAngleThres * 3),
                            [&d, &minAngleDist, &nearestBndPiece](const std::pair<Vec3, int> & bndPieceD) {
                            double angleDist = AngleBetweenDirections(bndPieceD.first, d);
                            if (angleDist < minAngleDist) {
                                minAngleDist = angleDist;
                                nearestBndPiece = bndPieceD.second;
                            }
                            return true;
                        });
                        nearestSeg = segs(ToPixel(view.camera.toScreen(d)));
                    }
                    bool neighborChanged = (nearestBndPiece != lastDetectedBndPiece || 
                        (nearestBndPiece == -1 && lastDetectedSeg != nearestSeg) ||
                        j == samples.size() - 1) && !(lastDetectedBndPiece == -1 && lastDetectedSeg == -1);
                    if (neighborChanged) {
                        double len = AngleBetweenDirections(collectedSamples.front(), collectedSamples.back());
                        if (collectedSamples.size() >= 2) {
                            if (lastDetectedBndPiece != -1) { // line piece bound to bnd piece
                                mg.linePiece2bndPiece.push_back(lastDetectedBndPiece);
                                int linePieceId = mg.linePiece2bndPiece.size() - 1;
                                mg.line2linePieces[i].push_back(linePieceId);
                                mg.linePiece2line.push_back(i);
                                mg.linePiece2samples.push_back(std::move(collectedSamples));
                                mg.linePiece2length.push_back(len);
                                mg.linePiece2seg.push_back(-1);
                                mg.linePiece2attachment.push_back(AttachmentRelation::Unknown);
                                auto & bndPieceDirs = mg.bndPiece2dirs[lastDetectedBndPiece];
                                assert(bndPieceDirs.size() > 1);
                                Vec3 bndRotNormal = bndPieceDirs.front().cross(bndPieceDirs.back());
                                mg.linePiece2bndPieceInSameDirection.push_back(lineRotNormal.dot(bndRotNormal) > 0);
                                mg.bndPiece2linePieces[lastDetectedBndPiece].push_back(linePieceId);
                            } else { // line piece bound to seg
                                mg.linePiece2bndPiece.push_back(-1);
                                int linePieceId = mg.linePiece2bndPiece.size() - 1;
                                mg.line2linePieces[i].push_back(linePieceId);
                                mg.linePiece2line.push_back(i);
                                mg.linePiece2samples.push_back(std::move(collectedSamples));
                                mg.linePiece2length.push_back(len);
                                mg.linePiece2seg.push_back(lastDetectedSeg);
                                mg.seg2linePieces[lastDetectedSeg].push_back(linePieceId);
                                mg.linePiece2attachment.push_back(AttachmentRelation::Unknown);
                                mg.linePiece2bndPieceInSameDirection.push_back(true);
                            }
                        }
                        collectedSamples.clear();
                    }
                    if (j < samples.size()) {
                        collectedSamples.push_back(samples[j]);
                    }

                    lastDetectedBndPiece = nearestBndPiece;
                    lastDetectedSeg = nearestSeg;
                }
            }

            // build line relations
            for (int i = 0; i < mg.lines.size(); i++) {                
                for (int j = i + 1; j < mg.lines.size(); j++) {
                    auto & linei = mg.lines[i].component;
                    int clazi = mg.lines[i].claz;
                    Vec3 ni = normalize(linei.first.cross(linei.second));
                    auto & linej = mg.lines[j].component;
                    int clazj = mg.lines[j].claz;
                    Vec3 nj = normalize(linej.first.cross(linej.second));

                    auto nearest = DistanceBetweenTwoLines(linei, linej);
                    double d = AngleBetweenDirections(nearest.second.first.position, nearest.second.second.position);

                    if (clazi == clazj && clazi >= 0) { // incidences for classified lines
                        auto conCenter = normalize(nearest.second.first.position + nearest.second.second.position);
                        auto conDir = (nearest.second.first.position - nearest.second.second.position);

                        auto & vp = vps[clazi];

                        if (AngleBetweenDirections(vp, conCenter) < intersectionAngleThreshold)
                            continue;

                        if (d < incidenceAngleAlongDirectionThreshold &&
                            AngleBetweenUndirectedVectors(ni, nj) < incidenceAngleVerticalDirectionThreshold) {
                            //LineRelationData lrd;
                            //lrd.type = LineRelationData::Type::Incidence;
                            //lrd.normalizedRelationCenter = conCenter;
                            //lrd.junctionWeight = 5.0;

                            if (HasValue(conCenter, IsInfOrNaN<double>))
                                continue;

                            mg.lineRelation2lines.emplace_back(i, j);
                            int lineRelationId = mg.lineRelation2lines.size() - 1;
                            mg.line2lineRelations[i].push_back(lineRelationId);
                            mg.line2lineRelations[j].push_back(lineRelationId);
                            mg.lineRelation2anchor.push_back(conCenter);
                            mg.lineRelation2weight.push_back(5.0);
                            mg.lineRelation2IsIncidence.push_back(true);

                            //mg.addConstraint(std::move(lrd), lhs[i], lhs[j]);
                        }

                    } else if (clazi != clazj && clazi >= 0 && clazj >= 0) { // intersections for classified lines
                        if (d < intersectionAngleThreshold) {
                            auto conCenter = normalize(ni.cross(nj));

                            if (DistanceFromPointToLine(conCenter, linei).first > intersectionAngleThreshold * 4 ||
                                DistanceFromPointToLine(conCenter, linej).first > intersectionAngleThreshold * 4)
                                continue;

                            //LineRelationData lrd;
                            //lrd.type = LineRelationData::Type::Intersection;
                            //lrd.normalizedRelationCenter = conCenter;
                            //lrd.junctionWeight = 3.0;

                            if (HasValue(conCenter, IsInfOrNaN<double>))
                                continue;

                            mg.lineRelation2lines.emplace_back(i, j);
                            int lineRelationId = mg.lineRelation2lines.size() - 1;
                            mg.line2lineRelations[i].push_back(lineRelationId);
                            mg.line2lineRelations[j].push_back(lineRelationId);
                            mg.lineRelation2anchor.push_back(conCenter);
                            mg.lineRelation2weight.push_back(3.0); 
                            mg.lineRelation2IsIncidence.push_back(false);

                            //mg.addConstraint(std::move(lrd), lhs[i], lhs[j]);
                        }
                    }
                }
            }

            // line relation weights
            // compute junction weights
            static const double angleThreshold = M_PI / 32;
            static const double sigma = 0.1;

            enum LineVotingDirection : int {
                TowardsVanishingPoint = 0,
                TowardsOppositeOfVanishingPoint = 1
            };
            enum class JunctionType : int {
                L, T, Y, W, X
            };
            for (int i = 0; i < mg.lineRelation2lines.size(); i++) {
                //std::cout << lr.topo.hd.id << std::endl;
                if (mg.lineRelation2IsIncidence[i]) {
                    mg.lineRelation2weight[i] = IncidenceJunctionWeight(false);
                } else {
                    Mat<float, 3, 2> votingData;
                    std::fill(std::begin(votingData), std::end(votingData), 0);

                    for (int lineid = 0; lineid < mg.lines.size(); lineid ++) {
                        auto & line = mg.lines[lineid].component;
                        int claz = mg.lines[lineid].claz;
                        if (claz == -1 || claz >= 3)
                            continue;

                        auto & vp = vps[claz];
                        Vec3 center = normalize(line.center());

                        Vec3 center2vp = normalize(center.cross(vp));
                        Vec3 center2pos = normalize(center.cross(mg.lineRelation2anchor[i]));

                        double angle = AngleBetweenUndirectedVectors(center2vp, center2pos);
                        double angleSmall = angle > M_PI_2 ? (M_PI - angle) : angle;
                        if (IsInfOrNaN(angleSmall))
                            continue;

                        assert(angleSmall >= 0 && angleSmall <= M_PI_2);

                        double angleScore =
                            exp(-(angleSmall / angleThreshold) * (angleSmall / angleThreshold) / sigma / sigma / 2);

                        auto proj = ProjectionOfPointOnLine(mg.lineRelation2anchor[i], line);
                        double projRatio = BoundBetween(proj.ratio, 0.0, 1.0);

                        Vec3 lined = line.first.cross(line.second);
                        double lineSpanAngle = AngleBetweenDirections(line.first, line.second);
                        if (AngleBetweenDirections(center2vp, lined) < M_PI_2) { // first-second-vp
                            votingData(claz, TowardsVanishingPoint) += angleScore * lineSpanAngle * (1 - projRatio);
                            votingData(claz, TowardsOppositeOfVanishingPoint) += angleScore * lineSpanAngle * projRatio;
                        } else { // vp-first-second
                            votingData(claz, TowardsOppositeOfVanishingPoint) += angleScore * lineSpanAngle * (1 - projRatio);
                            votingData(claz, TowardsVanishingPoint) += angleScore * lineSpanAngle * projRatio;
                        }
                    }
                    mg.lineRelation2weight[i] = ComputeIntersectionJunctionWeightWithLinesVotes(
                        votingData);
                }
            }
        
            return mg;
        }






        View<PartialPanoramicCamera, Imageub> PerfectSegMaskView(const PIGraph & mg, int seg, double focal) {
            auto & contours = mg.seg2contours[seg];
            double radiusAngle = 0.0;
            for (auto & cs : contours) {
                for (auto & c : cs) {
                    double angle = AngleBetweenDirections(mg.seg2center[seg], c);
                    if (angle > radiusAngle) {
                        radiusAngle = angle;
                    }
                }
            }
            int ppcSize = std::ceil(2 * radiusAngle * focal);
            Vec3 x;
            std::tie(x, std::ignore) = ProposeXYDirectionsFromZDirection(mg.seg2center[seg]);
            PartialPanoramicCamera ppc(ppcSize, ppcSize, focal, Point3(0, 0, 0), mg.seg2center[seg], x);
            Imageub mask = Imageub::zeros(ppc.screenSize());

            // project contours to ppc
            std::vector<std::vector<Point2i>> contourProjs(contours.size());
            for (int k = 0; k < contours.size(); k++) {
                auto & contourProj = contourProjs[k];
                contourProj.reserve(contours[k].size());
                for (auto & d : contours[k]) {
                    contourProj.push_back(ecast<int>(ppc.toScreen(d)));
                }
            }
            cv::fillPoly(mask, contourProjs, (uint8_t)1);
            return View<PartialPanoramicCamera, Imageub>{mask, ppc};
        }



    }

}