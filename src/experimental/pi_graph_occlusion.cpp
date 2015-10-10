#include "../core/utility.hpp"
#include "../core/containers.hpp"
#include "../ml/factor_graph.hpp"

#include "../gui/canvas.hpp"

#include "pi_graph_occlusion.hpp"
#include "pi_graph_vis.hpp"


namespace pano {
    namespace experimental {

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

        struct SegLabel {
            int orientationClaz, orientationNotClaz;
            inline int dof() const {
                if (orientationClaz != -1) {
                    return 1;
                }
                if (orientationNotClaz != -1) {
                    return 2;
                }
                return 3;
            }
            inline bool operator == (const SegLabel & sl) const {
                return std::tie(orientationClaz, orientationNotClaz) == std::tie(sl.orientationClaz, sl.orientationNotClaz);
            }
            inline bool operator < (const SegLabel & sl) const {
                return std::tie(orientationClaz, orientationNotClaz) < std::tie(sl.orientationClaz, sl.orientationNotClaz);
            }
        };

        void DetectOcclusions(PIGraph & mg) {

            ml::FactorGraph fg;

            std::cout << "building factor graph" << std::endl;

            //// vars

            // segs
            std::vector<ml::FactorGraph::VarHandle> seg2vh(mg.nsegs);
            std::vector<std::vector<SegLabel>> seg2allowedLabels(mg.nsegs);
            std::map<ml::FactorGraph::VarHandle, int> vh2seg;
            for (int i = 0; i < mg.nsegs; i++) {
                auto & c = mg.seg2control[i];
                if (!c.used) continue;

                std::vector<SegLabel> allowedLabels;
                double c_i = mg.seg2areaRatio[i] * 4 * M_PI;
                
                if (c.orientationClaz != -1) {
                    allowedLabels = { {c.orientationClaz, -1} };                   
                } else if (c.orientationNotClaz != -1) {
                    allowedLabels = { {-1, c.orientationNotClaz} };
                } else {
                    allowedLabels.push_back({ -1, -1 });
                }

                auto vh = fg.addVar(fg.addVarCategory(allowedLabels.size(), c_i));
                seg2vh[i] = vh;
                seg2allowedLabels[i] = std::move(allowedLabels);
                vh2seg[vh] = i;
            }


            // bndPieces
            enum BndPieceLabel : int { 
                //Coplanar = (int)SegRelation::Coplanar,
                Connected = (int)SegRelation::Connected, 
                LeftIsFront = (int)SegRelation::LeftIsFront, 
                RightIsFront = (int)SegRelation::RightIsFront
            };
            std::vector<ml::FactorGraph::VarHandle> bndPiece2vh(mg.nbndPieces());
            std::vector<std::vector<BndPieceLabel>> bndPiece2allowedLabels(mg.nbndPieces());
            std::map<ml::FactorGraph::VarHandle, int> vh2bndPiece;
            for (int i = 0; i < mg.nbndPieces(); i++) {
                std::vector<BndPieceLabel> allowedLabels;
                double c_i = mg.bndPiece2length[i];

                if (mg.bndPiece2segRelation[i] != SegRelation::Unknown) {
                    allowedLabels = { (BndPieceLabel)mg.bndPiece2segRelation[i] };
                } else if (mg.bndPiece2length[i] < DegreesToRadians(1)) { // always disconnect tiny bndPieces
                    allowedLabels = { LeftIsFront, RightIsFront };
                } else {
                    allowedLabels = { Connected, LeftIsFront, RightIsFront };
                }
                int vc = fg.addVarCategory(allowedLabels.size(), c_i);
                auto vh = fg.addVar(vc);
                bndPiece2vh[i] = vh;
                bndPiece2allowedLabels[i] = std::move(allowedLabels);
                vh2bndPiece[vh] = i;
            }

            // linePieces
            enum LinePieceLabel {
                Attached = (int)SegLineRelation::Attached, 
                Detached = (int)SegLineRelation::Detached
            };
            std::vector<ml::FactorGraph::VarHandle> linePiece2vh(mg.nlinePieces());
            std::vector<std::vector<LinePieceLabel>> linePiece2allowedLabels(mg.nlinePieces());
            std::map<ml::FactorGraph::VarHandle, int> vh2linePiece;
            for (int i = 0; i < mg.nlinePieces(); i++) {
                std::vector<LinePieceLabel> allowedLabels;
                if (mg.linePiece2segLineRelation[i] != SegLineRelation::Unknown) {
                    allowedLabels = { (LinePieceLabel)mg.linePiece2segLineRelation[i] };
                } else if(mg.linePiece2length[i] < DegreesToRadians(1)) {
                    allowedLabels = { Detached };
                } else {
                    allowedLabels = { Attached, Detached };
                }
                int vc = fg.addVarCategory(allowedLabels.size(), mg.linePiece2length[i]);
                auto vh = fg.addVar(vc);
                linePiece2vh[i] = vh;
                linePiece2allowedLabels[i] = std::move(allowedLabels);
                vh2linePiece[vh] = i;
            }


            //// factors

            // seg
            for (int seg = 0; seg < mg.nsegs; seg++) {
                if (seg2vh[seg].invalid()) {
                    continue;
                }
                int fc = fg.addFactorCategory([seg, &seg2allowedLabels, &mg](const int * varlabels, size_t nvar,
                    ml::FactorGraph::FactorCategoryId fcid, void * givenData)->double {
                    assert(nvar == 1);
                    return 0.0;
                }, 1.0);
                fg.addFactor({ seg2vh[seg] }, fc);
            }

            // bndPiece
            for (int bndPiece = 0; bndPiece < mg.nbndPieces(); bndPiece++) {
                if (bndPiece2vh[bndPiece].invalid()) {
                    continue;
                }
                int fc = fg.addFactorCategory([bndPiece, &bndPiece2allowedLabels, &mg](const int * varlabels, size_t nvar,
                    ml::FactorGraph::FactorCategoryId fcid, void * givenData)->double {
                    assert(nvar == 1);
                    BndPieceLabel label = bndPiece2allowedLabels[bndPiece][varlabels[0]];
                    double len = mg.bndPiece2length[bndPiece];
                    switch (label){
                    case Connected:
                        return 0.0;
                    case LeftIsFront:
                    case RightIsFront:
                    default:
                        return len;
                    }
                }, 1.0);
                fg.addFactor({ bndPiece2vh[bndPiece] }, fc);
            }

            // linePiece
            for (int linePiece = 0; linePiece < mg.nlinePieces(); linePiece++) {
                if (linePiece2vh[linePiece].invalid()) {
                    continue;
                }
                int fc = fg.addFactorCategory([linePiece, &linePiece2allowedLabels, &mg](const int * varlabels, size_t nvar,
                    ml::FactorGraph::FactorCategoryId fcid, void * givenData)->double {
                    assert(nvar == 1);
                    LinePieceLabel label = linePiece2allowedLabels[linePiece][varlabels[0]];
                    double len = mg.linePiece2length[linePiece];
                    switch (label){
                    case Attached:
                        return 0.0;
                    case Detached:
                    default:
                        return len;
                    }
                }, 1.0);
                fg.addFactor({ linePiece2vh[linePiece] }, fc);
            }






            // orientation consistency term between seg and line
            const auto costBetweenSegAndLine = [&mg](const SegLabel & seglabel, const LinePieceLabel & linePieceLabel, int seg, int linePiece) -> double {
                int lineClaz = mg.lines[mg.linePiece2line[linePiece]].claz;
                double len = mg.linePiece2length[linePiece];
                if (linePieceLabel == Detached) {
                    return 0.0;
                } else if (lineClaz == -1) {
                    return 0.0;
                } else if (seglabel.orientationNotClaz == lineClaz) {
                    return len * 20.0;
                } else {
                    return 0.0;
                }
            };
            // seg-linePiece-line, orientation consistency term between seg and line
            for (int linePiece = 0; linePiece < mg.nlinePieces(); linePiece++) {
                if (linePiece2vh[linePiece].invalid()) {
                    continue;
                }
                int seg = mg.linePiece2seg[linePiece];
                if (seg == -1 || seg2vh[seg].invalid()) {
                    continue;
                }
                int fc = fg.addFactorCategory([seg, linePiece, &seg2allowedLabels, &linePiece2allowedLabels, &mg, &costBetweenSegAndLine](
                    const int * varlabels, size_t nvar,
                    ml::FactorGraph::FactorCategoryId fcid, void * givenData)->double {
                    assert(nvar == 2);
                    const SegLabel & seglabel = seg2allowedLabels[seg][varlabels[0]];
                    LinePieceLabel lplabel = linePiece2allowedLabels[linePiece][varlabels[1]];
                    return costBetweenSegAndLine(seglabel, lplabel, seg, linePiece);
                }, 1.0);
                fg.addFactor({ seg2vh[seg], linePiece2vh[linePiece] }, fc);
            }

            // seg1/seg2-bndPiece-linePiece-line, orientation consistency term between seg and line
            for (int linePiece = 0; linePiece < mg.nlinePieces(); linePiece++) {
                if (linePiece2vh[linePiece].invalid()) {
                    continue;
                }
                int bndPiece = mg.linePiece2bndPiece[linePiece];
                if (bndPiece == -1 || bndPiece2vh[bndPiece].invalid()) {
                    continue;
                }
                int seg1 = -1, seg2 = -1;
                std::tie(seg1, seg2) = mg.bnd2segs[mg.bndPiece2bnd[bndPiece]];
                if (seg2vh[seg1].invalid() || seg2vh[seg2].invalid()) {
                    continue;
                }
                int line = mg.linePiece2line[linePiece];
                int fc = fg.addFactorCategory([seg1, seg2, bndPiece, linePiece, line, 
                    &seg2allowedLabels, &bndPiece2allowedLabels, &linePiece2allowedLabels, &mg, &costBetweenSegAndLine](
                    const int * varlabels, size_t nvar,
                    ml::FactorGraph::FactorCategoryId fcid, void * givenData)->double {
                    assert(nvar == 4); // 0-seg1 1-seg2 2-bndPiece 3-linePiece
                    int segs[] = { seg1, seg2 };
                    BndPieceLabel bndPieceLabel = bndPiece2allowedLabels[bndPiece][varlabels[2]];
                    bool segRelated[] = {
                        bndPieceLabel == LeftIsFront || bndPieceLabel == Connected,
                        bndPieceLabel == RightIsFront || bndPieceLabel == Connected
                    };
                    double cost = 0.0;
                    for (int i = 0; i < 2; i++) {
                        if (segRelated[i]) {
                            cost += costBetweenSegAndLine(seg2allowedLabels[segs[i]][varlabels[i]],
                                linePiece2allowedLabels[linePiece][varlabels[3]], segs[i], linePiece);
                        }
                    }
                    return cost;
                }, 1.0);
                fg.addFactor({ seg2vh[seg1], seg2vh[seg2], bndPiece2vh[bndPiece], linePiece2vh[linePiece] }, fc);
            }




            // {line-lineRelation-line} -> linePiece -> bndPiece, occlusion hints suggested by line's t-junction
            static const double angleThresForTJunction = DegreesToRadians(5);
            std::vector<std::pair<int, int>> tjunctionLines; // [hline, vline]
            std::vector<bool> tjunctionVLineLiesOnTheRightOfHLine;
            for (int lineRelation = 0; lineRelation < mg.nlineRelations(); lineRelation++) {
                if (mg.lineRelation2IsIncidence[lineRelation]) {
                    continue;
                }
                int line1 = -1, line2 = -1;
                std::tie(line1, line2) = mg.lineRelation2lines[lineRelation];
                // whether this forms a t-junction?
                if (mg.lines[line1].claz == mg.lines[line2].claz) {
                    continue;
                }
                auto l1 = normalize(mg.lines[line1].component);
                Vec3 n1 = normalize(l1.first.cross(l1.second));
                auto l2 = normalize(mg.lines[line2].component);
                Vec3 n2 = normalize(l2.first.cross(l2.second));
                if (AngleBetweenUndirectedVectors(n1, n2) < DegreesToRadians(3)) {
                    continue;
                }

                if (AngleBetweenDirections(l1.first, l1.second) < DegreesToRadians(3) ||
                    AngleBetweenDirections(l2.first, l2.second) < DegreesToRadians(3)) {
                    continue;
                }

                int hline = -1, vline = -1;
                double lambda1 = 0, lambda2 = 0;
                DistanceBetweenTwoLines(l1.ray(), l2.ray(), &lambda1, &lambda2);
                static const double shrinkRatio = 0.2;
                if (lambda1 < (1.0 - shrinkRatio) && lambda1 > shrinkRatio && (lambda2 < -shrinkRatio || lambda2 > (1.0 + shrinkRatio))) {
                    hline = line1; vline = line2;
                } else if (lambda2 < (1.0 - shrinkRatio) && lambda2 > shrinkRatio && (lambda1 < -shrinkRatio || lambda1 >(1.0 + shrinkRatio))) {
                    hline = line2; vline = line1;
                } else {
                    continue;
                }

                // Q: what if three lines form a T structure, and the center lies between the two hlines?
                // A: if the two hlines are close and are colinear, then they should have been merged after the line extraction step

                tjunctionLines.emplace_back(hline, vline);
                auto vlcenter = normalize(mg.lines[vline].component.center());
                auto & hl = mg.lines[hline].component;
                auto hlleft = hl.first.cross(hl.second);
                bool vlIsOnRightOfHl = (vlcenter - hl.first).dot(hlleft) < 0;
                tjunctionVLineLiesOnTheRightOfHLine.push_back(vlIsOnRightOfHl);

                bool leftOfLineIsFront = vlIsOnRightOfHl;
            }


            if (true) {
                auto pim = Print(mg, ConstantFunctor<gui::Color>(gui::White),
                    [&tjunctionLines, &mg](int lp) {                    
                    return gui::Transparent;
                }, ConstantFunctor<gui::Color>(gui::LightGray), 2, 3);
                const auto drawBndPiece = [&mg, &pim](int bp, const gui::Color & color) {
                    auto & e = mg.bndPiece2dirs[bp];
                    for (int i = 1; i < e.size(); i++) {
                        auto p1 = core::ToPixel(mg.view.camera.toScreen(e[i - 1]));
                        auto p2 = core::ToPixel(mg.view.camera.toScreen(e[i]));
                        if (Distance(p1, p2) >= pim.cols / 2) {
                            continue;
                        }
                        cv::clipLine(cv::Rect(0, 0, pim.cols, pim.rows), p1, p2);
                        cv::line(pim, p1, p2, (cv::Scalar)color / 255.0, 2);
                    }
                };
                const auto drawLine = [&mg, &pim](const Line3 & line, const gui::Color & color) {
                    double angle = AngleBetweenDirections(line.first, line.second);
                    std::vector<Pixel> ps;
                    for (double a = 0.0; a <= angle; a += 0.01) {
                        ps.push_back(ToPixel(mg.view.camera.toScreen(RotateDirection(line.first, line.second, a))));
                    }
                    for (int i = 1; i < ps.size(); i++) {
                        auto p1 = ps[i - 1];
                        auto p2 = ps[i];
                        if (Distance(p1, p2) >= pim.cols / 2) {
                            continue;
                        }
                        cv::clipLine(cv::Rect(0, 0, pim.cols, pim.rows), p1, p2);
                        cv::line(pim, p1, p2, (cv::Scalar)color / 255.0, 2);
                    }
                };
                auto randctable = gui::CreateRandomColorTableWithSize(tjunctionLines.size());
                int i = 0;
                for (auto & tj : tjunctionLines) {
                    drawLine(mg.lines[tj.first].component, randctable[i]);
                    //drawLine(mg.lines[tj.second].component, randctable[i]);
                }
                gui::AsCanvas(pim).show();
            }

            struct TJunctionHLineVotes {
                double rightCanNeverBeFrontVoteSum;
                double leftCanNeverBeFrontVoteSum;
                TJunctionHLineVotes() : rightCanNeverBeFrontVoteSum(0), leftCanNeverBeFrontVoteSum(0) {}
            };
            std::map<int, TJunctionHLineVotes> tjunctionHLine2votes;
            for (int i = 0; i < tjunctionLines.size(); i++) {
                auto & votes = tjunctionHLine2votes[tjunctionLines[i].first];
                auto & vl = mg.lines[tjunctionLines[i].second].component;
                double vllen = AngleBetweenDirections(vl.first, vl.second);
                if (tjunctionVLineLiesOnTheRightOfHLine[i]) {
                    votes.rightCanNeverBeFrontVoteSum += vllen;
                } else {
                    votes.leftCanNeverBeFrontVoteSum += vllen;
                }
            }
            for (auto & tjunctionHLineWithVotes : tjunctionHLine2votes) {
                int hline = tjunctionHLineWithVotes.first;
                auto & votes = tjunctionHLineWithVotes.second;
                for (int lpOnHLine : mg.line2linePieces[hline]) {
                    int bp = mg.linePiece2bndPiece[lpOnHLine];
                    if (bp == -1) {
                        continue;
                    }
                    double bpRightCanNeverBeFrontVoteSum = votes.rightCanNeverBeFrontVoteSum;
                    double bpLeftCanNeverBeFrontVoteSum = votes.leftCanNeverBeFrontVoteSum;
                    if (!mg.linePiece2bndPieceInSameDirection[lpOnHLine]) {
                        std::swap(bpRightCanNeverBeFrontVoteSum, bpLeftCanNeverBeFrontVoteSum);
                    }

                    // add a factor about this
                    int fc = fg.addFactorCategory([bp, bpRightCanNeverBeFrontVoteSum, bpLeftCanNeverBeFrontVoteSum, &mg, &bndPiece2allowedLabels](
                        const int * varlabels, size_t nvar,
                        ml::FactorGraph::FactorCategoryId fcid, void * givenData) -> double {
                        assert(nvar == 1);
                        BndPieceLabel label = bndPiece2allowedLabels[bp][varlabels[0]];
                        switch (label) {
                        case LeftIsFront:
                            return bpLeftCanNeverBeFrontVoteSum;
                        case RightIsFront:
                            return bpRightCanNeverBeFrontVoteSum;
                        case Connected:
                        default:
                            return 0.0;
                        }
                    }, 1.0);
                    fg.addFactor({ bndPiece2vh[bp] }, fc);
                }
            }




            // seg-bndPiece-seg, occlusion consistency term
            for (int bnd = 0; bnd < mg.nbnds(); bnd++) {
                int seg1 = mg.bnd2segs[bnd].first; // left
                int seg2 = mg.bnd2segs[bnd].second; // right
                auto sv1 = seg2vh[seg1], sv2 = seg2vh[seg2];            
                if (sv1.invalid() || sv2.invalid()) continue;

                auto & bndPieces = mg.bnd2bndPieces[bnd];
                for (int i = 0; i < bndPieces.size(); i++) {
                    int bndPiece = bndPieces[i];

                    int fc = fg.addFactorCategory([seg1, seg2, bndPiece, &mg, &seg2allowedLabels, &bndPiece2allowedLabels](
                        const int * varlabels, size_t nvar,
                        ml::FactorGraph::FactorCategoryId fcid, void * givenData)->double {
                        assert(nvar == 3);
                        
                        const auto & seglabel1 = seg2allowedLabels[seg1][varlabels[0]];
                        const auto & seglabel2 = seg2allowedLabels[seg2][varlabels[1]];
                        BndPieceLabel bndPieceLabel = bndPiece2allowedLabels[bndPiece][varlabels[2]];

                        int bndPieceClaz = mg.bndPiece2classes[bndPiece];
                        double bndPieceLen = mg.bndPiece2length[bndPiece];

                        if (seglabel1.dof() == 3 || seglabel2.dof() == 3) {
                            return 0.0;
                        }

                        if (seglabel1.dof() == 1 && seglabel2.dof() == 1 && 
                            seglabel1.orientationClaz != seglabel2.orientationClaz && 
                            bndPieceClaz == -1 && bndPieceLabel == Connected) {
                            return bndPieceLen * 0.1;
                        }
                        if (seglabel1.dof() == 1 && seglabel2.dof() == 1 &&
                            seglabel1.orientationClaz != seglabel2.orientationClaz &&
                            (bndPieceClaz == seglabel1.orientationClaz || bndPieceClaz == seglabel2.orientationClaz) && 
                            bndPieceLabel == Connected) {
                            return bndPieceLen * 20.0;
                        }


                        if (seglabel1.dof() == 1 && seglabel2.dof() == 1 &&
                            seglabel1.orientationClaz != seglabel2.orientationClaz &&
                            bndPieceClaz == seglabel1.orientationClaz && bndPieceLabel != RightIsFront) {
                            return bndPieceLen * 20.0;
                        }
                        if (seglabel1.dof() == 1 && seglabel2.dof() == 1 &&
                            seglabel1.orientationClaz != seglabel2.orientationClaz &&
                            bndPieceClaz == seglabel2.orientationClaz && bndPieceLabel != LeftIsFront) {
                            return bndPieceLen * 20.0;
                        }


                        if (seglabel1.dof() == 1 && seglabel2.dof() == 2 &&
                            seglabel1.orientationClaz == seglabel2.orientationNotClaz &&
                            bndPieceClaz == seglabel1.orientationClaz && bndPieceLabel != RightIsFront) {
                            return bndPieceLen * 20.0;
                        }
                        if (seglabel1.dof() == 2 && seglabel2.dof() == 1 &&
                            seglabel1.orientationNotClaz == seglabel2.orientationClaz &&
                            bndPieceClaz == seglabel2.orientationClaz && bndPieceLabel != LeftIsFront) {
                            return bndPieceLen * 20.0;
                        }                        

                        return 0.0;
                    }, 1.0);
                    fg.addFactor({ sv1, sv2, bndPiece2vh[bndPiece] }, fc);
                }    
            }

            //// seg supporting by bnd
            //for (int seg = 0; seg < mg.nsegs; seg++) {
            //    auto sv = seg2vh[seg];
            //    if (sv.invalid()) continue;

            //    int bndPiecesNum = 0;
            //    for (int bnd : mg.seg2bnds[seg]) {
            //        bndPiecesNum += mg.bnd2bndPieces[bnd].size();
            //    }

            //    for (int bnd : mg.seg2bnds[seg]) {
            //        for (int bndPiece : mg.bnd2bndPieces[bnd]) {
            //            int fc = fg.addFactorCategory([seg, bndPiece, bndPiecesNum, &mg, &seg2allowedLabels, &bndPiece2allowedLabels](
            //                const int * varlabels, size_t nvar,
            //                ml::FactorGraph::FactorCategoryId fcid, void * givenData) {
            //                assert(nvar == 2);
            //                auto seglabel = seg2allowedLabels[seg][varlabels[0]];
            //                BndPieceLabel bndPieceLabel = bndPiece2allowedLabels[bndPiece][varlabels[1]];
            //                
            //                int bndPieceClaz = mg.bndPiece2classes[bndPiece];

            //                if (bndPieceLabel == Connected) {
            //                    return 0.0;
            //                }
            //                return 1.0 / bndPiecesNum * 0.01;
            //                //return 0.0;
            //            }, 1.0);
            //            fg.addFactor({ sv, bndPiece2vh[bndPiece] }, fc);
            //        }
            //    }
            //}


            // seg hinted by bnd
            /*for (int seg = 0; seg < mg.nsegs; seg++) {
                auto sv = seg2vh[seg];
                if (sv.invalid()) continue;

                double surroundingBndPiecesLen = 0;
                for (int bnd : mg.seg2bnds[seg]) {
                    for (int bndPiece : mg.bnd2bndPieces[bnd]) {
                        surroundingBndPiecesLen += mg.bndPiece2length[bndPiece];
                    }
                }

                for (int bnd : mg.seg2bnds[seg]) {
                    for (int bndPiece : mg.bnd2bndPieces[bnd]) {
                        int fc = fg.addFactorCategory([seg, bndPiece, surroundingBndPiecesLen, &mg, &seg2allowedLabels, &bndPiece2allowedLabels](
                            const int * varlabels, size_t nvar,
                            ml::FactorGraph::FactorCategoryId fcid, void * givenData) {
                            assert(nvar == 2);
                            auto seglabel = seg2allowedLabels[seg][varlabels[0]];
                            BndPieceLabel bndPieceLabel = bndPiece2allowedLabels[bndPiece][varlabels[1]];

                            int bndPieceClaz = mg.bndPiece2classes[bndPiece];
                            double bndPieceLen = mg.bndPiece2length[bndPiece];

                            static const double bndPieceLenThres = DegreesToRadians(5);

                            if (seglabel.dof() == 1 && bndPieceClaz == seglabel.orientationClaz && bndPieceLen >= bndPieceLenThres) {
                                return bndPieceLen * 0.1;
                            }

                            return 0.0;

                        }, 1.0);
                        fg.addFactor({ sv, bndPiece2vh[bndPiece] }, fc);
                    }
                }
            }*/


            // line consistency
            for (int line = 0; line < mg.nlines(); line++) {
                auto & linePieces = mg.line2linePieces[line];
                if (linePieces.size() == 1) continue;
                // bndPiece - bndPiece
                for (int i = 0; i < linePieces.size(); i++) {
                    for (int j = i + 1; j < linePieces.size(); j++) {

                        int linePiece1 = linePieces[i];
                        int linePiece2 = linePieces[j];

                        int bndPiece1 = mg.linePiece2bndPiece[linePiece1];
                        int bndPiece2 = mg.linePiece2bndPiece[linePiece2];
                        if (bndPiece1 == -1 || bndPiece2 == -1) {
                            continue;
                        }

                        bool sameDirection = mg.linePiece2bndPieceInSameDirection[linePiece1] == mg.linePiece2bndPieceInSameDirection[linePiece2];

                        int fc = fg.addFactorCategory(
                            [linePiece1, linePiece2, bndPiece1, bndPiece2, sameDirection,
                            &mg, &bndPiece2allowedLabels](
                            const int * varlabels, size_t nvar, ml::FactorGraph::FactorCategoryId fcid, void * givenData) {
                            
                            assert(nvar == 2);
                            BndPieceLabel bndPieceLabel1 = bndPiece2allowedLabels[bndPiece1][varlabels[0]];
                            BndPieceLabel bndPieceLabel2 = bndPiece2allowedLabels[bndPiece2][varlabels[1]];
                            double bndPieceLen1 = mg.bndPiece2length[bndPiece1];
                            double bndPieceLen2 = mg.bndPiece2length[bndPiece2];
                            double weight = bndPieceLen1 + bndPieceLen2;

                            if (bndPieceLabel1 == Connected && bndPieceLabel2 == Connected) {
                                return 0.0;
                            }
                            if ((bndPieceLabel1 == Connected) || (bndPieceLabel2 == Connected)) {
                                return weight * 20;
                            }
                            if (sameDirection && bndPieceLabel1 != bndPieceLabel2) {
                                return weight * 20;
                            }
                            if (!sameDirection && bndPieceLabel1 == bndPieceLabel2) {
                                return weight * 20;
                            }

                            return 0.0;

                        }, 1.0);

                        fg.addFactor({ bndPiece2vh[bndPiece1], bndPiece2vh[bndPiece2] }, fc);

                    }
                }
                // bndPiece - linePiece
                for (int i = 0; i < linePieces.size(); i++) {
                    int bndPiece = mg.linePiece2bndPiece[linePieces[i]];
                    if (bndPiece == -1) {
                        continue;
                    }
                    for (int j = 0; j < linePieces.size(); j++) {
                        int linePiece = linePieces[j];
                        if (mg.linePiece2bndPiece[linePiece] != -1) {
                            continue;
                        }
                        int fc = fg.addFactorCategory(
                            [bndPiece, linePiece, &mg, &bndPiece2allowedLabels, &linePiece2allowedLabels](
                            const int * varlabels, size_t nvar, ml::FactorGraph::FactorCategoryId fcid, void * givenData) {

                            assert(nvar == 2);
                            BndPieceLabel bndPieceLabel = bndPiece2allowedLabels[bndPiece][varlabels[0]];
                            LinePieceLabel linePieceLabel = linePiece2allowedLabels[linePiece][varlabels[1]];
                            double bndPieceLen = mg.bndPiece2length[bndPiece];
                            double linePieceLen = mg.linePiece2length[linePiece];
                            double weight = bndPieceLen + linePieceLen;

                            if (bndPieceLabel != Connected) {
                                if (linePieceLabel == Attached) {
                                    return weight * 10;
                                } else {
                                    return 0.0;
                                }
                            }
                            return 0.0;                           
                        }, 1.0);
                        fg.addFactor({ bndPiece2vh[bndPiece], linePiece2vh[linePiece] }, fc);
                    }
                }
            }

            // bnd consistency
            for (int bnd = 0; bnd < mg.nbnds(); bnd++) {
                auto & bndPieces = mg.bnd2bndPieces[bnd];
                if (bndPieces.size() == 1) continue;
                for (int i = 1; i < bndPieces.size(); i++) {
                    int bndPiece1 = bndPieces[i - 1];
                    int bndPiece2 = bndPieces[i];

                    int fc = fg.addFactorCategory(
                        [bndPiece1, bndPiece2, &mg, &bndPiece2allowedLabels](
                        const int * varlabels, size_t nvar, ml::FactorGraph::FactorCategoryId fcid, void * givenData) {

                        assert(nvar == 2);
                        BndPieceLabel bndPieceLabel1 = bndPiece2allowedLabels[bndPiece1][varlabels[0]];
                        BndPieceLabel bndPieceLabel2 = bndPiece2allowedLabels[bndPiece2][varlabels[1]];
                        double bndPieceLen1 = mg.bndPiece2length[bndPiece1];
                        double bndPieceLen2 = mg.bndPiece2length[bndPiece2];
                        double weight = bndPieceLen1 + bndPieceLen2;

                        if (bndPieceLabel1 == Connected && bndPieceLabel2 == Connected) {
                            return 0.0;
                        }
                        if ((bndPieceLabel1 == Connected) || (bndPieceLabel2 == Connected)) {
                            return weight * 5;
                        }
                        if (bndPieceLabel1 != bndPieceLabel2) {
                            return weight * 5;
                        }
                        return 0.0;

                    }, 1.0);

                    fg.addFactor({ bndPiece2vh[bndPiece1], bndPiece2vh[bndPiece2] }, fc);
                }
            }


            // junction validity
            for (int junc = 0; junc < mg.njuncs(); junc++) {
                auto & bnds = mg.junc2bnds[junc];
                if (bnds.size() <= 2) {
                    continue;
                }
                
                std::vector<int> bndPieces(bnds.size());
                std::vector<bool> bndPiecesTowardJunc(bnds.size(), false);
                std::vector<ml::FactorGraph::VarHandle> bpvhs(bnds.size());

                double bndPieceLenSum = 0.0;

                for (int i = 0; i < bnds.size(); i++) {
                    int bnd = bnds[i];
                    auto & juncs = mg.bnd2juncs[bnd];
                    bool toThisJunc = juncs.second == junc;
                    int adjacentBndPiece = toThisJunc ? mg.bnd2bndPieces[bnd].back() : mg.bnd2bndPieces[bnd].front();
                    bndPieces[i] = adjacentBndPiece;
                    bndPiecesTowardJunc[i] = toThisJunc;
                    bpvhs[i] = bndPiece2vh[adjacentBndPiece];

                    bndPieceLenSum += mg.bndPiece2length[adjacentBndPiece];
                }

                int fc = fg.addFactorCategory([bndPieceLenSum, bndPieces, bndPiecesTowardJunc, &mg, &bndPiece2allowedLabels](
                    const int * varlabels, size_t nvar, ml::FactorGraph::FactorCategoryId fcid, void * givenData) {
                    
                    assert(nvar == bndPieces.size() && nvar == bndPiecesTowardJunc.size());
                    int clockwiseFrontCount = 0, counterClockwiseFrontCount = 0;
                    for (int i = 0; i < nvar; i++) {
                        int bndPiece = bndPieces[i];
                        bool towardJunc = bndPiecesTowardJunc[i];
                        BndPieceLabel bplabel = bndPiece2allowedLabels[bndPiece][varlabels[i]];
                        switch (bplabel) {
                        case Connected: break;
                        case LeftIsFront: towardJunc ? (clockwiseFrontCount+=1) : (counterClockwiseFrontCount+=1); break;
                        case RightIsFront: towardJunc ? (counterClockwiseFrontCount += 1) : (clockwiseFrontCount += 1); break;
                        }
                    }

                    if (clockwiseFrontCount + counterClockwiseFrontCount >= 3) {
                        return bndPieceLenSum * 3.0;
                    }
                    if (clockwiseFrontCount + counterClockwiseFrontCount == 1) {
                        return 0.0;
                    }
                    if (clockwiseFrontCount != counterClockwiseFrontCount) {
                        return bndPieceLenSum * 5.0;
                    }

                    return 0.0;

                }, 1.0);
                
                assert(bpvhs.size() >= 3);
                fg.addFactor(bpvhs.begin(), bpvhs.end(), fc);

            }

            

            // solve
            std::cout << "solving factor graph" << std::endl;

            ml::FactorGraph::ResultTable bestLabels;
            double minEnergy = std::numeric_limits<double>::infinity();
            fg.solve(1000, 10, [&bestLabels, &minEnergy](int epoch, double energy, double denergy, const ml::FactorGraph::ResultTable & results) -> bool {
                std::cout << "epoch: " << epoch << "\t energy: " << energy << std::endl;
                if (energy < minEnergy) {
                    bestLabels = results;
                    minEnergy = energy;
                }
                if (denergy / std::max(energy, 1.0) > 0.01) {
                    return false;
                }
                return true;
            });


            // fill back labels
            // segs
            for (int seg = 0; seg < mg.nsegs; seg++) {
                auto sv = seg2vh[seg];
                if (sv.invalid()) {
                    continue;
                }
                int varlabel = bestLabels[sv];
                const SegLabel & seglabel = seg2allowedLabels[seg][varlabel];
                auto & c = mg.seg2control[seg];
                c.orientationClaz = seglabel.orientationClaz;
                c.orientationNotClaz = seglabel.orientationNotClaz;
            }

            for (int bndPiece = 0; bndPiece < mg.nbndPieces(); bndPiece++) {
                auto bpv = bndPiece2vh[bndPiece];
                int varlabel = bestLabels[bpv];
                BndPieceLabel bplabel = bndPiece2allowedLabels[bndPiece][varlabel];
                mg.bndPiece2segRelation[bndPiece] = (SegRelation)bplabel;
            }

        }


        struct ComparePixel {
            inline bool operator ()(const Pixel & a, const Pixel & b) const {
                if (a.x != b.x)
                    return a.x < b.x;
                return a.y < b.y;
            }
        };


        // assume that all oclcusions are described by lines
        void DetectOcclusions2(PIGraph & mg, double minAngleSizeOfLineInTJunction,
            double lambdaShrinkForHLineDetectionInTJunction,
            double lambdaShrinkForVLineDetectionInTJunction,
            double angleSizeForPixelsNearLines) {

            int width = mg.segs.cols;
            int height = mg.segs.rows;

            RTreeMap<Vec3, std::pair<Line3, int>> lineSamplesTree;
            for (int i = 0; i < mg.nlines(); i++) {
                auto & line = mg.lines[i].component;
                double spanAngle = AngleBetweenDirections(line.first, line.second);
                for (double a = 0.0; a < spanAngle; a += angleSizeForPixelsNearLines / 3.0) {
                    Vec3 sample1 = normalize(RotateDirection(line.first, line.second, a));
                    Vec3 sample2 = normalize(RotateDirection(line.first, line.second, a + angleSizeForPixelsNearLines / 3.0));
                    lineSamplesTree.emplace(normalize(sample1 + sample2), std::make_pair(Line3(sample1, sample2), i));
                }
            }


            // collect lines' nearby pixels and segs
            std::vector<std::set<Pixel, ComparePixel>> line2nearbyPixels(mg.nlines());
            std::vector<std::map<int, Vec3>> line2nearbySegsWithLocalCenterDir(mg.nlines());
            std::vector<std::map<int, bool>> line2nearbySegsWithOnLeftFlag(mg.nlines());
            std::vector<std::map<int, double>> line2nearbySegsWithWeight(mg.nlines());

            for (auto it = mg.segs.begin(); it != mg.segs.end(); ++it) {
                Pixel p = it.pos();
                double weight = cos((p.y - (height-1) / 2.0) / (height-1) * M_PI);
                Vec3 dir = normalize(mg.view.camera.toSpace(p));
                int seg = *it;
                lineSamplesTree.search(BoundingBox(dir).expand(angleSizeForPixelsNearLines * 3), 
                    [&mg, &dir, &line2nearbyPixels, &line2nearbySegsWithLocalCenterDir, &line2nearbySegsWithWeight, 
                    angleSizeForPixelsNearLines, p, seg, weight](
                    const std::pair<Vec3, std::pair<Line3, int>> & lineSample) {
                    auto line = normalize(mg.lines[lineSample.second.second].component);
                    auto dirOnLine = DistanceFromPointToLine(dir, lineSample.second.first).second.position;
                    double angleDist = AngleBetweenDirections(dir, dirOnLine);
                    double lambda = ProjectionOfPointOnLine(dir, line).ratio; // the projected position on line
                    if (angleDist < angleSizeForPixelsNearLines && IsBetween(lambda, 0.0, 1.0)) {
                        line2nearbyPixels[lineSample.second.second].insert(p);
                        auto & localCenterDir = line2nearbySegsWithLocalCenterDir[lineSample.second.second][seg];
                        localCenterDir += dir * weight;
                        double & segWeight = line2nearbySegsWithWeight[lineSample.second.second][seg];
                        segWeight += weight * Gaussian(lambda - 0.5, 0.1); // the closer to the center, the more important it is!
                    }
                    return true;
                });
            }
            for (int i = 0; i < mg.nlines(); i++) {
                auto & nearbySegsWithLocalCenterDir = line2nearbySegsWithLocalCenterDir[i];
                auto & line = mg.lines[i].component;
                Vec3 lineRight = line.first.cross(line.second);
                for (auto & segWithDir : nearbySegsWithLocalCenterDir) {
                    Vec3 centerDir = normalize(segWithDir.second);
                    bool onLeft = (centerDir - line.first).dot(lineRight) < 0;
                    line2nearbySegsWithOnLeftFlag[i][segWithDir.first] = onLeft;
                }
            }

            auto paintLineWithNearbyPixels = [&mg, &line2nearbyPixels, &line2nearbySegsWithOnLeftFlag](
                Image3f & pim, int line, 
                const gui::Color & leftColor, const gui::Color & rightColor) {
                auto & pixels = line2nearbyPixels[line];
                auto & nearbySegsWithLeftFlags = line2nearbySegsWithOnLeftFlag[line];
                for (auto & p : pixels) {
                    int seg = mg.segs(p);
                    bool onLeft = nearbySegsWithLeftFlags.at(seg);
                    if (onLeft) {
                        gui::AsCanvas(pim).color(leftColor).add(p);
                    } else {
                        gui::AsCanvas(pim).color(rightColor).add(p);
                    }
                }
            };
            if (true) {
                auto pim = Print(mg, ConstantFunctor<gui::Color>(gui::White),
                    ConstantFunctor<gui::Color>(gui::Transparent),
                    ConstantFunctor<gui::Color>(gui::LightGray), 1, 2);
                for (int i = 0; i < mg.nlines(); i++) {
                    paintLineWithNearbyPixels(pim, i, gui::Green, gui::Blue);
                }
                gui::AsCanvas(pim).show();
            }




            ml::FactorGraph fg;

            std::cout << "building factor graph" << std::endl;

            std::vector<std::vector<LineLabel>> line2allowedLabels(mg.nlines());
            std::vector<ml::FactorGraph::VarHandle> line2vh(mg.nlines());
            std::map<ml::FactorGraph::VarHandle, int> vh2line;
            for (int line = 0; line < mg.nlines(); line++) {             
                if (mg.lines[line].claz == -1) {
                    continue;
                }
                auto & l = mg.lines[line].component;
                std::vector<LineLabel> allowedLabels = { { true, true }, { true, false }, { false, true }, { false, false } };
                auto vh = fg.addVar(fg.addVarCategory(allowedLabels.size(), AngleBetweenDirections(l.first, l.second)));
                line2vh[line] = vh;
                vh2line[vh] = line;
                line2allowedLabels[line] = std::move(allowedLabels);
            }

            // factors preparation

            // orientation consistency term between seg and line
           
            std::vector<LineLabelCost> line2labelCost(mg.nlines(), { 0, 1e-8, 1e-8, 1.0 });

            for (int line = 0; line < mg.nlines(); line++) {
                auto vh = line2vh[line];
                if (vh.invalid()) {
                    continue;
                }
                int claz = mg.lines[line].claz;
                auto & l = mg.lines[line].component;
                const auto & lps = mg.line2linePieces[line];
                auto & labelCosts = line2labelCost[line];

                // consider nearby segs
                for (auto & segAndWeight : line2nearbySegsWithWeight[line]) {
                    int seg = segAndWeight.first;
                    double weight = segAndWeight.second;
                    bool onLeft = line2nearbySegsWithOnLeftFlag[line].at(seg);
                    auto & segControl = mg.seg2control[seg];
                    if (segControl.orientationClaz != -1 && segControl.orientationClaz == claz) {
                        if (onLeft) {
                            labelCosts.connectLeftConnectRight += weight;
                            labelCosts.connectLeftDisconnectRight += weight;
                        } else {
                            labelCosts.connectLeftConnectRight += weight;
                            labelCosts.disconnectLeftConnectRight += weight;
                        }
                    }
                }

                //// consider linepiece-bndpiece connectivities
                //for (int lp : lps) {                 
                //    if (mg.linePiece2seg[lp] != -1) {
                //        double weight = mg.linePiece2length[lp] * 10;
                //        int seg = mg.linePiece2seg[lp];
                //        auto & segControl = mg.seg2control[seg];
                //        if (segControl.orientationClaz != -1 && segControl.orientationClaz == claz) {
                //            labelCosts.connectLeftConnectRight += weight;
                //            labelCosts.connectLeftDisconnectRight += weight;
                //            labelCosts.disconnectLeftConnectRight += weight;
                //        }
                //    } else {
                //        int bp = mg.linePiece2bndPiece[lp];
                //        double weight = mg.bndPiece2length[bp] * 10;
                //        int seg1, seg2;
                //        std::tie(seg1, seg2) = mg.bnd2segs[mg.bndPiece2bnd[bp]];
                //        if (!mg.linePiece2bndPieceInSameDirection[lp]) {
                //            std::swap(seg1, seg2);
                //        }
                //        auto & segControl1 = mg.seg2control[seg1];
                //        if (segControl1.orientationClaz != -1 && segControl1.orientationClaz == claz) {
                //            labelCosts.connectLeftConnectRight += weight;
                //            labelCosts.connectLeftDisconnectRight += weight;
                //        }
                //        auto & segControl2 = mg.seg2control[seg2];
                //        if (segControl2.orientationClaz != -1 && segControl2.orientationClaz == claz) {
                //            labelCosts.connectLeftConnectRight += weight;
                //            labelCosts.disconnectLeftConnectRight += weight;
                //        }
                //    }
                //}
            }


            // {line-lineRelation-line} -> linePiece -> bndPiece, occlusion hints suggested by line's t-junction
            std::vector<std::pair<int, int>> tjunctionLines; // [hline, vline]
            std::vector<bool> tjunctionVLineLiesOnTheRightOfHLine;
            for (int lineRelation = 0; lineRelation < mg.nlineRelations(); lineRelation++) {
                if (mg.lineRelation2IsIncidence[lineRelation]) {
                    continue;
                }
                int line1 = -1, line2 = -1;
                std::tie(line1, line2) = mg.lineRelation2lines[lineRelation];
                // whether this forms a t-junction?
                if (mg.lines[line1].claz == mg.lines[line2].claz) {
                    continue;
                }
                auto l1 = normalize(mg.lines[line1].component);
                Vec3 n1 = normalize(l1.first.cross(l1.second));
                auto l2 = normalize(mg.lines[line2].component);
                Vec3 n2 = normalize(l2.first.cross(l2.second));
                if (AngleBetweenUndirectedVectors(n1, n2) < DegreesToRadians(3)) {
                    continue;
                }

                if (AngleBetweenDirections(l1.first, l1.second) < minAngleSizeOfLineInTJunction ||
                    AngleBetweenDirections(l2.first, l2.second) < minAngleSizeOfLineInTJunction) {
                    continue;
                }

                int hline = -1, vline = -1;
                double lambda1 = 0, lambda2 = 0;
                DistanceBetweenTwoLines(l1.ray(), l2.ray(), &lambda1, &lambda2);
                const double shrinkRatio = lambdaShrinkForHLineDetectionInTJunction, shrinkRatio2 = lambdaShrinkForVLineDetectionInTJunction;
                if (lambda1 < (1.0 - shrinkRatio) && lambda1 > shrinkRatio && (lambda2 < -shrinkRatio2 || lambda2 >(1.0 + shrinkRatio2))) {
                    hline = line1; vline = line2;
                } else if (lambda2 < (1.0 - shrinkRatio) && lambda2 > shrinkRatio && (lambda1 < -shrinkRatio2 || lambda1 >(1.0 + shrinkRatio2))) {
                    hline = line2; vline = line1;
                } else {
                    continue;
                }

                // Q: what if three lines form a T structure, and the center lies between the two hlines?
                // A: if the two hlines are close and are colinear, then they should have been merged after the line extraction step

                tjunctionLines.emplace_back(hline, vline);
                auto vlcenter = normalize(mg.lines[vline].component.center());
                auto & hl = mg.lines[hline].component;
                auto hlright = hl.first.cross(hl.second);
                bool vlIsOnRightOfHl = (vlcenter - hl.first).dot(hlright) > 0;
                tjunctionVLineLiesOnTheRightOfHLine.push_back(vlIsOnRightOfHl);
            }

            for (int i = 0; i < tjunctionLines.size(); i++) {
                auto & labelCosts = line2labelCost[tjunctionLines[i].first];
                auto & vl = mg.lines[tjunctionLines[i].second].component;
                double vllen = AngleBetweenDirections(vl.first, vl.second);
                double weight = Gaussian(vllen, DegreesToRadians(10)) * DegreesToRadians(5);
                if (tjunctionVLineLiesOnTheRightOfHLine[i]) {
                    labelCosts.disconnectLeftConnectRight += vllen;
                } else {
                    labelCosts.connectLeftDisconnectRight += vllen;
                }
            }


            if (true) {
                auto pim = Print(mg, ConstantFunctor<gui::Color>(gui::White),
                    [&line2labelCost, &mg](int lp) {
                    return gui::Transparent;
                }, ConstantFunctor<gui::Color>(gui::LightGray), 2, 3);
                const auto drawLine = [&mg, &pim](const Line3 & line, const gui::Color & color, const std::string & text) {
                    double angle = AngleBetweenDirections(line.first, line.second);
                    std::vector<Pixel> ps;
                    for (double a = 0.0; a <= angle; a += 0.01) {
                        ps.push_back(ToPixel(mg.view.camera.toScreen(RotateDirection(line.first, line.second, a))));
                    }
                    for (int i = 1; i < ps.size(); i++) {
                        auto p1 = ps[i - 1];
                        auto p2 = ps[i];
                        if (Distance(p1, p2) >= pim.cols / 2) {
                            continue;
                        }
                        cv::clipLine(cv::Rect(0, 0, pim.cols, pim.rows), p1, p2);
                        cv::line(pim, p1, p2, (cv::Scalar)color / 255.0, 2);
                    }
                    cv::circle(pim, ps.back(), 3.0, (cv::Scalar)color / 255.0, 2);
                    cv::putText(pim, text, ps.back(), 1, 1.0, color);
                };
                auto randctable = gui::CreateRandomColorTableWithSize(tjunctionLines.size());
                for (int line = 0; line < line2labelCost.size(); line++) {
                    auto & labelCosts = line2labelCost[line];
                    double connectLeftCost = labelCosts.connectLeftConnectRight + labelCosts.connectLeftDisconnectRight;
                    double connectRightCost = labelCosts.connectLeftConnectRight + labelCosts.disconnectLeftConnectRight;
                    double allCost = labelCosts.connectLeftConnectRight + labelCosts.connectLeftDisconnectRight +
                        labelCosts.disconnectLeftConnectRight + labelCosts.disconnectAll;
                    gui::Color color = Vec3(0.0, connectRightCost, connectLeftCost) / allCost;
                    drawLine(mg.lines[line].component, color * 0.5, std::to_string(line));
                }
                gui::AsCanvas(pim).show();
            }


            // make factors
            // data term
            for (int line = 0; line < mg.nlines(); line++) {
                auto vh = line2vh[line];
                if (vh.invalid()) {
                    continue;
                }
                int fc = fg.addFactorCategory([&line2labelCost, &mg, &line2allowedLabels, line](
                    const int * varlabels, size_t nvar, ml::FactorGraph::FactorCategoryId fcid, void * givenData)->double {
                    assert(nvar == 1);
                    auto & labelCosts = line2labelCost[line];
                    
                    double costSum = labelCosts.connectLeftConnectRight + 
                        labelCosts.connectLeftDisconnectRight +
                        labelCosts.disconnectLeftConnectRight + 
                        labelCosts.disconnectAll/* + 1e-3*/;

                    const LineLabel & label = line2allowedLabels[line][varlabels[0]];
                    
                    double ret = 0.0;
                    if (label.connectLeft && label.connectRight) {
                        ret += labelCosts.connectLeftConnectRight;
                    } else if (label.connectLeft && !label.connectRight) {
                        ret += labelCosts.connectLeftDisconnectRight;
                    } else if (!label.connectLeft && label.connectRight) {
                        ret += labelCosts.disconnectLeftConnectRight;
                    } else if (!label.connectLeft && !label.connectRight) {
                        ret += labelCosts.disconnectAll;
                    }

                    return ret / costSum * 1.0;
                }, 1.0);
                fg.addFactor({ vh }, fc);
            }

            // smooth term
            for (int lineRelation = 0; lineRelation < mg.nlineRelations(); lineRelation++) {
                if (!mg.lineRelation2IsIncidence[lineRelation]) {
                    continue;
                }
                int line1, line2;
                std::tie(line1, line2) = mg.lineRelation2lines[lineRelation];
                auto vh1 = line2vh[line1], vh2 = line2vh[line2];
                if (vh1.invalid() || vh2.invalid()) {
                    continue;
                }
                auto & l1 = mg.lines[line1].component;
                auto & l2 = mg.lines[line2].component;
                auto nearest = DistanceBetweenTwoLines(l1, l2);
                double angleDist = AngleBetweenDirections(nearest.second.first.position, nearest.second.second.position);
                double weight = Gaussian(angleDist, DegreesToRadians(5));
                bool sameDirection = l1.first.cross(l1.second).dot(l2.first.cross(l2.second)) > 0;

                int fc = fg.addFactorCategory([&line2allowedLabels, line1, line2, sameDirection, weight](
                    const int * varlabels, size_t nvar, ml::FactorGraph::FactorCategoryId fcid, void * givenData)->double {
                    assert(nvar == 2);
                    const LineLabel & label1 = line2allowedLabels[line1][varlabels[0]];
                    const LineLabel & label2 = line2allowedLabels[line2][varlabels[1]];
                    if ((sameDirection && label1 != label2) || (!sameDirection && label1 != label2.leftRightSwapped())) {
                        return weight;
                    }
                    return 0.0;
                }, 1.0);
                fg.addFactor({ vh1, vh2 }, fc);
            }



            // solve
            std::cout << "solving factor graph" << std::endl;

            ml::FactorGraph::ResultTable bestLabels;
            double minEnergy = std::numeric_limits<double>::infinity();
            fg.solve(20, 10, [&bestLabels, &minEnergy](int epoch, double energy, double denergy, const ml::FactorGraph::ResultTable & results) -> bool {
                std::cout << "epoch: " << epoch << "\t energy: " << energy << std::endl;
                if (energy < minEnergy) {
                    bestLabels = results;
                    minEnergy = energy;
                }
                if (denergy / std::max(energy, 1.0) >= 1e-3) {
                    return false;
                }
                return true;
            });


            if (true) {
                auto pim = Print(mg, ConstantFunctor<gui::Color>(gui::White),
                    [&line2labelCost, &mg](int lp) {
                    return gui::Transparent;
                }, ConstantFunctor<gui::Color>(gui::LightGray), 2, 3);
                auto drawLine = [&mg, &pim](const Line3 & line, const std::string & text, gui::Color & color) {
                    double angle = AngleBetweenDirections(line.first, line.second);
                    std::vector<Pixel> ps;
                    for (double a = 0.0; a <= angle; a += 0.01) {
                        ps.push_back(ToPixel(mg.view.camera.toScreen(RotateDirection(line.first, line.second, a))));
                    }
                    for (int i = 1; i < ps.size(); i++) {
                        auto p1 = ps[i - 1];
                        auto p2 = ps[i];
                        if (Distance(p1, p2) >= pim.cols / 2) {
                            continue;
                        }
                        cv::clipLine(cv::Rect(0, 0, pim.cols, pim.rows), p1, p2);
                        cv::line(pim, p1, p2, (cv::Scalar)color / 255.0, 2);
                    }
                    cv::circle(pim, ps.back(), 2.0, (cv::Scalar)color / 255.0, 2);
                    cv::putText(pim, text, ps.back(), 1, 0.7, color);
                };
                for (int line = 0; line < line2labelCost.size(); line++) {
                    auto vh = line2vh[line];
                    if (vh.invalid()) {
                        continue;
                    }
                    const LineLabel & label = line2allowedLabels[line][bestLabels[vh]];
                    if (label.connectLeft && !label.connectRight) {
                        drawLine(mg.lines[line].component, std::to_string(line), gui::Color(gui::Green));
                    } else if (!label.connectLeft && label.connectRight) {
                        drawLine(mg.lines[line].component, std::to_string(line), gui::Color(gui::Blue));
                    } else if (!label.connectLeft && !label.connectRight) {
                        drawLine(mg.lines[line].component, std::to_string(line), gui::Color(gui::Red));
                    } else {
                        drawLine(mg.lines[line].component, std::to_string(line), gui::Color(gui::Black));
                    }
                }
                gui::AsCanvas(pim).show(0, "line labels");
            }





            // fill back labels
            // initialize all seg-seg, seg-line and line-line relations
            for (int bp = 0; bp < mg.nbndPieces(); bp++) {
                mg.bndPiece2segRelation[bp] = SegRelation::Connected;
            }
            for (int lp = 0; lp < mg.nlinePieces(); lp++) {
                mg.linePiece2segLineRelation[lp] = SegLineRelation::Attached;
            }
            for (int lr = 0; lr < mg.nlineRelations(); lr++) {
                mg.lineRelations[lr] = LineRelation::Attached;
            }
            // set relations according to the occluding lines
            for (int line = 0; line < mg.nlines(); line++) {
                auto vh = line2vh[line];
                if (vh.invalid()) {
                    continue;
                }
                const LineLabel & label = line2allowedLabels[line][bestLabels[vh]];
                if (label.connectLeft && label.connectRight) {
                    continue;
                }


                if (!label.connectLeft && !label.connectRight) {
                    for (int lp : mg.line2linePieces[line]) {
                        mg.linePiece2segLineRelation[lp] = SegLineRelation::Detached;
                    }
                    continue;
                } 
                
                std::set<int> leftSegs, rightSegs;
                for (auto & pp : line2nearbySegsWithOnLeftFlag[line]) {
                    int seg = pp.first;
                    bool onLeft = pp.second;
                    if (onLeft) {
                        leftSegs.insert(seg);
                    } else {
                        rightSegs.insert(seg);
                    }
                }

                if (label.connectLeft && !label.connectRight) { // connect left side only
                    for (int frontSeg : leftSegs) {
                        for (int bnd : mg.seg2bnds[frontSeg]) {
                            int anotherSeg = mg.bnd2segs[bnd].first;
                            if (anotherSeg == frontSeg) {
                                anotherSeg = mg.bnd2segs[bnd].second;
                            }
                            if (Contains(rightSegs, anotherSeg)) {
                                bool frontSegIsOnBndLeft = frontSeg == mg.bnd2segs[bnd].first;
                                if (frontSegIsOnBndLeft) {
                                    for (int bndPiece : mg.bnd2bndPieces[bnd]) {
                                        auto & d1 = mg.bndPiece2dirs[bndPiece].front();
                                        auto & d2 = mg.bndPiece2dirs[bndPiece].back();
                                        auto & l = mg.lines[line].component;
                                        double angleRadius = (AngleBetweenDirections(l.first, l.second) / 2.0);
                                        if (AngleBetweenDirections(d1, l.center()) <= angleRadius && 
                                            AngleBetweenDirections(d2, l.center()) <= angleRadius) {
                                            mg.bndPiece2segRelation[bndPiece] = SegRelation::LeftIsFront;
                                        }
                                    }
                                } else {
                                    for (int bndPiece : mg.bnd2bndPieces[bnd]) {
                                        auto & d1 = mg.bndPiece2dirs[bndPiece].front();
                                        auto & d2 = mg.bndPiece2dirs[bndPiece].back();
                                        auto & l = mg.lines[line].component;
                                        double angleRadius = (AngleBetweenDirections(l.first, l.second) / 2.0);
                                        if (AngleBetweenDirections(d1, l.center()) <= angleRadius &&
                                            AngleBetweenDirections(d2, l.center()) <= angleRadius) {
                                            mg.bndPiece2segRelation[bndPiece] = SegRelation::RightIsFront;
                                        }
                                    }
                                }
                            }
                        }
                    }
                } else if (!label.connectLeft && label.connectRight) { // connect right side only
                    for (int frontSeg : rightSegs) {
                        for (int bnd : mg.seg2bnds[frontSeg]) {
                            int anotherSeg = mg.bnd2segs[bnd].first;
                            if (anotherSeg == frontSeg) {
                                anotherSeg = mg.bnd2segs[bnd].second;
                            }
                            if (Contains(leftSegs, anotherSeg)) {
                                bool frontSegIsOnBndLeft = frontSeg == mg.bnd2segs[bnd].first;
                                if (frontSegIsOnBndLeft) {
                                    for (int bndPiece : mg.bnd2bndPieces[bnd]) {
                                        auto & d1 = mg.bndPiece2dirs[bndPiece].front();
                                        auto & d2 = mg.bndPiece2dirs[bndPiece].back();
                                        auto & l = mg.lines[line].component;
                                        double angleRadius = (AngleBetweenDirections(l.first, l.second) / 2.0);
                                        if (AngleBetweenDirections(d1, l.center()) <= angleRadius &&
                                            AngleBetweenDirections(d2, l.center()) <= angleRadius) {
                                            mg.bndPiece2segRelation[bndPiece] = SegRelation::LeftIsFront;
                                        }
                                    }
                                } else {
                                    for (int bndPiece : mg.bnd2bndPieces[bnd]) {
                                        auto & d1 = mg.bndPiece2dirs[bndPiece].front();
                                        auto & d2 = mg.bndPiece2dirs[bndPiece].back();
                                        auto & l = mg.lines[line].component;
                                        double angleRadius = (AngleBetweenDirections(l.first, l.second) / 2.0);
                                        if (AngleBetweenDirections(d1, l.center()) <= angleRadius &&
                                            AngleBetweenDirections(d2, l.center()) <= angleRadius) {
                                            mg.bndPiece2segRelation[bndPiece] = SegRelation::RightIsFront;
                                        }
                                    }
                                }
                            }
                        }
                    }
                }
                
                // line relations
                for (int lr : mg.line2lineRelations[line]) {
                    int anotherLine = mg.lineRelation2lines[lr].first;
                    if (anotherLine == line) {
                        anotherLine = mg.lineRelation2lines[lr].second;
                    }
                    auto anotherLineCenter = normalize(mg.lines[anotherLine].component.center());
                    auto & thisLine = mg.lines[line].component;
                    auto thisLineRightSide = thisLine.first.cross(thisLine.second);
                    bool anotherLineIsOnTheRightSide = (anotherLineCenter - thisLine.first).dot(thisLineRightSide) > 0;
                    if (!label.connectLeft && !anotherLineIsOnTheRightSide ||
                        !label.connectRight && anotherLineIsOnTheRightSide) {
                        mg.lineRelations[lr] = LineRelation::Detached;
                    }
                }              
            }

        }

    }
}