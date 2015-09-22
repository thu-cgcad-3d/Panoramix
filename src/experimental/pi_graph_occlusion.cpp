#include "../core/utility.hpp"
#include "../ml/factor_graph.hpp"
#include "pi_graph_occlusion.hpp"

namespace pano {
    namespace experimental {

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
                double c_i = mg.seg2area[i] * 10.0;
                
                if (c.orientationClaz != -1) {
                    allowedLabels = { {c.orientationClaz, -1} };                   
                } else if (c.orientationNotClaz != -1) {
                    allowedLabels = { {-1, c.orientationNotClaz} };
                    //for (int k = 0; k < mg.vps.size(); k++) {
                    //    if (k == c.orientationNotClaz) continue;
                    //    assert(IsFuzzyPerpendicular(mg.vps[k], mg.vps[c.orientationNotClaz]));
                    //    allowedLabels.push_back({ k, -1 });
                    //}
                } else {
                    allowedLabels.push_back({ -1, -1 });
                    //for (int k = 0; k < mg.vps.size(); k++) {
                    //    //allowedLabels.push_back({ k, -1 });
                    //    allowedLabels.push_back({ -1, k });
                    //}
                }

                auto vh = fg.addVar(fg.addVarCategory(allowedLabels.size(), c_i));
                seg2vh[i] = vh;
                seg2allowedLabels[i] = std::move(allowedLabels);
                vh2seg[vh] = i;
            }


            // bndPieces
            enum BndPieceLabel : int { 
                Connected = (int)OcclusionRelation::Connected, 
                LeftIsFront = (int)OcclusionRelation::LeftIsFront, 
                RightIsFront = (int)OcclusionRelation::RightIsFront
            };
            std::vector<ml::FactorGraph::VarHandle> bndPiece2vh(mg.nbndPieces());
            std::vector<std::vector<BndPieceLabel>> bndPiece2allowedLabels(mg.nbndPieces());
            std::map<ml::FactorGraph::VarHandle, int> vh2bndPiece;
            for (int i = 0; i < mg.nbndPieces(); i++) {
                std::vector<BndPieceLabel> allowedLabels;
                double c_i = mg.bndPiece2length[i];

                if (mg.bndPiece2occlusion[i] != OcclusionRelation::Unknown) {
                    allowedLabels = { (BndPieceLabel)mg.bndPiece2occlusion[i] };
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

            // linePieces (that are not bound to bndPieces)
            enum LinePieceLabel {
                Attached = (int)AttachmentRelation::Attached, 
                Detached = (int)AttachmentRelation::Detached
            };
            std::vector<ml::FactorGraph::VarHandle> linePiece2vh(mg.nlinePieces());
            std::vector<std::vector<LinePieceLabel>> linePiece2allowedLabels(mg.nlinePieces());
            std::map<ml::FactorGraph::VarHandle, int> vh2linePiece;
            for (int i = 0; i < mg.nlinePieces(); i++) {
                bool hasBndPiece = mg.linePiece2bndPiece[i] != -1;
                if (hasBndPiece) {
                    continue;
                }
                std::vector<LinePieceLabel> allowedLabels;
                if (mg.linePiece2attachment[i] != AttachmentRelation::Unknown) {
                    allowedLabels = { (LinePieceLabel)mg.linePiece2attachment[i] };
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
                    int dof = seg2allowedLabels[seg][varlabels[0]].dof();
                    //return mg.seg2area[seg] * (dof - 1) * 0.01;
                    return 0.0;
                }, 1.0);
                fg.addFactor({ seg2vh[seg] }, fc);
            }

            // bndPiece
            for (int bndPiece = 0; bndPiece < mg.nbndPieces(); bndPiece++) {
                int fc = fg.addFactorCategory([bndPiece, &bndPiece2allowedLabels, &mg](const int * varlabels, size_t nvar,
                    ml::FactorGraph::FactorCategoryId fcid, void * givenData)->double {
                    assert(nvar == 1);
                    BndPieceLabel label = bndPiece2allowedLabels[bndPiece][varlabels[0]];
                    double len = mg.bndPiece2length[bndPiece];
                    if (label != Connected) {
                        return mg.bndPiece2linePieces[bndPiece].empty() ? len * 0.01 : 0.0;
                    }
                    return 0.0;
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
                    if (label != Attached) {
                        return len * 0.005;
                    }
                    return 0.0;
                }, 1.0);
                fg.addFactor({ linePiece2vh[linePiece] }, fc);
            }

            // seg-bndPiece-seg
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
                        
                        auto seglabel1 = seg2allowedLabels[seg1][varlabels[0]];
                        auto seglabel2 = seg2allowedLabels[seg2][varlabels[1]];
                        BndPieceLabel bndPieceLabel = bndPiece2allowedLabels[bndPiece][varlabels[2]];

                        int bndPieceClaz = mg.bndPiece2classes[bndPiece];
                        double bndPieceLen = mg.bndPiece2length[bndPiece];

                        if (seglabel1.dof() == 3 || seglabel2.dof() == 3) {
                            return 0.0;
                        }

                        // 1 -- 1
                        //if (seglabel1.dof() == 1 && seglabel2.dof() == 1 &&
                        //    seglabel1 == seglabel2 &&
                        //    (bndPieceClaz == seglabel1.orientationClaz || bndPieceClaz == -1) && 
                        //    bndPieceLabel != Connected) {
                        //    return bndPieceLen * 3.0;
                        //}

                        if (seglabel1.dof() == 1 && seglabel2.dof() == 1 && 
                            seglabel1.orientationClaz != seglabel2.orientationClaz && 
                            bndPieceClaz == -1 && bndPieceLabel == Connected) {
                            return bndPieceLen * 0.1;
                        }
                        if (seglabel1.dof() == 1 && seglabel2.dof() == 1 &&
                            seglabel1.orientationClaz != seglabel2.orientationClaz &&
                            (bndPieceClaz == seglabel1.orientationClaz || bndPieceClaz == seglabel2.orientationClaz) && bndPieceLabel == Connected) {
                            return bndPieceLen * 5.0;
                        }


                        if (seglabel1.dof() == 1 && seglabel2.dof() == 1 &&
                            seglabel1.orientationClaz != seglabel2.orientationClaz &&
                            bndPieceClaz == seglabel1.orientationClaz && bndPieceLabel != RightIsFront) {
                            return bndPieceLen * 5.0;
                        }
                        if (seglabel1.dof() == 1 && seglabel2.dof() == 1 &&
                            seglabel1.orientationClaz != seglabel2.orientationClaz &&
                            bndPieceClaz == seglabel2.orientationClaz && bndPieceLabel != LeftIsFront) {
                            return bndPieceLen * 5.0;
                        }

                        // 1 -- 2
                        if (seglabel1.dof() == 1 && seglabel2.dof() == 2 && 
                            seglabel1.orientationClaz == seglabel2.orientationNotClaz &&
                            bndPieceClaz == -1 && bndPieceLabel != RightIsFront) {
                            return bndPieceLen * 1.0;
                        }
                        if (seglabel1.dof() == 2 && seglabel2.dof() == 1 &&
                            seglabel1.orientationNotClaz == seglabel2.orientationClaz &&
                            bndPieceClaz == -1 && bndPieceLabel != LeftIsFront) {
                            return bndPieceLen * 1.0;
                        }

                        if (seglabel1.dof() == 1 && seglabel2.dof() == 2 &&
                            seglabel1.orientationClaz == seglabel2.orientationNotClaz &&
                            bndPieceClaz == seglabel1.orientationClaz && bndPieceLabel != RightIsFront) {
                            return bndPieceLen * 100.0;
                        }
                        if (seglabel1.dof() == 2 && seglabel2.dof() == 1 &&
                            seglabel1.orientationNotClaz == seglabel2.orientationClaz &&
                            bndPieceClaz == seglabel2.orientationClaz && bndPieceLabel != LeftIsFront) {
                            return bndPieceLen * 100.0;
                        }                        

                        // 2 -- 2
                        //if (seglabel1.dof() == 2 && seglabel2.dof() == 2 &&
                        //    seglabel1.orientationNotClaz == seglabel2.orientationNotClaz &&
                        //    bndPieceClaz != -1 && bndPieceClaz != seglabel1.orientationNotClaz && bndPieceLabel != Connected) {
                        //    return bndPieceLen * 3.0;
                        //}

                        return 0.0;
                    }, 1.0);
                    fg.addFactor({ sv1, sv2, bndPiece2vh[bndPiece] }, fc);
                }    
            }

            // seg supporting by bnd
            for (int seg = 0; seg < mg.nsegs; seg++) {
                auto sv = seg2vh[seg];
                if (sv.invalid()) continue;

                int bndPiecesNum = 0;
                for (int bnd : mg.seg2bnds[seg]) {
                    bndPiecesNum += mg.bnd2bndPieces[bnd].size();
                }

                for (int bnd : mg.seg2bnds[seg]) {
                    for (int bndPiece : mg.bnd2bndPieces[bnd]) {
                        int fc = fg.addFactorCategory([seg, bndPiece, bndPiecesNum, &mg, &seg2allowedLabels, &bndPiece2allowedLabels](
                            const int * varlabels, size_t nvar,
                            ml::FactorGraph::FactorCategoryId fcid, void * givenData) {
                            assert(nvar == 2);
                            auto seglabel = seg2allowedLabels[seg][varlabels[0]];
                            BndPieceLabel bndPieceLabel = bndPiece2allowedLabels[bndPiece][varlabels[1]];
                            
                            int bndPieceClaz = mg.bndPiece2classes[bndPiece];

                            if (bndPieceLabel == Connected) {
                                return 0.0;
                            }
                            return 1.0 / bndPiecesNum * 0.01;
                            //return 0.0;
                        }, 1.0);
                        fg.addFactor({ sv, bndPiece2vh[bndPiece] }, fc);
                    }
                }
            }


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


            // line alignment
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
                                return weight * 100;
                            }
                            if (sameDirection && bndPieceLabel1 != bndPieceLabel2) {
                                return weight * 100;
                            }
                            if (!sameDirection && bndPieceLabel1 == bndPieceLabel2) {
                                return weight * 100;
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

                            if (bndPieceLabel == Connected && linePieceLabel == Attached ||
                                bndPieceLabel != Connected && linePieceLabel == Detached) {
                                return 0.0;
                            }
                            return weight * 10;
                        }, 1.0);
                        fg.addFactor({ bndPiece2vh[bndPiece], linePiece2vh[linePiece] }, fc);
                    }
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
                mg.bndPiece2occlusion[bndPiece] = (OcclusionRelation)bplabel;
            }

        }


        void AssumeThereAreNoOcclusions(PIGraph & mg) {
            for (int bp = 0; bp < mg.nbndPieces(); bp++) {
                mg.bndPiece2occlusion[bp] = OcclusionRelation::Connected;
            }
            for (int lp = 0; lp < mg.nlinePieces(); lp++) {
                mg.linePiece2attachment[lp] = AttachmentRelation::Attached;
            }
        }

    }
}