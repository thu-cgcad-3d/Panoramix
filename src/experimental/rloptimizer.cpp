#include "../core/utility.hpp"
#include "../core/containers.hpp"
#include "projective_solver.hpp"
#include "rloptimizer.hpp"

namespace panoramix {
    namespace experimental {

        std::vector<RLGraphPatch> MakePatches(const RLGraph & g, int n, double angleRadius) {
            auto rhVisited = g.createComponentTable<RData>(false);
            std::vector<RHandle> rhs;
            rhs.reserve(g.internalComponents<RData>().size());
            for (auto & r : g.components<RData>()) {
                rhs.push_back(r.topo.hd);
            }

            std::vector<RLGraphPatch> patches;
            patches.reserve(rhs.size());

            for (auto centerRh : rhs) {
                RLGraphPatch patch;
                MaxHeap<RHandle, double> heapr;
                MaxHeap<LHandle, double> heapl;
                heapr.push(centerRh, 1.0);

                auto & center = g.data(centerRh).normalizedCenter;

                while ((!heapr.empty() || !heapl.empty()) && patch.size() < n) {
                    RHandle curRh;
                    double curRhScore = -1.0;
                    if (!heapr.empty()) {
                        curRh = heapr.top();
                        curRhScore = heapr.topScore();
                    }
                    LHandle curLh;
                    double curLhScore = -1.0;
                    if (!heapl.empty()) {
                        curLh = heapl.top();
                        curLhScore = heapl.topScore();
                    }

                    std::vector<RHandle> rhcands;
                    std::vector<LHandle> lhcands;
                    if (curLhScore < curRhScore) {
                        assert(curRh.valid());
                        patch.insert(curRh);
                        heapr.pop();
                      
                        for (auto rrh : g.topo(curRh).constraints<RRData>()) {
                            auto anotherRh = g.topo(rrh).component<0>();
                            if (anotherRh == curRh) {
                                anotherRh = g.topo(rrh).component<1>();
                            }
                            rhcands.push_back(anotherRh);
                        }
                        for (auto rlh : g.topo(curRh).constraints<RLData>()) {
                            LHandle lh = g.topo(rlh).component<1>();
                            lhcands.push_back(lh);
                        }
                        for (auto rrlh : g.topo(curRh).constraints<RRLData>()) {
                            RHandle anotherRh = g.topo(rrlh).component<0>();
                            if (anotherRh == curRh) {
                                anotherRh = g.topo(rrlh).component<1>();
                            }
                            rhcands.push_back(anotherRh);                           
                            LHandle lh = g.topo(rrlh).component<2>();
                            lhcands.push_back(lh);
                        }
                    } else {
                        assert(curLh.valid());
                        patch.insert(curLh);
                        heapl.pop();

                        for (auto llh : g.topo(curLh).constraints<LLData>()) {
                            auto anotherLh = g.topo(llh).component<0>();
                            if (anotherLh == curLh) {
                                anotherLh = g.topo(llh).component<1>();
                            }
                            lhcands.push_back(anotherLh);
                        }
                    }

                    for (auto rhcand : rhcands) {
                        if (patch.contains(rhcand)) {
                            continue;
                        }
                        double angleDist = AngleBetweenDirections(center, g.data(rhcand).normalizedCenter);
                        if (angleDist > angleRadius) {
                            continue;
                        }
                        double score = 1.0 / angleDist;
                        heapr.push(rhcand, score);
                    }
                    for (auto lhcand : lhcands) {
                        if (patch.contains(lhcand)) {
                            continue;
                        }
                        auto nearest = DistanceFromPointToLine(center, g.data(lhcand).normalizedLine).second.position;
                        double angleDist = AngleBetweenDirections(center, nearest);
                        if (angleDist > angleRadius) {
                            continue;
                        }
                        double score = 1.0 / angleDist;
                        heapl.push(lhcand, score);
                    }
                }

                // fill constraints
                for (auto rh : patch.container<RHandle>()) {
                    for (auto h : g.topo(rh).constraints<RRData>()) {
                        if (patch.contains(g.topo(h).component<0>()) && patch.contains(g.topo(h).component<1>())) {
                            patch.insert(h);
                        }
                    }
                    for (auto h : g.topo(rh).constraints<RLData>()) {
                        if (patch.contains(g.topo(h).component<0>()) && patch.contains(g.topo(h).component<1>())) {
                            patch.insert(h);
                        }
                    }
                    for (auto h : g.topo(rh).constraints<RRLData>()) {
                        if (patch.contains(g.topo(h).component<0>()) && 
                            patch.contains(g.topo(h).component<1>()) &&
                            patch.contains(g.topo(h).component<2>())) {
                            patch.insert(h);
                        }
                    }
                }
                for (auto lh : patch.container<LHandle>()) {
                    for (auto h : g.topo(lh).constraints<LLData>()) {
                        if (patch.contains(g.topo(h).component<0>()) && patch.contains(g.topo(h).component<1>())) {
                            patch.insert(h);
                        }
                    }
                }

                patches.push_back(std::move(patch));
            }
            
            return patches;
        }


        namespace {
            bool NextSub(int * subs, const int * dims, int len) {
                if (len == 0)
                    return false;
                subs[0]++;
                int k = 0;
                while (k < len && subs[k] >= dims[k]) {
                    subs[k] = 0;
                    if (k + 1 < len) subs[k + 1] ++;
                    k++;
                }
                if (k == len) {
                    return false; 
                }
                return true;
            }
            bool NextSub(int * subs, int dim, int len) {
                if (len == 0)
                    return false;
                subs[0]++;
                int k = 0;
                while (k < len && subs[k] >= dim) {
                    subs[k] = 0;
                    if (k + 1 < len) subs[k + 1] ++;
                    k++;
                }
                if (k == len) {
                    return false;
                }
                return true;
            }
            int Sub2Ind(const int * subs, const int * dims, int len) {
                int ind = 0;
                for (int k = len - 1; k >= 0; k--) {
                    ind = ind * dims[k] + subs[k];
                }
                return ind;
            }
            int Sub2Ind(const int * subs, int dim, int len) {
                int ind = 0;
                for (int k = len - 1; k >= 0; k--) {
                    ind = ind * dim + subs[k];
                }
                return ind;
            }
        }


        // feature layout
#pragma pack(push)
#pragma pack(1)
#define Part(name, nbytes) struct {uint8_t content[(nbytes)];} name
        struct FeatureLayout {
            Part(R, 8);
            Part(L, 1);
            Part(RR, 6);
            Part(LL, LLData::ManhattanJunctionTypeNum * 2);
            Part(RRRorRRRR, 5);
            Part(Patch, 21);
        };
#undef Part
#pragma pack(pop)
        static const int FullFeatureLength = sizeof(FeatureLayout);

#define FeatureBeginPos(name) (offsetof(FeatureLayout, name))
#define FeatureEndPos(name) (offsetof(FeatureLayout, name) + sizeof(std::declval<FeatureLayout>().name))
#define FeatureLength(name) (sizeof(std::declval<FeatureLayout>().name))

        static const int KKK = sizeof(std::declval<FeatureLayout>().R);

        void RLOptimizer::preprocess() {

#pragma region Build vars for factor graph       
            _varhs = RLGraphVarHandleTable({
                _g.internalComponents<RData>().size(),
                _g.internalComponents<LData>().size(),
                _g.internalConstraints<RRData>().size(),
                _g.internalConstraints<LLData>().size()
            }, ml::FactorGraph::VarHandle());

            _fg.clear();
            // register (discrete) var cats
            // region label: free [tovp1 tovp2 tovp3] vertical notplanar
#define R_Label_Free 0
#define R_Label_ToVP(vpid) ((vpid)+1)
#define R_Label_Vertical (_vps.size()+1)
#define R_Label_NotPlanar (_vps.size()+2)
#define R_Label_Num (_vps.size() + 3)
#define VP_From_R_Label(label) ((label)-1)
            auto catR = _fg.addVarCategory(R_Label_Num, 1.0);
            // line label: free [tovp1 tovp2 tovp3]
#define L_Label_Free 0
#define L_Label_ToVP(vpid) ((vpid)+1)
#define L_Label_Num (_vps.size() + 1)
#define VP_From_L_Label(label) ((label)-1)
            auto catL = _fg.addVarCategory(L_Label_Num, 1.0);
            // rr label: connected, 1closer, 2closer
#define RR_Label_Connected 0
#define RR_Label_FirstIsCloser 1
#define RR_Label_SecondIsCloser 2
#define RR_Label_Num 3
            auto catRR = _fg.addVarCategory(RR_Label_Num, 1.0);
            // ll label: connected, disconnected
#define LL_Label_Connected 0
#define LL_Label_Disconnected 1
#define LL_Label_Num 2
            auto catLL = _fg.addVarCategory(LL_Label_Num, 1.0);

            
            // add vars
            // r
            for (auto && c : _g.components<RData>()) {
                _varhs[c.topo.hd] = _fg.addVar(catR);
            }
            // l
            for (auto && c : _g.components<LData>()) {
                _varhs[c.topo.hd] = _fg.addVar(catL);
            }
            // rr
            for (auto && c : _g.constraints<RRData>()) {
                _varhs[c.topo.hd] = _fg.addVar(catRR);
            }
            // ll
            for (auto && c : _g.constraints<LLData>()) {
                _varhs[c.topo.hd] = _fg.addVar(catLL);
            }

#pragma endregion Build vars for factor graph


#pragma region Build factors for factor graph  

            ///// build factor graph factors 

            // resize feature tables
            _rFeatureTable.resize(_g.internalComponents<RData>().size());
            _lFeatureTable.resize(_g.internalComponents<LData>().size());
            _rrFeatureTable.resize(_g.internalConstraints<RRData>().size());
            _llFeatureTable.resize(_g.internalConstraints<LLData>().size());
            _rrrFeatureTable.resize(_g.internalConstraints<RRRData>().size());
            _rrrrFeatureTable.resize(_g.internalConstraints<RRRRData>().size());

            _patchFeatureTable.resize(_patches.size());

            // r
            {
                _rAreaSum = 0.0;
                for (auto & r : _g.components<RData>()) {
                    _rAreaSum += r.data.area;
                }
            }
            findPeakyRhs(DegreesToRadians(6));
            findHorizonRhs(DegreesToRadians(5));
            for (auto & r : _g.components<RData>()) {
                // precompute feature table
                auto & feaTable = _rFeatureTable[r.topo.hd];
                feaTable.resize(R_Label_Num);
                feaTable[R_Label_Free] = featureInR(r.topo.hd, R_Label_Free);
                for (int vpid = 0; vpid < _vps.size(); vpid++) {
                    feaTable[R_Label_ToVP(vpid)] = featureInR(r.topo.hd, R_Label_ToVP(vpid));
                }
                feaTable[R_Label_Vertical] = featureInR(r.topo.hd, R_Label_Vertical);
                feaTable[R_Label_NotPlanar] = featureInR(r.topo.hd, R_Label_NotPlanar);

                // add factor
                ml::FactorGraph::VarHandle vh = _varhs[r.topo.hd];
                auto cost = [this, &r](const int * varlabels, size_t nvar,
                    ml::FactorGraph::FactorCategoryId fcid, void * givenData) -> double {
                    assert(nvar == 1);
                    int label = varlabels[0];
                    assert(label >= 0 && label < R_Label_Num);
                    auto weights = static_cast<const double*>(givenData);
                    auto & fea = _rFeatureTable[r.topo.hd][label];
                    assert(fea.size() == FeatureLength(R));
                    return std::inner_product(weights + FeatureBeginPos(R), weights + FeatureEndPos(R), fea.begin(), 0.0);
                };
                _fg.addFactor({ vh }, _fg.addFactorCategory(ml::FactorGraph::FactorCategory{ std::move(cost), 1.0 }));
            }

            // l
            {
                _lLengthSum = 0;
                for (auto & l : _g.components<LData>()) {
                    _lLengthSum += l.data.normalizedLine.length();
                }
            }
            for (auto & l : _g.components<LData>()) {
                // precompute feature table
                auto & feaTable = _lFeatureTable[l.topo.hd];
                feaTable.resize(L_Label_Num);
                feaTable[L_Label_Free] = featureInL(l.topo.hd, L_Label_Free);
                for (int vpid = 0; vpid < _vps.size(); vpid++) {
                    feaTable[L_Label_ToVP(vpid)] = featureInL(l.topo.hd, L_Label_ToVP(vpid));
                }                
                // add factor
                ml::FactorGraph::VarHandle vh = _varhs[l.topo.hd];
                auto cost = [this, &l](const int * varlabels, size_t nvar, 
                    ml::FactorGraph::FactorCategoryId fcid, void * givenData) -> double {
                    assert(nvar == 1);
                    int label = varlabels[0];
                    assert(label >= 0 && label < L_Label_Num);
                    auto weights = static_cast<const double*>(givenData);
                    auto & fea = _lFeatureTable[l.topo.hd][label];
                    assert(fea.size() == FeatureLength(L));
                    return std::inner_product(weights + FeatureBeginPos(L), weights + FeatureEndPos(L), fea.begin(), 0.0);
                };
                _fg.addFactor({ vh }, _fg.addFactorCategory(ml::FactorGraph::FactorCategory{ std::move(cost), 1.0 }));
            }

            // rr
            {
                _rrLengthSum = 0.0;
                for (auto & rr : _g.constraints<RRData>()) {
                    _rrLengthSum += rr.data.length;
                }
            }
            for (auto & rr : _g.constraints<RRData>()) {
                // precompute feature table
                auto & feaTable = _rrFeatureTable[rr.topo.hd];
                feaTable.resize(RR_Label_Num);
                feaTable[RR_Label_Connected] = featureInRR(rr.topo.hd, RR_Label_Connected);
                feaTable[RR_Label_FirstIsCloser] = featureInRR(rr.topo.hd, RR_Label_FirstIsCloser);
                feaTable[RR_Label_SecondIsCloser] = featureInRR(rr.topo.hd, RR_Label_SecondIsCloser);

                // add factor
                ml::FactorGraph::VarHandle vh = _varhs[rr.topo.hd];
                auto cost = [this, &rr](const int * varlabels, size_t nvar,
                    ml::FactorGraph::FactorCategoryId fcid, void * givenData) -> double {
                    assert(nvar == 1);
                    int label = varlabels[0];
                    assert(label >= 0 && label < RR_Label_Num);
                    auto weights = static_cast<const double*>(givenData);
                    auto & fea = _rrFeatureTable[rr.topo.hd][label];
                    assert(fea.size() == FeatureLength(RR));
                    return std::inner_product(weights + FeatureBeginPos(RR), weights + FeatureEndPos(RR), fea.begin(), 0.0);
                };
                _fg.addFactor({ vh }, _fg.addFactorCategory(ml::FactorGraph::FactorCategory{ std::move(cost), 1.0 }));
            }

            // ll 
            countIncidenceAndIntersectionLLhs();
            for (auto & ll : _g.constraints<LLData>()) {
                // precompute feature table
                auto & feaTable = _llFeatureTable[ll.topo.hd];
                feaTable.resize(LL_Label_Num);
                feaTable[LL_Label_Connected] = featureInLL(ll.topo.hd, LL_Label_Connected);
                feaTable[LL_Label_Disconnected] = featureInLL(ll.topo.hd, LL_Label_Disconnected);

                // add factor
                ml::FactorGraph::VarHandle vh = _varhs[ll.topo.hd];
                auto cost = [this, &ll](const int * varlabels, size_t nvar,
                    ml::FactorGraph::FactorCategoryId fcid, void * givenData) -> double {
                    assert(nvar == 1);
                    int label = varlabels[0];
                    assert(label >= 0 && label < LL_Label_Num);
                    auto weights = static_cast<const double*>(givenData);
                    auto & fea = _llFeatureTable[ll.topo.hd][label];
                    assert(fea.size() == FeatureLength(LL));
                    return std::inner_product(weights + FeatureBeginPos(LL), weights + FeatureEndPos(LL), fea.begin(), 0.0);
                };
                _fg.addFactor({ vh }, _fg.addFactorCategory(ml::FactorGraph::FactorCategory{ std::move(cost), 1.0 }));
            }

            // rrr
            for (auto & rrr : _g.constraints<RRRData>()) {
                std::vector<RHandle> rhs = {
                    rrr.topo.component<0>(),
                    rrr.topo.component<1>(),
                    rrr.topo.component<2>()
                };
                std::vector<RRHandle> rrhs(rhs.size());
                std::vector<int> rrhascend(rhs.size(), true);
                for (int i = 0; i < rhs.size(); i++) {
                    auto & nextrh = rhs[(i + 1) % 3];
                    for (const RRHandle & rrh : _g.topo(rhs[i]).constraints<RRData>()) {
                        if (_g.topo(rrh).component<0>() == nextrh || _g.topo(rrh).component<1>() == nextrh) {
                            rrhs[i] = rrh;
                            rrhascend[i] = _g.topo(rrh).component<1>() == nextrh;
                            break;
                        }
                    }
                }
                for (auto & rrh : rrhs) {
                    assert(rrh.valid());
                }
                std::vector<ml::FactorGraph::VarHandle> vhs(rrhs.size());
                for (int i = 0; i < vhs.size(); i++) {
                    vhs[i] = _varhs[rrhs[i]];
                }
                std::vector<int> rrlabels(rrhs.size(), 0);
                auto & feaTable = _rrrFeatureTable[rrr.topo.hd];
                feaTable.resize(std::pow(RR_Label_Num, rrlabels.size()));
                int caseId = 0;
                do {
                    feaTable[caseId] = featureInRRRorRRRR(rrlabels.data(), rrhascend.data(), rrlabels.size());
                    caseId++;
                } while (NextSub(rrlabels.data(), RR_Label_Num, rrlabels.size()));

                auto cost = [this, &rrr](const int * varlabels, size_t nvar,
                    ml::FactorGraph::FactorCategoryId fcid, void * givenData) -> double {
                    assert(nvar == 3);                    
                    int caseId = Sub2Ind(varlabels, RR_Label_Num, nvar);
                    auto weights = static_cast<const double*>(givenData);
                    auto & fea = _rrrFeatureTable[rrr.topo.hd][caseId];
                    assert(fea.size() == FeatureLength(RRRorRRRR));
                    return std::inner_product(weights + FeatureBeginPos(RRRorRRRR), weights + FeatureEndPos(RRRorRRRR), fea.begin(), 0.0);
                };
                ml::FactorGraph::FactorCategory fc;
                fc.costs = std::move(cost);
                fc.c_alpha = 1.0;
                _fg.addFactor(vhs.begin(), vhs.end(), _fg.addFactorCategory(std::move(fc)));
            }

            // rrrr
            for (auto & rrrr : _g.constraints<RRRRData>()) {
                std::vector<RHandle> rhs = {
                    rrrr.topo.component<0>(),
                    rrrr.topo.component<1>(),
                    rrrr.topo.component<2>(),
                    rrrr.topo.component<3>()
                };
                std::vector<RRHandle> rrhs(rhs.size());
                std::vector<int> rrhascend(rhs.size(), true);
                for (int i = 0; i < rhs.size(); i++) {
                    auto & nextrh = rhs[(i + 1) % 3];
                    for (const RRHandle & rrh : _g.topo(rhs[i]).constraints<RRData>()) {
                        if (_g.topo(rrh).component<0>() == nextrh || _g.topo(rrh).component<1>() == nextrh) {
                            rrhs[i] = rrh;
                            rrhascend[i] = _g.topo(rrh).component<1>() == nextrh;
                            break;
                        }
                    }
                }
                for (auto & rrh : rrhs) {
                    assert(rrh.valid());
                }
                std::vector<ml::FactorGraph::VarHandle> vhs(rrhs.size());
                for (int i = 0; i < vhs.size(); i++) {
                    vhs[i] = _varhs[rrhs[i]];
                }
                std::vector<int> rrlabels(rrhs.size(), 0);
                auto & feaTable = _rrrrFeatureTable[rrrr.topo.hd];
                feaTable.resize(std::pow(RR_Label_Num, rrlabels.size()));
                int caseId = 0;
                do {
                    feaTable[caseId] = featureInRRRorRRRR(rrlabels.data(), rrhascend.data(), rrlabels.size());
                    caseId++;
                } while (NextSub(rrlabels.data(), RR_Label_Num, rrlabels.size()));

                auto cost = [this, &rrrr](const int * varlabels, size_t nvar,
                    ml::FactorGraph::FactorCategoryId fcid, void * givenData) -> double {
                    assert(nvar == 4);
                    int caseId = Sub2Ind(varlabels, RR_Label_Num, nvar);
                    auto weights = static_cast<const double*>(givenData);
                    auto & fea = _rrrrFeatureTable[rrrr.topo.hd][caseId];
                    assert(fea.size() == FeatureLength(RRRorRRRR));
                    return std::inner_product(weights + FeatureBeginPos(RRRorRRRR), weights + FeatureEndPos(RRRorRRRR), fea.begin(), 0.0);
                };
                ml::FactorGraph::FactorCategory fc;
                fc.costs = std::move(cost);
                fc.c_alpha = 1.0;
                _fg.addFactor(vhs.begin(), vhs.end(), _fg.addFactorCategory(std::move(fc)));
            }


            // reconstructed patches
            for (int i = 0; i < _patches.size(); i++) {
                auto & patch = _patches[i];

                std::vector<ml::FactorGraph::VarHandle> vhs;
                vhs.reserve(
                    patch.container<RHandle>().size() + 
                    patch.container<LHandle>().size() +
                    patch.container<RRHandle>().size() +
                    patch.container<LLHandle>().size()
                );
                for (auto & h : patch.container<RHandle>()) {
                    vhs.push_back(_varhs[h]);
                }
                for (auto & h : patch.container<LHandle>()) {
                    vhs.push_back(_varhs[h]);
                }
                for (auto & h : patch.container<RRHandle>()) {
                    vhs.push_back(_varhs[h]);
                }
                for (auto & h : patch.container<LLHandle>()) {
                    vhs.push_back(_varhs[h]);
                }

                auto & feaTable = _patchFeatureTable[i];
                feaTable.reserve(
                    std::pow(R_Label_Num, patch.container<RHandle>().size()) *
                    std::pow(L_Label_Num, patch.container<LHandle>().size()) *
                    std::pow(RR_Label_Num, patch.container<RRHandle>().size()) *
                    std::pow(LL_Label_Num, patch.container<LLHandle>().size())
                );
                std::vector<int> featabledims(vhs.size(), -1);
                for (int i = 0; i < vhs.size(); i++) {
                    featabledims[i] = _fg.varCategory(vhs[i]).nlabels;
                }
                std::vector<int> alllabels(featabledims.size(), 0);
                int caseId = 0;
                std::unordered_map<RHandle, Plane3> planes;
                std::unordered_map<LHandle, Line3> lines;
                do {
                    feaTable[caseId] = featureInReconstructedPatch(i, alllabels.data(), planes, lines);
                    caseId++;
                } while (NextSub(alllabels.data(), featabledims.data(), vhs.size()));

                auto cost = [this, i, featabledims](const int * varlabels, size_t nvar,
                    ml::FactorGraph::FactorCategoryId fcid, void * givenData) -> double {
                    assert(nvar == featabledims.size());
                    int caseId = Sub2Ind(varlabels, featabledims.data(), nvar);
                    auto weights = static_cast<const double*>(givenData);
                    auto & fea = _patchFeatureTable[i][caseId];
                    assert(fea.size() == FeatureLength(Patch));
                    return std::inner_product(weights + FeatureBeginPos(Patch), weights + FeatureEndPos(Patch), fea.begin(), 0.0);
                };
                ml::FactorGraph::FactorCategory fc;
                fc.costs = std::move(cost);
                fc.c_alpha = 1.0;
                _fg.addFactor(vhs.begin(), vhs.end(), _fg.addFactorCategory(std::move(fc)));
            }

#pragma endregion Build factors for factor graph 
         
        }

        int RLOptimizer::featureLength() const {
            return FullFeatureLength;
        }

        std::vector<double> RLOptimizer::featureInR(RHandle rh, int rlabel) const {
            auto & rd = _g.data(rh);
            static_assert(FeatureLength(R) == 8, "");
            std::vector<double> fea(FeatureLength(R), 0.0);
            double areaRatio = rd.area / _rAreaSum;
            if (Contains(_peakyRhs[0], rh)) { // tend to be horizontal
                if (rlabel == R_Label_ToVP(0)) { // judging as horizontal
                    fea[0] = areaRatio;
                } else if(rlabel == R_Label_Vertical){ // judging as vertical
                    fea[1] = - areaRatio;
                } else if(rlabel == R_Label_NotPlanar){ // ... as not planar
                    fea[2] = - 0.5 * areaRatio;
                } else { // judging as free
                    fea[3] = 0.5 * areaRatio;
                }
            } else if (Contains(_horizonRhs, rh)) { // tend to be vertical
                if (rlabel == R_Label_Vertical) { // juding as vertical
                    fea[4] = areaRatio;
                } else if (rlabel == R_Label_ToVP(0)) { // judging as horizontal
                    fea[5] = -areaRatio;
                } else if (rlabel == R_Label_NotPlanar) { // .. as not planar
                    fea[6] = - 0.5 * areaRatio;
                } else { // judging as free
                    fea[7] = 0.5 * areaRatio;
                }
            } else {

            }
            return fea;
        }

        std::vector<double> RLOptimizer::featureInL(LHandle lh, int llabel) const {
            auto & ld = _g.data(lh);
            static_assert(FeatureLength(L) == 1, "");
            std::vector<double> fea(FeatureLength(L), 0.0);
            double lenRatio = ld.normalizedLine.length() / _lLengthSum;
            double vpScoreSum = std::accumulate(ld.vpScores.begin(), ld.vpScores.end(), 0.0);
            if (llabel == L_Label_Free) {
                return{ 1.0 / ld.vpScores.size() * lenRatio };
            } else {
                return{ (ld.vpScores[VP_From_L_Label(llabel)]) / vpScoreSum * lenRatio };
            }
        }


        std::vector<double> RLOptimizer::featureInRR(RRHandle rrh, int rrlabel) const {
            auto & rrd = _g.data(rrh);
            static_assert(FeatureLength(RR) == 6, "");
            std::vector<double> fea(FeatureLength(RR), 0.0);
            int resp = ToUnderlying(rrd.occDetectionResult);
            int id = resp * 3 + rrlabel;
            double lenRatio = rrd.length / _rrLengthSum;
            fea[id] = lenRatio;
            return fea;
        }

        std::vector<double> RLOptimizer::featureInLL(LLHandle llh, int lllabel) const {
            auto & lld = _g.data(llh);
            static_assert(FeatureLength(LL) == LLData::ManhattanJunctionTypeNum * 2, "");
            auto lh1 = _g.topo(llh).component<0>();
            auto lh2 = _g.topo(llh).component<1>();
            auto & line1 = _g.data(lh1).normalizedLine;
            auto & line2 = _g.data(lh2).normalizedLine;

            std::vector<double> fea(FeatureLength(LL), 0.0);
            fea[ToUnderlying(lld.mjType) * 2 + lllabel] = 1.0 / (_incidenceLLNum + _intersectionLLNum);
            return fea;
        }


        std::vector<double> RLOptimizer::featureInRRRorRRRR(const int * rrlabels, const int * rrascends, int rrnum) const {
            static_assert(FeatureLength(RRRorRRRR) == 5, "");
            assert(rrnum == 3 || rrnum == 4);
            for (int i = 0; i < rrnum; i++) {
                assert(rrlabels[i] >= 0 && rrlabels[i] < RR_Label_Num);
            }
            std::vector<int> labels;
            labels.reserve(rrnum);
            for (int i = 0; i < rrnum; i++) {
                int rrlabel = rrlabels[i];
                if (!rrascends[i]) {
                    if (rrlabel == RR_Label_FirstIsCloser) {
                        rrlabel = RR_Label_SecondIsCloser;
                    } else if (rrlabel == RR_Label_SecondIsCloser) {
                        rrlabel = RR_Label_FirstIsCloser;
                    }
                }
                if (rrlabel != RR_Label_Connected) {
                    labels.push_back(rrlabel);
                }
            }

            int nRRRorRRRR = _g.internalConstraints<RRRData>().size() +
                _g.internalConstraints<RRRRData>().size();
            std::vector<double> fea(FeatureLength(RRRorRRRR), 0.0);
            if (labels.size() == 0) {
                fea[0] = 1.0;
            } else if (labels.size() == 1) {
                fea[1] = 1.0;
            } else if (labels.size() == 2) {
                if (labels[0] == labels[1]) {
                    fea[2] = 1.0;
                } else {
                    fea[3] = 1.0;
                }
            } else if (labels.size() >= 3) {
                fea[4] = 1.0;
            }
            for (auto & f : fea) {
                f /= nRRRorRRRR;
            }
            return fea;
        }



        std::vector<double> RLOptimizer::featureInReconstructedPatch(int patchId, const int * allLabelsInPatch,
            std::unordered_map<RHandle, Plane3> & planes, std::unordered_map<LHandle, Line3> & lines) const {
            
            auto & patch = _patches[patchId];
            // fill labels
            RLGraphPatchDict<int> labels;
            int id = 0;
            for (auto & h : patch.container<RHandle>()) {
                labels[h] = allLabelsInPatch[id++];
            }
            for (auto & h : patch.container<LHandle>()) {
                labels[h] = allLabelsInPatch[id++];
            }
            for (auto & h : patch.container<RRHandle>()) {
                labels[h] = allLabelsInPatch[id++];
            }
            for (auto & h : patch.container<LLHandle>()) {
                labels[h] = allLabelsInPatch[id++];
            }

            int labelsNum = id;
            planes.clear();
            lines.clear();

            // setup solver
            // entities
            ProjectiveSolver solver;
            std::unordered_map<RHandle, int> rh2solverEntId;
            for (auto & c : labels.container<RHandle>()) {
                int label = c.second;
                auto & plane = planes[c.first];
                if (label == 0) { // free
                    rh2solverEntId[c.first] = solver.bindPlaneDoF3(plane);
                } else if (label >= 1 && label <= _vps.size()) { // to vp
                    plane = Plane3(_vps[label - 1], _vps[label - 1]);
                    rh2solverEntId[c.first] = solver.bindPlaneDoF1(plane);
                } else if (label == _vps.size() + 1) { // vertical
                    rh2solverEntId[c.first] = solver.bindPlaneDoF2(plane, up());
                } else if (label == _vps.size() + 2) { // nonplanar                    
                } else {
                    assert(false);
                }
            }
            std::unordered_map<LHandle, int> lh2solverEntId;
            for (auto & c : labels.container<LHandle>()) {
                int label = c.second;
                auto & line = lines[c.first];
                line = _g.data(c.first).normalizedLine;
                if (label == 0) { // free
                    lh2solverEntId[c.first] = solver.bindLineDoF2(line);
                } else if (label >= 1 && label <= _vps.size()) { // to vp
                    Ray3 ray(line.center(), _vps[label - 1]);
                    auto p1 = DistanceBetweenTwoLines(ray, Ray3(Origin(), line.first)).second.first;
                    auto p2 = DistanceBetweenTwoLines(ray, Ray3(Origin(), line.second)).second.second;
                    line = Line3(p1, p2);
                    lh2solverEntId[c.first] = solver.bindLineDoF1(line);
                } else {
                    assert(false);
                }
            }

            // anchors
            static const double depthLB = 1e-2;
            static const double depthUB = 1e2;

            // rr
            for (auto & c : labels.container<RRHandle>()) {
                int label = c.second;
                auto rrh = c.first;
                auto rh1 = _g.topo(rrh).component<0>();
                auto rh2 = _g.topo(rrh).component<1>();
                if (labels.at(rh1) == R_Label_NotPlanar || labels.at(rh2) == R_Label_NotPlanar) { // one of them is nonplanar
                    continue;
                }
                int ent1 = rh2solverEntId.at(rh1);
                int ent2 = rh2solverEntId.at(rh2);
                auto & points = _g.data(rrh).normalizedSampledPoints;
                // register anchors
                std::vector<int> anchorIds;
                for (auto & ps : points) {
                    for (auto & p : ps) {
                        anchorIds.push_back(solver.makeNormalAnchor(p));
                    }
                }
                if (label == 0) { // connected
                    for (int aid : anchorIds) {
                        solver.makeAEqualToBAt(ent1, ent2, aid);
                    }
                } else if (label == 1) { // 1 closer
                    for (int aid : anchorIds) {
                        solver.makeACloserThanBAt(ent1, ent2, aid);
                    }
                } else if (label == 2) { // 2 closer
                    for (int aid : anchorIds) {
                        solver.makeAFartherThanDepthAt(ent1, ent2, aid);
                    }
                } else {
                    assert(false);
                }

                // depth lb/ub
                for (int aid : anchorIds) {
                    solver.makeACloserThanDepthAt(ent1, depthUB, aid);
                    solver.makeAFartherThanDepthAt(ent1, depthLB, aid);
                    solver.makeACloserThanDepthAt(ent2, depthUB, aid);
                    solver.makeAFartherThanDepthAt(ent2, depthLB, aid);
                }
            }

            // ll
            for (auto & c : labels.container<LLHandle>()) {
                int label = c.second;
                auto llh = c.first;
                auto lh1 = _g.topo(llh).component<0>();
                auto lh2 = _g.topo(llh).component<1>();
                int ent1 = lh2solverEntId.at(lh1);
                int ent2 = lh2solverEntId.at(lh2);
                auto & point = _g.data(llh).normalizedRelationCenter;
                // register anchor
                int anchorId = solver.makeNormalAnchor(point);
                if (label == 0) { // connected
                    solver.makeAEqualToBAt(ent1, ent2, anchorId);
                } else if (label == 1) { // disconnected
                } else {
                    assert(false);
                }

                // depth lb/ub
                {
                    int aid = anchorId;
                    solver.makeACloserThanDepthAt(ent1, depthUB, aid);
                    solver.makeAFartherThanDepthAt(ent1, depthLB, aid);
                    solver.makeACloserThanDepthAt(ent2, depthUB, aid);
                    solver.makeAFartherThanDepthAt(ent2, depthLB, aid);
                }
            }

            // rl, rrl
            std::vector<std::pair<RHandle, LHandle>> rlcons;
            std::vector<const std::vector<Vec3>*> nanchorsPtrTable;
            appliedRLCons(patch, labels, rlcons, nanchorsPtrTable);

            assert(rlcons.size() == nanchorsPtrTable.size());

            for (int i = 0; i < rlcons.size(); i++) {
                auto rh = rlcons[i].first;
                auto lh = rlcons[i].second;
                int ent1 = rh2solverEntId.at(rh);
                int ent2 = lh2solverEntId.at(lh);
                auto & as = *nanchorsPtrTable[i];
                for (auto & a : as) {
                    int aid = solver.makeNormalAnchor(a);
                    solver.makeAEqualToBAt(ent1, ent2, aid);

                    solver.makeACloserThanDepthAt(ent1, depthUB, aid);
                    solver.makeAFartherThanDepthAt(ent1, depthLB, aid);
                    solver.makeACloserThanDepthAt(ent2, depthUB, aid);
                    solver.makeAFartherThanDepthAt(ent2, depthLB, aid);
                }
            }

            // how to put arbitrary depth anchors!!!?
            {
                // try longest line
                LHandle longestLh;
                double maxLen = 0.0;
                for (auto & lh : patch.container<LHandle>()) {
                    double len = _g.data(lh).normalizedLine.length();
                    if (len > maxLen) {
                        maxLen = len;
                        longestLh = lh;
                    }
                }
                if (longestLh.valid()) {
                    // put a fixed depth anchor to its center
                    auto center = normalize(_g.data(longestLh).normalizedLine.center());
                    int aid = solver.makeNormalAnchor(center);
                    solver.makeAEqualToDepthAt(lh2solverEntId.at(longestLh), 1.0, aid);
                } else {
                    // or try largest region
                    RHandle largestRh;
                    double maxArea = 0.0;
                    for (auto & rh : patch.container<RHandle>()) {
                        if (labels[rh] == R_Label_NotPlanar) {
                            continue;
                        }
                        double a = _g.data(rh).area;
                        if (a > maxArea) {
                            largestRh = rh;
                            maxArea = a;
                        }
                    }
                    // put a fixed depth anchor to its center
                    auto center = normalize(_g.data(largestRh).normalizedCenter);
                    int aid = solver.makeNormalAnchor(center);
                    solver.makeAEqualToDepthAt(rh2solverEntId.at(largestRh), 1.0, aid);
                }
            }

            // solve
            double nanOrInfRatio = 0.0;
            bool solvable = solver.solve(&nanOrInfRatio);

            std::vector<double> fea(FeatureLength(Patch), 0.0);
            if (!solvable) {
                assert(nanOrInfRatio > 0);
                fea[0] = nanOrInfRatio;
            } else {
                // normalize resulted planes and lines
                double centerDepthsSum = 0.0;
                int centerDepthsNum = 0;
                for (auto & rp : planes) {
                    int label = labels[rp.first];
                    if (label == R_Label_NotPlanar) { // non planar
                        // todo
                        continue;
                    }

                    auto & center = _g.data(rp.first).normalizedCenter;
                    const Plane3 & plane = rp.second;
                    double d = norm(IntersectionOfLineAndPlane(Ray3(Origin(), center), plane).position);
                    centerDepthsSum += d;
                    centerDepthsNum++;
                }
                double centerDepthMean = centerDepthsSum / centerDepthsNum;

                // manhattan region area ratio
                double manhattanRegionArea = 0.0;
                double allRegionArea = 0.0;
                for (auto & p : planes) {
                    int label = labels[p.first];
                    if (label == R_Label_NotPlanar) { // non planar
                        // todo
                        continue;
                    }

                    double area = _g.data(p.first).area;
                    allRegionArea += area;

                    auto & plane = p.second;
                    static const double angleThreshold = DegreesToRadians(5);
                    double angle = M_PI;
                    int nearestVPId = -1;
                    for (int k = 0; k < _vps.size(); k++) {
                        double a = AngleBetweenUndirectedVectors(plane.normal, _vps[k]);
                        if (a < angle) {
                            angle = a;
                            nearestVPId = k;
                        }
                    }
                    if (angle < angleThreshold) {
                        manhattanRegionArea += area;
                    }
                }

                // gc confusion mat
                Mat<double, 2, 3> gcWeightedMeanConfusionMat;
                for (auto & p : planes) {
                    int label = labels[p.first];
                    if (label == R_Label_NotPlanar) { // non planar
                        // todo
                        continue;
                    }

                    double area = _g.data(p.first).area;
                    auto & gcMean = _g.data(p.first).gcMean;

                    auto & plane = p.second;
                    double horizscore = abs(normalize(plane.normal).dot(normalize(up())));
                    double vertscore = 1.0 - horizscore;

                    double gcHoriz = gcMean[ToUnderlying(GeometricContextIndex::FloorOrGround)] +
                        gcMean[ToUnderlying(GeometricContextIndex::CeilingOrSky)];
                    double gcVert = gcMean[ToUnderlying(GeometricContextIndex::Vertical)];
                    double gcClutterOrPorous = gcMean[ToUnderlying(GeometricContextIndex::ClutterOrPorous)];

                    Mat<double, 2, 3> gcConfusionMat;

                    gcConfusionMat(0, 0) = horizscore * gcHoriz;
                    gcConfusionMat(0, 1) = horizscore * gcVert;
                    gcConfusionMat(0, 2) = horizscore * gcClutterOrPorous;

                    gcConfusionMat(1, 0) = vertscore * gcHoriz;
                    gcConfusionMat(1, 1) = vertscore * gcVert;
                    gcConfusionMat(1, 2) = vertscore * gcClutterOrPorous;

                    gcWeightedMeanConfusionMat += (area * gcConfusionMat);
                }
                gcWeightedMeanConfusionMat *= (1.0 / allRegionArea);

                // plane orientation/position confusion mat
                // corresponding to the (Note this!)THREE principle vanishing points
                Mat<double, 3, 3> planeWeightedMeanOPConfusionMat;
                assert(_vps.size() >= 3);
                for (auto & p : planes) {
                    int label = labels[p.first];
                    if (label == R_Label_NotPlanar) { // non planar
                        // todo
                        continue;
                    }

                    double area = _g.data(p.first).area;
                    auto & center = _g.data(p.first).normalizedCenter;

                    auto & plane = p.second;
                    double normalDotsToVPs[3];
                    for (int i = 0; i < 3; i++)
                        normalDotsToVPs[i] = abs(normalize(plane.normal).dot(normalize(_vps[i])));
                    double centerDotsToVPs[3];
                    for (int i = 0; i < 3; i++)
                        centerDotsToVPs[i] = abs(center.dot(normalize(_vps[i])));
                    Mat<double, 3, 3> confusionMat;
                    for (int i = 0; i < 3; i++) {
                        for (int j = 0; j < 3; j++) {
                            confusionMat(i, j) = normalDotsToVPs[i] * centerDotsToVPs[j];
                        }
                    }
                    planeWeightedMeanOPConfusionMat += (area * confusionMat);
                }
                planeWeightedMeanOPConfusionMat *= (1.0 / allRegionArea);


                // region connectivity
                double rrDistanceSum = 0;
                int appliedRRAnchorsNum = 0;
                for (auto & c : labels.container<RRHandle>()) {
                    int label = c.second;
                    auto rrh = c.first;
                    auto rh1 = _g.topo(rrh).component<0>();
                    auto rh2 = _g.topo(rrh).component<1>();
                    if (labels.at(rh1) == R_Label_NotPlanar || labels.at(rh2) == R_Label_NotPlanar) { // one of them is nonplanar
                        continue;
                    }

                    auto & plane1 = planes.at(rh1);
                    auto & plane2 = planes.at(rh2);
                    auto & points = _g.data(rrh).normalizedSampledPoints;
                    for (auto & ps : points) {
                        for (auto & p : ps) {
                            double d1 = norm(IntersectionOfLineAndPlane(Ray3(Origin(), p), plane1).position);
                            double d2 = norm(IntersectionOfLineAndPlane(Ray3(Origin(), p), plane2).position);
                            rrDistanceSum += Gaussian(abs(d1 - d2), 1.0);
                            appliedRRAnchorsNum++;
                        }
                    }
                }
                double rrDistanceMean = rrDistanceSum / appliedRRAnchorsNum;

                // region line connectivity
                double rlDistanceSum = 0.0;
                int appliedRLAnchorsNum = 0;
                for (int i = 0; i < rlcons.size(); i++) {
                    auto rh = rlcons[i].first;
                    auto lh = rlcons[i].second;
                    const Plane3 & plane = planes.at(rh);
                    const Line3 & line = lines.at(lh);
                    auto & anchors = *nanchorsPtrTable[i];
                    for (auto & a : anchors) {
                        double d1 = norm(IntersectionOfLineAndPlane(Ray3(Origin(), a), plane).position);
                        double d2 = norm(DistanceBetweenTwoLines(Ray3(Origin(), a), line.ray()).second.first);
                        rlDistanceSum += Gaussian(abs(d1 - d2), 1.0);
                        appliedRLAnchorsNum++;
                    }
                }
                double rlDistanceMean = rlDistanceSum / appliedRLAnchorsNum;

                int allAnchorsNum = 0;
                for (auto & c : patch.container<RRHandle>()) {
                    for (auto & ps : _g.data(c).normalizedSampledPoints) {
                        allAnchorsNum += ps.size();
                    }
                }
                for (auto & c : patch.container<RLHandle>()) {
                    allAnchorsNum += _g.data(c).normalizedAnchors.size();
                }
                for (auto & c : patch.container<RRLHandle>()) {
                    allAnchorsNum += _g.data(c).normalizedAnchors.size() * 2;
                }

                int appliedAnchorsNum = appliedRRAnchorsNum + appliedRLAnchorsNum;
                int notAppliedAnchors = allAnchorsNum - appliedAnchorsNum;
                assert(notAppliedAnchors > 0);

                fea[0] = 0.0; // 0
                fea[1] = manhattanRegionArea / allRegionArea; // 1
                std::copy(gcWeightedMeanConfusionMat.val, gcWeightedMeanConfusionMat.val + 6, fea.begin() + 2); // 2:7
                std::copy(planeWeightedMeanOPConfusionMat.val, planeWeightedMeanOPConfusionMat.val + 9, fea.begin() + 8); // 8:16
                fea[17] = double(appliedAnchorsNum) / allAnchorsNum; // 17
                fea[18] = double(notAppliedAnchors) / allAnchorsNum; // 18
                fea[19] = rrDistanceMean / centerDepthMean; // 19
                fea[20] = rlDistanceMean / centerDepthMean; // 20
            }

            // todo 
            // consider the size(area of rhs) of patch plz
            return fea;
        }

        void RLOptimizer::findPeakyRhs(double rangeAngle) {
            _peakyRhs.clear();
            _peakyRhs.resize(_vps.size());

            // find peaky regions
            for (auto & r : _g.components<RData>()) {
                auto h = r.topo.hd;

                auto & contours = r.data.normalizedContours;
                double radiusAngle = 0.0;
                for (auto & cs : r.data.normalizedContours) {
                    for (auto & c : cs) {
                        double angle = AngleBetweenDirections(r.data.normalizedCenter, c);
                        if (angle > radiusAngle) {
                            radiusAngle = angle;
                        }
                    }
                }

                bool mayCrossAnyVP = false;
                for (auto & vp : _vps) {
                    double angle = AngleBetweenUndirectedVectors(vp, r.data.normalizedCenter);
                    if (angle < radiusAngle) {
                        mayCrossAnyVP = true;
                        break;
                    }
                }

                if (!mayCrossAnyVP) {
                    continue;
                }

                static const double focal = 100.0;
                auto maskView = PerfectRegionMaskView(r.data.normalizedContours, r.data.normalizedCenter, focal);

                // intersection test
                for (int i = 0; i < _vps.size(); i++) {
                    auto p1 = ToPixel(maskView.camera.toScreen(_vps[i]));
                    auto p2 = ToPixel(maskView.camera.toScreen(-_vps[i]));

                    int dilateSize = focal * rangeAngle;
                    bool intersected = false;
                    for (int x = -dilateSize; x <= dilateSize; x++) {
                        if (intersected)
                            break;
                        for (int y = -dilateSize; y <= dilateSize; y++) {
                            if (intersected)
                                break;
                            auto pp1 = Pixel(p1.x + x, p1.y + y);
                            auto pp2 = Pixel(p2.x + x, p2.y + y);
                            if (Contains(maskView.image, pp1) && maskView.image(pp1)) {
                                _peakyRhs[i].insert(r.topo.hd);
                                intersected = true;
                            } else if (Contains(maskView.image, pp2) && maskView.image(pp2)) {
                                _peakyRhs[i].insert(r.topo.hd);
                                intersected = true;
                            }
                        }
                    }
                }
            }
        }

        void RLOptimizer::findHorizonRhs(double rangeAngle) {
            _horizonRhs.clear();
            for (auto & r : _g.components<RData>()) {
                auto h = r.topo.hd;
                auto & contours = r.data.normalizedContours;
                bool intersected = false;
                for (auto & cs : r.data.normalizedContours) {
                    if (intersected)
                        break;
                    for (auto & c : cs) {
                        double angle = M_PI_2 - AngleBetweenUndirectedVectors(c, up());
                        if (angle <= rangeAngle) {
                            intersected = true;
                            break;
                        }
                    }
                }
                if (intersected) {
                    _horizonRhs.insert(h);
                }
            }
        }


        void RLOptimizer::countIncidenceAndIntersectionLLhs() {
            _incidenceLLNum = 0;
            _intersectionLLNum = 0;
            for (auto & ll : _g.constraints<LLData>()) {
                if (ll.data.mjType == LLData::I || ll.data.mjType == LLData::IAcrossView) {
                    _incidenceLLNum++;
                } else {
                    _intersectionLLNum++;
                }
            }
        }

        void RLOptimizer::appliedRLCons(const RLGraphPatch & patch, const RLGraphPatchDict<int> & cl,
            std::vector<std::pair<RHandle, LHandle>> & rlcons, std::vector<const std::vector<Vec3>*> & nanchorsPtrTable) const {

            rlcons.clear();
            rlcons.reserve(patch.container<RLHandle>().size() + patch.container<RRLHandle>().size());

            nanchorsPtrTable.clear();
            nanchorsPtrTable.reserve(rlcons.capacity());

            for (auto & c : patch.container<RLHandle>()) {
                auto rlh = c;
                auto rh = _g.topo(rlh).component<0>();
                if (cl.at(rh) == R_Label_NotPlanar) { // rh is nonplanar
                    continue;
                }
                auto lh = _g.topo(rlh).component<1>();
                rlcons.emplace_back(rh, lh);
                nanchorsPtrTable.push_back(&_g.data(rlh).normalizedAnchors);
            }

            for (auto & c : patch.container<RRLHandle>()) {
                auto rrlh = c;
                auto rh1 = _g.topo(rrlh).component<0>();
                auto rh2 = _g.topo(rrlh).component<1>();
                auto lh = _g.topo(rrlh).component<2>();
                // find boundary between rh1 and rh2
                RRHandle rrh;
                for (const RRHandle & h : _g.topo(rh1).constraints<RRData>()) {
                    if (_g.topo(h).component<0>() == rh1 && _g.topo(h).component<1>() == rh2) {
                        rrh = h;
                        break;
                    }
                }
                assert(rrh.valid());
                int rlabel1 = cl.at(rh1);
                int rlabel2 = cl.at(rh2);
                bool r1isnonplanar = rlabel1 == R_Label_NotPlanar;
                bool r2isnonplanar = rlabel2 == R_Label_NotPlanar;
                if (r1isnonplanar && r2isnonplanar) { // all of them is nonplanar
                    continue;
                }

                int rrlabel = cl.at(rrh);
                bool connect_rh1_lh = (rrlabel == RR_Label_Connected || rrlabel == RR_Label_FirstIsCloser) && !r1isnonplanar;
                bool connect_rh2_lh = (rrlabel == RR_Label_Connected || rrlabel == RR_Label_SecondIsCloser) && !r2isnonplanar;

                if (connect_rh1_lh) {
                    rlcons.emplace_back(rh1, lh);
                    nanchorsPtrTable.push_back(&_g.data(rrlh).normalizedAnchors);
                }
                if (connect_rh2_lh) {
                    rlcons.emplace_back(rh2, lh);
                    nanchorsPtrTable.push_back(&_g.data(rrlh).normalizedAnchors);
                }
            }

        }


        std::vector<double> RLOptimizer::suggestedWeights() const {
            NOT_IMPLEMENTED_YET();
        }



        void RLOptimizer::inference(const std::vector<double> & weights,
            HandledTable<RHandle, Plane3> & planes, HandledTable<LHandle, Line3> & lines) const {
            assert(weights.size() == featureLength());
            auto predictedLabels = _fg.solveWithSimpleCallback(1e2, 10, nullptr, (void*)weights.data());

            
            // whole graph as a cluster

            // todo

            // todo

        }

        void RLOptimizer::inferenceWithLossFunction(const std::vector<double> & weights, 
            std::function<double(const HandledTable<RHandle, Plane3> & planes, const HandledTable<LHandle, Line3> & lines)> & lossFun) const {


        }

    }
}