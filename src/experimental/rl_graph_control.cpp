#include "../core/algorithms.hpp"
#include "../core/clock.hpp"
#include "../core/containers.hpp"

#include "rl_graph_control.hpp"

namespace pano {
    namespace experimental {



        int ConnectedComponents(const RLGraph & mg, RLGraphComponentTable<int> & ccIds) {
            std::vector<RLGraph> ccs;

            struct HComp {
                enum Type { Region, Line };
                Type type;
                int id;
            };
            struct HCons {
                enum Type { RegionBoundary, LineRelation, RegionLine };
                Type type;
                int id;
            };

            std::vector<HComp> comps;
            std::vector<HCons> conss;
            comps.reserve(mg.allComponentsNum());
            conss.reserve(mg.allConstraintsNum());

            auto comp2htable = MakeHandledTableForAllComponents<size_t>(mg);
            auto cons2htable = MakeHandledTableForAllConstraints<size_t>(mg);

            for (auto & r : mg.components<RegionData>()) {
                comp2htable[r.topo.hd] = comps.size();
                comps.push_back({ HComp::Region, r.topo.hd.id });
            }
            for (auto & l : mg.components<LineData>()) {
                comp2htable[l.topo.hd] = comps.size();
                comps.push_back({ HComp::Line, l.topo.hd.id });
            }

            for (auto & c : mg.constraints<RegionBoundaryData>()) {
                cons2htable[c.topo.hd] = conss.size();
                conss.push_back({ HCons::RegionBoundary, c.topo.hd.id });
            }
            for (auto & c : mg.constraints<LineRelationData>()) {
                cons2htable[c.topo.hd] = conss.size();
                conss.push_back({ HCons::LineRelation, c.topo.hd.id });
            }
            for (auto & c : mg.constraints<RegionLineConnectionData>()) {
                cons2htable[c.topo.hd] = conss.size();
                conss.push_back({ HCons::RegionLine, c.topo.hd.id });
            }

            std::vector<size_t> compIds(comps.size());
            std::iota(compIds.begin(), compIds.end(), 0ull);

            auto getNeighborComps = [&comps, &mg, &comp2htable](size_t compId) {
                auto & comp = comps[compId];
                std::vector<size_t> neighbors;
                if (comp.type == HComp::Region) {
                    RegionHandle rh(comp.id);
                    for (auto & h : mg.topo(rh).constraints<RegionBoundaryData>()) {
                        auto another = mg.topo(h).component<0>();
                        if (another == rh) {
                            another = mg.topo(h).component<1>();
                        }
                        neighbors.push_back(comp2htable[another]);
                    }
                    for (auto & h : mg.topo(rh).constraints<RegionLineConnectionData>()) {
                        LineHandle another = mg.topo(h).component<1>();
                        neighbors.push_back(comp2htable[another]);
                    }
                } else {
                    LineHandle lh(comp.id);
                    for (auto & h : mg.topo(lh).constraints<LineRelationData>()) {
                        auto another = mg.topo(h).component<0>();
                        if (another == lh) {
                            another = mg.topo(h).component<1>();
                        }
                        neighbors.push_back(comp2htable[another]);
                    }
                    for (auto & h : mg.topo(lh).constraints<RegionLineConnectionData>()) {
                        RegionHandle another = mg.topo(h).component<0>();
                        neighbors.push_back(comp2htable[another]);
                    }
                }
                return neighbors;
            };

            return core::ConnectedComponents(compIds.begin(), compIds.end(), std::move(getNeighborComps),
                [&ccIds, &comps](size_t v, int ccid) {
                const HComp & hcomp = comps[v];
                if (hcomp.type == HComp::Line) {
                    ccIds[LineHandle(hcomp.id)] = ccid;
                } else {
                    ccIds[RegionHandle(hcomp.id)] = ccid;
                }
            });
        }

        int ConnectedComponents(const RLGraph & mg, const RLGraphControls & controls,
            RLGraphComponentTable<int> & ccIds,
            const std::function<bool(const RLGraphConstraintControl &)> & constraintAsConnected) {

            assert(constraintAsConnected);

            std::vector<RLGraph> ccs;

            struct HComp {
                enum Type { Region, Line };
                Type type;
                int id;
            };
            struct HCons {
                enum Type { RegionBoundary, LineRelation, RegionLine };
                Type type;
                int id;
            };

            std::vector<HComp> comps;
            std::vector<HCons> conss;
            comps.reserve(mg.allComponentsNum());
            conss.reserve(mg.allConstraintsNum());

            auto comp2htable = MakeHandledTableForAllComponents<size_t>(mg);
            auto cons2htable = MakeHandledTableForAllConstraints<size_t>(mg);

            for (auto & r : mg.components<RegionData>()) {
                comp2htable[r.topo.hd] = comps.size();
                comps.push_back({ HComp::Region, r.topo.hd.id });
            }
            for (auto & l : mg.components<LineData>()) {
                comp2htable[l.topo.hd] = comps.size();
                comps.push_back({ HComp::Line, l.topo.hd.id });
            }

            for (auto & c : mg.constraints<RegionBoundaryData>()) {
                cons2htable[c.topo.hd] = conss.size();
                conss.push_back({ HCons::RegionBoundary, c.topo.hd.id });
            }
            for (auto & c : mg.constraints<LineRelationData>()) {
                cons2htable[c.topo.hd] = conss.size();
                conss.push_back({ HCons::LineRelation, c.topo.hd.id });
            }
            for (auto & c : mg.constraints<RegionLineConnectionData>()) {
                cons2htable[c.topo.hd] = conss.size();
                conss.push_back({ HCons::RegionLine, c.topo.hd.id });
            }

            std::vector<size_t> compIds(comps.size());
            std::iota(compIds.begin(), compIds.end(), 0ull);

            auto getNeighborComps = [&comps, &mg, &comp2htable, &constraintAsConnected, &controls](size_t compId) {
                auto & comp = comps[compId];
                std::vector<size_t> neighbors;
                if (comp.type == HComp::Region) {
                    RegionHandle rh(comp.id);
                    for (auto & h : mg.topo(rh).constraints<RegionBoundaryData>()) {
                        auto & consControl = controls[h];
                        if (!constraintAsConnected(consControl))
                            continue;
                        auto another = mg.topo(h).component<0>();
                        if (another == rh) {
                            another = mg.topo(h).component<1>();
                        }
                        neighbors.push_back(comp2htable[another]);
                    }
                    for (auto & h : mg.topo(rh).constraints<RegionLineConnectionData>()) {
                        auto & consControl = controls[h];
                        if (!constraintAsConnected(consControl))
                            continue;
                        LineHandle another = mg.topo(h).component<1>();
                        neighbors.push_back(comp2htable[another]);
                    }
                } else {
                    LineHandle lh(comp.id);
                    for (auto & h : mg.topo(lh).constraints<LineRelationData>()) {
                        auto & consControl = controls[h];
                        if (!constraintAsConnected(consControl))
                            continue;
                        auto another = mg.topo(h).component<0>();
                        if (another == lh) {
                            another = mg.topo(h).component<1>();
                        }
                        neighbors.push_back(comp2htable[another]);
                    }
                    for (auto & h : mg.topo(lh).constraints<RegionLineConnectionData>()) {
                        auto & consControl = controls[h];
                        if (!constraintAsConnected(consControl))
                            continue;
                        RegionHandle another = mg.topo(h).component<0>();
                        neighbors.push_back(comp2htable[another]);
                    }
                }
                return neighbors;
            };

            return core::ConnectedComponents(compIds.begin(), compIds.end(), std::move(getNeighborComps),
                [&ccIds, &comps](size_t v, int ccid) {
                const HComp & hcomp = comps[v];
                if (hcomp.type == HComp::Line) {
                    ccIds[LineHandle(hcomp.id)] = ccid;
                } else {
                    ccIds[RegionHandle(hcomp.id)] = ccid;
                }
            });
        }


        template <class T1, class T2>
        std::pair<T2, T1> Inversed(const std::pair<T1, T2> & p) {
            return std::make_pair(p.second, p.first);
        }

        std::vector<RLGraph> Decompose(const RLGraph & mg, const RLGraphComponentTable<int> & ccids, int ccnum,
            RLGraphOldToNew * old2new, RLGraphNewToOld * new2old) {
            std::vector<RLGraph> ccs(ccnum);

            RLGraphComponentTable<int> newCompIds = MakeHandledTableForAllComponents<int>(mg);
            for (auto & c : mg.components<RegionData>()) {
                int ccid = ccids[c.topo.hd];
                auto h = ccs[ccid].addComponent(c.data);
                newCompIds[c.topo.hd] = h.id;
                auto pair = std::make_pair(c.topo.hd, std::make_pair(ccid, h));
                if (old2new) {
                    (*old2new)[pair.first] = pair.second;
                }
                if (new2old) {
                    (*new2old)[pair.second] = pair.first;
                }
            }
            for (auto & c : mg.components<LineData>()) {
                int ccid = ccids[c.topo.hd];
                auto h = ccs[ccid].addComponent(c.data);
                newCompIds[c.topo.hd] = h.id;
                auto pair = std::make_pair(c.topo.hd, std::make_pair(ccid, h));
                if (old2new) {
                    (*old2new)[pair.first] = pair.second;
                }
                if (new2old) {
                    (*new2old)[pair.second] = pair.first;
                }
            }

            for (auto & c : mg.constraints<RegionBoundaryData>()) {
                int ccid1 = ccids[c.topo.component<0>()];
                int ccid2 = ccids[c.topo.component<1>()];
                if (ccid1 != ccid2) {
                    continue;
                }
                int ccid = ccid1;
                auto h = ccs[ccid].addConstraint(c.data,
                    RegionHandle(newCompIds[c.topo.component<0>()]),
                    RegionHandle(newCompIds[c.topo.component<1>()]));
                auto pair = std::make_pair(c.topo.hd, std::make_pair(ccid, h));
                if (old2new) {
                    (*old2new)[pair.first] = pair.second;
                }
                if (new2old) {
                    (*new2old)[pair.second] = pair.first;
                }
            }
            for (auto & c : mg.constraints<LineRelationData>()) {
                int ccid1 = ccids[c.topo.component<0>()];
                int ccid2 = ccids[c.topo.component<1>()];
                if (ccid1 != ccid2) {
                    continue;
                }
                int ccid = ccid1;
                auto h = ccs[ccid].addConstraint(c.data,
                    LineHandle(newCompIds[c.topo.component<0>()]),
                    LineHandle(newCompIds[c.topo.component<1>()]));
                auto pair = std::make_pair(c.topo.hd, std::make_pair(ccid, h));
                if (old2new) {
                    (*old2new)[pair.first] = pair.second;
                }
                if (new2old) {
                    (*new2old)[pair.second] = pair.first;
                }
            }
            for (auto & c : mg.constraints<RegionLineConnectionData>()) {
                int ccid1 = ccids[c.topo.component<0>()];
                int ccid2 = ccids[c.topo.component<1>()];
                if (ccid1 != ccid2) {
                    continue;
                }
                int ccid = ccid1;
                auto h = ccs[ccid].addConstraint(c.data,
                    RegionHandle(newCompIds[c.topo.component<0>()]),
                    LineHandle(newCompIds[c.topo.component<1>()]));
                auto pair = std::make_pair(c.topo.hd, std::make_pair(ccid, h));
                if (old2new) {
                    (*old2new)[pair.first] = pair.second;
                }
                if (new2old) {
                    (*new2old)[pair.second] = pair.first;
                }
            }

            return ccs;
        }





        namespace {

            //template <class FunctorT>
            //struct DataCostFunctorWrapper : GCoptimization::DataCostFunctor{
            //    inline DataCostFunctorWrapper(FunctorT && f) : fun(std::forward<FunctorT>(f)) {}
            //    virtual GCoptimization::EnergyTermType compute(GCoptimization::SiteID s, GCoptimization::LabelID l) override {
            //        return fun(s, l);
            //    }
            //    FunctorT fun;
            //};
            //template <class FunctorT>
            //inline DataCostFunctorWrapper<FunctorT> * AllocDataCostFunctor(FunctorT && f) {
            //    return new DataCostFunctorWrapper<FunctorT>(std::forward<FunctorT>(f));
            //}

            //template <class FunctorT>
            //struct SmoothCostFunctorWrapper : GCoptimization::SmoothCostFunctor {
            //    inline SmoothCostFunctorWrapper(FunctorT && f) : fun(std::forward<FunctorT>(f)){}
            //    virtual GCoptimization::EnergyTermType compute(
            //        GCoptimization::SiteID s1, GCoptimization::SiteID s2,
            //        GCoptimization::LabelID l1, GCoptimization::LabelID l2) override {
            //        return fun(s1, s2, l1, l2);
            //    }
            //    FunctorT fun;
            //};
            //template <class FunctorT>
            //inline SmoothCostFunctorWrapper<FunctorT> * AllocSmoothCostFunctor(FunctorT && f) {
            //    return new SmoothCostFunctorWrapper<FunctorT>(std::forward<FunctorT>(f));
            //}


            namespace {

                int SwappedComponent(const Vec3 & orientation) {
                    for (int i = 0; i < 2; i++) {
                        if (abs(orientation[i]) >= 1e-8) {
                            return i;
                        }
                    }
                    return 2;
                }

                struct InitializeVariablesForEachHandle {
                    const RLGraph & mg;
                    const RLGraphControls & controls;
                    RLGraphVars & vars;

                    void operator()(const LineHandle & lh) const {
                        auto & lp = controls[lh];
                        auto & v = vars[lh];
                        if (!lp.used) {
                            v.variables = {};
                            return;
                        }
                        if (lp.orientationClaz == -1) {
                            // (1/cornerDepth1, 1/cornerDepth2) for LineFree
                            v.variables = { 1.0, 1.0 };
                        } else {
                            // 1/centerDepth for LineOriented,
                            v.variables = { 1.0 };
                        }
                    }

                    void operator()(const RegionHandle & rh) const {
                        auto & rd = mg.data(rh);
                        auto & rp = controls[rh];
                        auto & v = vars[rh];
                        if (!rp.used) {
                            v.variables = {};
                            return;
                        }
                        if (rp.orientationClaz == -1 && rp.orientationNotClaz == -1) {
                            // (a, b, c) for RegionFree ax+by+c=1, 
                            v.variables = { rd.normalizedCenter[0], rd.normalizedCenter[1], rd.normalizedCenter[2] };
                        } else if (rp.orientationClaz == -1 && rp.orientationNotClaz >= 0) {
                            // (a, b), {or (b, c) or (a, c)} for RegionAlongFixedAxis  ax+by+c=1,
                            Vec3 anotherAxis = rd.normalizedCenter.cross(normalize(controls.vanishingPoints[rp.orientationNotClaz]));
                            Vec3 trueNormal = controls.vanishingPoints[rp.orientationNotClaz].cross(anotherAxis);
                            Plane3 plane(rd.normalizedCenter, trueNormal);
                            auto eq = Plane3ToEquation(plane);
                            int c = SwappedComponent(normalize(controls.vanishingPoints[rp.orientationNotClaz]));
                            std::swap(eq[c], eq[2]);
                            v.variables = { eq[0], eq[1] };
                        } else {
                            // 1/centerDepth for RegionWithFixedNormal
                            v.variables = { 1.0 };
                        }
                    }
                };




            }

        }





        void RLGraphControls::disable(RegionHandle h, const RLGraph & mg) {
            componentControls[h].used = false;
            for (auto ch : mg.topo(h).constraints<RegionBoundaryData>()) {
                constraintControls[ch].used = false;
            }
            for (auto ch : mg.topo(h).constraints<RegionLineConnectionData>()) {
                constraintControls[ch].used = false;
            }
        }

        void RLGraphControls::disable(LineHandle h, const RLGraph & mg) {
            componentControls[h].used = false;
            for (auto ch : mg.topo(h).constraints<LineRelationData>()) {
                constraintControls[ch].used = false;
            }
            for (auto ch : mg.topo(h).constraints<RegionLineConnectionData>()) {
                constraintControls[ch].used = false;
            }
        }

        void RLGraphControls::enable(RegionHandle h, const RLGraph & mg) {
            componentControls[h].used = true;
        }

        void RLGraphControls::enable(LineHandle h, const RLGraph & mg) {
            componentControls[h].used = true;
        }


        void RLGraphControls::disableAllInvalidConstraints(const RLGraph & mg) {
            for (auto & c : mg.constraints<RegionBoundaryData>()) {
                if (!componentControls[c.topo.component<0>()].used ||
                    !componentControls[c.topo.component<1>()].used) {
                    constraintControls[c.topo.hd].used = false;
                }
            }
            for (auto & c : mg.constraints<LineRelationData>()) {
                if (!componentControls[c.topo.component<0>()].used ||
                    !componentControls[c.topo.component<1>()].used) {
                    constraintControls[c.topo.hd].used = false;
                }
            }
            for (auto & c : mg.constraints<RegionLineConnectionData>()) {
                if (!componentControls[c.topo.component<0>()].used ||
                    !componentControls[c.topo.component<1>()].used) {
                    constraintControls[c.topo.hd].used = false;
                }
            }
        }


        void RLGraphControls::enableAll() {
            for (auto & ct : componentControls.data) {
                for (auto & c : ct) {
                    c.used = true;
                }
            }
            for (auto & ct : constraintControls.data) {
                for (auto & c : ct) {
                    c.used = true;
                }
            }
        }


        RLGraphControls::RLGraphControls(const RLGraph & mg, const std::vector<Vec3> & vps) : vanishingPoints(vps) {
            componentControls = MakeHandledTableForAllComponents<RLGraphComponentControl>(mg);
            constraintControls = MakeHandledTableForAllConstraints<RLGraphConstraintControl>(mg);
            for (auto & l : mg.internalComponents<LineData>()) {
                componentControls[l.topo.hd].orientationClaz = l.data.initialClaz;
                componentControls[l.topo.hd].orientationNotClaz = -1;
            }
            for (auto & r : mg.internalComponents<RegionData>()) {
                componentControls[r.topo.hd].orientationClaz = componentControls[r.topo.hd].orientationNotClaz = -1;
            }
            for (auto & dtable : componentControls.data) {
                for (auto & d : dtable) {
                    d.used = true;
                }
            }
            for (auto & dtable : constraintControls.data) {
                for (auto & d : dtable) {
                    d.used = true;
                }
            }
            SetNecessaryConstraintWeightedAnchors(mg, *this);
        }










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

            struct ExtractNecessaryAnchorsForBinary {
                std::vector<Vec3> operator()(const RLGraph & mg, RegionBoundaryHandle bh) const {
                    size_t n = ElementsNum(mg.data(bh).normalizedSampledPoints);

                    assert(n > 0);
                    const auto & points = mg.data(bh).normalizedSampledPoints;

                    if (n == 1) {
                        return{ points.front().front() };
                    }

                    const Vec3 * pp1 = nullptr;
                    const Vec3 * pp2 = nullptr;
                    double maxAngle = -1;
                    for (int i = 0; i < points.size(); i++) {
                        for (int ii = 0; ii < points[i].size(); ii++) {
                            for (int j = i; j < points.size(); j++) {
                                for (int jj = 0; jj < points[j].size(); jj++) {
                                    if (std::tie(i, ii) >= std::tie(j, jj))
                                        continue;
                                    double angle = AngleBetweenDirections(points[i][ii], points[j][jj]);
                                    if (angle > maxAngle) {
                                        maxAngle = angle;
                                        pp1 = &points[i][ii];
                                        pp2 = &points[j][jj];
                                    }
                                }
                            }
                        }
                    }
                    auto & p1 = *pp1;
                    auto & p2 = *pp2;

                    if (n == 2) {
                        return{ p1, p2 };
                    }

                    auto normal12 = p1.cross(p2);
                    auto pp3 = &p1;
                    for (auto & ps : points) {
                        for (auto & p : ps) {
                            if (abs(p.dot(normal12)) > abs(pp3->dot(normal12))) {
                                pp3 = &p;
                            }
                        }
                    }
                    auto & p3 = *pp3;
                    IMPROVABLE_HERE(? );
                    return{ p1, p2, p3 };
                }

                inline std::vector<Vec3> operator()(const RLGraph & mg, LineRelationHandle bh) const {
                    return{ mg.data(bh).normalizedRelationCenter };
                }

                inline std::vector<Vec3> operator()(const RLGraph & mg, RegionLineConnectionHandle bh) const {
                    return{ mg.data(bh).normalizedAnchors.front(), mg.data(bh).normalizedAnchors.back() };
                }
            };

            struct ExtractAllAnchorsForBinary {
                std::vector<Vec3> operator()(const RLGraph & mg, RegionBoundaryHandle bh) const {
                    std::vector<Vec3> anchors;
                    anchors.reserve(ElementsNum(mg.data(bh).normalizedSampledPoints));
                    for (auto & ps : mg.data(bh).normalizedSampledPoints) {
                        for (auto & p : ps) {
                            anchors.push_back(p);
                        }
                    }
                    return anchors;
                }
                inline std::vector<Vec3> operator()(const RLGraph & mg, LineRelationHandle bh) const {
                    return{ mg.data(bh).normalizedRelationCenter };
                }
                inline std::vector<Vec3> operator()(const RLGraph & mg, RegionLineConnectionHandle bh) const {
                    return mg.data(bh).normalizedAnchors;
                }
            };

            struct ExtractMoreAnchorsForBinary {
                std::vector<Vec3> operator()(const RLGraph & mg, RegionBoundaryHandle bh) const {
                    std::vector<Vec3> anchors;
                    anchors.reserve(ElementsNum(mg.data(bh).normalizedSampledPoints));
                    for (auto & ps : mg.data(bh).normalizedSampledPoints) {
                        if (ps.size() == 1) {
                            anchors.push_back(ps[0]);
                            continue;
                        }
                        anchors.push_back(ps[0]);
                        for (int i = 1; i < ps.size(); i++) {
                            anchors.push_back(ps[i]);
                            Vec3 pdir = normalize(ps[i].cross(ps[i - 1]));
                            auto newp = RotateDirection(ps[i], ps[i] + pdir, expandAngle);
                            anchors.push_back(normalize(newp));
                            newp = RotateDirection(ps[i], ps[i] + pdir, -expandAngle);
                            anchors.push_back(normalize(newp));
                        }
                    }
                    return anchors;
                }
                inline std::vector<Vec3> operator()(const RLGraph & mg, LineRelationHandle bh) const {
                    return{ mg.data(bh).normalizedRelationCenter };
                }
                inline std::vector<Vec3> operator()(const RLGraph & mg, RegionLineConnectionHandle bh) const {
                    std::vector<Vec3> anchors;
                    anchors.reserve(mg.data(bh).normalizedAnchors.size() * 3);
                    auto & ps = mg.data(bh).normalizedAnchors;
                    if (ps.size() == 1)
                        return ps;
                    Vec3 pdir = normalize(ps.front().cross(ps.back()));
                    for (auto & p : ps) {
                        anchors.push_back(p);
                        anchors.push_back(RotateDirection(p, p + pdir, expandAngle));
                        anchors.push_back(RotateDirection(p, p + pdir, -expandAngle));
                    }
                    return anchors;
                }
                double expandAngle;
            };

        }


        template <class ConstraintAnchorExtractorT>
        void SetConstraintWeightedAnchors(const RLGraph & mg, RLGraphControls & controls, ConstraintAnchorExtractorT && extractor) {
            static const bool constantWeightForEachAnchor = true;
            for (auto & c : mg.constraints<RegionBoundaryData>()) {
                auto & control = controls[c.topo.hd];
                if (!control.used) {
                    continue;
                }
                auto as = extractor(mg, c.topo.hd);
                assert(as.size() >= 1);
                double fullWeight = c.data.length / M_PI * 20.0;
                assert(fullWeight >= 0);
                control.weightedAnchors.clear();
                control.weightedAnchors.reserve(as.size());
                for (const Point3 & anchor : as) {
                    control.weightedAnchors.push_back(WeightAs(anchor, constantWeightForEachAnchor ? 1.0 : sqrt(fullWeight / std::max(as.size(), 1ull))));
                }
            }
            for (auto & c : mg.constraints<LineRelationData>()) {
                auto & control = controls[c.topo.hd];
                if (!control.used) {
                    continue;
                }
                auto as = extractor(mg, c.topo.hd);
                assert(as.size() >= 1);
                double fullWeight = c.data.junctionWeight;
                assert(fullWeight >= 0);
                control.weightedAnchors.clear();
                control.weightedAnchors.reserve(as.size());
                for (const Point3 & anchor : as) {
                    control.weightedAnchors.push_back(WeightAs(anchor, constantWeightForEachAnchor ? 1.0 : sqrt(fullWeight / std::max(as.size(), 1ull))));
                }
            }
            for (auto & c : mg.constraints<RegionLineConnectionData>()) {
                auto & control = controls[c.topo.hd];
                if (!control.used) {
                    continue;
                }
                auto as = extractor(mg, c.topo.hd);
                assert(as.size() >= 1);
                double fullWeight = c.data.length / M_PI * 20.0;
                assert(fullWeight >= 0);
                control.weightedAnchors.clear();
                control.weightedAnchors.reserve(as.size());
                for (const Point3 & anchor : as) {
                    control.weightedAnchors.push_back(WeightAs(anchor, constantWeightForEachAnchor ? 1.0 : sqrt(fullWeight / std::max(as.size(), 1ull))));
                }
            }
        }


        void SetNecessaryConstraintWeightedAnchors(const RLGraph & mg, RLGraphControls & controls) {
            SetConstraintWeightedAnchors(mg, controls, ExtractNecessaryAnchorsForBinary());
        }

        void SetFullConstraintWeightedAnchors(const RLGraph & mg, RLGraphControls & controls) {
            SetConstraintWeightedAnchors(mg, controls, ExtractAllAnchorsForBinary());
        }




        std::vector<RLGraphControls> Decompose(const RLGraph & mg,
            const RLGraphControls & controls, const RLGraphComponentTable<int> & ccids, int ccnum) {
            std::vector<RLGraphControls> ccs(ccnum);
            RLGraphComponentTable<int> newCompIds = MakeHandledTableForAllComponents<int>(mg);
            for (auto & cc : ccs) {
                cc.vanishingPoints = controls.vanishingPoints;
            }
            for (auto & c : mg.components<RegionData>()) {
                int ccid = ccids[c.topo.hd];
                auto & regionControls = ccs[ccid].componentControls.dataOfType<RegionHandle>();
                newCompIds[c.topo.hd] = regionControls.size();
                regionControls.push_back(controls[c.topo.hd]);
            }
            for (auto & c : mg.components<LineData>()) {
                int ccid = ccids[c.topo.hd];
                auto & lineControls = ccs[ccid].componentControls.dataOfType<LineHandle>();
                newCompIds[c.topo.hd] = lineControls.size();
                lineControls.push_back(controls[c.topo.hd]);
            }
            for (auto & c : mg.constraints<RegionBoundaryData>()) {
                int ccid1 = ccids[c.topo.component<0>()];
                int ccid2 = ccids[c.topo.component<1>()];
                if (ccid1 != ccid2) {
                    continue;
                }
                int ccid = ccid1;
                auto & ccontrols = ccs[ccid].constraintControls.dataOfType<RegionBoundaryHandle>();
                ccontrols.push_back(controls[c.topo.hd]);
            }
            for (auto & c : mg.constraints<LineRelationData>()) {
                int ccid1 = ccids[c.topo.component<0>()];
                int ccid2 = ccids[c.topo.component<1>()];
                if (ccid1 != ccid2) {
                    continue;
                }
                int ccid = ccid1;
                auto & ccontrols = ccs[ccid].constraintControls.dataOfType<LineRelationHandle>();
                ccontrols.push_back(controls[c.topo.hd]);
            }
            for (auto & c : mg.constraints<RegionLineConnectionData>()) {
                int ccid1 = ccids[c.topo.component<0>()];
                int ccid2 = ccids[c.topo.component<1>()];
                if (ccid1 != ccid2) {
                    continue;
                }
                int ccid = ccid1;
                auto & ccontrols = ccs[ccid].constraintControls.dataOfType<RegionLineConnectionHandle>();
                ccontrols.push_back(controls[c.topo.hd]);
            }
            return ccs;
        }





        template <class FunT>
        inline void ForeachRLGraphComponentHandle(const RLGraph & mg, FunT && fun) {
            for (auto & c : mg.components<LineData>()) {
                fun(c.topo.hd);
            }
            for (auto & c : mg.components<RegionData>()) {
                fun(c.topo.hd);
            }
        }

        template <class FunT>
        inline void ForeachRLGraphConstraintHandle(const RLGraph & mg, FunT && fun) {
            for (auto & c : mg.constraints<RegionBoundaryData>()) {
                fun(c.topo.hd);
            }
            for (auto & c : mg.constraints<LineRelationData>()) {
                fun(c.topo.hd);
            }
            for (auto & c : mg.constraints<RegionLineConnectionData>()) {
                fun(c.topo.hd);
            }
        }




        RLGraphVars MakeVariables(const RLGraph & mg, const RLGraphControls & controls, bool randomized) {
            RLGraphVars vars = MakeHandledTableForAllComponents<RLGraphVar>(mg);
            ForeachRLGraphComponentHandle(mg, InitializeVariablesForEachHandle{ mg, controls, vars });
            if (randomized) {
                for (auto & t : vars.data) {
                    for (auto & vv : t) {
                        for (auto & v : vv.variables) {
                            v = ((size_t)std::rand() % 1000 + 1) / 5000.0 + 1.0;
                        }
                    }
                }
            }
            return vars;
        }

        bool HasInfOrNaNValue(const RLGraphVars & v) {
            for (auto & vt : v.data) {
                for (auto & var : vt) {
                    for (auto & v : var.variables) {
                        if (IsInfOrNaN(v))
                            return true;
                    }
                }
            }
            return false;
        }



        namespace {

            inline double NonZeroize(double d) {
                return d == 0.0 ? 1e-6 : d;
            }

        }





        Line3 InstanceGivenVariables(const RLGraph & mg, const double * variables,
            const RLGraphControls & controls, const LineHandle & lh) {
            auto & ld = mg.data(lh);
            auto & c = controls[lh];
            //if (!lp.used)
            //    return Line3();
            if (c.orientationClaz >= 0) {
                Ray3 infLine(normalize(ld.line.center()) / NonZeroize(variables[0]), controls.vanishingPoints[c.orientationClaz]);
                return Line3(DistanceBetweenTwoLines(Ray3(Point3(0, 0, 0), normalize(ld.line.first)), infLine).second.second,
                    DistanceBetweenTwoLines(Ray3(Point3(0, 0, 0), normalize(ld.line.second)), infLine).second.second);
            } else /*if (line.type == MGUnary::LineFree)*/ {
                /*           | sin(theta) | | p | | q |
                len:---------------------------------   => 1/len: [(1/q)sin(phi) - (1/p)sin(phi-theta)] / sin(theta)
                | p sin(phi) - q sin(phi - theta) |
                */
                // variables[0] -> 1/p
                // variables[1] -> 1/q
                return Line3(normalize(ld.line.first) / NonZeroize(variables[0]), normalize(ld.line.second) / NonZeroize(variables[1]));
            }
        }


        Plane3 InstanceGivenVariables(const RLGraph & mg, const double * variables,
            const RLGraphControls & controls, const RegionHandle & rh) {
            auto & rd = mg.data(rh);
            auto & c = controls[rh];
            /*if (!rp.used)
            return Plane3();*/
            if (c.orientationClaz >= 0) {
                return Plane3(rd.normalizedCenter / NonZeroize(variables[0]), controls.vanishingPoints[c.orientationClaz]);
            } else if (c.orientationClaz == -1 && c.orientationNotClaz >= 0) {
                double vs[] = { variables[0], variables[1], 0.0 }; // fake vs
                // v1 * o1 + v2 * o2 + v3 * o3 = 0, since region.orientation is orthogonal to normal
                auto orientation = normalize(controls.vanishingPoints[c.orientationNotClaz]);
                int c = SwappedComponent(orientation);
                std::swap(orientation[c], orientation[2]); // now fake orientation
                vs[2] = (-vs[0] * orientation[0] - vs[1] * orientation[1])
                    / orientation[2];
                std::swap(vs[c], vs[2]); // now real vs
                return Plane3FromEquation(vs[0], vs[1], vs[2]);
            } else /*if (region.type == MGUnary::RegionWithFixedNormal)*/ {
                return Plane3FromEquation(variables[0], variables[1], variables[2]);
            }
        }


        Line3 Instance(const RLGraph & mg, const RLGraphControls & controls, const RLGraphVars & vars, const LineHandle & lh) {
            //auto & ld = mg.data(lh);
            //auto & c = controls[lh];
            //auto & v = vars[lh];
            ////if (!lp.used)
            ////    return Line3();
            //if (c.orientationClaz >= 0){
            //    assert(v.variables.size() == 1);
            //    Ray3 infLine(normalize(ld.line.center()) / NonZeroize(v.variables[0]), controls.vanishingPoints[c.orientationClaz]);
            //    return Line3(DistanceBetweenTwoLines(Ray3(Point3(0, 0, 0), normalize(ld.line.first)), infLine).second.second,
            //        DistanceBetweenTwoLines(Ray3(Point3(0, 0, 0), normalize(ld.line.second)), infLine).second.second);
            //}
            //else /*if (line.type == MGUnary::LineFree)*/{
            //    assert(v.variables.size() == 2);
            //    /*           | sin(theta) | | p | | q |
            //    len:---------------------------------   => 1/len: [(1/q)sin(phi) - (1/p)sin(phi-theta)] / sin(theta)
            //    | p sin(phi) - q sin(phi - theta) |
            //    */
            //    // variables[0] -> 1/p
            //    // variables[1] -> 1/q
            //    return Line3(normalize(ld.line.first) / NonZeroize(v.variables[0]), normalize(ld.line.second) / NonZeroize(v.variables[1]));
            //}
            return InstanceGivenVariables(mg, vars.at(lh).variables.data(), controls, lh);
        }



        Plane3 Instance(const RLGraph & mg, const RLGraphControls & controls, const RLGraphVars & vars, const RegionHandle & rh) {
            //auto & rd = mg.data(rh);
            //auto & c = controls[rh];
            //auto & v = vars[rh];
            ///*if (!rp.used)
            //    return Plane3();*/
            //if (c.orientationClaz >= 0){
            //    assert(v.variables.size() == 1);
            //    return Plane3(rd.normalizedCenter / NonZeroize(v.variables[0]), controls.vanishingPoints[c.orientationClaz]);
            //}
            //else if (c.orientationClaz == -1 && c.orientationNotClaz >= 0){
            //    assert(v.variables.size() == 2);
            //    double vs[] = { v.variables[0], v.variables[1], 0.0 }; // fake vs
            //    // v1 * o1 + v2 * o2 + v3 * o3 = 0, since region.orientation is orthogonal to normal
            //    auto orientation = normalize(controls.vanishingPoints[c.orientationNotClaz]);
            //    int c = SwappedComponent(orientation);
            //    std::swap(orientation[c], orientation[2]); // now fake orientation
            //    vs[2] = (-vs[0] * orientation[0] - vs[1] * orientation[1])
            //        / orientation[2];
            //    std::swap(vs[c], vs[2]); // now real vs
            //    return Plane3FromEquation(vs[0], vs[1], vs[2]);
            //}
            //else /*if (region.type == MGUnary::RegionWithFixedNormal)*/{
            //    assert(v.variables.size() == 3);
            //    return Plane3FromEquation(v.variables[0], v.variables[1], v.variables[2]);
            //}
            return InstanceGivenVariables(mg, vars.at(rh).variables.data(), controls, rh);
        }




        std::vector<Polygon3> RegionPolygon(const RLGraph & mg, const RLGraphControls & controls,
            const RLGraphVars & vars, RegionHandle rh) {
            if (!controls[rh].used)
                return std::vector<Polygon3>();
            std::vector<Polygon3> ps;
            Plane3 plane = Instance(mg, controls, vars, rh);
            for (auto & contour : mg.data(rh).normalizedContours) {
                Polygon3 polygon;
                polygon.corners.reserve(contour.size());
                for (auto & c : contour) {
                    polygon.corners.push_back(PointAt(c, plane));
                }
                polygon.normal = plane.normal;
                ps.push_back(std::move(polygon));
            }
            return ps;
        }

        HandledTable<RegionHandle, std::vector<Polygon3>> RegionPolygons(const RLGraph & mg,
            const RLGraphControls & controls,
            const RLGraphVars & vars) {

            auto polygons = mg.createComponentTable<RegionData, std::vector<Polygon3>>();
            for (auto & r : mg.components<RegionData>()) {
                if (!controls[r.topo.hd].used)
                    continue;
                polygons[r.topo.hd] = RegionPolygon(mg, controls, vars, r.topo.hd);
            }
            return polygons;
        }





        std::vector<double> VariableCoefficientsForInverseDepthAtDirection(const RLGraph & mg,
            const RLGraphControls & controls, const RLGraphVars & vars, const Vec3 & direction, const LineHandle & lh) {
            auto & ld = mg.data(lh);
            auto & c = controls[lh];
            auto & v = vars[lh];
            if (c.orientationClaz >= 0) {
                assert(v.variables.size() == 1);
                Ray3 infLine(normalize(ld.line.center()), controls.vanishingPoints[c.orientationClaz]);
                // variable is 1.0/centerDepth
                // corresponding coeff is 1.0/depthRatio
                // so that 1.0/depth = 1.0/centerDepth * 1.0/depthRatio -> depth = centerDepth * depthRatio
                double depthRatio = norm(DistanceBetweenTwoLines(Ray3(Point3(0, 0, 0), direction), infLine).second.first);
                return std::vector<double>{1.0 / depthRatio};
            } else /*if(u.type == MGUnary::LineFree)*/ {
                assert(v.variables.size() == 2);
                double theta = AngleBetweenDirections(normalize(ld.line.first), normalize(ld.line.second));
                double phi = AngleBetweenDirections(normalize(ld.line.first), direction);
                /*           | sin(theta) | | p | | q |
                len:---------------------------------   => 1/len: [(1/q)sin(phi) - (1/p)sin(phi-theta)] / sin(theta)
                | p sin(phi) - q sin(phi - theta) |
                */
                // variables[0] -> 1/p
                // variables[1] -> 1/q
                double coeffFor1_p = -sin(phi - theta) / sin(theta);
                double coeffFor1_q = sin(phi) / sin(theta);
                assert(!IsInfOrNaN(coeffFor1_p) && !IsInfOrNaN(coeffFor1_q));
                return std::vector<double>{coeffFor1_p, coeffFor1_q};
            }
        }


        std::vector<double> VariableCoefficientsForInverseDepthAtDirection(const RLGraph & mg,
            const RLGraphControls & controls, const RLGraphVars & vars, const Vec3 & direction, const RegionHandle & rh) {
            auto & rd = mg.data(rh);
            auto & c = controls[rh];
            auto & v = vars[rh];
            if (c.orientationClaz >= 0) {
                assert(v.variables.size() == 1);
                Plane3 plane(rd.normalizedCenter, controls.vanishingPoints[c.orientationClaz]);
                // variable is 1.0/centerDepth
                // corresponding coeff is 1.0/depthRatio
                // so that 1.0/depth = 1.0/centerDepth * 1.0/depthRatio -> depth = centerDepth * depthRatio
                double depthRatio = norm(IntersectionOfLineAndPlane(Ray3(Point3(0, 0, 0), direction), plane).position);
                return std::vector<double>{1.0 / depthRatio};
            } else if (c.orientationClaz == -1 && c.orientationNotClaz >= 0) {
                assert(v.variables.size() == 2);
                // depth = 1.0 / (ax + by + cz) where (x, y, z) = direction, (a, b, c) = variables
                // -> 1.0/depth = ax + by + cz
                // -> 1.0/depth = ax + by + (- o1a - o2b)/o3 * z = a(x-o1*z/o3) + b(y-o2*z/o3)
                auto orientation = normalize(controls.vanishingPoints[c.orientationNotClaz]);
                int sc = SwappedComponent(orientation);
                Vec3 forientation = orientation;
                std::swap(forientation[sc], forientation[2]);
                Vec3 fdirection = direction;
                std::swap(fdirection[sc], fdirection[2]);
                return std::vector<double>{
                    fdirection[0] - forientation[0] * fdirection[2] / forientation[2],
                        fdirection[1] - forientation[1] * fdirection[2] / forientation[2]
                };
            } else {
                assert(v.variables.size() == 3);
                // depth = 1.0 / (ax + by + cz) where (x, y, z) = direction, (a, b, c) = variables
                // -> 1.0/depth = ax + by + cz
                return std::vector<double>{direction[0], direction[1], direction[2]};
            }
        }





        double DepthAtDirectionGivenVariables(const RLGraph & mg, const double * variables,
            const RLGraphControls & controls, const Vec3 & direction, const LineHandle & lh) {
            auto & ld = mg.data(lh);
            auto & lp = controls[lh];
            if (lp.orientationClaz >= 0) {
                //assert(lp.variables.size() == 1);
                Ray3 infLine(normalize(ld.line.center()), controls.vanishingPoints[lp.orientationClaz]);
                // variable is 1.0/centerDepth
                // depths = depthRatio * centerDepth
                double depthRatio = norm(DistanceBetweenTwoLines(Ray3(Point3(0, 0, 0), direction), infLine).second.first);
                double inversedCenterDepth = variables[0];
                return depthRatio / inversedCenterDepth;
            } else /*if(u.type == MGUnary::LineFree)*/ {
                //assert(lp.variables.size() == 2);
                double theta = AngleBetweenDirections(normalize(ld.line.first), normalize(ld.line.second));
                double phi = AngleBetweenDirections(normalize(ld.line.first), direction);
                /*           | sin(theta) | | p | | q |
                len:---------------------------------   => 1/len: [(1/q)sin(phi) - (1/p)sin(phi-theta)] / sin(theta)
                | p sin(phi) - q sin(phi - theta) |
                */
                // variables[0] -> 1/p
                // variables[1] -> 1/q
                return sin(theta) / (variables[1] * sin(phi) - variables[0] * sin(phi - theta));
            }
        }


        double DepthAtDirectionGivenVariables(const RLGraph & mg, const double * variables,
            const RLGraphControls & controls, const Vec3 & direction, const RegionHandle & rh) {
            auto & rd = mg.data(rh);
            auto & rp = controls[rh];
            if (rp.orientationClaz >= 0) {
                //assert(rp.variables.size() == 1);
                Plane3 plane(rd.normalizedCenter, controls.vanishingPoints[rp.orientationClaz]);
                // variable is 1.0/centerDepth
                // depths = depthRatio * centerDepth
                double depthRatio = norm(IntersectionOfLineAndPlane(Ray3(Point3(0, 0, 0), direction), plane).position);
                double inversedCenterDepth = variables[0];
                return depthRatio / inversedCenterDepth;
            } else if (rp.orientationClaz == -1 && rp.orientationNotClaz >= 0) {
                //assert(rp.variables.size() == 2);
                // depth = 1.0 / (ax + by + cz) where (x, y, z) = direction, (a, b, c) = variables
                // -> 1.0/depth = ax + by + cz
                // -> 1.0/depth = ax + by + (- o1a - o2b)/o3 * z = a(x-o1*z/o3) + b(y-o2*z/o3)
                auto orientation = normalize(controls.vanishingPoints[rp.orientationNotClaz]);
                int c = SwappedComponent(orientation);
                Vec3 forientation = orientation;
                std::swap(forientation[c], forientation[2]);
                Vec3 fdirection = direction;
                std::swap(fdirection[c], fdirection[2]);
                double a = variables[0];
                double b = variables[1];

                //NOT_TESTED_YET();

                return forientation[2] / (
                    a * (fdirection[0] * forientation[2] - forientation[0] * fdirection[2]) +
                    b * (fdirection[1] * forientation[2] - forientation[1] * fdirection[2]));
            } else {
                // depth = 1.0 / (ax + by + cz) where (x, y, z) = direction, (a, b, c) = variables
                // -> 1.0/depth = ax + by + cz
                return 1.0 / (direction[0] * variables[0] + direction[1] * variables[1] + direction[2] * variables[2]);
            }
        }




        bool AttachAnchorToCenterOfLargestRegionIfNoAnchorExists(const RLGraph & mg,
            RLGraphControls & controls,
            double depth, double weight) {
            std::cout << "calling this feature may cause EYE SKEW!" << std::endl;
            RegionHandle largest;
            double maxArea = 0.0;
            for (auto & r : mg.components<RegionData>()) {
                if (!controls[r.topo.hd].used)
                    continue;
                if (controls[r.topo.hd].weightedAnchors.size() > 0)
                    return true;
                if (r.data.area > maxArea) {
                    largest = r.topo.hd;
                    maxArea = r.data.area;
                }
            }
            if (largest.valid())
                controls[largest].weightedAnchors.push_back(ScoreAs(mg.data(largest).normalizedCenter * depth, weight));
            return largest.valid();
        }


        int NumberOfComponentWeightedAnchors(const RLGraphControls & controls) {
            int nanchor = 0;
            for (auto & table : controls.componentControls.data) {
                for (auto & c : table) {
                    if (!c.used)
                        continue;
                    nanchor += c.weightedAnchors.size();
                }
            }
            return nanchor;
        }

        int NumberOfComponentBoundedAnchors(const RLGraphControls & controls) {
            int nanchor = 0;
            for (auto & table : controls.componentControls.data) {
                for (auto & c : table) {
                    if (!c.used)
                        continue;
                    nanchor += c.boundedAnchors.size();
                }
            }
            return nanchor;
        }

        int NumberOfConstraintWeightedAnchors(const RLGraphControls & controls) {
            int nanchor = 0;
            for (auto & table : controls.constraintControls.data) {
                for (auto & c : table) {
                    if (!c.used)
                        continue;
                    nanchor += c.weightedAnchors.size();
                }
            }
            return nanchor;
        }

        int NumberOfConstraintBoundedAnchors(const RLGraphControls & controls) {
            int nanchor = 0;
            for (auto & table : controls.constraintControls.data) {
                for (auto & c : table) {
                    if (!c.used)
                        continue;
                    nanchor += c.boundedAnchors.size();
                }
            }
            return nanchor;
        }


        bool AttachWeightedAnchorToCenterOfLargestLineIfNoExists(const RLGraph & mg,
            RLGraphControls & controls,
            double depth, double weight, bool orientedOnly) {
            LineHandle largest;
            double maxArea = 0.0;
            for (auto & r : mg.components<LineData>()) {
                if (!controls[r.topo.hd].used)
                    continue;
                if (controls[r.topo.hd].orientationClaz == -1 && orientedOnly)
                    continue;
                if (controls[r.topo.hd].weightedAnchors.size() > 0)
                    return true;
                if (r.data.line.length() > maxArea) {
                    largest = r.topo.hd;
                    maxArea = r.data.line.length();
                }
            }
            if (largest.valid())
                controls[largest].weightedAnchors.push_back(ScoreAs(normalize(mg.data(largest).line.center()) * depth, weight));
            return largest.valid();
        }

        bool AttachWeightedAnchorToCenterOfLargestRegionIfNoExists(const RLGraph & mg,
            RLGraphControls & controls,
            double depth, double weight, bool orientedOnly) {
            RegionHandle largest;
            double maxArea = 0.0;
            for (auto & r : mg.components<RegionData>()) {
                if (!controls[r.topo.hd].used)
                    continue;
                if (controls[r.topo.hd].orientationClaz == -1 && orientedOnly)
                    continue;
                if (controls[r.topo.hd].weightedAnchors.size() > 0)
                    return true;
                if (r.data.area > maxArea) {
                    largest = r.topo.hd;
                    maxArea = r.data.area;
                }
            }
            if (largest.valid())
                controls[largest].weightedAnchors.push_back(ScoreAs(normalize(mg.data(largest).normalizedCenter) * depth, weight));
            return largest.valid();
        }


        void ClearAllComponentAnchors(RLGraphControls & controls) {
            for (auto & table : controls.componentControls.data) {
                for (auto & c : table) {
                    c.weightedAnchors.clear();
                }
            }
        }

        void ClearAllConstraintAnchors(RLGraphControls & controls) {
            for (auto & table : controls.constraintControls.data) {
                for (auto & c : table) {
                    c.weightedAnchors.clear();
                }
            }
        }



        namespace {

            // returns false if confliction occurs
            bool MakeRegionPlaneUsable(RegionHandle rh, bool usable, RLGraphControls & controls) {
                auto & p = controls[rh];
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
            bool MakeRegionPlaneToward(RegionHandle rh, int normalVPId, RLGraphControls & controls) {
                auto & p = controls[rh];
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
                auto & dir = controls.vanishingPoints[p.orientationNotClaz];
                if (IsFuzzyPerpendicular(controls.vanishingPoints[normalVPId], dir)) {
                    p.orientationClaz = normalVPId;
                    p.orientationNotClaz = -1;
                    return true;
                }
                return false;
            }

            // returns false if confliction occurs
            bool MakeRegionPlaneAlsoAlong(RegionHandle rh, int alongVPId, RLGraphControls & controls) {
                auto & p = controls[rh];
                if (!p.used)
                    return true;
                assert(alongVPId != -1);
                auto & dir = controls.vanishingPoints[alongVPId];
                if (p.orientationClaz != -1) {
                    auto & normal = controls.vanishingPoints[p.orientationClaz];
                    return IsFuzzyPerpendicular(normal, dir);
                }
                if (p.orientationNotClaz == -1) {
                    p.orientationNotClaz = alongVPId;
                    return true;
                }
                if (p.orientationNotClaz == alongVPId)
                    return true;

                auto newNormal = dir.cross(controls.vanishingPoints[p.orientationNotClaz]);
                double minAngle = M_PI;
                for (int i = 0; i < controls.vanishingPoints.size(); i++) {
                    double angle = AngleBetweenUndirectedVectors(controls.vanishingPoints[i], newNormal);
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







        static const int regionLineRelatedThreshold = 2;

        void AttachPrincipleDirectionConstraints(const RLGraph & mg, RLGraphControls & controls,
            double rangeAngle, bool avoidLineConflictions) {

            SetClock();

            // find peaky regions
            std::vector<std::vector<RegionHandle>> peakyRegionHandles(controls.vanishingPoints.size());
            for (auto & r : mg.components<RegionData>()) {
                auto h = r.topo.hd;
                if (!controls[h].used)
                    continue;

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
                for (auto & vp : controls.vanishingPoints) {
                    double angle = AngleBetweenUndirectedVectors(vp, r.data.normalizedCenter);
                    if (angle <= radiusAngle) {
                        mayCrossAnyVP = true;
                        break;
                    }
                }

                if (!mayCrossAnyVP) {
                    continue;
                }

                float ppcFocal = 100.0f;
                int ppcSize = 2 * radiusAngle * ppcFocal + 10;
                Vec3 x;
                std::tie(x, std::ignore) = ProposeXYDirectionsFromZDirection(r.data.normalizedCenter);
                PartialPanoramicCamera ppc(ppcSize, ppcSize, ppcFocal, Point3(0, 0, 0), r.data.normalizedCenter, x);
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

                // intersection test
                for (int i = 0; i < controls.vanishingPoints.size(); i++) {
                    auto p1 = ToPixel(ppc.toScreen(controls.vanishingPoints[i]));
                    auto p2 = ToPixel(ppc.toScreen(-controls.vanishingPoints[i]));

                    int dilateSize = ppcFocal * rangeAngle;
                    bool intersected = false;
                    for (int x = -dilateSize; x <= dilateSize; x++) {
                        if (intersected)
                            break;
                        for (int y = -dilateSize; y <= dilateSize; y++) {
                            if (intersected)
                                break;
                            auto pp1 = Pixel(p1.x + x, p1.y + y);
                            auto pp2 = Pixel(p2.x + x, p2.y + y);
                            if (Contains(mask, pp1) && mask(pp1)) {
                                peakyRegionHandles[i].push_back(r.topo.hd);
                                intersected = true;
                            } else if (Contains(mask, pp2) && mask(pp2)) {
                                peakyRegionHandles[i].push_back(r.topo.hd);
                                intersected = true;
                            }
                        }
                    }
                }
            }

            for (int i = 0; i < controls.vanishingPoints.size(); i++) {
                auto & vp = controls.vanishingPoints[i];
                auto & rhs = peakyRegionHandles[i];
                for (auto rh : rhs) {
                    if (avoidLineConflictions) {
                        bool hasLineConflictions = false;
                        for (auto conh : mg.topo(rh).constraints<RegionLineConnectionData>()) {
                            if (mg.data(conh).normalizedAnchors.size() < regionLineRelatedThreshold)
                                continue;
                            LineHandle lh = mg.topo(conh).component<1>();
                            auto & lprop = controls[lh];
                            if (lprop.orientationClaz == i) {
                                hasLineConflictions = true;
                                break;
                            }
                        }
                        if (!hasLineConflictions)
                            MakeRegionPlaneToward(rh, i, controls);
                    } else {
                        MakeRegionPlaneToward(rh, i, controls);
                        for (auto conh : mg.topo(rh).constraints<RegionLineConnectionData>()) {
                            if (mg.data(conh).normalizedAnchors.size() < regionLineRelatedThreshold)
                                continue;
                            LineHandle lh = mg.topo(conh).component<1>();
                            auto & lprop = controls[lh];
                            if (lprop.orientationClaz == i) {
                                lprop.orientationClaz = -1;
                            }
                        }
                    }
                }
            }

            controls.disableAllInvalidConstraints(mg);

        }


        void AttachWallConstriants(const RLGraph & mg, RLGraphControls & controls,
            double rangeAngle, const Vec3 & verticalSeed) {

            SetClock();

            int vertVPId = -1;
            double minAngle = M_PI;
            for (int i = 0; i < controls.vanishingPoints.size(); i++) {
                double angle = AngleBetweenUndirectedVectors(verticalSeed, controls.vanishingPoints[i]);
                if (angle < minAngle) {
                    minAngle = angle;
                    vertVPId = i;
                }
            }
            assert(vertVPId != -1);
            auto & vertical = controls.vanishingPoints[vertVPId];

            std::vector<RegionHandle> horizontalRegionHandles;
            for (auto & r : mg.components<RegionData>()) {
                auto h = r.topo.hd;
                if (!controls[h].used)
                    continue;
                auto & contours = r.data.normalizedContours;
                bool intersected = false;
                for (auto & cs : r.data.normalizedContours) {
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
                    horizontalRegionHandles.push_back(h);
                }
            }

            for (auto h : horizontalRegionHandles) {
                MakeRegionPlaneAlsoAlong(h, vertVPId, controls);
            }

            controls.disableAllInvalidConstraints(mg);
        }



        void AttachFloorAndCeilingConstraints(const RLGraph & mg,
            RLGraphControls & controls, const RLGraphVars & vars,
            double eyeHeightRatioLowerBound, double eyeHeightRatioUpperBound,
            double angleThreshold, const Vec3 & verticalSeed) {

            SetClock();

            IMPROVABLE_HERE("we need a better way!");

            int vertVPId = -1;
            double minAngle = M_PI;
            for (int i = 0; i < controls.vanishingPoints.size(); i++) {
                double angle = AngleBetweenUndirectedVectors(verticalSeed, controls.vanishingPoints[i]);
                if (angle < minAngle) {
                    minAngle = angle;
                    vertVPId = i;
                }
            }
            assert(vertVPId != -1);
            auto & vertical = controls.vanishingPoints[vertVPId];

            int horizVPId = -1;
            double maxAngle = 0;
            for (int i = 0; i < controls.vanishingPoints.size(); i++) {
                double angle = AngleBetweenUndirectedVectors(vertical, controls.vanishingPoints[i]);
                if (angle > maxAngle) {
                    maxAngle = angle;
                    horizVPId = i;
                }
            }
            assert(horizVPId != -1);
            auto xDir = normalize(controls.vanishingPoints[horizVPId]);
            auto yDir = normalize(vertical.cross(xDir));

            // compute horizontal bounding range of may-be ceilings and floors
            double xmin = std::numeric_limits<double>::max(), ymin = std::numeric_limits<double>::max();
            double xmax = std::numeric_limits<double>::lowest(), ymax = std::numeric_limits<double>::lowest();

            std::vector<Scored<RegionHandle>> maybeCeilinsOrFloors[2];

            auto regionPlanes = MakeHandledTableForComponents<Plane3, RegionData>(mg);
            for (auto & r : mg.components<RegionData>()) {
                if (!controls[r.topo.hd].used)
                    continue;

                auto plane = Instance(mg, controls, vars, r.topo.hd);
                regionPlanes[r.topo.hd] = plane;
                double angleToVert = AngleBetweenUndirectedVectors(plane.normal, vertical);
                if (angleToVert < angleThreshold) {
                    int belonging = r.data.normalizedCenter.dot(vertical) < 0 ? 0 : 1;
                    double centerZ = IntersectionOfLineAndPlane(Ray3(Point3(0, 0, 0), r.data.normalizedCenter), plane).position.dot(vertical);
                    maybeCeilinsOrFloors[belonging].push_back(ScoreAs(r.topo.hd, centerZ));

                    for (auto & cs : r.data.normalizedContours) {
                        for (auto & c : cs) {
                            Point3 point = IntersectionOfLineAndPlane(Ray3(Point3(0, 0, 0), c), plane).position;
                            float x = point.dot(xDir), y = point.dot(yDir);
                            if (x < xmin) xmin = x;
                            if (x > xmax) xmax = x;
                            if (y < ymin) ymin = y;
                            if (y > ymax) ymax = y;
                        }
                    }
                }
            }

            if (maybeCeilinsOrFloors[0].empty() && maybeCeilinsOrFloors[1].empty())
                return;

            // estimate height or ceiling/floor, remove non-ceiling/floor horizontal regions
            double medianCenterDepth = MedianCenterDepth(mg, controls, vars);
            static const double heightAffectRange = medianCenterDepth * 0.01;
            for (int i = 0; i < 2; i++) {
                auto & rhsWithHeights = maybeCeilinsOrFloors[i];
                std::sort(rhsWithHeights.begin(), rhsWithHeights.end());
                std::vector<double> votes(rhsWithHeights.size(), 0.0);
                for (int a = 0; a < rhsWithHeights.size(); a++) {
                    for (int b = 0; b < rhsWithHeights.size(); b++) {
                        votes[a] += Gaussian(rhsWithHeights[a].score - rhsWithHeights[b].score,
                            heightAffectRange);
                    }
                }
                int maxId = std::distance(votes.begin(), std::max_element(votes.begin(), votes.end()));

                std::vector<Scored<RegionHandle>> hs;
                hs.reserve(rhsWithHeights.size());
                for (auto rhh : rhsWithHeights) {
                    if (Distance(rhh.score, rhsWithHeights[maxId].score) < heightAffectRange) {
                        hs.push_back(rhh);
                    }
                }
                rhsWithHeights = std::move(hs);
            }

            assert(xmax > xmin && ymax > ymin);
            int imSize = 300;
            double scale = imSize / ((xmax + ymax - xmin - ymin) / 2.0);
            int width = static_cast<int>((xmax - xmin) * scale + 1);
            int height = static_cast<int>((ymax - ymin) * scale + 1);
            Imageub ceilingOrFloorHorizontalRegionMasks[] = {
                Imageub::zeros(cv::Size(width, height)),
                Imageub::zeros(cv::Size(width, height))
            };
            for (int i = 0; i < 2; i++) {
                for (auto & rhh : maybeCeilinsOrFloors[i]) {
                    auto rh = rhh.component;
                    auto & contours = mg.data(rh).normalizedContours;
                    std::vector<std::vector<Point2i>> contourResized;
                    contourResized.reserve(contours.size());
                    for (auto & cs : contours) {
                        std::vector<Point2i> csresized;
                        csresized.reserve(cs.size());
                        for (auto & c : cs) {
                            Point3 point = IntersectionOfLineAndPlane(Ray3(Point3(0, 0, 0), c), regionPlanes[rh]).position;
                            float x = point.dot(xDir), y = point.dot(yDir);
                            csresized.emplace_back(static_cast<int>((x - xmin) * scale),
                                static_cast<int>((y - ymin) * scale));
                        }
                        contourResized.push_back(std::move(csresized));
                    }
                    cv::fillPoly(ceilingOrFloorHorizontalRegionMasks[i], contourResized, (uint8_t)1);
                }
            }

            Imagei ceilingOrFloorRegionMasks[] = {
                Imagei(cv::Size(width, height), -1),
                Imagei(cv::Size(width, height), -1)
            };
            for (auto & r : mg.components<RegionData>()) {
                const auto & prop = controls[r.topo.hd];
                if (!prop.used)
                    continue;
                if (prop.orientationClaz != -1) // already strongly constrained
                    continue;
                if (prop.orientationNotClaz == vertVPId) // constrained to be vertical
                    continue;

                auto & plane = regionPlanes[r.topo.hd];
                auto & contours = r.data.normalizedContours;

                std::vector<std::vector<Point2i>> contourResized;
                contourResized.reserve(contours.size());

                bool hasPosDot = false, hasNegDot = false;
                for (auto & cs : contours) {
                    std::vector<Point2i> csresized;
                    csresized.reserve(cs.size());
                    for (auto & c : cs) {
                        Point3 point = IntersectionOfLineAndPlane(Ray3(Point3(0, 0, 0), c), plane).position;
                        float x = point.dot(xDir), y = point.dot(yDir);
                        float z = point.dot(vertical);
                        hasPosDot |= z >= 0;
                        hasNegDot |= z <= 0;
                        csresized.emplace_back(static_cast<int>((x - xmin) * scale),
                            static_cast<int>((y - ymin) * scale));
                    }
                    contourResized.push_back(std::move(csresized));
                }

                if (hasPosDot && hasNegDot) { // this region crosses the horizon
                    continue;
                }
                if (contourResized.empty()) {
                    continue;
                }

                auto & mask = ceilingOrFloorRegionMasks[hasNegDot ? 0 : 1];
                cv::fillPoly(mask, contourResized, r.topo.hd.id);
            }

            // match masks
            for (auto & horizontalMask : ceilingOrFloorHorizontalRegionMasks) {
                static const int erosionSize = 5;
                cv::erode(horizontalMask, horizontalMask, cv::getStructuringElement(cv::MORPH_ELLIPSE,
                    cv::Size(2 * erosionSize + 1, 2 * erosionSize + 1),
                    cv::Point(erosionSize, erosionSize)));
            }

            //cv::imshow("0", ceilingOrFloorHorizontalRegionMasks[0] * 255);
            //cv::imshow("1", ceilingOrFloorHorizontalRegionMasks[1] * 255);
            //cv::waitKey();

            auto regionAffectedCounts = mg.createComponentTable<RegionData>(0);

            for (int i = 0; i < 2; i++) {
                auto & regionMasksHere = ceilingOrFloorRegionMasks[i];
                auto & horizontalMaskThere = ceilingOrFloorHorizontalRegionMasks[1 - i];
                for (auto it = regionMasksHere.begin(); it != regionMasksHere.end(); ++it) {
                    if (!horizontalMaskThere(it.pos()))
                        continue;
                    RegionHandle regionHandle(*it);
                    if (regionHandle.invalid())
                        continue;

                    regionAffectedCounts[regionHandle] ++;
                    // add horizontal constraint
                    if (regionAffectedCounts[regionHandle] > Square(imSize * 0.02)) {
                        MakeRegionPlaneToward(regionHandle, vertVPId, controls);
                    }
                }
            }

            controls.disableAllInvalidConstraints(mg);
        }



        //void AttachGeometricContextConstraints(const RLGraph & mg, RLGraphControls & controls,
        //    const std::vector<GeometricContextEstimator::Feature> & perspectiveGCs,
        //    const std::vector<PerspectiveCamera> & gcCameras,
        //    int shrinkRegionOrientationIteration, bool considerGCVerticalConstraint){

        //    auto orientationVotes = mg.createComponentTable<RegionData, std::unordered_map<OrientationHint, double>>(
        //        std::unordered_map<OrientationHint, double>((size_t)OrientationHint::Count));

        //    assert(perspectiveGCs.size() == gcCameras.size());

        //    // region views
        //    auto regionMaskViews = mg.createComponentTable<RegionData, View<PartialPanoramicCamera, Imageub>>();
        //    for (auto & r : mg.components<RegionData>()){
        //        auto h = r.topo.hd;
        //        auto & contours = r.data.normalizedContours;
        //        double radiusAngle = 0.0;
        //        for (auto & cs : r.data.normalizedContours){
        //            for (auto & c : cs){
        //                double angle = AngleBetweenDirections(r.data.normalizedCenter, c);
        //                if (angle > radiusAngle){
        //                    radiusAngle = angle;
        //                }
        //            }
        //        }
        //        float ppcFocal = 100.0f;
        //        int ppcSize = 2 * radiusAngle * ppcFocal;
        //        Vec3 x;
        //        std::tie(x, std::ignore) = ProposeXYDirectionsFromZDirection(r.data.normalizedCenter);
        //        PartialPanoramicCamera ppc(ppcSize, ppcSize, ppcFocal, Point3(0, 0, 0), r.data.normalizedCenter, x);
        //        Imageub mask = Imageub::zeros(ppc.screenSize());

        //        // project contours to ppc
        //        std::vector<std::vector<Point2i>> contourProjs(contours.size());
        //        for (int k = 0; k < contours.size(); k++){
        //            auto & contourProj = contourProjs[k];
        //            contourProj.reserve(contours[k].size());
        //            for (auto & d : contours[k]){
        //                contourProj.push_back(ecast<int>(ppc.toScreen(d)));
        //            }
        //        }
        //        cv::fillPoly(mask, contourProjs, (uint8_t)1);
        //        regionMaskViews[h].camera = ppc;
        //        regionMaskViews[h].image = mask;
        //    }

        //    auto regionAreas = mg.createComponentTable<RegionData, double>(0.0);
        //    for (auto & r : mg.components<RegionData>()){
        //        auto h = r.topo.hd;
        //        auto & ppMask = regionMaskViews[h];
        //        double area = cv::sum(ppMask.image).val[0];
        //        regionAreas[h] = area;
        //        for (int i = 0; i < perspectiveGCs.size(); i++){
        //            auto sampler = MakeCameraSampler(ppMask.camera, gcCameras[i]);
        //            auto occupation = sampler(Imageub::ones(gcCameras[i].screenSize()), cv::BORDER_CONSTANT, 0.0);
        //            double areaOccupied = cv::sum(ppMask.image & occupation).val[0];
        //            double ratio = areaOccupied / area;
        //            assert(ratio <= 1.01);
        //            GeometricContextLabel maxLabel = GeometricContextLabel::None;
        //            double maxScore = 0.0;
        //            for (auto & gcc : perspectiveGCs[i]){
        //                auto label = gcc.first;
        //                Imaged ppGC = sampler(gcc.second);
        //                double score = cv::mean(ppGC, ppMask.image).val[0];
        //                if (score > maxScore){
        //                    maxScore = score;
        //                    maxLabel = label;
        //                }
        //            }
        //            if (maxLabel != GeometricContextLabel::None){
        //                auto & v = orientationVotes[r.topo.hd][ToOrientationHint(maxLabel, considerGCVerticalConstraint)];
        //                v = std::max(v, ratio);
        //            }
        //        }
        //    }

        //    // optimize
        //    GCoptimizationGeneralGraph graph(mg.internalComponents<RegionData>().size(), (int)OrientationHint::Count);

        //    double maxBoundaryLength = 0;
        //    for (auto & b : mg.constraints<RegionBoundaryData>()){
        //        maxBoundaryLength = std::max(maxBoundaryLength, b.data.length);
        //    }
        //    for (auto & b : mg.constraints<RegionBoundaryData>()){
        //        graph.setNeighbors(b.topo.component<0>().id, b.topo.component<1>().id,
        //            100 * b.data.length / maxBoundaryLength);
        //    }
        //    for (auto & r : mg.components<RegionData>()){
        //        auto & orientationVote = orientationVotes[r.topo.hd];
        //        // normalize votes
        //        double votesSum = 0.0;
        //        for (auto & v : orientationVote){
        //            votesSum += v.second;
        //        }
        //        assert(votesSum >= 0.0);
        //        if (votesSum > 0.0){
        //            for (auto & v : orientationVote){
        //                assert(v.second >= 0.0);
        //                v.second /= votesSum;
        //            }
        //        }
        //        for (int label = 0; label < (int)OrientationHint::Count; label++){
        //            double vote = Contains(orientationVote, (OrientationHint)label) ?
        //                orientationVote.at((OrientationHint)label) : 0.0;
        //            graph.setDataCost(r.topo.hd.id, label,
        //                votesSum == 0.0 ? 0 : (10000 * (1.0 - vote)));
        //        }
        //    }
        //    for (int label1 = 0; label1 < (int)OrientationHint::Count; label1++){
        //        for (int label2 = 0; label2 < (int)OrientationHint::Count; label2++){
        //            if (label1 == label2){
        //                graph.setSmoothCost(label1, label2, 0);
        //            }
        //            else{
        //                graph.setSmoothCost(label1, label2, 1);
        //            }
        //        }
        //    }

        //    graph.expansion();
        //    graph.swap();

        //    // get the most vertical vp id
        //    int vVPId = -1;
        //    double angleToVert = std::numeric_limits<double>::max();
        //    for (int i = 0; i < controls.vanishingPoints.size(); i++){
        //        double a = AngleBetweenUndirectedVectors(controls.vanishingPoints[i], Vec3(0, 0, 1));
        //        if (a < angleToVert){
        //            vVPId = i;
        //            angleToVert = a;
        //        }
        //    }
        //    assert(vVPId != -1);

        //    // disorient outsided oriented gc labels
        //    auto regionOrientations = mg.createComponentTable<RegionData, OrientationHint>(OrientationHint::Void);
        //    for (auto & r : mg.components<RegionData>()){
        //        regionOrientations[r.topo.hd] = (OrientationHint)(graph.whatLabel(r.topo.hd.id));
        //    }

        //    for (int i = 0; i < shrinkRegionOrientationIteration; i++){
        //        auto shrinked = regionOrientations;
        //        for (auto & b : mg.constraints<RegionBoundaryData>()){
        //            auto l1 = regionOrientations[b.topo.component<0>()];
        //            auto l2 = regionOrientations[b.topo.component<1>()];
        //            if (l1 == OrientationHint::Horizontal && l2 != OrientationHint::Horizontal){
        //                shrinked[b.topo.component<0>()] = OrientationHint::OtherPlanar;
        //            }
        //            else if (l1 != OrientationHint::Horizontal && l2 == OrientationHint::Horizontal){
        //                shrinked[b.topo.component<1>()] = OrientationHint::OtherPlanar;
        //            }
        //            if (l1 == OrientationHint::Vertical && l2 != OrientationHint::Vertical){
        //                shrinked[b.topo.component<0>()] = OrientationHint::OtherPlanar;
        //            }
        //            else if (l1 != OrientationHint::Vertical && l2 == OrientationHint::Vertical){
        //                shrinked[b.topo.component<1>()] = OrientationHint::OtherPlanar;
        //            }
        //        }
        //        regionOrientations = std::move(shrinked);
        //    }

        //    // install gc labels
        //    for (auto & r : mg.components<RegionData>()){
        //        OrientationHint oh = regionOrientations[r.topo.hd];
        //        switch (oh) {
        //        case OrientationHint::Void:
        //            /*controls[r.topo.hd].used = false;
        //            controls[r.topo.hd].orientationClaz = -1;
        //            controls[r.topo.hd].orientationNotClaz = -1;*/
        //            MakeRegionPlaneUsable(r.topo.hd, false, controls);
        //            break;
        //        case OrientationHint::Horizontal:
        //            /*controls[r.topo.hd].used = true;
        //            controls[r.topo.hd].orientationClaz = vVPId;
        //            controls[r.topo.hd].orientationNotClaz = -1;*/
        //            MakeRegionPlaneToward(r.topo.hd, vVPId, controls);
        //            break;
        //        case OrientationHint::Vertical:
        //            //controls[r.topo.hd].used = true;
        //            //controls[r.topo.hd].orientationClaz = -1;
        //            //controls[r.topo.hd].orientationNotClaz = vVPId;
        //            MakeRegionPlaneAlsoAlong(r.topo.hd, vVPId, controls);
        //            break;
        //        case OrientationHint::OtherPlanar:
        //            /*controls[r.topo.hd].used = true;
        //            controls[r.topo.hd].orientationClaz = -1;
        //            controls[r.topo.hd].orientationNotClaz = -1;*/
        //            //MakeRegionPlaneUsable(r.topo.hd, true, controls);
        //            break;
        //        case OrientationHint::NonPlanar:
        //            /*controls[r.topo.hd].used = false;
        //            controls[r.topo.hd].orientationClaz = -1;
        //            controls[r.topo.hd].orientationNotClaz = -1;*/
        //            MakeRegionPlaneUsable(r.topo.hd, false, controls);
        //            break;
        //        default:
        //            assert(0);
        //            break;
        //        }
        //    }

        //    /*
        //    // detect big aligned rectangular regions and orient them
        //    double regionAreaSum = std::accumulate(regionAreas.begin(), regionAreas.end(), 0.0);
        //    assert(regionAreaSum > 0);

        //    static const double dotThreshold = 0.001;
        //    static const double spanAngleThreshold = DegreesToRadians(10);
        //    static const double wholeDotThreshold = 0.01;
        //    static const double wholeSpanAngleThreshold = DegreesToRadians(10);

        //    HandledTable<RegionBoundaryHandle, std::vector<double>> boundaryMaxSpanAnglesForVPs(mg.internalConstraints<RegionBoundaryData>().size());
        //    for (auto & b : mg.constraints<RegionBoundaryData>()){
        //    auto & edges = b.data.normalizedEdges;
        //    std::vector<double> maxSpanAngleForVPs(props.vanishingPoints.size(), 0.0);
        //    if (b.topo.hd.id == 849){
        //    std::cout << std::endl;
        //    }

        //    for (auto & edge : edges){
        //    std::vector<std::vector<bool>> edgeVPFlags(props.vanishingPoints.size(),
        //    std::vector<bool>(edge.size() - 1, false));
        //    for (int i = 0; i < edge.size() - 1; i++){
        //    Vec3 n = normalize(edge[i].cross(edge[i + 1]));
        //    for (int k = 0; k < props.vanishingPoints.size(); k++){
        //    auto vp = normalize(props.vanishingPoints[k]);
        //    double dotv = abs(vp.dot(n));
        //    if (dotv <= dotThreshold){
        //    edgeVPFlags[k][i] = true;
        //    }
        //    }
        //    }
        //    // detect aligned longest line spanAngles from edge
        //    for (int k = 0; k < props.vanishingPoints.size(); k++){
        //    auto vp = normalize(props.vanishingPoints[k]);
        //    auto & edgeFlags = edgeVPFlags[k];
        //    int lastHead = -1, lastTail = -1;
        //    bool inChain = false;

        //    double & maxSpanAngle = maxSpanAngleForVPs[k];
        //    for (int j = 0; j <= edgeFlags.size(); j++){
        //    if (!inChain && j < edgeFlags.size() && edgeFlags[j]){
        //    lastHead = j;
        //    inChain = true;
        //    }
        //    else if (inChain && (j == edgeFlags.size() || !edgeFlags[j])){
        //    lastTail = j;
        //    inChain = false;
        //    // examine current chain
        //    assert(lastHead != -1);
        //    // compute full span angle
        //    double spanAngle = 0.0;
        //    for (int i = lastHead; i < lastTail; i++){
        //    spanAngle += AngleBetweenDirections(edge[i], edge[i + 1]);
        //    }
        //    if (spanAngle < spanAngleThreshold){
        //    continue;
        //    }
        //    // fit line
        //    const Vec3 & midCorner = edge[(lastHead + lastTail) / 2];
        //    Vec3 commonNormal = normalize(midCorner.cross(vp));
        //    std::vector<double> dotsToCommonNormal(lastTail - lastHead + 1, 0.0);
        //    for (int i = lastHead; i <= lastTail; i++){
        //    dotsToCommonNormal[i - lastHead] = abs(edge[i].dot(commonNormal));
        //    }
        //    if (*std::max_element(dotsToCommonNormal.begin(), dotsToCommonNormal.end()) <= wholeDotThreshold){
        //    // acceptable!
        //    if (spanAngle > maxSpanAngle){
        //    maxSpanAngle = spanAngle;
        //    }
        //    }
        //    }
        //    }
        //    }
        //    }
        //    if ((b.topo.component<0>().id == 355 || b.topo.component<1>().id == 355) && std::accumulate(maxSpanAngleForVPs.begin(), maxSpanAngleForVPs.end(), 0.0) > 0){
        //    std::cout << std::endl;
        //    }
        //    boundaryMaxSpanAnglesForVPs[b.topo.hd] = std::move(maxSpanAngleForVPs);
        //    }
        //    for (auto & r : mg.components<RegionData>()){
        //    if (!props[r.topo.hd].used)
        //    continue;
        //    if (props[r.topo.hd].orientationClaz != -1 || props[r.topo.hd].orientationNotClaz != -1)
        //    continue;

        //    if (r.topo.hd.id == 355){
        //    std::cout << std::endl;
        //    }

        //    double area = regionAreas[r.topo.hd];
        //    if (area / regionAreaSum < 0.005){
        //    continue;
        //    }
        //    std::vector<double> spanAngleSumsForVPs(props.vanishingPoints.size(), 0.0);
        //    for (auto & h : r.topo.constraints<RegionBoundaryData>()){
        //    for (int i = 0; i < props.vanishingPoints.size(); i++){
        //    spanAngleSumsForVPs[i] += boundaryMaxSpanAnglesForVPs[h][i];
        //    }
        //    }
        //    auto maxIter = std::max_element(spanAngleSumsForVPs.begin(), spanAngleSumsForVPs.end());
        //    int firstMaxVPId = std::distance(spanAngleSumsForVPs.begin(), maxIter);
        //    double firstMaxSpanAngle = *maxIter;
        //    *maxIter = -1.0;
        //    maxIter = std::max_element(spanAngleSumsForVPs.begin(), spanAngleSumsForVPs.end());
        //    int secondMaxVPId = std::distance(spanAngleSumsForVPs.begin(), maxIter);
        //    double secondMaxSpanAngle = *maxIter;
        //    if (secondMaxSpanAngle >= wholeSpanAngleThreshold){
        //    assert(firstMaxVPId != secondMaxVPId);
        //    Vec3 normal = props.vanishingPoints[firstMaxVPId].cross(props.vanishingPoints[secondMaxVPId]);
        //    double minAngle = 0.1;
        //    int bestMatchedVPId = -1;
        //    for (int i = 0; i < props.vanishingPoints.size(); i++){
        //    double angle = AngleBetweenUndirectedVectors(normal, props.vanishingPoints[i]);
        //    if (angle < minAngle){
        //    minAngle = angle;
        //    bestMatchedVPId = i;
        //    }
        //    }
        //    if (bestMatchedVPId != -1){
        //    std::cout << "vpid:::: " << bestMatchedVPId << std::endl;
        //    props[r.topo.hd].orientationClaz = bestMatchedVPId;
        //    props[r.topo.hd].orientationNotClaz = -1;
        //    }
        //    }
        //    else if (firstMaxSpanAngle >= wholeSpanAngleThreshold){
        //    std::cout << "vpid-along:::: " << firstMaxVPId << std::endl;
        //    props[r.topo.hd].orientationClaz = -1;
        //    props[r.topo.hd].orientationNotClaz = firstMaxVPId;
        //    }
        //    }*/

        //    //for (auto & l : mg.internalComponents<LineData>()){
        //    //    props[l.topo.hd].used = true;
        //    //    props[l.topo.hd].orientationClaz = l.data.initialClaz;
        //    //    props[l.topo.hd].orientationNotClaz = -1;
        //    //}
        //    
        //    controls.disableAllInvalidConstraints(mg);

        //}







        void LooseOrientationConstraintsOnComponents(const RLGraph & mg,
            RLGraphControls & controls, const RLGraphVars & vars,
            double linesLoosableRatio, double regionsLoosableRatio, double distThresRatio) {

            SetClock();

            if (linesLoosableRatio <= 0.0 && regionsLoosableRatio <= 0.0)
                return;

            double medianCenterDepth = MedianCenterDepth(mg, controls, vars);
            double distThres = medianCenterDepth * distThresRatio;
            double nnSearchRange = distThres;

            bool directApproach = true;
            if (directApproach) {

                // collect current instances
                HandledTable<LineHandle, Line3> lines(mg.internalComponents<LineData>().size());
                HandledTable<RegionHandle, Plane3> planes(mg.internalComponents<RegionData>().size());

                for (auto & r : mg.components<RegionData>()) {
                    if (!controls[r.topo.hd].used)
                        continue;
                    planes[r.topo.hd] = Instance(mg, controls, vars, r.topo.hd);
                }
                for (auto & l : mg.components<LineData>()) {
                    if (!controls[l.topo.hd].used)
                        continue;
                    lines[l.topo.hd] = Instance(mg, controls, vars, l.topo.hd);
                }

                if (linesLoosableRatio > 0.0) {
                    std::vector<LineHandle> usableLineHandles;
                    usableLineHandles.reserve(mg.internalComponents<LineData>().size());
                    auto distanceToNearestRegions = mg.createComponentTable<LineData>(0.0);
                    for (auto & l : mg.components<LineData>()) {
                        if (!controls[l.topo.hd].used)
                            continue;
                        usableLineHandles.push_back(l.topo.hd);

                        Point3 points[] = { lines[l.topo.hd].first, lines[l.topo.hd].second };
                        double & ndist = distanceToNearestRegions[l.topo.hd];
                        for (auto & linePoint : points) {
                            double distToThis = std::numeric_limits<double>::max();
                            for (auto rlch : l.topo.constraints<RegionLineConnectionData>()) {
                                auto & anchors = mg.data(rlch).normalizedAnchors;
                                RegionHandle rh = mg.topo(rlch).component<0>();
                                Plane3 & plane = planes[rh];
                                for (auto & anchor : anchors) {
                                    auto planePoint = IntersectionOfLineAndPlane(Ray3(Point3(0, 0, 0), anchor), plane).position;
                                    double distance = Distance(planePoint, linePoint);
                                    if (distance < distToThis) {
                                        distToThis = distance;
                                    }
                                }
                            }
                            if (distToThis > ndist) {
                                ndist = distToThis;
                            }
                        }
                    }

                    std::sort(usableLineHandles.begin(), usableLineHandles.end(), [&distanceToNearestRegions](LineHandle l1, LineHandle l2) {
                        return distanceToNearestRegions[l1] > distanceToNearestRegions[l2];
                    });

                    // loose orientation constraints on these lines
                    for (int i = 0; i < usableLineHandles.size() * linesLoosableRatio; i++) {
                        double maxCornerDist = distanceToNearestRegions[usableLineHandles[i]];
                        if (maxCornerDist < distThres)
                            continue;
                        auto & lp = controls[usableLineHandles[i]];
                        assert(lp.used);
                        if (lp.orientationClaz >= 0) {
                            lp.orientationClaz = -1;
                        }
                        assert(lp.orientationNotClaz == -1);
                        if (lp.orientationClaz == -1) {
                            lp.used = false; // disable this line!!!!!!!!!
                        }
                    }
                }

                if (regionsLoosableRatio > 0.0) {
                    std::vector<RegionHandle> usableRegionHandles;
                    usableRegionHandles.reserve(mg.internalComponents<RegionData>().size());
                    auto maxDistanceToNearestComponents = mg.createComponentTable<RegionData>(0.0);

                    for (auto & r : mg.components<RegionData>()) {
                        if (!controls[r.topo.hd].used)
                            continue;
                        usableRegionHandles.push_back(r.topo.hd);
                        Plane3 & thisPlane = planes[r.topo.hd];

                        double & ndist = maxDistanceToNearestComponents[r.topo.hd];

                        // region region
                        for (auto rrch : r.topo.constraints<RegionBoundaryData>()) {
                            RegionHandle thatRh = mg.topo(rrch).component<0>();
                            if (thatRh == r.topo.hd) {
                                thatRh = mg.topo(rrch).component<1>();
                            }
                            Plane3 & thatPlane = planes[thatRh];
                            for (auto & cs : mg.data(rrch).normalizedSampledPoints) {
                                for (auto & c : cs) {
                                    Point3 pointHere = IntersectionOfLineAndPlane(Ray3(Point3(0, 0, 0), c), thisPlane).position;
                                    Point3 pointThere = IntersectionOfLineAndPlane(Ray3(Point3(0, 0, 0), c), thatPlane).position;
                                    double dist = Distance(pointHere, pointThere);
                                    if (dist > ndist) {
                                        ndist = dist;
                                    }
                                }
                            }
                        }

                        // region line
                        for (auto & rlch : r.topo.constraints<RegionLineConnectionData>()) {
                            LineHandle lh = mg.topo(rlch).component<1>();
                            Line3 & thatLine = lines[lh];
                            for (auto & a : mg.data(rlch).normalizedAnchors) {
                                Point3 pointHere = IntersectionOfLineAndPlane(Ray3(Point3(0, 0, 0), a), thisPlane).position;
                                Point3 pointThere = DistanceBetweenTwoLines(Ray3(Point3(0, 0, 0), a), thatLine.ray()).second.first;
                                double dist = Distance(pointHere, pointThere);
                                if (dist > ndist) {
                                    ndist = dist;
                                }
                            }
                        }
                    }

                    std::sort(usableRegionHandles.begin(), usableRegionHandles.end(), [&maxDistanceToNearestComponents](RegionHandle l1, RegionHandle l2) {
                        return maxDistanceToNearestComponents[l1] > maxDistanceToNearestComponents[l2];
                    });

                    // loose orientation constraints on these lines
                    for (int i = 0; i < usableRegionHandles.size() * regionsLoosableRatio; i++) {
                        double maxCornerDist = maxDistanceToNearestComponents[usableRegionHandles[i]];
                        if (maxCornerDist < distThres)
                            continue;
                        auto & rp = controls[usableRegionHandles[i]];
                        assert(rp.used);
                        if (rp.orientationClaz == -1 && rp.orientationNotClaz == -1) {
                            rp.used = false;
                        } else {
                            rp.orientationClaz = rp.orientationNotClaz = -1;
                        }
                    }
                }
            } else {

                // collect keypoints
                struct ComponentKeyPoint {
                    Point3 point;
                    int handleId;
                    bool isOnRegion;
                    inline ComponentKeyPoint() {}
                    inline ComponentKeyPoint(const Point3 & p, RegionHandle rh) : point(p), isOnRegion(true), handleId(rh.id) {}
                    inline ComponentKeyPoint(const Point3 & p, LineHandle lh) : point(p), isOnRegion(false), handleId(lh.id) {}
                    inline bool matches(RegionHandle rh) const { return isOnRegion && handleId == rh.id; }
                    inline bool matches(LineHandle lh) const { return !isOnRegion && handleId == lh.id; }
                };
                auto getBB = [distThres](const ComponentKeyPoint & p) {return BoundingBox(p.point).expand(distThres); };
                RTreeWrapper<ComponentKeyPoint, decltype(getBB)> keyPointsRTree(getBB);

                // collect current instances
                HandledTable<LineHandle, Line3> lines(mg.internalComponents<LineData>().size());
                HandledTable<RegionHandle, Plane3> planes(mg.internalComponents<RegionData>().size());

                // keypoints on regions
                for (auto & r : mg.components<RegionData>()) {
                    if (!controls[r.topo.hd].used)
                        continue;
                    auto plane = Instance(mg, controls, vars, r.topo.hd);
                    planes[r.topo.hd] = plane;
                    for (auto & cs : r.data.normalizedContours) {
                        for (auto & c : cs) {
                            double depth = DepthAt(c, plane);
                            keyPointsRTree.insert(ComponentKeyPoint(c * depth, r.topo.hd));
                        }
                    }
                }
                // keypoints on lines
                for (auto & l : mg.components<LineData>()) {
                    if (!controls[l.topo.hd].used)
                        continue;
                    auto line = Instance(mg, controls, vars, l.topo.hd);
                    lines[l.topo.hd] = line;
                    keyPointsRTree.insert(ComponentKeyPoint(line.first, l.topo.hd));
                    keyPointsRTree.insert(ComponentKeyPoint(line.second, l.topo.hd));
                }

                // find most isolated lines
                if (linesLoosableRatio > 0.0) {

                    auto lineMaxCornerNNDistances = mg.createComponentTable<LineData, double>();
                    std::vector<LineHandle> usableLineHandles;
                    for (auto & l : mg.components<LineData>()) {
                        if (!controls[l.topo.hd].used)
                            continue;
                        usableLineHandles.push_back(l.topo.hd);
                        auto & line = lines[l.topo.hd];
                        double nnDistances[] = { std::numeric_limits<double>::max(), std::numeric_limits<double>::max() };
                        keyPointsRTree.search(BoundingBox(line.first).expand(nnSearchRange), [&l, &line, &nnDistances](const ComponentKeyPoint & p) {
                            if (p.matches(l.topo.hd))
                                return true;
                            if (!p.isOnRegion)
                                return true;
                            nnDistances[0] = std::min(Distance(line.first, p.point), nnDistances[0]);
                            //nnDistances[1] = std::min(Distance(line.second, p.point), nnDistances[1]);
                            return true;
                        });
                        keyPointsRTree.search(BoundingBox(line.second).expand(nnSearchRange), [&l, &line, &nnDistances](const ComponentKeyPoint & p) {
                            if (p.matches(l.topo.hd))
                                return true;
                            //nnDistances[0] = std::min(Distance(line.first, p.point), nnDistances[0]);
                            nnDistances[1] = std::min(Distance(line.second, p.point), nnDistances[1]);
                            return true;
                        });
                        lineMaxCornerNNDistances[l.topo.hd] = std::max(nnDistances[0], nnDistances[1]);
                    }

                    std::sort(usableLineHandles.begin(), usableLineHandles.end(), [&lineMaxCornerNNDistances](LineHandle l1, LineHandle l2) {
                        return lineMaxCornerNNDistances[l1] > lineMaxCornerNNDistances[l2];
                    });

                    // loose orientation constraints on these lines
                    for (int i = 0; i < usableLineHandles.size() * linesLoosableRatio; i++) {
                        double maxCornerDist = lineMaxCornerNNDistances[usableLineHandles[i]];
                        if (maxCornerDist < distThres)
                            continue;
                        auto & lp = controls[usableLineHandles[i]];
                        assert(lp.used);
                        if (lp.orientationClaz >= 0) {
                            lp.orientationClaz = -1;
                        }
                        assert(lp.orientationNotClaz == -1);
                        if (lp.orientationClaz == -1) {
                            lp.used = false; // disable this line!!!!!!!!!
                        }
                    }
                }

                // find most isolated regions
                if (regionsLoosableRatio > 0.0) {

                    auto regionMaxCornerNNDistances = mg.createComponentTable<RegionData, double>();
                    std::vector<RegionHandle> usableRegionHandles;
                    for (auto & r : mg.components<RegionData>()) {
                        if (!controls[r.topo.hd].used)
                            continue;
                        usableRegionHandles.push_back(r.topo.hd);
                        auto & plane = planes[r.topo.hd];
                        double maxCornerNNDistance = 0.0;
                        for (int i = 0; i < r.data.normalizedContours.size(); i++) {
                            auto & contour = r.data.normalizedContours[i];
                            for (int j = 0; j < contour.size(); j++) {
                                auto c = contour[j];
                                double depth = DepthAt(c, plane);
                                c *= depth; // the corner position
                                auto id = std::make_pair(i, j); // the corner id
                                double cornerNNDistance = std::numeric_limits<double>::max(); // find nearest neighbor distance to this corner
                                keyPointsRTree.search(BoundingBox(c).expand(nnSearchRange), [&r, &plane, &cornerNNDistance, id, &c](const ComponentKeyPoint & p) {
                                    if (p.matches(r.topo.hd))
                                        return true;
                                    double d = Distance(c, p.point);
                                    cornerNNDistance = std::min(cornerNNDistance, d);
                                    return true;
                                });
                                if (cornerNNDistance > maxCornerNNDistance) {
                                    maxCornerNNDistance = cornerNNDistance; // get the maximum nn distance of all corners on this region
                                }
                            }
                        }
                        regionMaxCornerNNDistances[r.topo.hd] = maxCornerNNDistance;
                    }

                    std::sort(usableRegionHandles.begin(), usableRegionHandles.end(), [&regionMaxCornerNNDistances](RegionHandle l1, RegionHandle l2) {
                        return regionMaxCornerNNDistances[l1] > regionMaxCornerNNDistances[l2];
                    });

                    // loose orientation constraints on these lines
                    for (int i = 0; i < usableRegionHandles.size() * regionsLoosableRatio; i++) {
                        double maxCornerDist = regionMaxCornerNNDistances[usableRegionHandles[i]];
                        if (maxCornerDist < distThres)
                            continue;
                        auto & rp = controls[usableRegionHandles[i]];
                        assert(rp.used);
                        rp.orientationClaz = rp.orientationNotClaz = -1;
                        //if (lp.orientationClaz == -1){
                        //    lp.used = false; // disable this line!!!!!!!!!
                        //}
                    }
                }


            }

            controls.disableAllInvalidConstraints(mg);

        }









        double MedianCenterDepth(const RLGraph & mg, const RLGraphControls & controls,
            const RLGraphVars & vars) {

            std::vector<double> centerDepths;
            centerDepths.reserve(mg.internalComponents<RegionData>().size() + mg.internalComponents<LineData>().size());
            for (auto & c : mg.components<RegionData>()) {
                if (!controls[c.topo.hd].used)
                    continue;
                double d = DepthAt(c.data.normalizedCenter, Instance(mg, controls, vars, c.topo.hd));
                if (!IsInfOrNaN(d)) {
                    centerDepths.push_back(d);
                }
                //assert(!IsInfOrNaN(centerDepths.back()));
            }
            for (auto & c : mg.components<LineData>()) {
                if (!controls[c.topo.hd].used)
                    continue;
                double d = DepthAt(normalize(c.data.line.center()), Instance(mg, controls, vars, c.topo.hd));
                if (!IsInfOrNaN(d)) {
                    centerDepths.push_back(d);
                }
                //assert(!IsInfOrNaN(centerDepths.back()));
            }
            std::nth_element(centerDepths.begin(), centerDepths.begin() + centerDepths.size() / 2, centerDepths.end());
            double medianCenterDepth = centerDepths[centerDepths.size() / 2];
            return medianCenterDepth;
        }


        double Score(const RLGraph & mg, const RLGraphControls & controls,
            const RLGraphVars & vars) {

            // manhattan fitness
            double sumOfComponentWeightedFitness = 0.0;
            double sumOfComponentWeights = 0.0;

            // lines
            double maxLineSpanAngle = 0.0;
            for (auto & l : mg.components<LineData>()) {
                if (!controls[l.topo.hd].used)
                    continue;
                double lineSpanAngle = AngleBetweenDirections(l.data.line.first, l.data.line.second);
                if (lineSpanAngle > maxLineSpanAngle)
                    maxLineSpanAngle = lineSpanAngle;
            }

            HandledTable<LineHandle, Line3> lines(mg.internalComponents<LineData>().size());
            HandledTable<RegionHandle, Plane3> planes(mg.internalComponents<RegionData>().size());

            for (auto & l : mg.components<LineData>()) {
                if (!controls[l.topo.hd].used)
                    continue;
                lines[l.topo.hd] = Instance(mg, controls, vars, l.topo.hd);

                double fitness = 0.0;
                double weight = 0.0;
                double lineSpanAngle = AngleBetweenDirections(l.data.line.first, l.data.line.second);
                auto & prop = controls[l.topo.hd];
                if (prop.orientationClaz >= 0) {
                    fitness = 1.0;
                    weight = lineSpanAngle / maxLineSpanAngle;
                } else {
                    static const double angleThreshold = DegreesToRadians(10);
                    double minAngle = angleThreshold;
                    int bestVPId = -1;
                    for (int i = 0; i < controls.vanishingPoints.size(); i++) {
                        double angle = AngleBetweenUndirectedVectors(controls.vanishingPoints[i], lines[l.topo.hd].direction());
                        if (angle < minAngle) {
                            minAngle = angle;
                            bestVPId = i;
                        }
                    }
                    fitness = Gaussian(minAngle / angleThreshold, 0.1);
                    weight = lineSpanAngle / maxLineSpanAngle;
                }
                sumOfComponentWeightedFitness += (fitness * weight);
                sumOfComponentWeights += weight;
            }

            // regions
            double maxRegionArea = 0.0;
            for (auto & r : mg.components<RegionData>()) {
                if (!controls[r.topo.hd].used)
                    continue;
                if (r.data.area > maxRegionArea) {
                    maxRegionArea = r.data.area;
                }
            }

            for (auto & r : mg.components<RegionData>()) {
                if (!controls[r.topo.hd].used)
                    continue;
                planes[r.topo.hd] = Instance(mg, controls, vars, r.topo.hd);

                double fitness = 0.0;
                double weight = 1.0;
                auto & prop = controls[r.topo.hd];
                if (prop.orientationClaz >= 0) {
                    fitness = 1.0;
                    weight = r.data.area / maxRegionArea;
                } else {
                    static const double angleThreshold = DegreesToRadians(10);
                    double minAngle = angleThreshold;
                    int bestVPId = -1;
                    for (int i = 0; i < controls.vanishingPoints.size(); i++) {
                        double angle = AngleBetweenUndirectedVectors(controls.vanishingPoints[i], planes[r.topo.hd].normal);
                        if (angle < minAngle) {
                            minAngle = angle;
                            bestVPId = i;
                        }
                    }
                    fitness = Gaussian(minAngle / angleThreshold, 0.1);
                    weight = r.data.area / maxRegionArea;
                }
                sumOfComponentWeightedFitness += (fitness * weight);
                sumOfComponentWeights += weight;
            }


            // constraint fitness
            double sumOfConstraintWeightedFitness = 0.0;
            double sumOfConstraintWeights = 0.0;
            double sumOfNotUsedConstraintWeights = 0.0;
            for (auto & c : mg.constraints<RegionBoundaryData>()) {
                static const double typeWeight = 1.0;
                if (!controls[c.topo.hd].used) {
                    sumOfNotUsedConstraintWeights += ElementsNum(c.data.normalizedSampledPoints) * typeWeight;
                    continue;
                }

                auto & plane1 = planes[c.topo.component<0>()];
                auto & plane2 = planes[c.topo.component<1>()];
                for (auto & ss : c.data.normalizedSampledPoints) {
                    for (auto & s : ss) {
                        double d1 = DepthAt(s, plane1);
                        double d2 = DepthAt(s, plane2);
                        if (IsInfOrNaN(d1) || IsInfOrNaN(d2)) {
                            std::cout << "nan/inf plane!" << std::endl;
                            continue;
                        }
                        assert(d1 >= 0 && d2 >= 0);

                        double fitness = Gaussian(abs(d1 - d2) / (d1 + d2), 0.1);
                        double weight = typeWeight;

                        sumOfConstraintWeightedFitness += (fitness * weight);
                        sumOfConstraintWeights += weight;
                    }
                }
            }
            //double maxJunctionWeight = 0.0;
            //for (auto & c : mg.constraints<LineRelationData>()){
            //    if (c.data.junctionWeight > maxJunctionWeight){
            //        maxJunctionWeight = c.data.junctionWeight;
            //    }
            //}
            for (auto & c : mg.constraints<LineRelationData>()) {
                static const double typeWeight = 8.0;
                if (!controls[c.topo.hd].used) {
                    sumOfNotUsedConstraintWeights += /* c.data.junctionWeight / maxJunctionWeight **/ typeWeight;
                    continue;
                }

                auto & line1 = lines[c.topo.component<0>()];
                auto & line2 = lines[c.topo.component<1>()];
                double d1 = DepthAt(c.data.normalizedRelationCenter, line1);
                double d2 = DepthAt(c.data.normalizedRelationCenter, line2);
                if (IsInfOrNaN(d1) || IsInfOrNaN(d2)) {
                    std::cout << "nan/inf line!" << std::endl;
                    continue;
                }

                double fitness = Gaussian(abs(d1 - d2) / (d1 + d2), 0.1);
                double weight = typeWeight;

                sumOfConstraintWeightedFitness += (fitness * weight);
                sumOfConstraintWeights += weight;
            }
            for (auto & c : mg.constraints<RegionLineConnectionData>()) {
                static const double typeWeight = 1.0;
                if (!controls[c.topo.hd].used) {
                    sumOfNotUsedConstraintWeights += ElementsNum(c.data.normalizedAnchors) * typeWeight;
                    continue;
                }

                auto & plane = planes[c.topo.component<0>()];
                auto & line = lines[c.topo.component<1>()];
                for (auto & a : c.data.normalizedAnchors) {
                    double d1 = DepthAt(a, plane);
                    double d2 = DepthAt(a, line);
                    if (IsInfOrNaN(d1) || IsInfOrNaN(d2)) {
                        std::cout << "nan/inf regionline!" << std::endl;
                        continue;
                    }
                    assert(d1 >= 0 && d2 >= 0);

                    double fitness = Gaussian(abs(d1 - d2) / (d1 + d2), 0.1);
                    double weight = typeWeight;

                    sumOfConstraintWeightedFitness += (fitness * weight);
                    sumOfConstraintWeights += weight;
                }
            }

            double score = sumOfComponentWeightedFitness / sumOfComponentWeights * 50
                + sumOfConstraintWeightedFitness / sumOfConstraintWeights * 1000
                - sumOfNotUsedConstraintWeights / (sumOfConstraintWeights + sumOfNotUsedConstraintWeights) * 5;

            return score;
        }





        void NormalizeVariables(const RLGraph & mg, const RLGraphControls & controls,
            RLGraphVars & vars) {
            double medianCenterDepth = MedianCenterDepth(mg, controls, vars);

            SetClock();

            std::cout << "median center depth: " << medianCenterDepth << std::endl;

            assert(!IsInfOrNaN(medianCenterDepth));

            // normalize variables
            for (auto & c : mg.components<RegionData>()) {
                if (!controls[c.topo.hd].used)
                    continue;
                for (int i = 0; i < vars[c.topo.hd].variables.size(); i++) {
                    vars[c.topo.hd].variables[i] *= medianCenterDepth;
                }
            }
            for (auto & c : mg.components<LineData>()) {
                if (!controls[c.topo.hd].used)
                    continue;
                for (int i = 0; i < vars[c.topo.hd].variables.size(); i++) {
                    vars[c.topo.hd].variables[i] *= medianCenterDepth;
                }
            }
        }


















    }
}