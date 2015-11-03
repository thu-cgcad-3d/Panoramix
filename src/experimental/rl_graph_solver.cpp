#include "../misc/eigen.hpp"
#include "../misc/clock.hpp"

#include "rl_graph_solver.hpp"


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




            int RegisterVariablePositions(const RLGraph & mg,
                const RLGraphControls & controls, const RLGraphVars & vars,
                RLGraphComponentTable<int> & uh2varStartPosition) {
                int varNum = 0;
                for (auto & c : mg.components<LineData>()) {
                    if (!controls[c.topo.hd].used)
                        continue;
                    uh2varStartPosition[c.topo.hd] = varNum;
                    varNum += vars[c.topo.hd].variables.size();
                }
                for (auto & c : mg.components<RegionData>()) {
                    if (!controls[c.topo.hd].used)
                        continue;
                    uh2varStartPosition[c.topo.hd] = varNum;
                    varNum += vars[c.topo.hd].variables.size();
                }
                return varNum;
            }


            void FillVariableValues(const RLGraph & mg,
                const RLGraphControls & controls, const RLGraphVars & vars,
                const RLGraphComponentTable<int> & uh2varStartPosition,
                std::vector<double> & X) {
                for (auto & c : mg.components<LineData>()) {
                    if (!controls[c.topo.hd].used)
                        continue;
                    for (int i = 0; i < vars[c.topo.hd].variables.size(); i++) {
                        X[uh2varStartPosition.at(c.topo.hd) + i] = vars[c.topo.hd].variables.at(i);
                    }
                }
                for (auto & c : mg.components<RegionData>()) {
                    if (!controls[c.topo.hd].used)
                        continue;
                    for (int i = 0; i < vars[c.topo.hd].variables.size(); i++) {
                        X[uh2varStartPosition.at(c.topo.hd) + i] = vars[c.topo.hd].variables.at(i);
                    }
                }
            }






            enum EquationType {
                OnComponentWeightedAnchor = 0,
                OnComponentBoundedAnchor = 1,
                OnConstraintWeightedAnchor = 2,
                OnConstraintBoundedAnchor = 3
            };

            // register component positions in equations
            // register constraint positions in equations 
            // collect applied constraint anchors
            // layout:
            //  [component1 weighted anchors][component1 bounded anchors]  [component2 ...] ...
            //  [constraint1 weighted anchors][constraint1 bounded anchors]  [constraint2 ...] ...
            // type order:
            //  Region, Line
            //  RegionBoundary, LineRelation, RegionLineConnection
            std::vector<EquationType> RegisterEquationPositions(const RLGraph & mg,
                const RLGraphControls & controls,
                RLGraphComponentTable<int> & uh2eqStartPosition,
                RLGraphConstraintTable<int> & bh2eqStartPosition) {

                std::vector<EquationType> equationTypes;

                int consNum = 0;

                // register component constraint positions (weighted anchors, bounded anchors)
                for (auto & c : mg.components<RegionData>()) {
                    if (!controls[c.topo.hd].used)
                        continue;
                    uh2eqStartPosition[c.topo.hd] = consNum;
                    auto & was = controls[c.topo.hd].weightedAnchors;
                    consNum += was.size();
                    equationTypes.insert(equationTypes.end(), was.size(), OnComponentWeightedAnchor);
                    auto & bas = controls[c.topo.hd].boundedAnchors;
                    consNum += bas.size();
                    equationTypes.insert(equationTypes.end(), bas.size(), OnComponentBoundedAnchor);
                }
                for (auto & c : mg.components<LineData>()) {
                    if (!controls[c.topo.hd].used)
                        continue;
                    uh2eqStartPosition[c.topo.hd] = consNum;
                    auto & was = controls[c.topo.hd].weightedAnchors;
                    consNum += was.size();
                    equationTypes.insert(equationTypes.end(), was.size(), OnComponentWeightedAnchor);
                    auto & bas = controls[c.topo.hd].boundedAnchors;
                    consNum += bas.size();
                    equationTypes.insert(equationTypes.end(), bas.size(), OnComponentBoundedAnchor);
                }

                // register constraint positions (weighted anchors, bounded anchors)
                for (auto & c : mg.constraints<RegionBoundaryData>()) {
                    if (!controls[c.topo.hd].used)
                        continue;
                    if (!(controls[c.topo.component<0>()].used && controls[c.topo.component<1>()].used))
                        continue;
                    bh2eqStartPosition[c.topo.hd] = consNum;
                    auto & was = controls[c.topo.hd].weightedAnchors;
                    consNum += was.size();
                    equationTypes.insert(equationTypes.end(), was.size(), OnConstraintWeightedAnchor);
                    auto & bas = controls[c.topo.hd].boundedAnchors;
                    consNum += bas.size();
                    equationTypes.insert(equationTypes.end(), bas.size(), OnConstraintBoundedAnchor);
                }
                for (auto & c : mg.constraints<LineRelationData>()) {
                    if (!controls[c.topo.hd].used)
                        continue;
                    if (!(controls[c.topo.component<0>()].used && controls[c.topo.component<1>()].used))
                        continue;
                    bh2eqStartPosition[c.topo.hd] = consNum;
                    auto & was = controls[c.topo.hd].weightedAnchors;
                    consNum += was.size();
                    equationTypes.insert(equationTypes.end(), was.size(), OnConstraintWeightedAnchor);
                    auto & bas = controls[c.topo.hd].boundedAnchors;
                    consNum += bas.size();
                    equationTypes.insert(equationTypes.end(), bas.size(), OnConstraintBoundedAnchor);
                }
                for (auto & c : mg.constraints<RegionLineConnectionData>()) {
                    if (!controls[c.topo.hd].used)
                        continue;
                    if (!(controls[c.topo.component<0>()].used && controls[c.topo.component<1>()].used))
                        continue;
                    bh2eqStartPosition[c.topo.hd] = consNum;
                    auto & was = controls[c.topo.hd].weightedAnchors;
                    consNum += was.size();
                    equationTypes.insert(equationTypes.end(), was.size(), OnConstraintWeightedAnchor);
                    auto & bas = controls[c.topo.hd].boundedAnchors;
                    consNum += bas.size();
                    equationTypes.insert(equationTypes.end(), bas.size(), OnConstraintBoundedAnchor);
                }
                assert(consNum == equationTypes.size());
                return equationTypes;
            }


            // fill A, W, B for component weighted anchors
            // fill A, B1, B2 for component bounded anchors
            template <class ComponentDataT, class SparseMatElementT>
            void FillComponentAnchorEquationMatrices(const RLGraph & mg,
                const RLGraphControls & controls,
                const RLGraphVars & vars,
                const RLGraphComponentTable<int> & uh2varStartPosition,
                const RLGraphComponentTable<int> & uh2eqStartPosition,
                std::vector<SparseMatElementT> & Atriplets,
                std::vector<SparseMatElementT> & Wtriplets,
                std::vector<double> & B,
                std::vector<double> & B1,
                std::vector<double> & B2) {
                for (auto & c : mg.components<ComponentDataT>()) {
                    if (!controls[c.topo.hd].used)
                        continue;
                    auto uh = c.topo.hd;
                    int uhVarNum = vars[uh].variables.size();
                    int uhVarStartPosition = uh2varStartPosition.at(uh);
                    int eid = uh2eqStartPosition[uh];
                    // fill weighted anchors
                    for (auto & wa : controls[c.topo.hd].weightedAnchors) {
                        const Point3 & anchor = wa.component;
                        double weight = wa.score;
                        auto uhVarCoeffsAtAnchorDirection = VariableCoefficientsForInverseDepthAtDirection(mg, controls, vars, anchor, uh);
                        assert(uhVarCoeffsAtAnchorDirection.size() == uhVarNum);
                        for (int i = 0; i < uhVarCoeffsAtAnchorDirection.size(); i++) {
                            //A.insert(eid, uhVarStartPosition + i) = uhVarCoeffsAtCenter[i];
                            Atriplets.emplace_back(eid, uhVarStartPosition + i, uhVarCoeffsAtAnchorDirection[i]);
                        }
                        B[eid] = 1.0 / norm(anchor);
                        B1[eid] = 0.0;
                        B2[eid] = 0.0;
                        Wtriplets.emplace_back(eid, eid, weight);
                        eid++;
                    }
                    // fill bounded anchors
                    for (auto & ba : controls[c.topo.hd].boundedAnchors) {
                        const Point3 & anchor = ba.component;
                        double lb = ba.lowerBound();
                        double ub = ba.upperBound();
                        assert(lb > 0 && ub > 0);
                        auto uhVarCoeffsAtAnchorDirection = VariableCoefficientsForInverseDepthAtDirection(mg, controls, vars, anchor, uh);
                        assert(uhVarCoeffsAtAnchorDirection.size() == uhVarNum);
                        for (int i = 0; i < uhVarCoeffsAtAnchorDirection.size(); i++) {
                            //A.insert(eid, uhVarStartPosition + i) = uhVarCoeffsAtCenter[i];
                            Atriplets.emplace_back(eid, uhVarStartPosition + i, uhVarCoeffsAtAnchorDirection[i]);
                        }
                        B[eid] = 0.0;
                        B1[eid] = 1.0 / ub;
                        B2[eid] = 1.0 / lb;
                        Wtriplets.emplace_back(eid, eid, 1.0);
                        eid++;
                    }
                }
            }


            // fill A, W, B for constraint weighted anchors
            // fill A, B1, B2 for constraint bounded anchors
            template <class ConstraintDataT, class SparseMatElementT>
            void FillConstraintAnchorEquationMatrices(const RLGraph & mg,
                const RLGraphControls & controls,
                const RLGraphVars & vars,
                const RLGraphComponentTable<int> & uh2varStartPosition,
                const RLGraphConstraintTable<int> & bh2eqStartPosition,
                std::vector<SparseMatElementT> & Atriplets,
                std::vector<SparseMatElementT> & Wtriplets,
                std::vector<double> & B,
                std::vector<double> & B1,
                std::vector<double> & B2) {

                for (auto & c : mg.constraints<ConstraintDataT>()) {
                    if (!controls[c.topo.hd].used)
                        continue;
                    if (!(controls[c.topo.component<0>()].used && controls[c.topo.component<1>()].used))
                        continue;
                    auto bh = c.topo.hd;
                    auto uh1 = mg.topo(bh).component<0>();
                    auto uh2 = mg.topo(bh).component<1>();
                    assert(controls[uh1].used || controls[uh2].used);
                    auto & u1 = mg.data(uh1);
                    auto & u2 = mg.data(uh2);

                    int u1VarStartPosition = uh2varStartPosition.at(uh1);
                    int u1VarNum = vars[uh1].variables.size();

                    int u2VarStartPosition = uh2varStartPosition.at(uh2);
                    int u2VarNum = vars[uh2].variables.size();

                    int eid = bh2eqStartPosition[bh];

                    auto & was = controls[bh].weightedAnchors;
                    auto & bas = controls[bh].boundedAnchors;

                    // weighted anchors
                    for (auto & wa : was) {
                        double weight = wa.weight();
                        const Vec3 & a = wa.component;
                        B[eid] = 0.0;
                        B1[eid] = 0.0;
                        B2[eid] = 0.0;
                        Wtriplets.emplace_back(eid, eid, weight);
                        {
                            auto u1VarCoeffs = VariableCoefficientsForInverseDepthAtDirection(mg, controls, vars, a, uh1);
                            assert(u1VarCoeffs.size() == u1VarNum);
                            for (int i = 0; i < u1VarCoeffs.size(); i++) {
                                //A.insert(eid, u1VarStartPosition + i) = u1VarCoeffs[i]; // pos
                                Atriplets.emplace_back(eid, u1VarStartPosition + i, u1VarCoeffs[i]);
                            }
                        }
                        {
                            auto u2VarCoeffs = VariableCoefficientsForInverseDepthAtDirection(mg, controls, vars, a, uh2);
                            assert(u2VarCoeffs.size() == u2VarNum);
                            for (int i = 0; i < u2VarCoeffs.size(); i++) {
                                //A.insert(eid, u2VarStartPosition + i) = -u2VarCoeffs[i]; // neg
                                Atriplets.emplace_back(eid, u2VarStartPosition + i, -u2VarCoeffs[i]);
                            }
                        }
                        eid++;
                    }

                    // bounded anchors
                    for (auto & ba : bas) {
                        double lb = ba.lowerBound();
                        double ub = ba.upperBound();
                        const Vec3 & a = ba.component;
                        B[eid] = 0.0;
                        B1[eid] = -ub;
                        B2[eid] = -lb;
                        Wtriplets.emplace_back(eid, eid, 1.0);
                        {
                            auto u1VarCoeffs = VariableCoefficientsForInverseDepthAtDirection(mg, controls, vars, a, uh1);
                            assert(u1VarCoeffs.size() == u1VarNum);
                            for (int i = 0; i < u1VarCoeffs.size(); i++) {
                                //A.insert(eid, u1VarStartPosition + i) = u1VarCoeffs[i]; // pos
                                Atriplets.emplace_back(eid, u1VarStartPosition + i, u1VarCoeffs[i]);
                            }
                        }
                        {
                            auto u2VarCoeffs = VariableCoefficientsForInverseDepthAtDirection(mg, controls, vars, a, uh2);
                            assert(u2VarCoeffs.size() == u2VarNum);
                            for (int i = 0; i < u2VarCoeffs.size(); i++) {
                                //A.insert(eid, u2VarStartPosition + i) = -u2VarCoeffs[i]; // neg
                                Atriplets.emplace_back(eid, u2VarStartPosition + i, -u2VarCoeffs[i]);
                            }
                        }
                        eid++;
                    }
                }
            }



            template <class SparseMatElementT>
            std::pair<int, int> PrepareEverything(const RLGraph & mg,
                const RLGraphControls & controls,
                const RLGraphVars & vars,
                RLGraphComponentTable<int> & uh2varStartPosition,
                RLGraphComponentTable<int> & uh2eqStartPosition,
                RLGraphConstraintTable<int> & bh2eqStartPosition,
                std::vector<double> & Xdata,
                std::vector<SparseMatElementT> & Atriplets,
                std::vector<SparseMatElementT> & Wtriplets,
                std::vector<double> & Bdata,
                std::vector<double> & B1data,
                std::vector<double> & B2data,
                std::vector<EquationType> & equationTypes) {

                int varNum = RegisterVariablePositions(mg, controls, vars, uh2varStartPosition);
                Xdata.resize(varNum, 1.0);
                FillVariableValues(mg, controls, vars, uh2varStartPosition, Xdata);

                // register cons
                equationTypes = RegisterEquationPositions(mg, controls, uh2eqStartPosition, bh2eqStartPosition);
                int consNum = equationTypes.size();

                Bdata.resize(consNum, 0.0);
                B1data.resize(consNum, 0.0);
                B2data.resize(consNum, 0.0);
                Atriplets.reserve(consNum * 6);
                Wtriplets.reserve(consNum);

                FillComponentAnchorEquationMatrices<RegionData>(mg, controls, vars,
                    uh2varStartPosition, uh2eqStartPosition, Atriplets, Wtriplets, Bdata, B1data, B2data);
                FillComponentAnchorEquationMatrices<LineData>(mg, controls, vars,
                    uh2varStartPosition, uh2eqStartPosition, Atriplets, Wtriplets, Bdata, B1data, B2data);

                FillConstraintAnchorEquationMatrices<RegionBoundaryData>(mg, controls, vars,
                    uh2varStartPosition, bh2eqStartPosition, Atriplets, Wtriplets, Bdata, B1data, B2data);
                FillConstraintAnchorEquationMatrices<LineRelationData>(mg, controls, vars,
                    uh2varStartPosition, bh2eqStartPosition, Atriplets, Wtriplets, Bdata, B1data, B2data);
                FillConstraintAnchorEquationMatrices<RegionLineConnectionData>(mg, controls, vars,
                    uh2varStartPosition, bh2eqStartPosition, Atriplets, Wtriplets, Bdata, B1data, B2data);

                return std::make_pair(varNum, consNum);
            }



            template <class VectorT>
            void InstallVariables(const RLGraph & mg, const RLGraphControls & controls,
                const RLGraphComponentTable<int> & uh2varStartPosition, const VectorT & X,
                RLGraphVars & vars) {
                for (auto & c : mg.components<RegionData>()) {
                    if (!controls[c.topo.hd].used)
                        continue;
                    int uhStartPosition = uh2varStartPosition.at(c.topo.hd);
                    for (int i = 0; i < vars[c.topo.hd].variables.size(); i++) {
                        vars[c.topo.hd].variables[i] = X(uhStartPosition + i);
                    }
                }
                for (auto & c : mg.components<LineData>()) {
                    if (!controls[c.topo.hd].used)
                        continue;
                    int uhStartPosition = uh2varStartPosition.at(c.topo.hd);
                    for (int i = 0; i < vars[c.topo.hd].variables.size(); i++) {
                        vars[c.topo.hd].variables[i] = X(uhStartPosition + i);
                    }
                }
            }

            void InstallVariables(const RLGraph & mg, const RLGraphControls & controls,
                const RLGraphComponentTable<int> & uh2varStartPosition, const std::vector<double> & X,
                RLGraphVars & vars) {
                for (auto & c : mg.components<RegionData>()) {
                    if (!controls[c.topo.hd].used)
                        continue;
                    int uhStartPosition = uh2varStartPosition.at(c.topo.hd);
                    for (int i = 0; i < vars[c.topo.hd].variables.size(); i++) {
                        vars[c.topo.hd].variables[i] = X[uhStartPosition + i];
                    }
                }
                for (auto & c : mg.components<LineData>()) {
                    if (!controls[c.topo.hd].used)
                        continue;
                    int uhStartPosition = uh2varStartPosition.at(c.topo.hd);
                    for (int i = 0; i < vars[c.topo.hd].variables.size(); i++) {
                        vars[c.topo.hd].variables[i] = X[uhStartPosition + i];
                    }
                }
            }

        }






        //// inverse depth optimization
        //RLGraphVars SolveVariablesWithoutBoundedAnchors(const RLGraph & mg,
        //    const RLGraphControls & controls,
        //    bool useWeights) {

        //    assert(NumberOfComponentWeightedAnchors(controls) +
        //        NumberOfComponentBoundedAnchors(controls) > 0);

        //    RLGraphVars vars = MakeVariables(mg, controls);

        //    SetClock();

        //    auto uh2varStartPosition = MakeHandledTableForAllComponents<int>(mg);
        //    auto uh2eqStartPosition = MakeHandledTableForAllComponents<int>(mg);
        //    auto bh2eqStartPosition = MakeHandledTableForAllConstraints<int>(mg);

        //    int varNum, consNum;

        //    std::vector<EquationType> equationTypes;
        //    std::vector<double> Xdata;
        //    std::vector<double> Bdata;
        //    std::vector<double> B1data, B2data;
        //    std::vector<Eigen::Triplet<double>> Atriplets;
        //    std::vector<Eigen::Triplet<double>> Wtriplets;

        //    std::tie(varNum, consNum) = PrepareEverything(mg, controls, vars, uh2varStartPosition,
        //        uh2eqStartPosition, bh2eqStartPosition, Xdata, Atriplets, Wtriplets, Bdata, B1data, B2data, equationTypes);


        //    // matrices
        //    Eigen::SparseMatrix<double> A;
        //    {
        //        Clock clock("form matrix A");
        //        A.resize(consNum, varNum);
        //        A.setFromTriplets(Atriplets.begin(), Atriplets.end());
        //    }
        //    Eigen::Map<const Eigen::VectorXd> B(Bdata.data(), Bdata.size());


        //    Eigen::SparseMatrix<double> WA;
        //    Eigen::VectorXd WB;
        //    if (useWeights) {
        //        Clock clock("form matrix WA");
        //        Eigen::SparseMatrix<double> W;
        //        W.resize(consNum, consNum);
        //        W.setFromTriplets(Wtriplets.begin(), Wtriplets.end());
        //        WA = W * A;
        //        WB = W * B;
        //    }

        //    static const bool useSPQR = true;

        //    Eigen::VectorXd X;
        //    if (!useSPQR) {
        //        Clock clock("solve equations using Eigen::SparseQR");

        //        A.makeCompressed();
        //        WA.makeCompressed();

        //        Eigen::SparseQR<Eigen::SparseMatrix<double>, Eigen::COLAMDOrdering<int>> solver;
        //        static_assert(!(Eigen::SparseMatrix<double>::IsRowMajor), "COLAMDOrdering only supports column major");

        //        solver.compute(useWeights ? WA : A);

        //        if (solver.info() != Eigen::Success) {
        //            assert(0);
        //            std::cout << "computation error" << std::endl;
        //            SHOULD_NEVER_BE_CALLED();
        //        }

        //        X = solver.solve(useWeights ? WB : B);
        //        if (solver.info() != Eigen::Success) {
        //            assert(0);
        //            std::cout << "solving error" << std::endl;
        //            SHOULD_NEVER_BE_CALLED();
        //        }
        //    } else {
        //        Clock clock("solve equations using Eigen::SPQR");

        //        Eigen::SPQR<Eigen::SparseMatrix<double, Eigen::ColMajor>> solver;
        //        solver.compute(useWeights ? WA : A);

        //        if (solver.info() != Eigen::Success) {
        //            assert(0);
        //            std::cout << "computation error" << std::endl;
        //            SHOULD_NEVER_BE_CALLED();
        //        }

        //        X = solver.solve(useWeights ? WB : B);
        //        if (solver.info() != Eigen::Success) {
        //            assert(0);
        //            std::cout << "solving error" << std::endl;
        //            SHOULD_NEVER_BE_CALLED();
        //        }
        //    }

        //    {
        //        Clock clock("install solved variables");
        //        InstallVariables(mg, controls, uh2varStartPosition, X, vars);
        //    }
        //    NormalizeVariables(mg, controls, vars);
        //    return vars;
        //}




        void ResetToFullArmorAnchors(const RLGraph & mg, RLGraphControls & controls) {

            for (auto & ct : controls.componentControls.data) {
                for (auto & c : ct) {
                    c.boundedAnchors.clear();
                }
            }

            for (auto & c : mg.components<RegionData>()) {
                for (auto & ps : c.data.normalizedContours) {
                    for (auto & p : ps) {
                        controls[c.topo.hd].boundedAnchors.push_back(BoundAs(p, 1.0, 10.0));
                    }
                }
            }

            for (auto & c : mg.components<LineData>()) {
                controls[c.topo.hd].boundedAnchors.push_back(BoundAs(c.data.line.first, 1.0, 10.0));
                controls[c.topo.hd].boundedAnchors.push_back(BoundAs(c.data.line.second, 1.0, 10.0));
            }

        }

        void ResetToSampledArmorAnchors(const RLGraph & mg, RLGraphControls & controls, double sampleStepAngle) {
            for (auto & ct : controls.componentControls.data) {
                for (auto & c : ct) {
                    c.boundedAnchors.clear();
                }
            }
            for (auto & c : mg.components<RegionData>()) {
                for (auto & ps : c.data.normalizedContours) {
                    for (auto & p : ps) {
                        auto & ba = controls[c.topo.hd].boundedAnchors;
                        if (!ba.empty()) {
                            if (AngleBetweenDirections(ba.back().component, p) >= sampleStepAngle) {
                                ba.push_back(BoundAs(p, 1.0, 10.0));
                            } else {
                                continue;
                            }
                        } else {
                            ba.push_back(BoundAs(p, 1.0, 10.0));
                        }
                    }
                }
            }
            for (auto & c : mg.components<LineData>()) {
                controls[c.topo.hd].boundedAnchors.push_back(BoundAs(c.data.line.first, 1.0, 10.0));
                controls[c.topo.hd].boundedAnchors.push_back(BoundAs(c.data.line.second, 1.0, 10.0));
            }
        }



        RLGraphVars SolveVariablesWithBoundedAnchors(misc::Matlab & matlab, const RLGraph & mg, const RLGraphControls & controls,
            bool useWeights, int tryNum) {

            assert(tryNum > 0);

            int n1 = NumberOfComponentWeightedAnchors(controls);
            int n2 = NumberOfComponentBoundedAnchors(controls);
            int n3 = NumberOfConstraintWeightedAnchors(controls);
            int n4 = NumberOfConstraintBoundedAnchors(controls);

            assert(NumberOfComponentWeightedAnchors(controls) +
                NumberOfComponentBoundedAnchors(controls) > 0);

            RLGraphVars vars = MakeVariables(mg, controls);

            SetClock();

            auto uh2varStartPosition = MakeHandledTableForAllComponents<int>(mg);
            auto uh2eqStartPosition = MakeHandledTableForAllComponents<int>(mg);
            auto bh2eqStartPosition = MakeHandledTableForAllConstraints<int>(mg);

            int varNum, consNum;

            std::vector<EquationType> equationTypes;
            std::vector<double> Xdata;
            std::vector<double> Bdata;
            std::vector<double> B1data, B2data;
            std::vector<SparseMatElement<double>> Atriplets;
            std::vector<SparseMatElement<double>> Wtriplets;

            std::tie(varNum, consNum) = PrepareEverything(mg, controls, vars, uh2varStartPosition,
                uh2eqStartPosition, bh2eqStartPosition, Xdata, Atriplets, Wtriplets, Bdata, B1data, B2data,
                equationTypes);

            SparseMatd A = MakeSparseMatFromElements(consNum, varNum, Atriplets.begin(), Atriplets.end());
            if (!useWeights) {
                for (int i = 0; i < Wtriplets.size(); i++) {
                    Wtriplets[i].col = Wtriplets[i].row = i;
                    Wtriplets[i].value = 1.0;
                }
            }
            SparseMatd W = MakeSparseMatFromElements(consNum, consNum, Wtriplets.begin(), Wtriplets.end());

            bool succeed = true;
            matlab << "clear;";
            succeed = matlab.setVar("A", A);
            succeed = matlab.setVar("W", W);
            succeed = matlab.setVar("B", cv::Mat(Bdata));
            assert(succeed);
            succeed = matlab.setVar("B1", cv::Mat(B1data));
            assert(succeed);
            succeed = matlab.setVar("B2", cv::Mat(B2data));
            assert(succeed);
            //matlab << "B = B'; B1 = B1'; B2 = B2';";

            DenseMatd equationTypesData(equationTypes.size(), 1, 0.0);
            for (int i = 0; i < equationTypes.size(); i++) {
                equationTypesData(i) = ToUnderlying(equationTypes[i]);
            }
            succeed = matlab.setVar("equationTypes", equationTypesData);

            double scale = 1.0;
            matlab.setVar("scale", scale);

            // now optimize!
            matlab
                << "consNum = size(A, 1);"
                << "varNum = size(A, 2);"

                << "was = equationTypes(:) == 0 | equationTypes(:) == 2;"
                << "bas = equationTypes(:) == 1 | equationTypes(:) == 3;"

                << "A_was = A(was, :);"
                << "W_was = W(was, was);"
                << "B_was = B(was, :);"

                << "A_bas = A(bas, :);"
                << "B1_bas = B1(bas, :);"
                << "B2_bas = B2(bas, :);";

            DenseMatd X;
            for (int t = 0; t < tryNum; t++) {
                matlab
                    << "cvx_begin"
                    << "variable X(varNum);"
                    << "minimize norm(W_was * (A_was * X - B_was))"
                    << "subject to "
                    << "    B1_bas / scale <= A_bas * X <= B2_bas * scale;"
                    << "cvx_end"
                    << "scale = scale * 2.0;";

                X = (DenseMatd)matlab.var("X");
                bool hasInfOrNaN = false;
                for (int i = 0; i < varNum; i++) {
                    if (IsInfOrNaN(X(i))) {
                        hasInfOrNaN = true;
                        break;
                    }
                }
                if (!hasInfOrNaN) {
                    break;
                }
            }

            InstallVariables(mg, controls, uh2varStartPosition, X, vars);
            NormalizeVariables(mg, controls, vars);
            return vars;

        }



        RLGraphVars SolveVariablesWithBoundedAnchors2(misc::Matlab & matlab, const RLGraph & mg, const RLGraphControls & controls,
            bool useWeights, int tryNum) {

            assert(tryNum > 0);

            int n1 = NumberOfComponentWeightedAnchors(controls);
            int n2 = NumberOfComponentBoundedAnchors(controls);
            int n3 = NumberOfConstraintWeightedAnchors(controls);
            int n4 = NumberOfConstraintBoundedAnchors(controls);

            assert(NumberOfComponentWeightedAnchors(controls) +
                NumberOfComponentBoundedAnchors(controls) > 0);

            RLGraphVars vars = MakeVariables(mg, controls);

            SetClock();

            auto uh2varStartPosition = MakeHandledTableForAllComponents<int>(mg);
            auto uh2eqStartPosition = MakeHandledTableForAllComponents<int>(mg);
            auto bh2eqStartPosition = MakeHandledTableForAllConstraints<int>(mg);

            int varNum, consNum;

            std::vector<EquationType> equationTypes;
            std::vector<double> Xdata;
            std::vector<double> Bdata;
            std::vector<double> B1data, B2data;
            std::vector<SparseMatElement<double>> Atriplets;
            std::vector<SparseMatElement<double>> Wtriplets;

            std::tie(varNum, consNum) = PrepareEverything(mg, controls, vars, uh2varStartPosition,
                uh2eqStartPosition, bh2eqStartPosition, Xdata, Atriplets, Wtriplets, Bdata, B1data, B2data,
                equationTypes);

            SparseMatd A = MakeSparseMatFromElements(consNum, varNum, Atriplets.begin(), Atriplets.end());
            if (!useWeights) {
                for (int i = 0; i < Wtriplets.size(); i++) {
                    Wtriplets[i].col = Wtriplets[i].row = i;
                    Wtriplets[i].value = 1.0;
                }
            }
            SparseMatd W = MakeSparseMatFromElements(consNum, consNum, Wtriplets.begin(), Wtriplets.end());

            bool succeed = true;
            matlab << "clear;";
            succeed = matlab.setVar("A", A);
            succeed = matlab.setVar("W", W);
            succeed = matlab.setVar("B", cv::Mat(Bdata));
            assert(succeed);
            succeed = matlab.setVar("B1", cv::Mat(B1data));
            assert(succeed);
            succeed = matlab.setVar("B2", cv::Mat(B2data));
            assert(succeed);
            //matlab << "B = B'; B1 = B1'; B2 = B2';";

            DenseMatd equationTypesData(equationTypes.size(), 1, 0.0);
            for (int i = 0; i < equationTypes.size(); i++) {
                equationTypesData(i) = ToUnderlying(equationTypes[i]);
            }
            succeed = matlab.setVar("equationTypes", equationTypesData);

            double scale = 1.0;
            matlab.setVar("scale", scale);

            // now optimize!
            matlab
                << "consNum = size(A, 1);"
                << "varNum = size(A, 2);"

                << "was = equationTypes(:) == 0 | equationTypes(:) == 2;"
                << "bas = equationTypes(:) == 1 | equationTypes(:) == 3;"

                << "A_was = A(was, :);"
                << "W_was = W(was, was);"
                << "B_was = B(was, :);"

                << "A_bas = A(bas, :);"
                << "B1_bas = B1(bas, :);"
                << "B2_bas = B2(bas, :);";

            DenseMatd X;
            double minE = std::numeric_limits<double>::infinity();

            for (int t = 0; t < tryNum; t++) {
                matlab
                    << "cvx_begin"
                    << "variable X(varNum);"
                    << "minimize norm(W_was * (A_was * X - B_was))"
                    << "subject to "
                    << "    B1_bas / scale <= A_bas * X <= B2_bas * scale;"
                    << "cvx_end"
                    << "scale = scale * 2.0;";

                matlab << "e = norm(W_was * (A_was * X - B_was)) * scale'";
                double e = matlab.var("e").scalar();
                if (IsInfOrNaN(e)) {
                    continue;
                }
                if (e >= minE) {
                    continue;
                }

                minE = e;
                X = (DenseMatd)matlab.var("X");
            }

            InstallVariables(mg, controls, uh2varStartPosition, X, vars);
            NormalizeVariables(mg, controls, vars);
            return vars;

        }


        namespace {

            template <class ComponentDataT, class SparseMatElementT>
            void FillCertainComponentAnchorDepthMatrix(const RLGraph & mg,
                const RLGraphControls & controls,
                const double * Xdata,
                const RLGraphComponentTable<int> & uh2varStartPosition,
                const RLGraphComponentTable<int> & uh2eqStartPosition,
                std::vector<SparseMatElementT> & compDtriplets) {

                for (auto & c : mg.components<ComponentDataT>()) {
                    if (!controls[c.topo.hd].used)
                        continue;
                    auto uh = c.topo.hd;
                    int uhVarStartPosition = uh2varStartPosition.at(uh);
                    int eid = uh2eqStartPosition[uh];
                    // fill weighted anchors
                    for (auto & wa : controls[c.topo.hd].weightedAnchors) {
                        double d = DepthAtDirectionGivenVariables(mg, Xdata + uh2varStartPosition.at(uh), controls, wa.component, uh);
                        compDtriplets.emplace_back(eid, eid, d);
                        eid++;
                    }
                    // fill bounded anchors
                    for (auto & ba : controls[c.topo.hd].boundedAnchors) {
                        double d = DepthAtDirectionGivenVariables(mg, Xdata + uh2varStartPosition.at(uh), controls, ba.component, uh);
                        compDtriplets.emplace_back(eid, eid, d);
                        eid++;
                    }
                }
            }

            template <class SparseMatElementT>
            void FillComponentAnchorDepthMatrix(const RLGraph & mg,
                const RLGraphControls & controls,
                const double * Xdata,
                const RLGraphComponentTable<int> & uh2varStartPosition,
                const RLGraphComponentTable<int> & uh2eqStartPosition,
                std::vector<SparseMatElementT> & compDtriplets) {
                FillCertainComponentAnchorDepthMatrix<RegionData>(mg, controls, Xdata,
                    uh2varStartPosition, uh2eqStartPosition, compDtriplets);
                FillCertainComponentAnchorDepthMatrix<LineData>(mg, controls, Xdata,
                    uh2varStartPosition, uh2eqStartPosition, compDtriplets);
            }


            template <class ConstraintDataT, class SparseMatElementT>
            void FillCertainConstriantAnchorDepthMatrix(const RLGraph & mg,
                const RLGraphControls & controls,
                const double * Xdata,
                const RLGraphComponentTable<int> & uh2varStartPosition,
                const RLGraphConstraintTable<int> & bh2consStartPosition,
                std::vector<SparseMatElementT> & D1triplets,
                std::vector<SparseMatElementT> & D2triplets) {

                for (auto & c : mg.constraints<ConstraintDataT>()) {
                    if (!controls[c.topo.hd].used)
                        continue;
                    if (!(controls[c.topo.component<0>()].used && controls[c.topo.component<1>()].used))
                        continue;

                    auto bh = c.topo.hd;
                    auto uh1 = mg.topo(bh).component<0>();
                    auto uh2 = mg.topo(bh).component<1>();
                    assert(controls[uh1].used || controls[uh2].used);
                    int eid = bh2consStartPosition[bh];

                    auto & was = controls[bh].weightedAnchors;
                    auto & bas = controls[bh].boundedAnchors;

                    for (auto & wa : was) {
                        double d1 = DepthAtDirectionGivenVariables(mg, Xdata + uh2varStartPosition.at(uh1), controls, wa.component, uh1);
                        double d2 = DepthAtDirectionGivenVariables(mg, Xdata + uh2varStartPosition.at(uh2), controls, wa.component, uh2);
                        D1triplets.emplace_back(eid, eid, d1);
                        D2triplets.emplace_back(eid, eid, d2);
                        eid++;
                    }

                    for (auto & ba : bas) {
                        double d1 = DepthAtDirectionGivenVariables(mg, Xdata + uh2varStartPosition.at(uh1), controls, ba.component, uh1);
                        double d2 = DepthAtDirectionGivenVariables(mg, Xdata + uh2varStartPosition.at(uh2), controls, ba.component, uh2);
                        D1triplets.emplace_back(eid, eid, d1);
                        D2triplets.emplace_back(eid, eid, d2);
                        eid++;
                    }
                }
            }

            template <class SparseMatElementT>
            void FillConstriantAnchorDepthMatrix(const RLGraph & mg,
                const RLGraphControls & controls,
                const double * Xdata,
                const RLGraphComponentTable<int> & uh2varStartPosition,
                const RLGraphConstraintTable<int> & bh2consStartPosition,
                std::vector<SparseMatElementT> & D1triplets,
                std::vector<SparseMatElementT> & D2triplets) {
                FillCertainConstriantAnchorDepthMatrix<RegionBoundaryData>(mg, controls, Xdata,
                    uh2varStartPosition, bh2consStartPosition, D1triplets, D2triplets);
                FillCertainConstriantAnchorDepthMatrix<LineRelationData>(mg, controls, Xdata,
                    uh2varStartPosition, bh2consStartPosition, D1triplets, D2triplets);
                FillCertainConstriantAnchorDepthMatrix<RegionLineConnectionData>(mg, controls, Xdata,
                    uh2varStartPosition, bh2consStartPosition, D1triplets, D2triplets);
            }

        }


        void OptimizeVariablesWithBoundedAnchors(misc::Matlab & matlab, const RLGraph & mg, const RLGraphControls & controls, RLGraphVars & vars,
            bool useWeights, int tryNum, int maxOptimizeNum,
            const std::function<bool(const RLGraphVars &)> & callback) {

            SetClock();

            auto uh2varStartPosition = MakeHandledTableForAllComponents<int>(mg);
            auto uh2eqStartPosition = MakeHandledTableForAllComponents<int>(mg);
            auto bh2eqStartPosition = MakeHandledTableForAllConstraints<int>(mg);

            int varNum, consNum;

            std::vector<EquationType> equationTypes;
            std::vector<double> Xdata;
            std::vector<double> Bdata;
            std::vector<double> B1data, B2data;
            std::vector<SparseMatElement<double>> Atriplets;
            std::vector<SparseMatElement<double>> Wtriplets;

            std::tie(varNum, consNum) = PrepareEverything(mg, controls, vars, uh2varStartPosition,
                uh2eqStartPosition, bh2eqStartPosition, Xdata, Atriplets, Wtriplets, Bdata, B1data, B2data,
                equationTypes);

            SparseMatd A = MakeSparseMatFromElements(consNum, varNum, Atriplets.begin(), Atriplets.end());
            if (!useWeights) {
                for (int i = 0; i < Wtriplets.size(); i++) {
                    Wtriplets[i].col = Wtriplets[i].row = i;
                    Wtriplets[i].value = 1.0;
                }
            }
            SparseMatd W = MakeSparseMatFromElements(consNum, consNum, Wtriplets.begin(), Wtriplets.end());

            bool succeed = true;
            matlab << "clear;";
            succeed = matlab.setVar("A", A);
            succeed = matlab.setVar("W", W);
            succeed = matlab.setVar("B", cv::Mat(Bdata));
            succeed = matlab.setVar("B1", cv::Mat(B1data));
            succeed = matlab.setVar("B2", cv::Mat(B2data));
            matlab << "B = B'; B1 = B1'; B2 = B2';";

            DenseMatd equationTypesData(equationTypes.size(), 1, 0.0);
            for (int i = 0; i < equationTypes.size(); i++) {
                equationTypesData(i) = ToUnderlying(equationTypes[i]);
            }
            succeed = matlab.setVar("equationTypes", equationTypesData);

            double scale = 1.0;
            matlab.setVar("scale", scale);

            // now optimize!
            matlab
                << "consNum = size(A, 1);"
                << "varNum = size(A, 2);"

                << "was = equationTypes(:) == 0 | equationTypes(:) == 2;"
                << "bas = equationTypes(:) == 1 | equationTypes(:) == 3;"

                << "A_was = A(was, :);"
                << "W_was = W(was, was);"
                << "B_was = B(was, :);"

                << "A_bas = A(bas, :);"
                << "B1_bas = B1(bas, :);"
                << "B2_bas = B2(bas, :);";

            std::vector<SparseMatElement<double>> compDtriplets, consD1triplets, consD2triplets;
            compDtriplets.reserve(consNum);
            consD1triplets.reserve(consNum);
            consD2triplets.reserve(consNum);

            for (int i = 0; i < maxOptimizeNum; i++) {

                // compute D matrix using Xdata
                compDtriplets.clear();
                consD1triplets.clear();
                consD2triplets.clear();
                FillComponentAnchorDepthMatrix(mg, controls, Xdata.data(),
                    uh2varStartPosition, uh2eqStartPosition, compDtriplets);
                FillConstriantAnchorDepthMatrix(mg, controls, Xdata.data(),
                    uh2varStartPosition, bh2eqStartPosition, consD1triplets, consD2triplets);

                matlab.setVar("compD",
                    MakeSparseMatFromElements(consNum, consNum,
                    compDtriplets.begin(), compDtriplets.end()));
                matlab.setVar("consD1",
                    MakeSparseMatFromElements(consNum, consNum,
                    consD1triplets.begin(), consD1triplets.end()));
                matlab.setVar("consD2",
                    MakeSparseMatFromElements(consNum, consNum,
                    consD2triplets.begin(), consD2triplets.end()));

                for (int ii = 0; ii < tryNum; ii++) {
                    misc::Clock clock("Solve LP");


                    matlab
                        << "cvx_begin"
                        << "variable X(varNum);"
                        << "minimize norm(consD1(was, was) * consD2(was, was) * W_was * (A_was * X - B_was))"
                        << "subject to "
                        << "    B1_bas / scale <= A_bas * X <= B2_bas * scale;"
                        << "cvx_end"
                        << "scale = scale * 2.0;";

                    Xdata.clear();
                    matlab << "X = X';";
                    Xdata = (DenseMatd)matlab.var("X");
                    if (!HasValue(Xdata, IsInfOrNaN<double>)) {
                        break;
                    }
                }

                if (HasValue(Xdata, IsInfOrNaN<double>)) {
                    break;
                }

                double n = 0;
                for (auto & v : Xdata) {
                    n += v * v;
                }
                n = sqrt(n);
                for (auto & v : Xdata) {
                    v /= n;
                }

                if (callback) {
                    InstallVariables(mg, controls, uh2varStartPosition, Xdata, vars);
                    NormalizeVariables(mg, controls, vars);
                    callback(vars);
                }
            }

            InstallVariables(mg, controls, uh2varStartPosition, Xdata, vars);
            NormalizeVariables(mg, controls, vars);
        }






        namespace {


            template <class T>
            struct InstanceTableByHandle_ {};
            template <>
            struct InstanceTableByHandle_<RegionHandle> {
                using type = HandledTable<RegionHandle, Plane3>;
            };
            template <>
            struct InstanceTableByHandle_<LineHandle> {
                using type = HandledTable<LineHandle, Line3>;
            };
            template <class T>
            using InstanceTableByHandle = typename InstanceTableByHandle_<T>::type;

            using RLGraphInstanceTable = MetaBind<InstanceTableByHandle, RegionHandle, LineHandle>;


            double DistanceToInstance(const Point3 & p, const Plane3 & plane) {
                return plane.distanceTo(p);
            }
            double DistanceToInstance(const Point3 & p, const Line3 & line) {
                return DistanceFromPointToLine(p, line).first;
            }


            template <class ComponentDataT>
            void SetComponentAnchorCosts(const Eigen::VectorXd & variables, Eigen::VectorXd & costs,
                const RLGraph & mg, const RLGraphControls & controls,
                const RLGraphComponentTable<int> & compAnchorConsStartPosition,
                const RLGraphInstanceTable & insts,
                bool useWeights) {
                for (auto & c : mg.components<ComponentDataT>()) {
                    if (!controls[c.topo.hd].used)
                        continue;
                    int anchorStartPosition = compAnchorConsStartPosition[c.topo.hd];
                    const auto & inst = insts[c.topo.hd];
                    int i = 0;
                    for (auto & wa : controls[c.topo.hd].weightedAnchors) {
                        const Point3 & anchor = wa.component;
                        double weight = wa.score;
                        double dist = DistanceToInstance(anchor, inst);
                        costs[anchorStartPosition + i] = dist * (useWeights ? weight : 1.0);
                        i++;
                    }
                }
            }

            template <class ConstarintDataT>
            void SetConstraintAnchorCosts(const Eigen::VectorXd & variables, Eigen::VectorXd & costs,
                const RLGraph & mg, const RLGraphControls & controls,
                const RLGraphConstraintTable<int> & consAnchorConsStartPosition,
                const RLGraphConstraintTable<std::vector<Vec3>> & appliedConsAnchors,
                const RLGraphInstanceTable & insts,
                bool useWeights) {
                for (auto & c : mg.constraints<ConstarintDataT>()) {
                    if (!controls[c.topo.hd].used)
                        continue;
                    // constraint anchors
                    int anchorStartPosition = consAnchorConsStartPosition[c.topo.hd];
                    auto & anchors = appliedConsAnchors.at(c.topo.hd);
                    auto uh1 = c.topo.component<0>();
                    auto uh2 = c.topo.component<1>();
                    auto & inst1 = insts[uh1];
                    auto & inst2 = insts[uh2];
                    for (int i = 0; i < anchors.size(); i++) {
                        double depth1 = DepthAt(anchors[i], inst1);
                        double depth2 = DepthAt(anchors[i], inst2);
                        costs(anchorStartPosition + i) = abs(depth1 - depth2) *
                            (useWeights ? (controls[c.topo.hd].weight / sqrt(anchors.size())) : 1.0);
                    }
                }
            }

        }



        //void OptimizeVariablesLevenbergMarquardt(const RLGraph & mg, const RLGraphControls & controls, RLGraphVars & vars,
        //    bool useWeights, bool useAllAnchors, 
        //    const std::function<bool(const RLGraphVars &)> & callback){

        //    SetClock();

        //    auto uh2varStartPosition = MakeHandledTableForAllComponents<int>(mg);

        //    // register vars
        //    int varNum = RegisterVariablePositions(mg, controls, vars, uh2varStartPosition);
        //    std::vector<double> Xdata(varNum, 1.0);
        //    RegisterVariableValues(mg, controls, vars, uh2varStartPosition, Xdata);


        //    auto compAnchorConsStartPosition = MakeHandledTableForAllComponents<int>(mg);

        //    auto appliedConsAnchors = MakeHandledTableForAllConstraints<std::vector<Vec3>>(mg);
        //    auto weightsForEachAppliedConsAnchor = MakeHandledTableForAllConstraints<double>(mg);
        //    auto consAnchorConsStartPosition = MakeHandledTableForAllConstraints<int>(mg);
        //    
        //    //auto regionBoundaryNormalConsStartPosition = MakeHandledTableForConstraints<int, RegionBoundaryData>(mg);


        //    // register cons
        //    int consNum = 0;
        //    std::vector<ConstraintEquationType> consTypes;
        //    // register comp anchors and cons anchors
        //    consNum += useAllAnchors ?
        //        RegisterConstraintPositions(mg, controls, uh2varStartPosition, compAnchorConsStartPosition,
        //        consAnchorConsStartPosition, appliedConsAnchors, weightsForEachAppliedConsAnchor,
        //        consTypes, ExtractAllAnchorsForBinary())
        //        :
        //        RegisterConstraintPositions(mg, controls, uh2varStartPosition, compAnchorConsStartPosition,
        //        consAnchorConsStartPosition, appliedConsAnchors, weightsForEachAppliedConsAnchor,
        //        consTypes, ExtractNecessaryAnchorsForBinary());
        //    //// register region boundary normal consistencies
        //    //for (auto & b : mg.constraints<RegionBoundaryData>()){
        //    //    regionBoundaryNormalConsStartPosition[b.topo.hd] = consNum;
        //    //    consNum++;
        //    //}

        //    auto costFunctor = misc::MakeGenericNumericDiffFunctor<double>(
        //        [&](const Eigen::VectorXd & v, Eigen::VectorXd & costs){

        //        Eigen::VectorXd variables = v;

        //        // compute current planes and lines
        //        RLGraphInstanceTable insts;
        //        insts.container<RegionHandle>() = mg.createComponentTable<RegionData, Plane3>();
        //        insts.container<LineHandle>() = mg.createComponentTable<LineData, Line3>();
        //        for (auto & c : mg.components<RegionData>()){
        //            insts[c.topo.hd] = InstanceGivenVariables(mg,
        //                variables.data() + uh2varStartPosition.at(c.topo.hd), controls, c.topo.hd);
        //        }
        //        for (auto & c : mg.components<LineData>()){
        //            insts[c.topo.hd] = InstanceGivenVariables(mg,
        //                variables.data() + uh2varStartPosition.at(c.topo.hd), controls, c.topo.hd);
        //        }
        //        
        //        SetComponentAnchorCosts<RegionData>(variables, costs, mg, controls, compAnchorConsStartPosition, insts, useWeights);
        //        SetComponentAnchorCosts<LineData>(variables, costs, mg, controls, compAnchorConsStartPosition, insts, useWeights);

        //        SetConstraintAnchorCosts<RegionBoundaryData>(variables, costs, mg, controls, 
        //            consAnchorConsStartPosition, appliedConsAnchors, insts, useWeights);
        //        SetConstraintAnchorCosts<LineRelationData>(variables, costs, mg, controls,
        //            consAnchorConsStartPosition, appliedConsAnchors, insts, useWeights);
        //        SetConstraintAnchorCosts<RegionLineConnectionData>(variables, costs, mg, controls,
        //            consAnchorConsStartPosition, appliedConsAnchors, insts, useWeights);

        //        //for (auto & c : mg.constraints<RegionBoundaryData>()){
        //        //    if (!controls[c.topo.hd].used)
        //        //        continue;
        //        //    // constraint anchors
        //        //    int anchorStartPosition = consAnchorConsStartPosition[c.topo.hd];
        //        //    auto & anchors = appliedConsAnchors.at(c.topo.hd);
        //        //    auto uh1 = c.topo.component<0>();
        //        //    auto uh2 = c.topo.component<1>();
        //        //    const Plane3 & inst1 = insts[uh1];
        //        //    const Plane3 & inst2 = insts[uh2];
        //        //    // region boundary normal consistency
        //        //    costs(regionBoundaryNormalConsStartPosition[c.topo.hd]) = 
        //        //        AngleBetweenUndirectedVectors(inst1.normal, inst2.normal) / 10.0 * c.data.length;
        //        //}
        //    }, varNum, consNum);


        //    Eigen::VectorXd X = Eigen::Map<const Eigen::VectorXd>(Xdata.data(), Xdata.size());
        //    Eigen::LevenbergMarquardt<decltype(costFunctor)> lm(costFunctor);
        //    lm.parameters.maxfev = 5000;
        //    //auto status = lm.minimize(X);
        //    Eigen::LevenbergMarquardtSpace::Status status = lm.minimizeInit(X);
        //    assert(status != Eigen::LevenbergMarquardtSpace::ImproperInputParameters);
        //    do {
        //        status = lm.minimizeOneStep(X);                
        //        std::cout << "iter: " << lm.iter << "\t\t lm.fnorm = " << lm.fnorm << std::endl;

        //        // test
        //        if(0){
        //            Eigen::VectorXd costs = Eigen::VectorXd::Zero(consNum);
        //            costFunctor(X, costs);
        //            std::cout << "norm(costFunctor(X)) = " << costs.norm() << std::endl;
        //            costFunctor(X * 2, costs);
        //            std::cout << "norm(costFunctor(X * 2)) = " << costs.norm() << std::endl;
        //            costFunctor(X * 3, costs);
        //            std::cout << "norm(costFunctor(X * 3)) = " << costs.norm() << std::endl;
        //            costFunctor(X * 4, costs);
        //            std::cout << "norm(costFunctor(X * 4)) = " << costs.norm() << std::endl;
        //        }

        //        X.normalize();
        //        if (callback){
        //            InstallVariables(mg, controls, uh2varStartPosition, X, vars);
        //            if (!callback(vars))
        //                break;
        //        }
        //    } while (status == Eigen::LevenbergMarquardtSpace::Running);


        //    // install X
        //    X.normalize();
        //    InstallVariables(mg, controls, uh2varStartPosition, X, vars);

        //}











    }
}