#ifndef PANORAMIX_CORE_OPTIMIZATION_HPP
#define PANORAMIX_CORE_OPTIMIZATION_HPP

#include "meta.hpp"
#include "basic_types.hpp"
#include "utilities.hpp"
#include "cons_graph.hpp"

namespace panoramix {
    namespace core {

        template <class ValueT>
        using CoeffMatrix = ImageOfType<ValueT>;

        struct RegisteredPosition {
            int64_t start;
            size_t size;
        };

        // used to register all things into one vector
        // 1. store all Ts
        // 2. store all T names
        // 3. map from Identity to VarIDs
        template <class T, class ... IdentityTs>
        struct Registration {
            static_assert(sizeof...(IdentityTs) > 0, "No Identity Type Assigned!");
            using IdentityTupleType = std::tuple<IdentityTs...>;

            using Position = RegisteredPosition;

            std::vector<T> values;
            std::vector<std::string> valueNames;
            std::tuple<std::unordered_map<IdentityTs, Position>...> positions;

            template <class IdentityT>
            inline void append(IdentityT id, size_t sz, std::string name){
                enum { _idx = TypeFirstLocationInTuple<IdentityT, IdentityTupleType>::value };
                std::get<_idx>(positions)[id] = { values.size(), sz };
                values.insert(values.end(), sz, T());
                valueNames.insert(valueNames.end(), sz, name);
            }

            template <class IdentityT>
            inline void append(IdentityT id, const std::vector<T> & vs, const std::vector<std::string> & ns) {
                enum { _idx = TypeFirstLocationInTuple<IdentityT, IdentityTupleType>::value };
                std::get<_idx>(positions)[id] = { values.size(), vs.size() };
                values.insert(values.end(), vs.begin(), vs.end());
                valueNames.insert(valueNames.end(), ns.begin(), ns.end());
            }
          

            // startPosition from identity
            template <class IdentityT>
            inline Position position(const IdentityT & id) const {
                enum { _idx = TypeFirstLocationInTuple<IdentityT, IdentityTupleType>::value };
                return std::get<_idx>(positions).at(id); 
            }
            template <class IdentityT>
            inline T operator()(const IdentityT & id, size_t nthValue) const {
                return values[position(id).start + nthValue];
            }
            template <class IdentityT>
            inline T & operator()(const IdentityT & id, size_t nthValue) {
                return values[position(id).start + nthValue];
            }

            // ...
            void reserve(size_t sz) {
                values.reserve(sz);
                valueNames.reserve(sz);
            }
            void clear() {
                values.clear();
                valueNames.clear();
                positions.clear();
            }

            template <class Archiver>
            inline void serialize(Archiver & ar) { ar(values, valueNames, positions); }


            // concat
            Registration & concat(const Registration & v) {
                size_t oldSize = values.size();
                values.insert(values.end(), v.values.begin(), v.values.end());
                valueNames.insert(valueNames.end(), v.valueNames.begin(), v.valueNames.end());
                mergePositionsUsingSequence(v, oldSize, typename SequenceGenerator<sizeof...(IdentityTs)>::type());
                return *this;
            }

        private:
            template <int I>
            inline bool mergePositions(const Registration & v, size_t oldSize) {
                for (auto && idvar : std::get<I>(v.positions)){
                    assert(!core::Contains(std::get<I>(positions), idvar.first));
                    std::get<I>(positions)[idvar.first] = { idvar.second.start + oldSize, idvar.second.size };
                }
                return true;
            }
            template <int ... Is>
            inline void mergePositionsUsingSequence(const Registration & v, size_t oldSize, Sequence<Is...>){
                bool[] = { mergePositions<Is>(v, oldSize)... };
            }


        };





        // vars <-> constraint graph

        template <class ValueT, class ConstraintGraphT>
        struct ComponentsRegistrationFromConstraintGraph {};

        template <class ValueT, class ... ComponentDataTs, class ConstraintConfigTupleT>
        struct ComponentsRegistrationFromConstraintGraph<ValueT, ConstraintGraph<std::tuple<ComponentDataTs...>, ConstraintConfigTupleT>> {
            using type = Registration<ValueT, ComponentHandle<ComponentDataTs>...>;
        };

        namespace {
            template <class ValueT, class VarRegT, class EnvironmentT>
            struct RegisterVarsViaTripletTopoHandle {
                inline RegisterVarsViaTripletTopoHandle(VarRegT & vr, const EnvironmentT & e) : varReg(vr), env(e) {}
                template <class DataT, class TopoT>
                inline void operator()(const Triplet<TopoT, DataT> & t) const {
                    if (!t.exists)
                        return;
                    std::vector<ValueT> values;
                    std::vector<std::string> names;
                    RegisterValuesAndNames(t.data, values, names, env);
                    varReg.append(t.topo.hd, values, names);
                }
                VarRegT & varReg;
                const EnvironmentT & env;
            };

            template <class ValueT, class VarRegT, class EnvironmentT>
            struct UpdateVarsViaTripletTopoHandle {
                inline UpdateVarsViaTripletTopoHandle(const VarRegT & vr, const EnvironmentT & e)
                : varReg(vr), env(e) {}
                template <class DataT, class TopoT>
                inline void operator()(Triplet<TopoT, DataT> & t) const {
                    if (!t.exists)
                        return;                    
                    auto getValue = [&t, this](size_t nthValue) -> ValueT {
                        return varReg(t.topo.hd, nthValue);
                    };
                    UpdateValues(t.data, getValue, env);
                }
                const VarRegT & varReg;
                const EnvironmentT & env;
            };
        }

        // register all components
        template <class ValueT, class ConstraintGraphT, class EnvironmentT>
        inline typename ComponentsRegistrationFromConstraintGraph<ValueT, ConstraintGraphT>::type RegisterComponents(
            const ConstraintGraphT & cg, 
            const EnvironmentT & env) {
            typename ComponentsRegistrationFromConstraintGraph<ValueT, ConstraintGraphT>::type varReg;
            IterateOver(cg.allComponents(), 
                RegisterVarsViaTripletTopoHandle<ValueT, decltype(varReg), EnvironmentT>(varReg, env));
            return varReg;
        }

        // update all components
        template <class ValueT, class ConstraintGraphT, class ... IdentityTs, class EnvironmentT>
        inline void UpdateComponents(ConstraintGraphT & cg, const Registration<ValueT, IdentityTs...> & varReg, 
            const EnvironmentT & env){
            IterateOver(cg.allComponents(),
                UpdateVarsViaTripletTopoHandle<ValueT, Registration<ValueT, IdentityTs...>, EnvironmentT>(varReg, env));
        }





        // constraint representations
        enum ConstraintRelation {
            SmallerThan,
            CloseTo
        };
        template <class ValueT>
        struct SparseMatrixElement {
            std::vector<int> indices; // -1 -> as constant in this dimension (for this component)
            ValueT coefficient;
        };
        template <class ValueT>
        struct Constraint {
            std::vector<RegisteredPosition> varPositions;
            std::vector<SparseMatrixElement<ValueT>> coefficients;
            ConstraintRelation relation;
            ValueT rhs;
            double weight;
            int priority; // satisfy bigger ones first
        };

        namespace {

            template <class ValueT, class ConstraintGraphT>
            struct ConstraintsRegistrationFromConstraintDataTuple {};

            template <class ValueT, class ... ConstraintDataTs>
            struct ConstraintsRegistrationFromConstraintDataTuple<ValueT, std::tuple<ConstraintDataTs...>> {
                using type = Registration<Constraint<ValueT>, ConstraintHandle<ConstraintDataTs>...>;
            };

        }

        template <class ValueT, class ConstraintGraphT>
        struct ConstraintsRegistrationFromConstraintGraph {};

        template <class ValueT, class ComponentDataTupleT, class ConstraintConfigTupleT>
        struct ConstraintsRegistrationFromConstraintGraph<ValueT, ConstraintGraph<ComponentDataTupleT, ConstraintConfigTupleT>> {
            using type = typename ConstraintsRegistrationFromConstraintDataTuple<ValueT, 
                typename ConstraintGraph<ComponentDataTupleT, ConstraintConfigTupleT>::ConstraintDataTupleType>::type;
        };

        namespace {

            template <class VarRegT>
            struct CollectComponentValuePosition {
                inline CollectComponentValuePosition(const VarRegT & vr, std::vector<RegisteredPosition> & poss)
                    : varReg(vr), positions(poss) {}
                template <class DataT>
                inline void operator()(ComponentHandle<DataT> h) const {
                    positions.push_back(varReg.position(h));
                }
                const VarRegT & varReg;
                std::vector<RegisteredPosition> & positions;
            };

            template <class ValueT, class VarRegT, class ConsRegT, class CGT, class EnvironmentT>
            struct RegisterConstraintsViaTripletTopoHandle {
                inline RegisterConstraintsViaTripletTopoHandle(const VarRegT & vr, ConsRegT & r, const CGT & c, const EnvironmentT & e)
                    : varReg(vr), consReg(r), cg(c), env(e) {}
                template <class DataT, class TopoT>
                inline void operator()(const Triplet<TopoT, DataT> & t) const {
                    if (!t.exists)
                        return;

                    auto & components = t.topo.allComponents;
                    std::vector<RegisteredPosition> compValuePositions;
                    IterateOver(components,
                        CollectComponentValuePosition<VarRegT>(varReg, compValuePositions));

                    std::vector<Constraint<ValueT>> ces;
                    std::vector<std::string> names;

                    auto addConstraint =
                        [&ces, &names, &compValuePositions](
                        const std::vector<SparseMatrixElement<ValueT>> & mat,
                        ConstraintRelation rel, ValueT rhs, double weight, int priority, std::string n){
                        if (!mat.empty()) {
                            for (auto & e : mat){
                                assert(e.indices.size() == compValuePositions.size());
                                for (int i = 0; i < compValuePositions.size(); i++){
                                    assert(e.indices[i] < compValuePositions[i].size);
                                }
                            }
                        }
                        ces.push_back(Constraint<ValueT>{compValuePositions, mat, rel, rhs, weight, priority});
                        names.push_back(n);
                    };

                    RegisterConstraint(t.data, t.topo, cg, addConstraint, env);
                    consReg.append(t.topo.hd, ces, names);
                }
                const VarRegT & varReg;
                ConsRegT & consReg;
                const CGT & cg;
                const EnvironmentT & env;
            };

        }


        // register all constraints
        template <class ValueT, class ConstraintGraphT, class ... IdentityTs, class EnvironmentT>
        inline typename ConstraintsRegistrationFromConstraintGraph<ValueT, ConstraintGraphT>::type RegisterConstraints(
            const ConstraintGraphT & cg, 
            const Registration<ValueT, IdentityTs...> & varReg, 
            const EnvironmentT & env) {
            typename ConstraintsRegistrationFromConstraintGraph<ValueT, ConstraintGraphT>::type consReg;
            IterateOver(cg.allConstraints(),
                RegisterConstraintsViaTripletTopoHandle<
                ValueT, 
                Registration<ValueT, IdentityTs...>, decltype(consReg), 
                ConstraintGraphT, EnvironmentT>(varReg, consReg, cg, env));
            return consReg;
        }


        void Optimize(const std::vector<Constraint<double>> & constraints, std::vector<double> & values);

    }
}


namespace std {

    //template <class ValueT>
    //struct hash<panoramix::core::Var<ValueT>> {
    //    inline size_t operator()(panoramix::core::Var<ValueT> var) const {
    //        return var.id;
    //    }
    //};

}
 
#endif