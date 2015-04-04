#ifndef PANORAMIX_ML_FACTOR_GRAPH_HPP
#define PANORAMIX_ML_FACTOR_GRAPH_HPP

#include <type_traits>

#include "../core/iterators.hpp"
#include "../core/basic_types.hpp"
#include "../core/homo_graph.hpp"
 
namespace panoramix {
    namespace ml {

        class FactorGraph {
        public:
            using CostFunction = std::function<double(const int * varlabels, size_t nvar)>;
            using CallbackFunction = std::function<bool(int epoch, double energy, double denergy)>;

            struct FactorCategory {
                CostFunction costs;
                double c_alpha;
            };
            struct VarCategory {
                size_t nlabels;
                double c_i;
            };
            using FactorCategoryId = int;
            using VarCategoryId = int;

            using Topology = core::HomogeneousGraph0x<VarCategoryId, FactorCategoryId>;
            using VarHandle = core::HandleOfTypeAtLevel<Topology, 0>;
            using FactorHandle = core::HandleOfTypeAtLevel<Topology, 1>;
            using ResultTable = core::HandledTable<VarHandle, int>;

        public:
            void reserveVarCategories(size_t cap) { _varCategories.reserve(cap); }
            void reserveFactorCategories(size_t cap) { _factorCategories.reserve(cap); }

            VarCategoryId addVarCategory(VarCategory && vc) { 
                _varCategories.push_back(std::move(vc)); return _varCategories.size() - 1;
            }
            VarCategoryId addVarCategory(const VarCategory & vc) { 
                _varCategories.push_back(vc); return _varCategories.size() - 1; 
            }

            FactorCategoryId addFactorCategory(FactorCategory && fc) {
                _factorCategories.push_back(std::move(fc)); return _factorCategories.size() - 1; 
            }
            FactorCategoryId addFactorCategory(const FactorCategory & fc) {
                _factorCategories.push_back(fc); return _factorCategories.size() - 1; 
            }

            void reserveVars(size_t cap) { _graph.internalElements<0>().reserve(cap); }
            void reserveFactors(size_t cap) { _graph.internalElements<1>().reserve(cap); }

            VarHandle addVar(VarCategoryId vc) { return _graph.add(vc); }
            FactorHandle addFactor(std::initializer_list<VarHandle> vhs, FactorCategoryId fc) { 
                return _graph.add<1>(vhs, fc); 
            }
            template <class IteratorT>
            FactorHandle addFactor(IteratorT vhsBegin, IteratorT vhsEnd, FactorCategoryId fc) { 
                return _graph.add<1>(vhsBegin, vhsEnd, fc); 
            }

            bool valid() const;
            double energy(const ResultTable & labels) const;

            // convex belief propagation
            ResultTable solve(int maxEpoch = 100, int innerLoopNum = 10, 
                const CallbackFunction & callback = nullptr) const;

        private:
            std::vector<VarCategory> _varCategories;
            std::vector<FactorCategory> _factorCategories;
            Topology _graph;
        };


    }
}
 
#endif