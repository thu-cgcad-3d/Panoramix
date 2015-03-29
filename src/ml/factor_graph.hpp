#ifndef PANORAMIX_ML_FACTOR_GRAPH_HPP
#define PANORAMIX_ML_FACTOR_GRAPH_HPP

#include "../core/basic_types.hpp"
#include "../core/homo_graph.hpp"
 
namespace panoramix {
    namespace ml {

        using CostFunction = std::function<double(const int*)>;
        using CallbackFunction = std::function<bool(int epoch, double energy, double denergy)>;
        struct FactorGraph {
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

            std::vector<VarCategory> varCategories;
            std::vector<FactorCategory> factorCategories;
            Topology graph;

            inline VarCategoryId addVarCategory(VarCategory && vc) { varCategories.push_back(std::move(vc)); return varCategories.size() - 1; }
            inline VarCategoryId addVarCategory(const VarCategory & vc) { varCategories.push_back(vc); return varCategories.size() - 1; }

            inline FactorCategoryId addFactorCategory(FactorCategory && fc) { factorCategories.push_back(std::move(fc)); return factorCategories.size() - 1; }
            inline FactorCategoryId addFactorCategory(const FactorCategory & fc) { factorCategories.push_back(fc); return factorCategories.size() - 1; }

            inline VarHandle addVar(VarCategoryId vc) { return graph.add(vc); }
            inline FactorHandle addFactor(std::initializer_list<VarHandle> vhs, FactorCategoryId fc) { return graph.add<1>(vhs, fc); }

            bool valid() const;
            double energy(const ResultTable & labels) const;
            ResultTable solve(int maxEpoch = 100, int innerLoopNum = 10, 
                const CallbackFunction & callback = nullptr) const;
        };


    }
}
 
#endif