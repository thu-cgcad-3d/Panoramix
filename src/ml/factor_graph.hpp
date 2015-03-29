#ifndef PANORAMIX_ML_FACTOR_GRAPH_HPP
#define PANORAMIX_ML_FACTOR_GRAPH_HPP

#include "../core/basic_types.hpp"
#include "../core/homo_graph.hpp"
 
namespace panoramix {
    namespace ml {

        using CostFunction = std::function<double(const int*)>;
        using CallbackFunction = std::function<bool(int epoch, double energy, double denergy)>;
        struct FactorGraph {
            struct FactorData {
                CostFunction costs;
                double c_alpha;
            };
            struct VarData {
                size_t nlabels;
                double c_i;
            };

            using Topology = core::HomogeneousGraph0x<int, int>;
            using VarHandle = core::HandleOfTypeAtLevel<Topology, 0>;
            using FactorHandle = core::HandleOfTypeAtLevel<Topology, 1>;
            using ResultTable = core::HandledTable<VarHandle, int>;

            std::vector<VarData> vars;
            std::vector<FactorData> factors;
            Topology graph;

            bool valid() const;
            double energy(const ResultTable & labels) const;
            ResultTable solve(int maxEpoch = 100, int innerLoopNum = 10, 
                const CallbackFunction & callback = nullptr) const;
        };


    }
}
 
#endif