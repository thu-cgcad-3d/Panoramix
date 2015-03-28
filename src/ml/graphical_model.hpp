#ifndef PANORAMIX_ML_GRPHICAL_MODEL_HPP
#define PANORAMIX_ML_GRPHICAL_MODEL_HPP

#include "../core/basic_types.hpp"
#include "../core/homo_graph.hpp"
 
namespace panoramix {
    namespace ml {

        struct FactorGraph {
            struct FactorData {
                core::DenseMatd costs;
                double c_alpha;
            };
            struct VarData {
                size_t nlabels;
                double c_i;
            };

            using Graph = core::HomogeneousGraph0x<int, int>;
            using VarHandle = core::HandleOfTypeAtLevel<Graph, 0>;
            using FactorHandle = core::HandleOfTypeAtLevel<Graph, 1>;
            using ResultTable = core::HandledTable<VarHandle, int>;

            std::vector<VarData> vars;
            std::vector<FactorData> factors;
            Graph graph;

            bool isValid() const;
            ResultTable solve() const;
        };




    }
}
 
#endif