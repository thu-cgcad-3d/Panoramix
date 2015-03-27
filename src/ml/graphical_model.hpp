#ifndef PANORAMIX_ML_GRPHICAL_MODEL_HPP
#define PANORAMIX_ML_GRPHICAL_MODEL_HPP

#include "../core/basic_types.hpp"
#include "../core/homo_graph.hpp"
 
namespace panoramix {
    namespace ml {

        struct FactorGraph {
            using labelnum_t = size_t;
            using costindex_t = int;

            using Graph = core::HomogeneousGraph0x<labelnum_t, costindex_t>;
            using VarHandle = core::HandleOfTypeAtLevel<FactorGraph, 0>;
            using FactorHandle = core::HandleOfTypeAtLevel<FactorGraph, 1>;
            using ResultTable = core::HandledTable<VarHandle, int>;

            Graph graph;
            std::vector<core::DenseMatd> costs;

            bool isValid() const;
            ResultTable solve() const;
        };

    }
}
 
#endif