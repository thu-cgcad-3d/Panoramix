#ifndef PANORAMIX_ML_GRPHICAL_MODEL_HPP
#define PANORAMIX_ML_GRPHICAL_MODEL_HPP

#include "../core/basic_types.hpp"
#include "../core/homo_graph.hpp"
 
namespace panoramix {
    namespace ml {

        struct NodeLabel {
            int labelNum;
            int finalLabel;
        };

        template <class T>
        using FactorGraph = core::HomogeneousGraph0x<NodeLabel, core::DenseMat<T>>;
        

        bool IsValid(const FactorGraph<double> & fg);
        void LoopyBeliefPropagation(FactorGraph<double> & fg);

    }
}
 
#endif