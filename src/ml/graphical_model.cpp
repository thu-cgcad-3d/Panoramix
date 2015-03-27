#include "graphical_model.hpp"
 
namespace panoramix {
    namespace ml {

        bool IsValid(const FactorGraph<double> & fg){
            for (auto & f : fg.elements<1>()){
                auto & mat = f.data;
                auto & nhs = f.topo.lowers;
                for (int i = 0; i < mat.dims; i++){
                    int sz = mat.size[i];
                    if (fg.data(nhs[i]).labelNum != sz)
                        return false;
                }
            }
            return true;
        }

        void LoopyBeliefPropagation(FactorGraph<double> & fg){
            
        }


    }
}
 