#include <Eigen/Dense>
#include <Eigen/StdVector>
#include "graphical_model.hpp"
 
namespace panoramix {
    namespace ml {

        using namespace core;

        bool FactorGraph::isValid() const {
            for (auto & f : graph.elements<1>()){
                auto & matId = f.data;
                auto & nhs = f.topo.lowers;
                if (matId < 0 || matId >= costs.size())
                    return false;
                auto & mat = costs[matId];
                for (int i = 0; i < mat.dims; i++){
                    int sz = mat.size[i];
                    if (graph.data(nhs[i]) != sz)
                        return false;
                }
            }
            return true;
        }

        FactorGraph::ResultTable FactorGraph::solve() const {
            assert(isValid());
            using namespace Eigen;
            HandledTable<VarHandle, VectorXd> vmsgs(graph.internalElements<0>().size());
            HandledTable<FactorHandle, VectorXd> fmsgs(graph.internalElements<1>().size());
            NOT_IMPLEMENTED_YET();
            while (true){

            }
        }


    }
}
 