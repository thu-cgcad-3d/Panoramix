#include "expression.hpp"

namespace panoramix {
    namespace core {

        bool ExpressionGraph::isForwardConnection(CHandle h) const {
            auto & topo = _g.topo(h);
            return topo.from().id < topo.to().id;
        }

        EHandle ExpressionGraph::addNode(std::shared_ptr<Op> op, const std::vector<EHandle>& inputs) {
            auto h = _g.addVertex(op);
            for (auto & ih : inputs){
                _g.addEdge(ih, h);
            }
            return h;
        }

        void ExpressionGraph::evaluate(EHandle h) const {
            // TODO: only traverse related handles
            // from old to new
            for (int i = 0; i <= h.id; i++){
                EHandle self(i);
                if (_g.removed(self))
                    continue;
                std::vector<EHandle> inputs; // input ops
                inputs.reserve(_g.topo(self).halfedges.size());
                for (auto hh : _g.topo(self).halfedges){
                    if (!_g.removed(hh) && isBackwardConnection(hh)){
                        inputs.push_back(_g.topo(hh).to());
                    }
                }
                _g.data(self)->eval(this, std::move(inputs)); 
                // IMPORTANT:Op does not know the structure at all!!!
                // and the structure does not know the values and operations at all!!!
                // all data types, data shapes and data operations are managed in Op, the graph does not know the details!
            }
        }


        std::vector<EHandle> ExpressionGraph::createDerivatives(EHandle cost, const std::vector<EHandle> & vars) {
            //assert(_g.data(cost).dshape.rank() == 0); // cost MUST be a scalar!

            // TODO: only traverse related handles

            // the id of expression handles record the topological order
            // idToDeriv[exprHandle.id] stores the derivative handle of exprHandle
            
            std::vector<std::vector<EHandle>> idToDerivs(cost.id + 1);
            idToDerivs[cost.id].push_back(addNode(_g.data(cost)->makeOne()));

            std::vector<EHandle> idToDerivsSumTable(cost.id + 1);
            for (int i = cost.id; i >= 0; i--){
                EHandle self(i);
                if (_g.removed(self))
                    continue;

                // this expression is not used by cost at all
                if (idToDerivs[self.id].empty())
                    continue;

                // sum all derivatives of self outputs
                auto idToDerivsSum = idToDerivs[self.id].size() == 1 ? idToDerivs[self.id].front() :
                    addPlus(idToDerivs[self.id]);

                // get all inputs
                std::vector<EHandle> inputs;
                for (auto conh : _g.topo(self).halfedges){
                    if (!_g.removed(conh) && isBackwardConnection(conh)){
                        auto input = _g.topo(conh).to();
                        inputs.push_back(input);
                    }
                }

                // compute input derivatives
                std::vector<EHandle> inputDerivs =
                    _g.data(self)->makeInputDerivatives(this, self,
                    std::move(inputs), idToDerivsSum);

                assert(inputDerivs.size() == inputs.size());

                // store the input derivatives
                for (int k = 0; k < inputs.size(); k++){
                    idToDerivs[inputs[k].id].push_back(inputDerivs[k]);
                }

                idToDerivsSumTable[self.id] = idToDerivsSum;
            }

            std::vector<EHandle> targetedDerivSums;
            for (auto v : vars){
                targetedDerivSums.push_back(idToDerivsSumTable[v.id]);
            }
            return targetedDerivSums;
            
        }

    }

}