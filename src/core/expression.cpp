#include "expression.hpp"

namespace panoramix {
    namespace core {

        bool ExpressionGraph::isForwardConnection(CHandle h) const {
            auto & topo = _g.topo(h);
            return topo.from().id < topo.to().id;
        }

        std::vector<EHandle> ExpressionGraph::inputs(EHandle self) const {
            std::vector<EHandle> inputs; // input ops
            inputs.reserve(_g.topo(self).halfedges.size());
            for (auto hh : _g.topo(self).halfedges){
                if (!_g.removed(hh) && isBackwardConnection(hh)){
                    inputs.push_back(_g.topo(hh).to());
                }
            }
            return inputs;
        }

        EHandle ExpressionGraph::addNode(std::shared_ptr<Op> op, const std::vector<EHandle>& inputs) {
            auto h = _g.addVertex(op);
            for (auto & ih : inputs){
                _g.addEdge(ih, h);
            }
            return h;
        }

        void ExpressionGraph::evaluate(EHandle result, const std::set<EHandle, EHandleComp> & vars) const {
            // collect related nodes
            std::vector<EHandle> related;
            related.reserve(result.id);
            related.push_back(result);
            int idx = 0;
            while (idx < related.size()){
                auto & curToCheck = related[idx];
                if (vars.find(curToCheck) == vars.end()){
                    // add all inputs
                    for (auto hh : _g.topo(curToCheck).halfedges){
                        if (!_g.removed(hh) && isBackwardConnection(hh)){
                            related.push_back(_g.topo(hh).to());
                        }
                    }
                }
                idx++;
            }

            // the id of expression handles record the topological order
            // sort from old to new
            std::sort(related.begin(), related.end(), EHandleComp());

            // from old to new
            for (const auto & self : related){
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
                // Op does not know the structure at all!!!
                // the structure does not know the values and operations at all!!!
                // all data types, data shapes and data operations are managed in Op, the graph does not know the details!
            }
        }


        std::vector<EHandle> ExpressionGraph::createDerivatives(EHandle cost, const std::vector<EHandle> & vars) {
            //assert(_g.data(cost).dshape.rank() == 0); // cost MUST be a scalar!

            // collect related nodes
            std::vector<EHandle> related;
            related.reserve(cost.id);
            related.push_back(cost);
            int idx = 0;
            while (idx < related.size()){
                auto & curToCheck = related[idx];
                if (std::find(vars.begin(), vars.end(), curToCheck) == vars.end()){
                    // add all inputs
                    for (auto hh : _g.topo(curToCheck).halfedges){
                        if (!_g.removed(hh) && isBackwardConnection(hh)){
                            related.push_back(_g.topo(hh).to());
                        }
                    }
                }
                idx++;
            }

            // the id of expression handles record the topological order
            // sort from new to old
            std::sort(related.begin(), related.end(), [](EHandle a, EHandle b){return a.id > b.id; });

            // idToDeriv[exprHandle.id] stores the derivative handle of exprHandle
            std::map<EHandle, std::vector<EHandle>, EHandleComp> idToDerivs;

            // an expression of cost, make a data same shape with cost and filled with ones
            idToDerivs[cost].push_back(addNode(_g.data(cost)->makeOnes(), {cost}));

            // stores the sum of the derivative handles
            std::map<EHandle, EHandle, EHandleComp> idToDerivsSumTable;

            for (const auto & self : related){
                if (_g.removed(self))
                    continue;

                // this expression is not used by cost at all
                if (idToDerivs[self].empty())
                    continue;

                // sum all derivatives of self outputs
                auto idToDerivsSum = idToDerivs[self].size() == 1 ? idToDerivs[self].front() :
                    addNode(op(idToDerivs[self].front()).makePlus(), idToDerivs[self]);

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
                    // ignore disconnected derivatives represented by invalid handles
                    // meaning that values in the corresponding inputs[k] does not affect at all the value of self
                    if (inputDerivs[k].isInValid())
                        continue;
                    idToDerivs[inputs[k]].push_back(inputDerivs[k]);
                }

                idToDerivsSumTable[self] = idToDerivsSum;
            }

            std::vector<EHandle> targetedDerivSums;
            for (auto v : vars){
                targetedDerivSums.push_back(idToDerivsSumTable[v]);
            }
            return targetedDerivSums;
            
        }


        std::ostream & ExpressionGraph::toString(std::ostream & ss, EHandle expr) const {
            op(expr).toString(ss) << "(";
            
            // get all inputs
            std::vector<EHandle> inputs;
            for (auto conh : _g.topo(expr).halfedges){
                if (!_g.removed(conh) && isBackwardConnection(conh)){
                    auto input = _g.topo(conh).to();
                    inputs.push_back(input);
                }
            }

            if (inputs.empty()){
                ss << ")";
                return ss;
            }

            for (int i = 0; i < inputs.size() - 1; i++){
                toString(ss, inputs[i]) << ",";
            }
            toString(ss, inputs.back());

            ss << ")";
            return ss;
        }

    }

}