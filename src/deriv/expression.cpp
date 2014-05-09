#include <map>

#include "expression.hpp"

namespace panoramix {
    namespace deriv {

        bool ExpressionGraph::isForwardConnection(CHandle h) const {
            auto & topo = _g.topo(h);
            return topo.from().id < topo.to().id;
        }

        void ExpressionGraph::reserve(size_t sz) {
            _g.internalVertices().reserve(sz);
            _g.internalHalfEdges().reserve(sz * 4);
        }

        void ExpressionGraph::invalidateAll() {
            _g.clear();
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
            _g.data(h)->graph = this;
            _g.data(h)->self = h;
            for (auto & ih : inputs){
                _g.addEdge(ih, h, Dummy(), Dummy(), false); // do not merge duplicate edges
            }
            return h;
        }

        void ExpressionGraph::forwardPropagateExecute(EHandle result, const std::set<EHandle, EHandleComp> & vars) const {
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
                _g.data(self)->forwardPropagateExecute();
            }
        }


        std::vector<EHandle> ExpressionGraph::backPropagateGradient(EHandle cost, const std::vector<EHandle> & vars) {

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
                        if (!_g.removed(hh) && isBackwardConnection(hh)) {
                            related.push_back(_g.topo(hh).to());
                        }
                    }
                }
                idx++;
            }

            // the id of expression handles record the topological order
            // sort from new to old
            std::sort(related.begin(), related.end(), [](EHandle a, EHandle b){return a.id > b.id; });
            auto relatedUniqueEnd = std::unique(related.begin(), related.end());

            // idToDeriv[exprHandle.id] stores the derivative handle of exprHandle
            std::map<EHandle, std::vector<EHandle>, EHandleComp> idToDerivs;

            // make a scalar one representing the derivative of the cost on itsself
            idToDerivs[cost].push_back(addNode(_g.data(cost)->one()));

            // stores the sum of the derivative handles
            std::map<EHandle, EHandle, EHandleComp> idToDerivsSumTable;

            for (auto i = related.begin(); i != relatedUniqueEnd; ++i){
                auto & self = *i;
                if (_g.removed(self))
                    continue;

                // this expression is not used by cost at all
                if (idToDerivs[self].empty())
                    continue;

                // sum all derivatives of self outputs
                auto idToDerivsSum = op(self).homomorphicSum(idToDerivs[self]);
                //toString(std::cout << "sum:", idToDerivsSum) << std::endl;

                // compute input derivatives
                std::vector<EHandle> inputDerivs = op(self).backPropagateGradient(idToDerivsSum);

                auto inp = inputs(self);
                assert(inputDerivs.size() == inp.size());

                // store the input derivatives
                for (int k = 0; k < inp.size(); k++){
                    // ignore disconnected derivatives represented by invalid handles
                    // meaning that values in the corresponding inputs[k] does not affect at all the value of self
                    if (inputDerivs[k].isInValid())
                        continue;
                    idToDerivs[inp[k]].push_back(inputDerivs[k]);
                }

                idToDerivsSumTable[self] = idToDerivsSum;

                //for (auto p : idToDerivs){
                //    std::cout << "[" << p.first.id << "]:" << std::endl;
                //    for (auto h : p.second){
                //        toString(std::cout << "\t", h) << std::endl;
                //    }
                //    std::cout << std::endl;
                //}
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
            std::vector<EHandle> inp = inputs(expr);

            if (inp.empty()){
                ss << ")";
                return ss;
            }

            for (int i = 0; i < inp.size() - 1; i++){
                toString(ss, inp[i]) << ",";
            }
            toString(ss, inp.back());

            ss << ")";
            return ss;
        }

    }

}