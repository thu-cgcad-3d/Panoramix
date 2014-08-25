#include <map>

#include "../core/utilities.hpp"

#include "expression.hpp"

namespace panoramix {
    namespace deriv {

        using namespace core;

        bool ExpressionGraph::isForwardConnection(CHandle h) const {
            return _g.data(h);
        }

        void ExpressionGraph::reserve(size_t sz) {
            _g.internalVertices().reserve(sz);
            _g.internalHalfEdges().reserve(sz * 4);
        }

        void ExpressionGraph::invalidateAll() {
            _g.clear();
        }

        //void ExpressionGraph::exchangeWhenUsedAsInputs(EHandle a, EHandle b){
        //    if (a == b)
        //        return;

        //    // for all expressions who accepts a as input, make it accepts b as input with the same position of a
        //    // for all expressions who accepts b as input, make it accepts a as input with the same position of b
        //    std::vector<CHandle> a2os, o2as;
        //    std::vector<CHandle> b2os, o2bs;
        //    for (auto hh : _g.topo(a).halfedges){
        //        if (!_g.removed(hh) && isForwardConnection(hh)){
        //            a2os.push_back(hh);
        //            o2as.push_back(_g.topo(hh).opposite);
        //        }
        //    }
        //    for (auto hh : _g.topo(b).halfedges){
        //        if (!_g.removed(hh) && isForwardConnection(hh)){
        //            b2os.push_back(hh);
        //            o2bs.push_back(_g.topo(hh).opposite);
        //        }
        //    }

        //    for (auto & a2o : a2os){
        //        _g.topo(a2o).from() = b;
        //    }
        //    for (auto & o2a : o2as){
        //        _g.topo(o2a).to() = b;
        //    }
        //    for (auto & b2o : b2os){
        //        _g.topo(b2o).from() = a;
        //    }
        //    for (auto & o2b : o2bs){
        //        _g.topo(o2b).to() = a;
        //    }

        //    _g.topo(a).halfedges.erase(
        //        std::remove_if(_g.topo(a).halfedges.begin(), _g.topo(a).halfedges.end(),
        //        [this](const CHandle & hh){return !_g.removed(hh) && isForwardConnection(hh); }), 
        //        _g.topo(a).halfedges.end());
        //    _g.topo(a).halfedges.insert(_g.topo(a).halfedges.end(), b2os.begin(), b2os.end());

        //    _g.topo(b).halfedges.erase(
        //        std::remove_if(_g.topo(b).halfedges.begin(), _g.topo(b).halfedges.end(),
        //        [this](const CHandle & hh){return !_g.removed(hh) && isForwardConnection(hh); }),
        //        _g.topo(b).halfedges.end());
        //    _g.topo(b).halfedges.insert(_g.topo(b).halfedges.end(), a2os.begin(), a2os.end());
        //}

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
                _g.addEdge(ih, h, true, false, false); // do not merge duplicate edges
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
            /*std::cout << "topological sort in backPropagateGradient with " << related.size() << " related variables started" << std::endl;
            std::vector<EHandle> temp;
            temp.reserve(related.size());            
            core::TopologicalSort(related.begin(), related.end(), std::back_inserter(temp), [this](EHandle h){
                std::vector<EHandle> forwardHandles;
                for (auto & con : forwardConnections(h)){
                    assert(_g.topo(con).from() == h);
                    forwardHandles.push_back(_g.topo(con).to());
                }
                return forwardHandles;
            });
            related = temp;
            std::cout << "topological sort in backPropagateGradient with " << related.size() << " related variables ended" << std::endl;*/

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
                    if (inputDerivs[k].isInvalid())
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