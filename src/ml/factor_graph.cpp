#include <Eigen/Dense>
#include <Eigen/StdVector>

#include "../core/cons_graph.hpp"
#include "factor_graph.hpp"
 
namespace panoramix {
    namespace ml {

        using namespace core;

        bool FactorGraph::isValid() const {
            for (auto & v : graph.elements<0>()){
                auto & vid = v.data;
                if (vid < 0 || vid >= vars.size())
                    return false;
            }
            for (auto & f : graph.elements<1>()){
                auto & fid = f.data;
                auto & nhs = f.topo.lowers;
                if (fid < 0 || fid >= factors.size())
                    return false;
                const DenseMatd & mat = factors[fid].costs;
                for (int i = 0; i < mat.dims; i++){
                    size_t sz = mat.size[i];
                    if (vars[graph.data(nhs[i])].nlabels != sz)
                        return false;
                }
            }
            return true;
        }

        FactorGraph::ResultTable FactorGraph::solve() const {

            // convex belief propagation

            assert(isValid());
            using namespace Eigen;
            
            struct F2VMessage {
                VectorXd values;
            };
            struct V2FMessage {
                VectorXd values;
            }; 

            ConstraintGraph<
                std::tuple<VarHandle, FactorHandle>,
                std::tuple<
                ConstraintConfig<V2FMessage, VarHandle, FactorHandle>,
                ConstraintConfig<F2VMessage, FactorHandle, VarHandle>
                >
            > messages;
            {
                messages.internalComponents<VarHandle>().reserve(graph.internalElements<0>().size());
                messages.internalComponents<FactorHandle>().reserve(graph.internalElements<1>().size());
                size_t vfConnectionNum = 0;
                for (auto & f : graph.elements<1>()){
                    vfConnectionNum += f.topo.lowers.size();
                }
                messages.internalConstraints<V2FMessage>().reserve(vfConnectionNum);
                messages.internalConstraints<F2VMessage>().reserve(vfConnectionNum);
            }

            using MGVHandle = ComponentHandle<VarHandle>;
            using MGFHandle = ComponentHandle<FactorHandle>;
            using MGV2FHandle = ConstraintHandle<V2FMessage>;
            using MGF2VHandle = ConstraintHandle<F2VMessage>;

            // install graph to message graph and initialize messages
            HandledTable<VarHandle, MGVHandle> vhToMGVH(graph.internalElements<0>().size());
            for (auto & v : graph.elements<0>()){
                auto vh = messages.addComponent(v.topo.hd);
                vhToMGVH[v.topo.hd] = vh;
            }
            for (auto & f : graph.elements<1>()){
                MGFHandle fh = messages.addComponent(f.topo.hd);
                for (auto & lh : f.topo.lowers){
                    MGVHandle vh = vhToMGVH[lh];
                    messages.addConstraint(V2FMessage{ VectorXd::Zero(vars[graph.data(lh)].nlabels) }, vh, fh);
                    messages.addConstraint(F2VMessage{ VectorXd::Zero(vars[graph.data(lh)].nlabels) }, fh, vh);
                }
            }

            HandledTable<VarHandle, double> c_i_hats(graph.internalElements<0>().size(), 0.0);
            for (auto & v : graph.elements<0>()){
                double & c_i_hat = c_i_hats[v.topo.hd];
                c_i_hat = vars[v.data].c_i;
                for (auto fh : v.topo.uppers){
                    c_i_hat += factors[graph.data(fh)].c_alpha;
                }
            }

            core::HandledTable<VarHandle, VectorXd> varMarginals(graph.internalElements<0>().size());
            for (auto & v : graph.elements<0>()){
                varMarginals[v.topo.hd] = VectorXd::Zero(vars[v.data].nlabels);
            }

            static const int innerLoopNum = 3;
            static const int nepoch = 100;

            for (int epoch = 0; epoch < nepoch; epoch ++){
                for (int l = 0; l < innerLoopNum; l++){
                    // var -> factor
                    for (auto & v2f : messages.constraints<V2FMessage>()){
                        MGVHandle vh = v2f.topo.component<0>();
                        MGFHandle fh = v2f.topo.component<1>();

                        auto & fdata = factors[graph.data(messages.data(fh))];

                        v2f.data.values.setZero();
                        MGF2VHandle oppose;
                        for (MGF2VHandle f2v : messages.topo(vh).constraints<F2VMessage>()){
                            v2f.data.values += messages.data(f2v).values;
                            if (messages.topo(f2v).component<0>() == fh){
                                oppose = f2v;
                            }
                        }
                        assert(oppose.valid());
                        v2f.data.values = v2f.data.values * fdata.c_alpha / c_i_hats[messages.data(vh)]
                            - messages.data(oppose).values;
                    }

                    // factor -> var
                    for (auto & f : messages.components<FactorHandle>()){
                        MGFHandle mgfh = f.topo.hd;
                        FactorHandle fh = f.data;
                        auto & vhs = graph.topo(fh).lowers;

                        auto & v2fmsghsSet = messages.topo(mgfh).constraints<V2FMessage>();
                        std::vector<MGV2FHandle> v2fmsghs(v2fmsghsSet.size());
                        auto & f2vmsghsSet = messages.topo(mgfh).constraints<F2VMessage>();
                        std::vector<MGF2VHandle> f2vmsghs(f2vmsghsSet.size());
                        {
                            for (int i = 0; i < vhs.size(); i++){
                                for (auto v2fmsgh : v2fmsghsSet){
                                    if (messages.topo(v2fmsgh).component<0>() == vhToMGVH[vhs[i]]){
                                        v2fmsghs[i] = v2fmsgh;
                                        break;
                                    }
                                }
                            }
                            for (int i = 0; i < vhs.size(); i++){
                                for (auto f2vmsgh : f2vmsghs){
                                    if (messages.topo(f2vmsgh).component<1>() == vhToMGVH[vhs[i]]){
                                        f2vmsghs[i] = f2vmsgh;
                                        break;
                                    }
                                }
                            }
                        }

                        // reset f2v messages
                        for (auto & v2f : v2fmsghs){
                            messages.data(v2f).values.setConstant(std::numeric_limits<double>::max());
                        }

                        // dispatch messages from this factor
                        int fid = graph.data(fh);
                        const DenseMatd & mat = factors[fid].costs;
                        assert(mat.dims > 0);
                        std::vector<int> idx(mat.dims, 0);
                        assert(vhs.size() == mat.dims);
                        size_t matTotal = mat.total();
                        while (true){
                            for (int i = 0; i < mat.dims; i++){
                                // i: output var
                                // others: input vars
                                double & outValue = messages.data(f2vmsghs[i]).values[idx[i]];
                                // theta
                                double theta = mat(idx.data());
                                // sum of other input vars
                                double sumOfOtherVars = 0.0;
                                for (int j = 0; j < mat.dims; j++){
                                    if (j == i) continue;
                                    sumOfOtherVars += messages.data(v2fmsghs[j]).values[idx[j]];
                                }

                                double score = theta + sumOfOtherVars;
                                if (score < outValue){
                                    outValue = score;
                                }
                            }

                            // add idx
                            idx[0]++;
                            int k = 0;
                            while (k < idx.size() && idx[k] >= mat.size[k]){
                                idx[k] = 0;
                                if (k + 1 < idx.size()) idx[k + 1] ++;
                                k++;
                            }
                            if (k == idx.size()){
                                break;
                            }
                        }
                    }
                }


                // marginalize on variables
                // reset marginals
                for (auto & v : graph.elements<0>()){
                    varMarginals[v.topo.hd].setZero();
                }
                for (auto & f2v : messages.constraints<F2VMessage>()){
                    MGVHandle vh = f2v.topo.component<1>();
                    varMarginals[messages.data(vh)] += f2v.data.values;                    
                }

                // when to stop ?
                // TODO
            }

            ResultTable results(graph.internalElements<0>().size());
            for (auto & v : graph.elements<0>()){
                int label = -1;
                double curCost = std::numeric_limits<double>::max();
                const VectorXd & marginal = varMarginals[v.topo.hd];
                for (int i = 0; i < marginal.size(); i++){
                    if (marginal[i] < curCost){
                        label = i;
                        curCost = marginal[i];
                    }
                }
                results[v.topo.hd] = label;
            }

            return results;
        }


    }
}
 