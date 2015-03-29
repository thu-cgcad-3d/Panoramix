#include <Eigen/Dense>
#include <Eigen/StdVector>

#include "../core/cons_graph.hpp"
#include "factor_graph.hpp"
 
namespace panoramix {
    namespace ml {

        using namespace core;

        bool FactorGraph::valid() const {
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
                auto & costFun = factors[fid].costs;
                if (!costFun)
                    return false;
            }
            return true;
        }

        double FactorGraph::energy(const ResultTable & labels) const{
            assert(valid());            
            double e = 0.0;
            for (auto & f : graph.elements<1>()){
                auto & vhs = f.topo.lowers;
                std::vector<int> idx(vhs.size());
                for (int i = 0; i < idx.size(); i++){
                    idx[i] = labels[vhs[i]];
                }
                double factorCost = factors[f.data].costs(idx.data());
                e += factorCost;
            }
            return e;
        }


        FactorGraph::ResultTable FactorGraph::solve(int maxEpoch, int innerLoopNum, const CallbackFunction & callback) const {

            // convex belief propagation

            assert(valid());
            using namespace Eigen;
            
            struct VHandleWrapper { VarHandle h; };
            struct FHandleWrapper { FactorHandle h; };
            struct F2VMessage { VectorXd values; };
            struct V2FMessage { VectorXd values; }; 

            ConstraintGraph<
                std::tuple<VHandleWrapper, FHandleWrapper>,
                std::tuple<
                ConstraintConfig<V2FMessage, VHandleWrapper, FHandleWrapper>,
                ConstraintConfig<F2VMessage, FHandleWrapper, VHandleWrapper>
                >
            > messages;
            {
                messages.internalComponents<VHandleWrapper>().reserve(graph.internalElements<0>().size());
                messages.internalComponents<FHandleWrapper>().reserve(graph.internalElements<1>().size());
                size_t vfConnectionNum = 0;
                for (auto & f : graph.elements<1>()){
                    vfConnectionNum += f.topo.lowers.size();
                }
                messages.internalConstraints<V2FMessage>().reserve(vfConnectionNum);
                messages.internalConstraints<F2VMessage>().reserve(vfConnectionNum);
            }

            using MGVHandle = ComponentHandle<VHandleWrapper>;
            using MGFHandle = ComponentHandle<FHandleWrapper>;
            using MGV2FHandle = ConstraintHandle<V2FMessage>;
            using MGF2VHandle = ConstraintHandle<F2VMessage>;

            // install graph to message graph and initialize messages
            HandledTable<VarHandle, MGVHandle> vhToMGVH(graph.internalElements<0>().size());
            for (auto & v : graph.elements<0>()){
                MGVHandle vh = messages.addComponent(VHandleWrapper{ v.topo.hd });
                vhToMGVH[v.topo.hd] = vh;
            }
            for (auto & f : graph.elements<1>()){
                MGFHandle fh = messages.addComponent(FHandleWrapper{ f.topo.hd });
                for (auto & lh : f.topo.lowers){
                    MGVHandle vh = vhToMGVH[lh];
                    messages.addConstraint(V2FMessage{ VectorXd::Zero(vars[graph.data(lh)].nlabels) }, vh, fh);
                    messages.addConstraint(F2VMessage{ VectorXd::Zero(vars[graph.data(lh)].nlabels) }, fh, vh);
                }
            }

            // precompute c_i_hats for each var
            HandledTable<VarHandle, double> c_i_hats(graph.internalElements<0>().size(), 0.0);
            for (auto & v : graph.elements<0>()){
                double & c_i_hat = c_i_hats[v.topo.hd];
                c_i_hat = vars[v.data].c_i;
                for (auto fh : v.topo.uppers){
                    c_i_hat += factors[graph.data(fh)].c_alpha;
                }
            }

            // store ordered f2v and v2f handles for each factor
            core::HandledTable<MGFHandle, std::vector<MGF2VHandle>> orderedF2Vmsghs(graph.internalElements<1>().size());
            core::HandledTable<MGFHandle, std::vector<MGV2FHandle>> orderedV2Fmsghs(graph.internalElements<1>().size());
            // store var dimensions of each factor
            core::HandledTable<MGFHandle, std::vector<size_t>> orderedVDims(graph.internalElements<1>().size());
            for (auto & f : messages.components<FHandleWrapper>()){
                MGFHandle mgfh = f.topo.hd;
                FactorHandle fh = f.data.h;
                auto & vhs = graph.topo(fh).lowers;
                auto & v2fmsghsSet = messages.topo(mgfh).constraints<V2FMessage>();
                auto & v2fmsghs = orderedV2Fmsghs[mgfh];
                v2fmsghs.resize(vhs.size());
                auto & f2vmsghsSet = messages.topo(mgfh).constraints<F2VMessage>();
                auto & f2vmsghs = orderedF2Vmsghs[mgfh];
                f2vmsghs.resize(vhs.size());
                for (int i = 0; i < vhs.size(); i++){
                    for (MGV2FHandle v2fmsgh : v2fmsghsSet){
                        if (messages.topo(v2fmsgh).component<0>() == vhToMGVH[vhs[i]]){
                            v2fmsghs[i] = v2fmsgh;
                            break;
                        }
                    }
                    for (MGF2VHandle f2vmsgh : f2vmsghsSet){
                        if (messages.topo(f2vmsgh).component<1>() == vhToMGVH[vhs[i]]){
                            f2vmsghs[i] = f2vmsgh;
                            break;
                        }
                    }
                }
                orderedVDims[mgfh].resize(vhs.size());
                for (int i = 0; i < vhs.size(); i++){
                    orderedVDims[mgfh][i] = vars[graph.data(vhs[i])].nlabels;
                }
            }

           


            // initialize marginals
            core::HandledTable<VarHandle, VectorXd> varMarginals(graph.internalElements<0>().size());
            for (auto & v : graph.elements<0>()){
                varMarginals[v.topo.hd] = VectorXd::Zero(vars[v.data].nlabels);
            }

            double lastE = std::numeric_limits<double>::max();
            ResultTable results(graph.internalElements<0>().size());

            for (int epoch = 0; epoch < maxEpoch; epoch ++){
                for (int l = 0; l < innerLoopNum; l++){
                    // var -> factor
                    for (auto & v2f : messages.constraints<V2FMessage>()){
                        MGVHandle vh = v2f.topo.component<0>();
                        MGFHandle fh = v2f.topo.component<1>();

                        auto & fdata = factors[graph.data(messages.data(fh).h)];

                        v2f.data.values.setZero();
                        MGF2VHandle oppose;
                        for (MGF2VHandle f2v : messages.topo(vh).constraints<F2VMessage>()){
                            v2f.data.values += messages.data(f2v).values;
                            if (messages.topo(f2v).component<0>() == fh){
                                oppose = f2v;
                            }
                        }
                        assert(oppose.valid());
                        v2f.data.values = v2f.data.values * fdata.c_alpha / c_i_hats[messages.data(vh).h]
                            - messages.data(oppose).values;
                    }

                    // factor -> var
                    for (auto & f : messages.components<FHandleWrapper>()){
                        MGFHandle mgfh = f.topo.hd;
                        FactorHandle fh = f.data.h;
                        auto & vhs = graph.topo(fh).lowers;

                        const std::vector<MGV2FHandle> & v2fmsghs = orderedV2Fmsghs[mgfh];
                        const std::vector<MGF2VHandle> & f2vmsghs = orderedF2Vmsghs[mgfh];
                        assert(vhs.size() == v2fmsghs.size());
                        assert(vhs.size() == f2vmsghs.size());

                        // reset f2v messages
                        for (auto & f2v : f2vmsghs){
                            messages.data(f2v).values.setConstant(std::numeric_limits<double>::max());
                        }

                        // dispatch messages from this factor
                        int fid = graph.data(fh);
                        auto & costFun = factors[fid].costs;
                        auto & dims = orderedVDims[mgfh];
                        assert(dims.size() > 0);
                        std::vector<int> idx(dims.size(), 0);
                        assert(vhs.size() == dims.size());
                        while (true){
                            for (int i = 0; i < dims.size(); i++){
                                // i: output var
                                // others: input vars
                                double & outValue = messages.data(f2vmsghs[i]).values[idx[i]];
                                // theta
                                double theta = costFun(idx.data());
                                // sum of other input vars
                                double sumOfOtherVars = 0.0;
                                for (int j = 0; j < dims.size(); j++){
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
                            while (k < idx.size() && idx[k] >= dims[k]){
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
                    varMarginals[messages.data(vh).h] += f2v.data.values;                    
                }
                // get result                
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

                double e = energy(results);
                if (callback && !callback(epoch, e, e - lastE)){
                    break;
                }
                lastE = e;
            }            

            return results;
        }


    }
}
 