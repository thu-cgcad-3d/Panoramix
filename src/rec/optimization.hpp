#ifndef PANORAMIX_REC_OPTIMIZATION_HPP
#define PANORAMIX_REC_OPTIMIZATION_HPP
 
#include <Eigen/Core>
#include <Eigen/Dense>
#include <Eigen/Sparse>

#include <unsupported/Eigen/NonLinearOptimization>

#include <glpk.h>
#include <setjmp.h>

#include "../deriv/derivative.hpp"

#include "../core/basic_types.hpp"
#include "../core/graphical_model.hpp"
#include "../core/utilities.hpp"


namespace panoramix {
    namespace rec {

        // for glpk
        struct sinfo {
            char * text;
            jmp_buf * env;
        };

        void glpErrorHook(void * in){
            sinfo * info = (sinfo*)in;
            glp_free_env();
            longjmp(*(info->env), 1);
        }
       


        


        //// optimize constraint graph
        //template <class ComponentDataT, class ConstraintDataT, class ComponentUpdaterT,
        //class ComponentExpressionMakerT, class ConstraintExpressionMakerT, class CostExpressionMakerT>
        //int OptimizeConstraintGraph(
        //    core::ConstraintGraph<ComponentDataT, ConstraintDataT> & consGraph, 
        //    ComponentUpdaterT compUpdator = ComponentUpdaterT(),
        //    ComponentExpressionMakerT compExprMaker = ComponentExpressionMakerT(),
        //    ConstraintExpressionMakerT consExprMaker = ConstraintExpressionMakerT(),
        //    CostExpressionMakerT costExprMaker = CostExpressionMakerT(),
        //    OptimizeUsingGradientTag tag = OptimizeUsingGradientTag()) {
        //    
        //    using ConsGraph = core::ConstraintGraph<ComponentDataT, ConstraintDataT>;
        //    using CompHandle = typename ConsGraph::ComponentHandle;
        //    using ConsHandle = typename ConsGraph::ConstraintHandle;

        //    deriv::ExpressionGraph graph;
        //    
        //    using ComponentExprType = decltype(compExprMaker(graph, std::declval<ComponentDataT>()));
        //    using ComponentExprContentType = typename ComponentExprType::Type;            

        //    std::map<CompHandle, ComponentExprType, typename CompareHandle<core::ComponentTopo>> compExprTable;
        //    std::vector<ComponentExprType> compExprs;
        //    std::vector<CompHandle> compHandles;
        //    compExprs.reserve(consGraph.internalComponents().size());
        //    compHandles.reserve(consGraph.internalComponents().size());
        //    for (auto & comp : consGraph.components()) {
        //        compExprs.push_back(compExprTable[comp.topo.hd] = compExprMaker(graph, comp.data).handle());
        //        comphandles.push_back(comp.topo.hd);
        //    }
        //    
        //    using ConstraintExprType = decltype(consExprMaker(graph, std::declval<ConstraintDataT>(), std::vector<ComponentExprType>()));
        //    using ConstraintExprContentType = typename ConstraintExprType::Type;

        //    std::map<ConsHandle, ConstraintExprType, typename CompareHandle<core::ConstraintTopo>> consExprTable;
        //    std::vector<ConstraintExprType> consExprs;
        //    consExprs.reserve(consGraph.internalConstraints().size());
        //    for (auto & cons : consGraph.constraints()) {
        //        std::vector<ComponentExprType> inputs;
        //        inputs.reserve(cons.topo.components.size());
        //        for (auto & comp : cons.topo.components) {
        //            inputs.push_back(compExprTable[comp]);
        //        }
        //        consExprs.push_back(consExprTable[cons.topo.hd] = consExprMaker(graph, cons.data, inputs));
        //    }

        //    using CostExprType = decltype(costExprMaker(graph, std::vector<ConstraintExprType>()));
        //    using CostExprContentType = typename CostExprType::Type;

        //    // get the final cost function
        //    CostExprType costExpr = costExprMaker(graph, consExprs);

        //    // compute gradients
        //    using DerivedComponentExprType = deriv::DerivativeExpression<ComponentExprContentType>;
        //    std::vector<DerivedComponentExprType> derivs;
        //    derivs.reserve(compExprs.size());
        //    costExpr.derivativesRange(compExprs.begin(), compExprs.end(), std::back_inserter(derivs));
        //    assert(derivs.size() == compExprs.size());

        //    // iterate
        //    int epoches = 0;
        //    while (true) {
        //        bool shouldStop = false;
        //        for (int i = 0; i < derivs.size(); i++){
        //            if (compUpdator(consGraph.data(compHandles[i]), derivs[i].execute())){
        //                shouldStop = true;
        //            }
        //        }
        //        if (shouldStop)
        //            break;
        //        epoches++;
        //    }

        //    return epoches;
        //}


        

    }
}
 
#endif