#pragma once

#include "../misc/matlab_api.hpp"

#include "rl_graph.hpp"
#include "rl_graph_control.hpp"


namespace pano {
    namespace experimental {


        RLGraphVars MakeVariables(const RLGraph & mg, const RLGraphControls & controls, bool randomized = false);
        bool HasInfOrNaNValue(const RLGraphVars & v);




        // solve equations
        RLGraphVars SolveVariablesWithoutBoundedAnchors(const RLGraph & mg,
            const RLGraphControls & controls,
            bool useWeights = false);

        void ResetToFullArmorAnchors(const RLGraph & mg, RLGraphControls & controls);
        void ResetToSampledArmorAnchors(const RLGraph & mg, RLGraphControls & controls, double sampleStepAngle);

        RLGraphVars SolveVariablesWithBoundedAnchors(misc::Matlab & matlab, const RLGraph & mg, const RLGraphControls & controls,
            bool useWeights = false, int tryNum = 10);

        void OptimizeVariablesWithBoundedAnchors(misc::Matlab & matlab, const RLGraph & mg, const RLGraphControls & controls, RLGraphVars & vars,
            bool useWeights = true, int tryNum = 100, int maxOptimizeNum = 10,
            const std::function<bool(const RLGraphVars &)> & callback = nullptr);      





    }
}


