#ifndef PANORAMIX_EXPERIMENTAL_RL_GRAPH_OCCLUSION_HPP
#define PANORAMIX_EXPERIMENTAL_RL_GRAPH_OCCLUSION_HPP

#include"rl_graph.hpp"

namespace panoramix {
    namespace experimental {

        // depth ordering
        enum class DepthRelationGuess {
            FirstMaybeCloser,
            SecondMaybeCloser,
            MaybeConnected,
            Unknown
        };

        HandledTable<LineRelationHandle, DepthRelationGuess> GuessLineDepthRelation(
            const RLGraph & mg, double maxDistance);

        

    }
}

#endif