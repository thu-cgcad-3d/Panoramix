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

        void AttachGeometricContextConstraints(const RLGraph & mg, RLGraphControls & controls,
            const PanoramicCamera & pcam, const ImageOfType<Vec<double, 5>> & gc, const Imagei & gcVotes,
            const std::function<void(RLGraphComponentControl &, const Vec<double, 5> &, double significancy)> & fun = nullptr);

        void AttachGeometricContextConstraints(const RLGraph & mg, RLGraphControls & controls,
            const PerspectiveCamera & pcam, const ImageOfType<Vec<double, 5>> & gc,
            const std::function<void(RLGraphComponentControl &, const Vec<double, 5> &, double significancy)> & fun = nullptr);

        void AttachGeometricContextConstraintsSmarter(const RLGraph & mg, RLGraphControls & controls,
            const PerspectiveCamera & pcam, const ImageOfType<Vec<double, 5>> & gc);

    }
}

#endif