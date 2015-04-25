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

        struct LineDetachStatus {
            bool lineRightMayDetach;
            bool lineLeftMayDetach;
        };
        HandledTable<LineHandle, LineDetachStatus> GuessLineDetachStatus(
            const RLGraph & mg, const HandledTable<LineRelationHandle, DepthRelationGuess> & ldr);
        inline HandledTable<LineHandle, LineDetachStatus> GuessLineDetachStatus(
            const RLGraph & mg, double maxDistance) {
            return GuessLineDetachStatus(mg, GuessLineDepthRelation(mg, maxDistance));
        }
        
        bool MayOccludes(const RLGraph & mg, const HandledTable<LineHandle, LineDetachStatus> & lrds,
            LineHandle lh,
            LineRelationHandle lrh);
        bool MayOccludes(const RLGraph & mg, const HandledTable<LineHandle, LineDetachStatus> & lrds,
            LineHandle lh,
            RegionLineConnectionHandle rlh);




        using RLGraphGCTable = HandledTable<RegionHandle, Vec<double, 5>>;


        void AttachGeometricContextConstraints(const RLGraph & mg, RLGraphControls & controls,
            const PanoramicCamera & pcam, const ImageOfType<Vec<double, 5>> & gc, const Imagei & gcVotes,
            const std::function<void(RLGraphComponentControl &, const Vec<double, 5> &, double significancy)> & fun = nullptr);

        void AttachGeometricContextConstraints(const RLGraph & mg, RLGraphControls & controls,
            const PerspectiveCamera & pcam, const ImageOfType<Vec<double, 5>> & gc,
            const std::function<void(RLGraphComponentControl &, const Vec<double, 5> &, double significancy)> & fun = nullptr);


    }
}

#endif