#pragma once

#include "../gui/canvas.hpp"

#include "rl_graph.hpp"
#include "rl_graph_control.hpp"

namespace pano {
    namespace experimental {

        // depth ordering
        enum class DepthRelationGuess {
            FirstMaybeCloser,
            SecondMaybeCloser,
            MaybeConnected,
            Unknown
        };


        // guess line depth relations
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






        std::vector<std::vector<Vec3>> SamplesOnBoundaries(const SegmentationTopo & segtopo,
            const PanoramicCamera & cam, double sampleStepAngle = DegreesToRadians(1));


        std::vector<int> ClassifyBoundaries(const std::vector<std::vector<Vec3>> & bndSamples,
            const std::vector<Vec3> & vps, double angleThres = DegreesToRadians(1.5));



        struct TStructure {
            int centerJunctionId;
            int longEndJunctionIds[2];
            int shortEndJunctionId;
            std::vector<int> longBndIds[2];
            std::vector<int> shortBndIds;
            int longVPId;
            int shortVPId;

            template <class Archiver>
            void serialize(Archiver & ar) {
                ar(centerJunctionId, longEndJunctionIds, shortEndJunctionId, 
                    longBndIds, shortBndIds, longVPId, shortVPId);
            }
        };


        // find TStructures in graph
        std::vector<TStructure> FindTStructures(const SegmentationTopo & segtopo,
            const PanoramicCamera & cam, const std::vector<Vec3> & vps, 
            double minSpanAngle = DegreesToRadians(20), double angleThres = DegreesToRadians(1.5));

        std::vector<TStructure> FindTStructuresFuzzy(const SegmentationTopo & segtopo,
            const std::vector<std::vector<Vec3>> & bndsamples,
            const std::vector<int> & bndclasses,
            const PanoramicCamera & cam, const std::vector<Vec3> & vps,
            double minSpanAngle, 
            double angleThres /*= DegreesToRadians(1.5)*/);


        // make mask from TStructure
        void TStructureMaskView(const PanoramicCamera & cam, 
            const SegmentationTopo & segtopo, const TStructure & ts, 
            std::array<View<PartialPanoramicCamera, Imageub>, 2> & smallParts,
            View<PartialPanoramicCamera, Imageub> & largePart,
            double focal = 50.0, double largePartWidthAngle = DegreesToRadians(5));


        // GuessOccludedTStructures
        std::vector<int> GuessOccludedTStructures(const std::vector<Vec3> & vps, 
            const SegmentationTopo & segtopo, const std::vector<TStructure> & tstructs, const PanoramicCamera & cam,
            const Image5d & gc);
    
        // ApplyOccludedTStructure
        void ApplyOccludedTStructure(const RLGraph & mg, RLGraphControls & controls, 
            const SegmentationTopo & segtopo, 
            const std::vector<RegionBoundaryHandle> & bhs,
            const PanoramicCamera & camera, const std::vector<TStructure> & occtstructs);

        // ShowTStructures
        void ShowTStructures(gui::Canvas3ub &canvas,
            const SegmentationTopo & segtopo, const std::vector<TStructure> & occtstructs);




    }
}

