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



        enum class DepthRelation {
            FirstIsFront,
            SecondIsFront,
            Connected,
            Disconnected,
            MaybeFolder,
            Unknown
        };

        std::vector<std::vector<Vec3>> SamplesOnBoundaries(const SegmentationTopo & segtopo,
            const PanoramicCamera & cam, double sampleStepAngle = DegreesToRadians(1));


        std::vector<int> ClassifyBoundaries(const std::vector<std::vector<Vec3>> & bndSamples,
            const std::vector<Vec3> & vps, double angleThres = DegreesToRadians(1.5));





        HandledTable<RegionBoundaryHandle, DepthRelation> DetectOcclusions(
            const RLGraph & mg, const RLGraphControls & controls, 
            const SegmentationTopo & segtopo, const std::vector<int> & bndclasses,
            const std::vector<RegionHandle> & rhs, const std::vector<RegionBoundaryHandle> & bhs,
            const std::vector<Vec3> & vps);

        HandledTable<RegionBoundaryHandle, DepthRelation> DetectOcclusions2(
            const RLGraph & mg, const RLGraphControls & controls,
            const Imagei & segs,
            const SegmentationTopo & segtopo, const std::vector<std::vector<Vec3>> & bndsamples, const std::vector<int> & bndclasses,
            const std::vector<RegionHandle> & rhs, const std::vector<RegionBoundaryHandle> & bhs,
            const std::vector<Vec3> & vps, double angleDistThres = DegreesToRadians(1));


        HandledTable<RegionBoundaryHandle, DepthRelation> DetectOcclusions3(
            const RLGraph & mg, const RLGraphControls & controls,
            const Imagei & segs,
            const SegmentationTopo & segtopo, const std::vector<std::vector<Vec3>> & bndsamples, const std::vector<int> & bndclasses,
            const std::vector<RegionHandle> & rhs, const std::vector<RegionBoundaryHandle> & bhs,
            const std::vector<Vec3> & vps, double angleDistThres = DegreesToRadians(1), double angleSampleStepOnLine = DegreesToRadians(1));

        
        struct OrientationControl {
            bool used;
            int orientationClaz;
            int orientationNotClaz;
            OrientationControl() : used(true), orientationClaz(-1), orientationNotClaz(-1) {}
        };
        HandledTable<RegionBoundaryHandle, DepthRelation> DetectOcclusions4(
            const RLGraph & mg, const HandledTable<RegionHandle, OrientationControl> & ocontrols,
            const Imagei & segs,
            const SegmentationTopo & segtopo, const std::vector<std::vector<Vec3>> & bndsamples, const std::vector<int> & bndclasses,
            const std::vector<RegionHandle> & rhs, const std::vector<RegionBoundaryHandle> & bhs,
            const std::vector<Vec3> & vps, double angleDistThres = DegreesToRadians(1), double angleSampleStepOnLine = DegreesToRadians(1));


        HandledTable<RegionBoundaryHandle, DepthRelation> DetectOcclusions5(
            const RLGraph & mg, const HandledTable<RegionHandle, OrientationControl> & ocontrols,
            const Imagei & segs,
            const SegmentationTopo & segtopo, const std::vector<std::vector<Vec3>> & bndsamples, const std::vector<int> & bndclasses,
            const std::vector<RegionHandle> & rhs, const std::vector<RegionBoundaryHandle> & bhs,
            const std::vector<Vec3> & vps, 
            double angleDistThres = DegreesToRadians(1), double angleSampleStepOnLine = DegreesToRadians(1));


        HandledTable<RegionBoundaryHandle, DepthRelation> DetectOcclusions6(
            const RLGraph & mg, const HandledTable<RegionHandle, OrientationControl> & ocontrols,
            const Imagei & segs,
            const SegmentationTopo & segtopo, const std::vector<std::vector<Vec3>> & bndsamples, const std::vector<int> & bndclasses,
            const std::vector<RegionHandle> & rhs, const std::vector<RegionBoundaryHandle> & bhs,
            const std::vector<Vec3> & vps,
            double angleDistThres = DegreesToRadians(1), double angleSampleStepOnLine = DegreesToRadians(1));


        void ApplyOcclusions(const RLGraph & mg, RLGraphControls & controls,
            const HandledTable<RegionBoundaryHandle, DepthRelation> & occlusions, bool spreadOnLineEachTime = true);

        void ApplyOcclusions2(const RLGraph & mg, RLGraphControls & controls,
            const HandledTable<RegionBoundaryHandle, DepthRelation> & occlusions);


        void DisableTJunctionsInLineRelations(const RLGraph & mg, RLGraphControls & controls, double tjuncRatioThres = 0.1);

        

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

