#include "../gui/visualize2d.hpp"

#include "rl_graph_occlusion.hpp"

namespace panoramix {
    namespace experimental {

        using namespace core;



        namespace {

            void ShrinkLine(Line3 & line, double ratio){
                double angle = AngleBetweenDirections(line.first, line.second);
                double dangle = angle * ratio;
                line = {
                    RotateDirection(line.first, line.second, dangle),
                    RotateDirection(line.second, line.first, dangle)
                };
            }

        }


        HandledTable<LineRelationHandle, DepthRelationGuess> GuessLineDepthRelation(const RLGraph & mg, double maxDistance){

            auto dr = mg.createConstraintTable<LineRelationData>(DepthRelationGuess::Unknown);

            // find T-junctions in line relations
            for (auto & lr : mg.constraints<LineRelationData>()){
                auto & lh1 = lr.topo.component<0>();
                auto & lh2 = lr.topo.component<1>();
                auto & ld1 = mg.data(lh1);
                auto & ld2 = mg.data(lh2);
                if (ld1.initialClaz != -1 && ld1.initialClaz == ld2.initialClaz){
                    // incidence // connected?
                    // todo
                    dr[lr.topo.hd] = DepthRelationGuess::MaybeConnected;
                    continue;
                }
                if (ld1.initialClaz == -1 && ld2.initialClaz == -1){
                    continue;
                }
                if (DistanceBetweenTwoLines(normalize(ld1.line), normalize(ld2.line)).first > maxDistance)
                    continue;

                Line3 sline1 = ld1.line; ShrinkLine(sline1, 0.3);
                Line3 sline2 = ld2.line; ShrinkLine(sline2, 0.3);

                Vec3 eq1 = ld1.line.first.cross(ld1.line.second);
                Vec3 eq2 = ld2.line.first.cross(ld2.line.second);
                Vec3 inter = normalize(eq1.cross(eq2));
                bool interInLine1 = (inter.cross(ld1.line.first)).dot(ld1.line.second.cross(inter)) >= 0;
                bool interInLine2 = (inter.cross(ld2.line.first)).dot(ld2.line.second.cross(inter)) >= 0;
                if (interInLine1 && interInLine2){
                    // todo
                    dr[lr.topo.hd] = DepthRelationGuess::MaybeConnected;
                    continue;
                }
                if (interInLine1 && !interInLine2){
                    // |-
                    // todo
                    dr[lr.topo.hd] = DepthRelationGuess::FirstMaybeCloser;
                    continue;
                }
                if (interInLine2 && !interInLine1){
                    // -|
                    // todo
                    dr[lr.topo.hd] = DepthRelationGuess::SecondMaybeCloser;
                    continue;
                }
            }

            return dr;

        }



        HandledTable<LineHandle, LineDetachStatus> GuessLineDetachStatus(
            const RLGraph & mg, const HandledTable<LineRelationHandle, DepthRelationGuess> & ldr){

            auto result = mg.createComponentTable<LineData>(LineDetachStatus{ false, false });
            for (auto & l : mg.components<LineData>()){
                if (l.data.initialClaz == -1){
                    continue;
                }
                Vec3 rightdir = normalize(l.data.line.first.cross(l.data.line.second));

                int rightOccCandNum = 0;
                int rightConNum = 0;
                int rightAll = 0;
                int leftOccCandNum = 0;
                int leftConNum = 0;
                int leftAll = 0;

                for (auto & lr : l.topo.constraints<LineRelationData>()){
                    const DepthRelationGuess & dr = ldr[lr];
                    auto another = mg.topo(lr).component<0>();
                    if (another == l.topo.hd){
                        another = mg.topo(lr).component<1>();
                    }
                    auto & anotherld = mg.data(another);                    

                    bool isOnRight = (anotherld.line.center() - l.data.line.center()).dot(rightdir) > 0;
                    if (mg.topo(lr).component<0>() == l.topo.hd && dr == DepthRelationGuess::FirstMaybeCloser ||
                        mg.topo(lr).component<1>() == l.topo.hd && dr == DepthRelationGuess::SecondMaybeCloser){
                        if (isOnRight){
                            rightOccCandNum++;
                        }
                        else{
                            leftOccCandNum++;
                        }
                    }

                    if (isOnRight){
                        rightAll++;
                    }
                    else{
                        leftAll++;
                    }

                    if (dr == DepthRelationGuess::MaybeConnected && mg.data(lr).type == LineRelationData::Intersection){
                        if (isOnRight){
                            rightConNum++;
                        }
                        else{
                            leftConNum++;
                        }

                    }
                }

                for (auto & rl : l.topo.constraints<RegionLineConnectionData>()){
                    bool isOnRight = (mg.data(mg.topo(rl).component<0>()).normalizedCenter - l.data.line.center()).dot(rightdir) > 0;
                    if (isOnRight){
                        rightAll++;
                    }
                    else{
                        leftAll++;
                    }
                }

                result[l.topo.hd].lineRightMayDetach = rightOccCandNum > 0 && rightConNum == 0 && 
                    leftAll > 0 &&
                    leftOccCandNum < rightOccCandNum;
                result[l.topo.hd].lineLeftMayDetach = leftOccCandNum > 0 && leftConNum == 0 &&
                    rightAll > 0 &&
                    rightOccCandNum < leftOccCandNum;
            }

            return result;

        }


        bool MayOccludes(const RLGraph & mg, const HandledTable<LineHandle, LineDetachStatus> & lrds,
            LineHandle lh,
            LineRelationHandle lrh){
            auto another = mg.topo(lrh).component<0>();
            if (another == lh){
                another = mg.topo(lrh).component<1>();
            }
            else{
                assert(mg.topo(lrh).component<1>() == lh);
            }
            auto & anotherld = mg.data(another);
            Vec3 rightdir = normalize(mg.data(lh).line.first.cross(mg.data(lh).line.second));
            bool isOnRight = (anotherld.line.center() - mg.data(lh).line.center()).dot(rightdir) > 0;
            if (isOnRight && lrds[lh].lineRightMayDetach)
                return true;
            if (!isOnRight && lrds[lh].lineLeftMayDetach)
                return true;
            return false;
        }

        bool MayOccludes(const RLGraph & mg, const HandledTable<LineHandle, LineDetachStatus> & lrds,
            LineHandle lh,
            RegionLineConnectionHandle rlh){
            auto another = mg.topo(rlh).component<0>();
            auto & anotherld = mg.data(another);
            Vec3 rightdir = normalize(mg.data(lh).line.first.cross(mg.data(lh).line.second));
            bool isOnRight = (anotherld.normalizedCenter - mg.data(lh).line.center()).dot(rightdir) > 0;
            if (isOnRight && lrds[lh].lineRightMayDetach)
                return true;
            if (!isOnRight && lrds[lh].lineLeftMayDetach)
                return true;
            return false;
        }



        template <class CameraT, class CallbackFunT>
        void AttachGeometricContextConstraintsTemplated(const RLGraph & mg, RLGraphControls & controls,
            const CameraT & pcam, const ImageOfType<Vec<double, 5>> & gc, const Imagei & gcVotes,
            CallbackFunT && fun){


            // set component orientations
            // use geometric context!
            SetComponentControl<RegionData>(controls,
                [&mg, &gc, &gcVotes, &pcam, &fun](RegionHandle rh, RLGraphComponentControl & c){
                auto regionMaskView = PerfectRegionMaskView(mg, rh);
                auto sampler = core::MakeCameraSampler(regionMaskView.camera, pcam);
                auto gcOnRegion = sampler(gc);
                auto gcVotesOnRegion = sampler(gcVotes);
                int votes = 0;
                core::Vec<double, 5> gcSum;
                for (auto it = regionMaskView.image.begin(); it != regionMaskView.image.end(); ++it){
                    if (!*it){
                        continue;
                    }
                    gcSum += gcOnRegion(it.pos());
                    votes += gcVotesOnRegion(it.pos());
                }
                auto gcMean = gcSum / std::max(votes, 1);
                auto gcMean2 = gcMean;
                std::sort(std::begin(gcMean2.val), std::end(gcMean2.val), std::greater<>());
                double significancy = gcMean2[0] / std::max(gcMean2[1], 1e-5);
                if (!fun){
                    if (significancy > 2){
                        int label = std::max_element(std::begin(gcMean.val), std::end(gcMean.val)) - gcMean.val;
                        // label: vp1 vp2 vp3 clutter other
                        if (label < 3){ // vp
                            c.orientationClaz = label;
                            c.orientationNotClaz = -1;
                        }
                        else{
                            c.used = false;
                            //c.orientationClaz = c.orientationNotClaz = -1;
                        }
                    }
                }
                else{
                    fun(c, gcMean, significancy);
                }
            });

            controls.disableAllInvalidConstraints(mg);
        }



        void AttachGeometricContextConstraints(const RLGraph & mg, RLGraphControls & controls,
            const PanoramicCamera & pcam, const ImageOfType<Vec<double, 5>> & gc, const Imagei & gcVotes,
            const std::function<void(RLGraphComponentControl &, const Vec<double, 5> &, double significancy)> & fun){
            AttachGeometricContextConstraintsTemplated(mg, controls, pcam, gc, gcVotes, fun);
        }
        void AttachGeometricContextConstraints(const RLGraph & mg, RLGraphControls & controls,
            const PerspectiveCamera & pcam, const ImageOfType<Vec<double, 5>> & gc,
            const std::function<void(RLGraphComponentControl &, const Vec<double, 5> &, double significancy)> & fun){
            AttachGeometricContextConstraintsTemplated(mg, controls, pcam, gc, Imagei::ones(gc.size()), fun);
        }




        void AttachGeometricContextConstraintsSmarter(const RLGraph & mg, RLGraphControls & controls,
            const PerspectiveCamera & pcam, const ImageOfType<Vec<double, 5>> & gc){



        }

    }
}