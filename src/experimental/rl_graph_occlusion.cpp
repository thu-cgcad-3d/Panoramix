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



    }
}