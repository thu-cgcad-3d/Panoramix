#pragma once

#include "../misc/matlab_api.hpp"
#include "pi_graph.hpp"
#include "pi_graph_annotation.hpp"

namespace pano {
    namespace experimental {
        
        void ReconstructLayoutAnnotation(PILayoutAnnotation & anno, misc::Matlab & matlab);

        struct PIConstraintGraph {
            struct Entity {
                enum Type {
                    IsSeg,
                    IsLine
                } type;
                bool isSeg() const { return type == IsSeg; }
                bool isLine() const { return type == IsLine; }
                int id;
                //int dof;
                //int orientationNotClaz; // for isSeg() && dof == 2
                double size;
                template <class Archiver>
                void serialize(Archiver & ar) { ar(type, id/*, dof, orientationNotClaz*/, size); }
            };
            struct Constraint {
                enum Type {
                    SegCoplanarity,
                    Connection
                } type;
                bool isSegCoplanarity() const { return type == SegCoplanarity; }
                bool isConnection() const { return type == Connection; }
                int ent1, ent2;
                double weight;
                std::vector<Vec3> anchors;
                //int anchorOrientation; // for anchors.size() >= 2
                template <class Archiver>
                void serialize(Archiver & ar) { ar(type, ent1, ent2, weight, anchors/*, anchorOrientation*/); }
            };

            std::vector<Entity> entities;
            std::vector<Constraint> constraints;

            std::vector<int> seg2ent;
            std::vector<int> line2ent;

            std::set<int> determinableEnts;

            template <class Archiver>
            inline void serialize(Archiver & ar) {
                ar(entities, constraints, seg2ent, line2ent, determinableEnts);
            }
        };


        PIConstraintGraph BuildPIConstraintGraph(const PIGraph & mg, 
            double minAngleThresForAWideEdge, double angleThresForColinearity);

        double Solve(PIGraph & mg, const PIConstraintGraph & cg, misc::Matlab & matlab);

    }
}