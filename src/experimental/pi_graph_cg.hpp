#pragma once

#include "pi_graph_control.hpp"
#include "pi_graph_annotation.hpp"

namespace pano {
    namespace experimental {


        struct PIConstraintGraph {
            struct Entity {
                enum Type {
                    IsSeg,
                    IsLine
                } type;
                bool isSeg() const { return type == IsSeg; }
                bool isLine() const { return type == IsLine; }
                int id;
                struct SupportingPlane {
                    int dof;
                    Vec3 center;
                    Vec3 toward;
                    Vec3 along;
                    Plane3 reconstructed;
                    DenseMatd matFromVarsToPlaneCoeffs() const;
                    template <class Archiver>
                    void serialize(Archiver & ar) { ar(dof, center, toward, along, reconstructed); }
                } supportingPlane;
                double size;
                template <class Archiver>
                void serialize(Archiver & ar) { ar(type, id, supportingPlane, size); }
            };
            struct Constraint {
                enum Type {
                    Coplanarity,
                    Connection
                } type;
                bool isCoplanarity() const { return type == Coplanarity; }
                bool isConnection() const { return type == Connection; }
                int ent1, ent2;
                double weight;
                std::vector<Vec3> anchors;
                template <class Archiver>
                void serialize(Archiver & ar) { ar(type, ent1, ent2, weight, anchors); }
            };

            std::vector<Entity> entities;
            std::vector<Constraint> constraints;
            std::vector<std::vector<int>> ent2cons;
            std::vector<bool> cons2enabled;

            std::vector<int> seg2ent;
            std::vector<int> line2ent;

            template <class Archiver>
            inline void serialize(Archiver & ar) {
                ar(entities, constraints, ent2cons, cons2enabled, seg2ent, line2ent);
            }
        };

        struct PICGDeterminablePart {
            int rootEnt;
            std::set<int> determinableEnts;
            std::set<int> consBetweenDeterminableEnts;
            bool empty() const { return rootEnt == -1; }
            template <class Archiver>
            inline void serialize(Archiver & ar) {
                ar(rootEnt, determinableEnts, consBetweenDeterminableEnts);
            }
        };


        PIConstraintGraph BuildPIConstraintGraph(const PIGraph & mg,
            double minAngleThresForAWideEdge);

        PICGDeterminablePart LocateDeterminablePart(const PIConstraintGraph & cg, double angleThres);




        PIConstraintGraph BuildPIConstraintGraph(const PILayoutAnnotation & anno,
            double minAngleThresForAWideEdge);

        PIConstraintGraph BuildPIConstraintGraphWithLines(const PILayoutAnnotation & anno,
            double minAngleThresForAWideEdge);




    }
}