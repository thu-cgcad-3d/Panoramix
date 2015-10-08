#pragma once

#include "../misc/matlab_api.hpp"
#include "pi_graph.hpp"
#include "pi_graph_annotation.hpp"

namespace pano {
    namespace experimental {

#if 0
        void BuildConstraintGraph(PIGraph & mg);


        std::vector<double> InverseDepthCoefficientsOfSegAtDirection(const PIGraph & mg,
            int seg, const Vec3 & direction);
        std::vector<double> InverseDepthCoefficientsOfLineAtDirection(const PIGraph & mg,
            int line, const Vec3 & direction);


        void SolvePIGraph(int ccid, PIGraph & mg, misc::Matlab & matlab, int tryNum);
#endif

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
                template <class Archiver>
                void serialize(Archiver & ar) { ar(type, id); }
            };
            struct Constraint {
                enum Type {
                    SegCoplanarity,
                    Connection
                } type;
                int ent1, ent2;
                double weight;
                std::vector<Vec3> anchors;
                template <class Archiver>
                void serialize(Archiver & ar) { ar(type, ent1, ent2, weight, anchors); }
            };

            std::vector<Entity> entities;
            std::vector<Constraint> constraints;

            std::vector<int> seg2ent;
            std::vector<int> line2ent;

            std::vector<int> ent2cc;
            std::map<int, std::vector<int>> cc2ents;
            int nccs;
            std::vector<int> ccsBigToSmall;

            PIGraph & mg;
            PIConstraintGraph(PIGraph & g);
            void reconstructLargestCC() const;
        };

    }
}