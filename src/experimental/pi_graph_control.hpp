#pragma once

#include "pi_graph.hpp"

namespace pano {
    namespace experimental {

        void AttachPrincipleDirectionConstraints(PIGraph & mg);
        void AttachWallConstraints(PIGraph & mg, double angle);

        void DisableTopSeg(PIGraph & mg);
        void DisableBottomSeg(PIGraph & mg);

        void AttachGCConstraints(PIGraph & mg, const Image5d & gc);

        // annotation
        struct AnnotedPolygon {
            Polygon3 polygon;
            SegControl control;
            template <class Archiver>
            inline void serialize(Archiver & ar) {
                ar(polygon, control);
            }
        };
        struct AnnotedOcclusion {
            Chain3 chain;
            bool leftIsFront;
            template <class Archiver>
            inline void serialize(Archiver & ar) {
                ar(chain, leftIsFront);
            }
        };

        struct PIAnnotation {
            std::vector<AnnotedPolygon> polygons;
            std::vector<AnnotedOcclusion> occlusions;
            template <class Archiver>
            inline void serialize(Archiver & ar) {
                ar(polygons, occlusions);
            }
        };

        void AttachAnnotations(PIGraph & mg, const PIAnnotation & anno);
 
    }
}