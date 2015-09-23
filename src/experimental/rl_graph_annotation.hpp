#include "../core/basic_types.hpp"
#include "../core/geometry.hpp"
#include "../core/cameras.hpp"

#include "rl_graph_control.hpp"

namespace pano {

    namespace experimental {

        using namespace core;

        struct OrientedPolygon {
            Polygon3 polygon;
            int towardVPId;
            int alongVPId;
            bool isClutter;
            template <class Archiver>
            inline void serialize(Archiver & ar) {
                ar(polygon, towardVPId, alongVPId, isClutter);
            }
        };

        struct PanoIndoorAnnotation {
            PanoramicView view;

            std::vector<Vec3> vps;
            int vertVPId;
            std::vector<Classified<Line3>> lines;

            std::vector<OrientedPolygon> polygons;
            std::vector<Chain3> occlusions;

            void reset() { 
                view.image = Image();
                vps.clear();
                vertVPId = 0;
                lines.clear();
                polygons.clear();
                occlusions.clear();
            }

            template <class Archiver>
            inline void serialize(Archiver & ar) {
                ar(view, vps, vertVPId, lines,
                    polygons, occlusions);
            }
        };

        void ProjectOn(const Polygon3 & polygon, const PanoramicCamera & cam, Imagei & canvas, int val);

        void AttachAnnotationConstraints(const RLGraph & mg,
            const PanoramicCamera & cam,
            const SegmentationTopo & segtopo,
            const std::vector<RegionHandle> & rhs,
            const std::vector<RegionBoundaryHandle> & bhs,
            RLGraphControls & controls,
            const std::vector<OrientedPolygon> & polygons,
            const std::vector<Chain3> & occlusions);




        // quantitative evaluations 
        struct QuantEvaluation {
            double surfaceOrientationEval;
            double depthEval;
        };
        QuantEvaluation Evaluate(const LayeredShape3 & shape, const PanoIndoorAnnotation & anno);


    }

}