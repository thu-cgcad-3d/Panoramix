#include "../core/basic_types.hpp"
#include "../core/geometry.hpp"
#include "../core/cameras.hpp"

namespace pano {

    namespace experimental {

        using namespace core;

        struct PanoIndoorAnnotation {
            PanoramicView view;

            std::vector<Vec3> vps;
            int vertVPId;
            std::vector<Classified<Line3>> lines;

            std::vector<Polygon3> polygons;
            std::vector<int> polygonTowardVPIds;
            std::vector<int> polygonAlongVPIds;
            std::vector<bool> polygonAreClutters;
            std::vector<Chain3> occlusions;

            void reset() { 
                view.image = Image();
                vps.clear();
                vertVPId = 0;
                lines.clear();
                polygons.clear();
                polygonTowardVPIds.clear();
                polygonAlongVPIds.clear();
                polygonAreClutters.clear();
                occlusions.clear();
            }

            template <class Archiver>
            inline void serialize(Archiver & ar) {
                ar(view, vps, vertVPId, lines,
                    polygons, polygonTowardVPIds, polygonAlongVPIds, polygonAreClutters, occlusions);
            }
        };


        // quantitative evaluations 
        struct QuantEvaluation {
            double surfaceOrientationEval;
            double depthEval;
        };
        QuantEvaluation Evaluate(const LayeredShape3 & shape, const PanoIndoorAnnotation & anno);


    }

}