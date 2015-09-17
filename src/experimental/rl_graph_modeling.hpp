#pragma once

#include "../ml/annotation.hpp"
#include "rl_graph_solver.hpp"

namespace pano {

    namespace experimental {

        // cut loop
        struct SectionalPiece {
            RegionHandle rh;
            std::pair<Point3, Point3> range;
        };
        std::vector<SectionalPiece> MakeSectionalPieces(const HandledTable<RegionHandle, std::vector<Polygon3>> & polygons,
            const Plane3 & cutplane, bool removeIllPosedPolygons = false);

        Chain3 MakeChain(const std::vector<SectionalPiece> & pieces, bool closed = true);
        inline double Area(const std::vector<SectionalPiece> & loop) { return Area(Polygon3(MakeChain(loop))); }



        // optimize loops
        LayeredShape3 OptimizedModel(const std::vector<SectionalPiece> & pieces, 
            const RLGraph & mg, const RLGraphControls & controls);



        // estimate 
        std::pair<double, double> EstimateEffectiveRangeAlongDirection(
            const HandledTable<RegionHandle, std::vector<Polygon3>> & polygons,
            const Vec3 & direction, double stepLen, double minEffectiveAreaRatio = 0.6,
            double gamma1 = 0.05, double gamma2 = 0.05);



    }

}