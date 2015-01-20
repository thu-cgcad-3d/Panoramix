#ifndef PANORAMIX_CORE_P3GRAPH_HPP
#define PANORAMIX_CORE_P3GRAPH_HPP
 
#include "optimization.hpp"
#include "cons_graph.hpp"
#include "view.hpp"

// perspective projection planer recovery graph

namespace panoramix {
    namespace core {

        // environment configurations
        struct P3Environment {
            std::vector<Vec3> vanishingPoints;
            template <class Archive> void serialize(Archive & ar) {
                ar(vanishingPoints);
            }
        };

        // unary variables
        // free line
        struct P3FreeLine {
            Line3 line;
            template <class Archive> 
            void serialize(Archive & ar) {
                ar(line);
            }
        };
        template <class ValuesT, class NamesT, class EnvironmentT>
        inline void RegisterValuesAndNames(const P3FreeLine & u, 
            ValuesT & values, NamesT & names, const EnvironmentT & env){
            values = { 1.0 / norm(u.line.first), 1.0 / norm(u.line.second) };
            names = { "inverseDepth1, inverseDepth2" };
        }
        template <class GetValuesT, class EnvironmentT >
        inline void UpdateValues(P3FreeLine & u, GetValuesT getValues, const EnvironmentT & env){
            u.line.first = normalize(u.line.first) / getValues(0);
            u.line.second = normalize(u.line.second) / getValues(1);
        }
        double InverseDepthAtDirection(const Vec3 & direction, double v1, double v2,
            const P3FreeLine & u);



        // classified line
        struct P3OrientedLine {
            Line3 line;
            int claz;
            template <class Archive>
            void serialize(Archive & ar) {
                ar(line, claz);
            }
        };
        template <class ValuesT, class NamesT, class EnvironmentT>
        inline void RegisterValuesAndNames(const P3OrientedLine & u, 
            ValuesT & values, NamesT & names, const EnvironmentT & env){
            values = { 1.0 / norm(u.line.center()) };
            names = { "inverseDepthOfCenter" };
        }
        template <class GetValuesT, class EnvironmentT >
        inline void UpdateValues(P3OrientedLine & u, GetValuesT getValues, const EnvironmentT & env){
            InfiniteLine3 infLine(u.line.center() / getValues(0), env.vanishingPoints[u.claz]);
            u.line = Line3(DistanceBetweenTwoLines(InfiniteLine3(Point3(0, 0, 0), u.line.first), infLine).second.second,
                DistanceBetweenTwoLines(InfiniteLine3(Point3(0, 0, 0), u.line.second), infLine).second.second);
        }
       /* double InverseDepthAtDirection(const Vec3 & direction, double v,
            const P3OrientedLine & u, const std::vector<Vec3> & vps);*/





        // region
        // free region
        struct P3FreeRegion {
            Plane3 plane;
            template <class Archive> 
            void serialize(Archive & ar) {
                ar(plane);
            }
        };
        template <class ValuesT, class NamesT, class EnvironmentT>
        inline void RegisterValuesAndNames(const P3FreeRegion & u,
            ValuesT & values, NamesT & names, const EnvironmentT & env){
            Vec3 eq = Plane3ToEquation(u.plane);
            values = { eq[0], eq[1], eq[2] };
            names = { "a", "b", "c" };
        }
        template <class GetValuesT, class EnvironmentT>
        inline void UpdateValues(P3FreeRegion & u, GetValuesT getValues, const EnvironmentT & env){
            u.plane = Plane3FromEquation(getValues(0), getValues(1), getValues(2));
        }


        
        // classified region





        // binary variable       
        struct P3RegionLineConnection {
            std::vector<Vec3> anchors;
            template <class Archive> 
            void serialize(Archive & ar) {
                ar(anchors);
            }
        };
        using P3RegionLineConnectionConfig = core::ConstraintConfig<P3RegionLineConnection,
            core::ComponentOccupation<P3FreeRegion, 1>,
            core::ComponentOccupation<P3FreeLine, 1>
        >;
        template <class TopoT, class CGT, class AddConsT, class EnvironmentT>
        void RegisterConstraint(const P3RegionLineConnection & b, const TopoT & topo, const CGT & cg,
            AddConsT & addCons, const EnvironmentT & env){
            const P3FreeRegion & region = cg.data(topo.components<P3FreeRegion>()[0]);
            const P3FreeLine & line = cg.data(topo.components<P3FreeLine>()[0]);            
            for (auto & a : b.anchors){
                //addCons({ {}, {}, {} })
            }
        }






        struct P3RegionRegionConnection {
            std::vector<Vec3> anchors;
            template <class Archive> void serialize(Archive & ar) {
                ar(anchors);
            }
        };
        using P3RegionRegionConnectionConfig = core::ConstraintConfig<P3RegionRegionConnection,
            core::ComponentOccupation<P3FreeRegion, 2>
        >;
        template <class TopoT, class CGT, class AddConsT, class EnvironmentT>
        void RegisterConstraint(const P3RegionRegionConnection & b, const TopoT & topo, const CGT & cg,
            AddConsT & addCons, const EnvironmentT & env){
            for (auto & a : b.anchors){
                //addCons({ {}, {}, {} })
            }
        }


        struct P3RegionRegionOverlap {            
            template <class Archive> void serialize(Archive & ar) {
                //ar(anchors);
            }
        };
        using P3RegionRegionOverlapConfig = core::ConstraintConfig<P3RegionRegionOverlap,
            core::ComponentOccupation<P3FreeRegion, Dynamic>
        >;
        template <class TopoT, class CGT, class AddConsT, class EnvironmentT>
        void RegisterConstraint(const P3RegionRegionOverlap & b, const TopoT & topo, const CGT & cg,
            AddConsT & addCons, const EnvironmentT & env){
            
        }


        struct P3LineLineIntersection {
            Vec3 anchor;
            template <class Archive> void serialize(Archive & ar) {
                ar(anchor);
            }
        };
        using P3LineLineIntersectionConfig = core::ConstraintConfig<P3LineLineIntersection,
            core::ComponentOccupation<P3FreeLine, 2>
        >;
        template <class TopoT, class CGT, class AddConsT, class EnvironmentT>
        void RegisterConstraint(const P3LineLineIntersection & b, const TopoT & topo, const CGT & cg,
            AddConsT & addCons, const EnvironmentT & env){

        }


        struct P3LineLineIncidence {
            Vec3 anchor;
            template <class Archive> void serialize(Archive & ar) {
                ar(anchor);
            }
        };
        using P3LineLineIncidenceConfig = core::ConstraintConfig<P3LineLineIncidence,
            core::ComponentOccupation<P3FreeLine, Dynamic>
        >;
        template <class TopoT, class CGT, class AddConsT, class EnvironmentT>
        void RegisterConstraint(const P3LineLineIncidence & b, const TopoT & topo, const CGT & cg, 
            AddConsT & addCons, const EnvironmentT & env){

        }


        struct P3LineLineParallelism {

        };


        struct P3LineLineOrthogonality {

        };



        using P3Graph = ConstraintGraph <
            std::tuple<P3FreeLine, P3FreeRegion>,
            std::tuple<
                P3LineLineIntersectionConfig,
                P3LineLineIncidenceConfig,
                P3RegionRegionConnectionConfig,
                P3RegionRegionOverlapConfig,
                P3RegionLineConnectionConfig
            >
        >;


        //P3Graph BuildP3Graph();

        //void EncodeMatrices();

    }
}
 
#endif