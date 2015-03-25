#ifndef PANORAMIX_MEX_GRAPH_OBJ_HPP
#define PANORAMIX_MEX_GRAPH_OBJ_HPP 

#include "../../src/misc/matlab.hpp"
#include "../../src/core/generic_topo.hpp"
#include "../../src/core/feature.hpp"
#include "../../src/gui/basic_types.hpp"

namespace mex {

    using namespace panoramix;
    using namespace panoramix::core;
    using namespace panoramix::misc;

    enum OcclusionType {
        LeftOcclusion,
        RightOcclusion,
        SameObject,
        ConnectedObjects
    };

    struct Region {
        Vec3 meanColor;
        double area;
        double perimeter;
        Box2 bbox;
    };
    struct Boundary {
        std::vector<std::vector<PixelLoc>> pixels;
        double strength; //
        double length;
        double smoothness;
        Vec2 orientation;
        double continuity[2];
        int longRange;

        OcclusionType ot;
    };

    using Graph = HomogeneousGraph02<Region, Boundary>;
    using RegionHandle = HandleOfTypeAtLevel<Graph, 0>;
    using BoundaryHandle = HandleOfTypeAtLevel<Graph, 1>;


        template <class To, class From>
    inline ImageOfType<To> im_static_cast(const ImageOfType<From> & im){
        ImageOfType<To> cim(im.size());
        for (auto it = im.begin(); it != im.end(); ++it){
            cim(it.pos()) = static_cast<To>(*it);
        }
        return cim;
    }

}

 
#endif