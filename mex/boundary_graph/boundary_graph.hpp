#ifndef PANORAMIX_MEX_BOUNDARY_GRAPH_HPP
#define PANORAMIX_MEX_BOUNDARY_GRAPH_HPP 

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
        std::vector<std::vector<PixelLoc>> contours;

        Vec3 meanColor;
        Imaged colorHist;
        double area;
        Point2 center;
        double perimeter;
        Vec<double, 7> meanGC;

        Plane3 fittedPlane;
    };

    struct Boundary {
        std::vector<std::vector<PixelLoc>> pixels;

        double strength; // pb?
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

    using Imaged7 = ImageOfType<Vec<double, 7>>;

    struct BoundaryGraph {
        Graph graph;

        Image image;
        Imagei segmentedImage;
        int segmentsNum;
        Imaged7 geometricContext;
        Imaged depth;

        BoundaryGraph(){}
        BoundaryGraph(const Image & im, const Imaged7 & gc, const Imaged & d = Imaged());

        void print() const;

        Imaged computeFeaturesMat() const;
        Imaged computeLabelsVec() const;
        void installLabelsVec(const Imaged & labels);
    };

}

 
#endif