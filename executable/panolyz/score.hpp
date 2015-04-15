#pragma once

#include "routines.hpp"
#include "../../src/experimental/rl_graph.hpp"
#include "../../src/experimental/rl_graph_occlusion.hpp"
#include "../../src/experimental/tools.hpp"


namespace panolyz {

    using namespace panoramix;
    using namespace core;
    using namespace experimental;


    struct Cache {
        RLGraphVars vars;
        InstanceTable<RegionData> planes;
        InstanceTable<LineData> lines;
    };

    struct Features {
        HandledTable<RegionHandle, Vec<double, 5>> gc;
        HandledTable<RegionBoundaryHandle, int> occ;
    };

    struct Cost {
        double manhattan;
        double connection;
        double occlusion;
        double nanPunishment;
        double lineIsolationPunishment;
        double gcViolation;
        double junctionPunishment;
    };

    double RegionManhattanCost(const RLGraph & mg, const RLGraphControls & controls, 
        const Cache & cache, const Features & features);

    double LineManhattanCost(const RLGraph & mg, const RLGraphControls & controls, 
        const Cache & cache, const Features & features);

    double RegionBoundaryConnectionCost(const RLGraph & mg, const RLGraphControls & controls,
        const Cache & cache, const Features & features);

    double RegionLineConnectionCost(const RLGraph & mg, const RLGraphControls & controls,
        const Cache & cache, const Features & features);

    double LineConnectionCost(const RLGraph & mg, const RLGraphControls & controls,
        const Cache & cache, const Features & features);







}