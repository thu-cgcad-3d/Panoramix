#pragma once

#include "../../src/misc/cache.hpp"

#include "../../src/gui/singleton.hpp"
#include "../../src/gui/utility.hpp"
#include "../../src/gui/canvas.hpp"

#include "../../src/experimental/pi_graph_annotation.hpp"
#include "../../src/experimental/pi_graph_solve.hpp"
#include "../../src/experimental/pi_graph_vis.hpp"

namespace panolyz {
    using namespace pano;
    using namespace pano::core;
    using namespace pano::experimental;

    struct ReconstructedModel {
        virtual double depthAt(const Vec3 & direction) const = 0;
        virtual bool isValid(const Vec3 & direction) const { return true; }
        virtual void visualize(const std::vector<Vec3> & directions = std::vector<Vec3>()) const = 0;
    };

    std::vector<Vec3> FibonacciDirections(int N);
    double AverageVariation(const ReconstructedModel & m1, const ReconstructedModel & m2, 
        const std::vector<Vec3> & directions);


    std::unique_ptr<ReconstructedModel> PredictionOfGT(const std::string & impath);

    struct PredictOptions {
        bool useGroundTruthOcclusions;
    };
    std::unique_ptr<ReconstructedModel> PredictionOfPanoramix(const std::string & impath, 
        const PredictOptions & options, misc::Matlab & matlab);
    
    std::unique_ptr<ReconstructedModel> PredictionOfPanoContext(const std::string & impath,
        misc::Matlab & matlab);

}
