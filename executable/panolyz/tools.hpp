#pragma once

#include "routines.hpp"
#include "../../src/experimental/rl_graph.hpp"

namespace panolyz {

    using namespace panoramix;
    using namespace core;
    using namespace experimental;

    struct Labels{
        std::vector<int> regionLabels;
        std::map<std::pair<int, int>, int> boundaryLabels;
        template <class Archive>
        inline void serialize(Archive & ar) { ar(regionLabels, boundaryLabels); }
    };

    void LabelIt(Labels & labels, 
        const Image & im, 
        const Imagei & segments, 
        const std::map<std::pair<int, int>, std::vector<std::vector<PixelLoc>>> & boundaryPixels,
        const std::vector<std::string> & regionLabelNames, 
        const std::vector<std::string> & boundaryLabelNames,
        const std::vector<gui::Color> & regionLabelColors,
        const std::vector<gui::Color> & boundaryLabelColors,
        bool doModal = true);


}