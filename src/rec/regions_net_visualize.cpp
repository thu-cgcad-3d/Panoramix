#include "regions_net_visualize.hpp"

namespace panoramix {
    namespace rec {

        using namespace core;

        Visualizer2D operator << (Visualizer2D viz, const RegionsNet & net) {
            viz.setImage(net.image());
            
            int height = net.image().rows;
            int width = net.image().cols;
            
            Image coloredOutput(net.segmentedRegions().size(), CV_8UC3);
            std::vector<cv::Vec<uint8_t, 3>> colors(net.regions().internalElements<0>().size());
            std::generate(colors.begin(), colors.end(), [](){
                return cv::Vec<uint8_t, 3>(uint8_t(std::rand() % 256),
                    uint8_t(std::rand() % 256),
                    uint8_t(std::rand() % 256));
            });
            for (int y = 0; y < height; y++) {
                for (int x = 0; x < width; x++) {
                    coloredOutput.at<cv::Vec<uint8_t, 3>>(cv::Point(x, y)) =
                        colors[net.segmentedRegions().at<int32_t>(cv::Point(x, y))];
                }
            }

            viz.params.alphaForNewImage = 0.3f;
            viz << coloredOutput;

            viz.params.thickness = 1;
            viz.params.color = ColorFromTag(ColorTag::Black);
            for (auto & h : net.regions().elements<1>()){
                auto fromCenter = net.regions().data(h.topo.lowers[0]).center;
                auto toCenter = net.regions().data(h.topo.lowers[1]).center;
                Line2 line;
                line.first = fromCenter;
                line.second = toCenter;
                viz << line;
            }

            return viz;
        }

    }
}