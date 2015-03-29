#include "../src/core/utilities.hpp"
#include "../src/gui/basic_types.hpp"
#include "../src/ml/factor_graph.hpp"
#include "config.hpp"

using namespace panoramix;
using namespace test;


TEST(FactorGraph, Simple){

    ml::FactorGraph fg;
    fg.varCategories = { ml::FactorGraph::VarCategory{ 2, 1.0 } };
    fg.factorCategories = { ml::FactorGraph::FactorCategory{
        [](const int * labels) -> double {
            return labels[0] == 0 ? 1.0 : 0.0;
        }, 1.0
    } };

    auto vh = fg.graph.add(0);
    auto fh = fg.graph.add<1>({ vh }, 0);

    auto results = fg.solve(50, 10, [](int epoch, double energy, double denergy){
        std::cout << "energy: " << energy << std::endl; 
        return denergy > - 1e-5; 
    });

    ASSERT(results[vh] == 1);

}


TEST(FactorGraph, Denoise){

    auto im = cv::imread(ProjectDataDirStrings::BPTests + "/horse.jpg");
    core::ResizeToMakeHeightUnder(im, 200);
    core::Imaged3 noised(im.size(), core::Vec3());
    for (auto it = noised.begin(); it != noised.end(); ++it){
        gui::Color color = gui::ColorFromImage(im, it.pos());
        core::Vec3 noise;
        for (double & n : noise.val){
            n = ((std::rand() % 1000) - 500) / 500.0;
        }
        *it = (core::Vec3)color + noise;
    }
    cv::imshow("original", im);
    cv::imshow("noised", noised);
    cv::waitKey();

    const core::Vec3 background = gui::Color(gui::ColorTag::White);
    const core::Vec3 foreground = gui::Color(gui::ColorTag::Black);

    ml::FactorGraph fg;
    fg.varCategories = { ml::FactorGraph::VarCategory{ 2, 1.0 } };
    fg.factorCategories.reserve(im.cols * im.rows + 2);
    fg.graph.internalElements<0>().reserve(im.cols * im.rows);
    fg.graph.internalElements<1>().reserve(im.cols * im.rows * 4);

    std::vector<ml::FactorGraph::VarHandle> vhs(im.cols * im.rows);
    
    // add varCategories and data costs
    for (auto it = noised.begin(); it != noised.end(); ++it){
        // add var node
        ml::FactorGraph::VarHandle vh = fg.graph.add(0);
        int id = core::EncodeSubscriptToIndex(it.pos(), noised.size());
        vhs[id] = vh;

        // append new factor type
        ml::FactorGraph::FactorCategory fd;
        double distToBackground = core::Distance(background, *it);
        double distToForeground = core::Distance(foreground, *it);
        fd.costs = [distToBackground, distToForeground](const int * labels) -> double {
            int label = labels[0];
            if (label == 0){ // judge as background
                return distToBackground;
            }
            else{
                return distToForeground;
            }
        };
        fd.c_alpha = 1.0;
        fg.factorCategories.push_back(std::move(fd));

        // add factor node
        fg.graph.add<1>({ vh }, fg.factorCategories.size() - 1);
    }

    // append smoothness factor types
    fg.factorCategories.push_back(ml::FactorGraph::FactorCategory{
       [](const int * labels) -> double {
            return labels[0] == labels[1] ? 0.0 : 0.5;
        }, 1.0
    });
    int smoothnessFid1 = fg.factorCategories.size() - 1;
    fg.factorCategories.push_back(ml::FactorGraph::FactorCategory{
        [](const int * labels) -> double {
            return labels[0] == labels[1] ? 0.0 : 0.3;
        }, 1.0
    });
    int smoothnessFid2 = fg.factorCategories.size() - 1;

    // add smoothness factor nodes
    for (int i = 0; i < noised.rows - 1; i++){
        for (int j = 0; j < noised.cols - 1; j++){
            auto vh = vhs[core::EncodeSubscriptToIndex(core::PixelLoc(j, i), noised.size())];
            fg.graph.add<1>({ vh, vhs[core::EncodeSubscriptToIndex(core::PixelLoc(j + 1, i), noised.size())] }, smoothnessFid1);
            fg.graph.add<1>({ vh, vhs[core::EncodeSubscriptToIndex(core::PixelLoc(j, i + 1), noised.size())] }, smoothnessFid1);
            fg.graph.add<1>({ vh, vhs[core::EncodeSubscriptToIndex(core::PixelLoc(j + 1, i + 1), noised.size())] }, smoothnessFid2);
            if (i > 0){
                fg.graph.add<1>({ vh, vhs[core::EncodeSubscriptToIndex(core::PixelLoc(j + 1, i - 1), noised.size())] }, smoothnessFid2);
            }
        }
    }

    auto results = fg.solve(100, 3, [](int epoch, double e, double de){
        std::cout << "#" << epoch << "  energy: " << e << std::endl;
        return true;
    });
    core::Imaged3 recovered(noised.size());
    for (auto it = recovered.begin(); it != recovered.end(); ++it){
        auto vh = vhs[core::EncodeSubscriptToIndex(it.pos(), noised.size())];
        int label = results[vh];
        *it = label == 0 ? background : foreground;
    }

    cv::imshow("recovered", recovered);
    cv::waitKey();

}