#include "../src/core/utilities.hpp"
#include "../src/gui/basic_types.hpp"
#include "../src/ml/factor_graph.hpp"
#include "config.hpp"

using namespace panoramix;
using namespace test;


TEST(FactorGraph, Simple){

    ml::FactorGraph fg;
    auto vcid = fg.addVarCategory(ml::FactorGraph::VarCategory{ 2, 1.0 });
    auto fcid = fg.addFactorCategory(ml::FactorGraph::FactorCategory{
        [](const int * labels) -> double {
            return labels[0] == 0 ? 1.0 : 0.0;
        }, 1.0
    });

    auto vh = fg.addVar(vcid);
    auto fh = fg.addFactor({ vh }, fcid);

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

    const core::Vec3 background = gui::Color(gui::White);
    const core::Vec3 foreground = gui::Color(gui::Black);

    ml::FactorGraph fg;
    auto vcid = fg.addVarCategory(ml::FactorGraph::VarCategory{ 2, 1.0 });

    fg.reserveFactorCategories(im.cols * im.rows + 2);
    fg.reserveVars(im.cols * im.rows);
    fg.reserveFactors(im.cols * im.rows * 4);

    std::vector<ml::FactorGraph::VarHandle> vhs(im.cols * im.rows);
    
    // add varCategories and data costs
    for (auto it = noised.begin(); it != noised.end(); ++it){
        // add var node
        ml::FactorGraph::VarHandle vh = fg.addVar(vcid);
        vhs[core::EncodeSubscriptToIndex(it.pos(), noised.size())] = vh;

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
        auto fcid = fg.addFactorCategory(std::move(fd));

        // add factor node
        fg.addFactor({ vh }, fcid);
    }

    // append smoothness factor types
    auto smoothnessfcid1 = fg.addFactorCategory(ml::FactorGraph::FactorCategory{
       [](const int * labels) -> double {
            return labels[0] == labels[1] ? 0.0 : 0.5;
        }, 1.0
    });
    auto smoothnessfcid2 = fg.addFactorCategory(ml::FactorGraph::FactorCategory{
        [](const int * labels) -> double {
            return labels[0] == labels[1] ? 0.0 : 0.3;
        }, 1.0
    });

    // add smoothness factor nodes
    for (int i = 0; i < noised.rows - 1; i++){
        for (int j = 0; j < noised.cols - 1; j++){
            auto vh = vhs[core::EncodeSubscriptToIndex(core::PixelLoc(j, i), noised.size())];
            fg.addFactor({ vh, vhs[core::EncodeSubscriptToIndex(core::PixelLoc(j + 1, i), noised.size())] }, smoothnessfcid1);
            fg.addFactor({ vh, vhs[core::EncodeSubscriptToIndex(core::PixelLoc(j, i + 1), noised.size())] }, smoothnessfcid1);
            fg.addFactor({ vh, vhs[core::EncodeSubscriptToIndex(core::PixelLoc(j + 1, i + 1), noised.size())] }, smoothnessfcid2);
            if (i > 0){
                fg.addFactor({ vh, vhs[core::EncodeSubscriptToIndex(core::PixelLoc(j + 1, i - 1), noised.size())] }, smoothnessfcid2);
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