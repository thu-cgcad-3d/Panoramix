#include "../src/core/utility.hpp"
#include "../src/gui/basic_types.hpp"
#include "../src/ml/factor_graph.hpp"
#include "config.hpp"

using namespace pano;
using namespace test;


TEST(FactorGraph, Simple){

    ml::FactorGraph fg;
    auto vcid = fg.addVarCategory(ml::FactorGraph::VarCategory{ 2, 1.0 });
    auto fcid = fg.addFactorCategory(ml::FactorGraph::FactorCategory{
        [](const int * labels, size_t nvars, ml::FactorGraph::FactorCategoryId fcid, void *) -> double {
            return labels[0] == 0 ? 1.0 : 0.0;
        }, 1.0
    });

    auto vh = fg.addVar(vcid);
    auto fh = fg.addFactor({ vh }, fcid);

    auto results = fg.solveWithSimpleCallback(50, 10, [](int epoch, double energy){
        std::cout << "energy: " << energy << std::endl; 
        return true; 
    });

    ASSERT(results[vh] == 1);

}


TEST(FactorGraph, Denoise){

    auto im = core::ImageRead(ProjectDataDirStrings::BPTests + "/horse.jpg");
    core::ResizeToMakeHeightUnder(im, 200);
    core::Image3d noised(im.size(), core::Vec3());
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
        fd.costs = [distToBackground, distToForeground](const int * labels, size_t nvars, 
            ml::FactorGraph::FactorCategoryId fcid, void *) -> double {
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
        [](const int * labels, size_t nvars, ml::FactorGraph::FactorCategoryId fcid, void *) -> double {
            return labels[0] == labels[1] ? 0.0 : 5;
        }, 1.0
    });
    auto smoothnessfcid2 = fg.addFactorCategory(ml::FactorGraph::FactorCategory{
        [](const int * labels, size_t nvars, ml::FactorGraph::FactorCategoryId fcid, void *) -> double {
            return labels[0] == labels[1] ? 0.0 : 3;
        }, 1.0
    });

    auto smoothnessfcid3 = fg.addFactorCategory(ml::FactorGraph::FactorCategory{
        [](const int * labels, size_t nvars, ml::FactorGraph::FactorCategoryId fcid, void *) -> double {
            if (labels[4] == 0 && std::accumulate(labels, labels + nvars, 0) == 8)
                return 1.0; // return std::numeric_limits<double>::infinity();
            if (labels[4] == 1 && std::accumulate(labels, labels + nvars, 0) == 1)
                return 1.0; // return std::numeric_limits<double>::infinity();
            return 0.0;
        }, 1.0
    });

    // add smoothness factor nodes
    for (int i = 0; i < noised.rows - 1; i++){
        for (int j = 0; j < noised.cols - 1; j++){
            auto vh = vhs[core::EncodeSubscriptToIndex(core::Pixel(j, i), noised.size())];
            fg.addFactor({ vh, vhs[core::EncodeSubscriptToIndex(core::Pixel(j + 1, i), noised.size())] }, smoothnessfcid1);
            fg.addFactor({ vh, vhs[core::EncodeSubscriptToIndex(core::Pixel(j, i + 1), noised.size())] }, smoothnessfcid1);
            fg.addFactor({ vh, vhs[core::EncodeSubscriptToIndex(core::Pixel(j + 1, i + 1), noised.size())] }, smoothnessfcid2);
            if (i > 0){
                fg.addFactor({ vh, vhs[core::EncodeSubscriptToIndex(core::Pixel(j + 1, i - 1), noised.size())] }, smoothnessfcid2);
            }
          /*std::vector<ml::FactorGraph::VarHandle> localVhs;
            for (int x = -1; x <= 1; x++) {
                for (int y = -1; y <= 1; y++) {
                    auto p = core::Pixel(j + x, i + y);
                    if (!core::Contains(noised, p))
                        continue;
                    localVhs.push_back(vhs[core::EncodeSubscriptToIndex(p, noised.size())]);
                }
            }
            fg.addFactor(localVhs.begin(), localVhs.end(), smoothnessfcid3);*/
        }
    }

    auto results = fg.solveWithSimpleCallback(100, 3, [](int epoch, double e){
        std::cout << "#" << epoch << "  energy: " << e << std::endl;
        return true;
    });
    core::Image3d recovered(noised.size());
    for (auto it = recovered.begin(); it != recovered.end(); ++it){
        auto vh = vhs[core::EncodeSubscriptToIndex(it.pos(), noised.size())];
        int label = results[vh];
        *it = label == 0 ? background : foreground;
    }

    cv::imshow("recovered", recovered);
    cv::waitKey();

}