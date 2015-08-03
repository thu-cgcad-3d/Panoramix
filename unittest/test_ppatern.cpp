#include "../src/core/cameras.hpp"
#include "../src/core/utility.hpp"
#include "../src/core/feature.hpp"
#include "../src/core/containers.hpp"
#include "../src/core/clock.hpp"
#include "../src/gui/canvas.hpp"
#include "../src/gui/utility.hpp"

#include "config.hpp"

using namespace pano;
using namespace test;

std::string ppatternFolder = "F://PPattern//image";

void ShowDir(const core::Image3ub & im, std::initializer_list<double> angles){
    gui::ColorTable ctable = gui::ColorTableDescriptor::RGB;
    ctable.appendRandomizedGreyColors(angles.size());
    int id = 0;
    auto canvas = gui::MakeCanvas(im);
    for (auto angle : angles){
        core::Vec2 maxVec = core::Direction<double>(angle);
        maxVec = core::normalize(maxVec) * 0.25 * sqrt(core::Area(im));
        core::Line2 visline(core::Center(im) - maxVec, core::Center(im) + maxVec);
        canvas.color(ctable[id++]).thickness(10).add(visline);
    }
    canvas.show();
}

//core::Vec2 EstimateGeneralGrad(const core::Imaged & dense){
//
//    using namespace core;
//    int winhsize = 50;
//    for (int x = winhsize; x < dense.cols - winhsize; x++){
//        for (int y = winhsize; y < dense.rows - winhsize; y++){
//
//        }
//    }
//
//}


TEST(PPattern, A){

    core::Image3ub image = gui::PickAnImage(ppatternFolder);
    core::ClipToSquare(image);
    assert(image.cols == image.rows);

    core::ResizeToHeight(image, 800);
    auto mask = core::Rotate(image, core::DegreesToRadians(0));

    auto kps = core::SIFT(image);

    gui::MakeCanvas(image).add(kps).show();

    core::RTreeWrapper<cv::KeyPoint, core::DefaultInfluenceBoxFunctor<double>> featureRtree(kps.begin(), kps.end(), 
        core::DefaultInfluenceBoxFunctor<double>(0.5));

    int winhsize = 100;
    int step = 3;
    core::Imaged dense(core::Sizei((image.cols - 2 * winhsize) / step + 1, (image.rows - 2 * winhsize) / step + 1), -1);
    core::RTree<core::Box2, core::Point3> denseSamples;
    for (int x = winhsize; x < image.cols - winhsize; x += step){
        for (int y = winhsize; y < image.rows - winhsize; y += step){
            if (!mask(y, x)){
                continue;
            }
            core::Box2 box = core::BoundingBox(core::Pixel(x, y)).expand(winhsize);
            int count = featureRtree.count(box);
            int i = (y - winhsize) / step;
            int j = (x - winhsize) / step;
            dense(i, j) = count;
            denseSamples.insert(core::BoundingBox(core::Pixel(x, y)).expand(1), core::Point3(x, y, count));
        }
    }
    double maxDense = core::MinMaxValOfImage(dense).second;
    gui::AsCanvas<double>(dense / maxDense).show();

    double scale = core::norm(core::Point2(image.cols, image.rows));
    std::vector<double> thetas;
    for (double theta = 0; theta < M_PI * 2; theta += M_PI / 10.0){
        thetas.push_back(theta);
    }
    std::vector<std::vector<double>> votes(thetas.size(), std::vector<double>(1000, 0.0));
    std::vector<double> entropies(thetas.size(), 0.0);
    std::vector<double> grads(thetas.size(), 0.0);
    double validRadius = std::min(image.cols, image.rows) / 2.0;
    for (int i = 0; i < votes.size(); i++){
        auto & vtable = votes[i];
        auto theta = thetas[i];
        auto thetaDir = core::Direction<double>(theta);
        std::vector<int> voteNums(vtable.size(), 0);
        for (auto it = dense.begin(); it != dense.end(); ++it){
            auto p = core::ecast<double>(it.pos());
            auto dir = p - core::Center(image);
            double angle = core::AngleBetweenDirections(thetaDir, dir);
            double projPos = cos(angle) * core::norm(dir);
            if (abs(projPos >= validRadius))
                continue;
            int id = core::BoundBetween(int((projPos / validRadius / 2.0 + 0.5) * (vtable.size()-1)), 0, vtable.size());
            vtable[id] += *it;
            voteNums[id] ++;
        }
        for (int j = 0; j < vtable.size(); j++){
            vtable[j] /= std::max(voteNums[j], 1);
        }
        double entropy = core::EntropyOfContainer(vtable);
        entropies[i] = entropy;
        grads[i] = abs(std::accumulate(vtable.begin(), vtable.begin() + 500, 0)
            - std::accumulate(vtable.begin() + 500, vtable.begin() + 1000, 0));
        //grads[i] /= std::accumulate(vtable.begin(), vtable.begin() + 1000, 0);
    }

    int maxEntropyId = std::max_element(entropies.begin(), entropies.end()) - entropies.begin();
    int minEntropyId = std::min_element(entropies.begin(), entropies.end()) - entropies.begin();

    int maxgrad = std::min_element(grads.begin(), grads.end()) - grads.begin();

    double thetaA = thetas[maxEntropyId];
    double thetaB = thetas[minEntropyId];

    ShowDir(image, { thetaA + M_PI_2});

    //// grad directions
    //int gradsize = 320 / step;

    //// collect geographical contour directions
    //// only collect directions sampled in the central region
    //std::vector<double> contourDirs;
    //{
    //    core::SetClock();
    //    for (int x = gradsize; x < dense.cols - gradsize; x += 5){
    //        for (int y = gradsize; y < dense.rows - gradsize; y += 5){

    //            core::Pixel pos(x * step + winhsize, y * step + winhsize);
    //            if (!mask(pos)){
    //                continue;
    //            }

    //            if (core::Distance(core::ecast<double>(pos), core::Center(image)) > 400)
    //                continue;

    //            for (int xx = x - gradsize; xx <= x + gradsize; xx++){
    //                for (int yy = y - gradsize; yy <= y + gradsize; yy++){

    //                    core::Pixel pos(xx * step + winhsize, yy * step + winhsize);
    //                    if (!mask(pos)){
    //                        continue;
    //                    }

    //                    if (core::Distance(core::ecast<double>(pos), core::Center(image)) > 400)
    //                        continue;

    //                    if (xx == x && yy == y) {
    //                        continue;
    //                    }
    //                    int count1 = dense(y, x);
    //                    int count2 = dense(yy, xx);
    //                    if (abs(count1 - count2) > 0.2 * maxDense)
    //                        continue;
    //                    contourDirs.push_back(core::Angle(core::Vec2(xx - x, yy - y)));
    //                    if (contourDirs.back() < 0){
    //                        contourDirs.back() += M_PI;
    //                    }
    //                }
    //            }
    //        }
    //    }
    //}

    //// vote
    //int vsize = 100;
    //std::vector<double> votes(vsize, 0.0);
    //{
    //    core::SetClock();
    //    for (double d : contourDirs){
    //        for (int i = 0; i < votes.size(); i++){
    //            double p = i / double(vsize) * M_PI;
    //            double dist = std::abs(d - p);
    //            if (dist > M_PI_2){
    //                dist = M_PI - dist;
    //            }
    //            double score = core::Gaussian(dist, 0.1);
    //            votes[i] += score;
    //        }
    //    }
    //}

    //double contourDirAngle = (std::max_element(votes.begin(), votes.end()) - votes.begin()) /  double(vsize) * M_PI;
    //core::Vec2 maxVec = core::Direction<double>(contourDirAngle);
    //maxVec = core::normalize(maxVec) * 0.25 * sqrt(core::Area(dense));

    //core::Line2 visline(core::Point2(image.cols, image.rows) / 2 - maxVec, core::Point2(image.cols, image.rows) / 2 + maxVec);

    //gui::MakeCanvas(image).color(gui::Red).thickness(10).add(visline).maxWidth(800).show();

}