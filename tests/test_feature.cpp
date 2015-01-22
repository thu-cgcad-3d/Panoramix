#include "../src/core/cameras.hpp"
#include "../src/core/utilities.hpp"
#include "../src/core/feature.hpp"
#include "../src/core/view.hpp"
#include "../src/vis/visualize2d.hpp"

#include "config.hpp"

using namespace panoramix;
using namespace test;

TEST(Feature, SegmentationExtractor) {
    {
        core::SegmentationExtractor::Params p;
        p.c = 5;
        p.minSize = 400;
        p.sigma = 1;
        core::SegmentationExtractor seg(p);
        core::Image im = cv::imread(ProjectDataDirStrings::Normal + "/75.jpg");
        vis::Visualizer2D(im) << vis::manip2d::Show();
        auto segs = seg(im);
        vis::Visualizer2D(segs.first)
            << vis::manip2d::SetColorTable(vis::CreateRandomColorTableWithSize(segs.second))
            << vis::manip2d::Show();
    }
    {
        core::SegmentationExtractor::Params p;
        p.c = 5;
        p.minSize = 400;
        p.sigma = 1;
        core::SegmentationExtractor seg(p);
        core::Image im = cv::imread(ProjectDataDirStrings::Normal + "/75.jpg");
        auto segs = seg(im, { core::Line2({ 0.0, 0.0 }, core::Point2(im.cols, im.rows)), core::Line2(core::Point2(im.cols, 0), core::Point2(0, im.rows)) });
        vis::Visualizer2D(segs.first)
            << vis::manip2d::SetColorTable(vis::CreateRandomColorTableWithSize(segs.second))
            << vis::manip2d::Show();
    }
    {
        core::SegmentationExtractor::Params p;
        p.algorithm = core::SegmentationExtractor::SLIC;
        p.superpixelSizeSuggestion = 3000;
        core::SegmentationExtractor seg(p);
        core::Image im = cv::imread(ProjectDataDirStrings::Normal + "/75.jpg");
        vis::Visualizer2D(im) << vis::manip2d::Show();
        auto segs = seg(im);
        vis::Visualizer2D(segs.first)
            << vis::manip2d::SetColorTable(vis::CreateRandomColorTableWithSize(segs.second))
            << vis::manip2d::Show();
    }
    {
        core::SegmentationExtractor::Params p;
        p.algorithm = core::SegmentationExtractor::QuickShiftCPU;
        core::SegmentationExtractor seg(p);
        core::Image im = cv::imread(ProjectDataDirStrings::Normal + "/75.jpg");
        vis::Visualizer2D(im) << vis::manip2d::Show();
        auto segs = seg(im);
        vis::Visualizer2D(segs.first)
            << vis::manip2d::SetColorTable(vis::CreateRandomColorTableWithSize(segs.second))
            << vis::manip2d::Show();
    }
}

TEST(Feature, LineSegmentExtractor) {
    core::LineSegmentExtractor lineseg;
    core::Image im = cv::imread(ProjectDataDirStrings::LocalManhattan + "/buildings2.jpg");
    //vis::Visualizer2D(im) 
    //    << vis::manip2d::SetColor(vis::ColorTag::Yellow) 
    //    << vis::manip2d::SetThickness(2) << 
    //    lineseg(im) 
    //    << vis::manip2d::Show();

    core::LineSegmentExtractor::Params params;
    params.algorithm = core::LineSegmentExtractor::LSD;
    core::LineSegmentExtractor lineseg2(params);
    vis::Visualizer2D(im)
        << vis::manip2d::SetColor(vis::ColorTag::Yellow)
        << vis::manip2d::SetThickness(2) <<
        lineseg2(im)
        << vis::manip2d::Show(0);
}

TEST(Feature, VanishingPointsDetector) {
    core::LineSegmentExtractor::Params lsParams;
    lsParams.minLength = 20;
    lsParams.xBorderWidth = lsParams.yBorderWidth = 20;
    core::LineSegmentExtractor lineseg(lsParams);
    core::VanishingPointsDetector vpdetector;

    core::Image im = cv::imread(ProjectDataDirStrings::Normal + "/room.png");
    core::ResizeToMakeWidthUnder(im, 1000);
    auto lines = lineseg(im);
    std::vector<int> lineClasses;
    std::array<core::HPoint2, 3> vps;
    double focalLength;
    std::tie(vps, focalLength, lineClasses) = vpdetector(lines, core::Point2(im.cols/2, im.rows/2));
    std::vector<core::Classified<core::Line2>> classifiedLines;
    for (int i = 0; i < lines.size(); i++) {
        classifiedLines.push_back({ lineClasses[i], lines[i]});
    }
    vis::Visualizer2D(im)
        << vis::manip2d::SetColorTable(vis::ColorTableDescriptor::RGB)
        << vis::manip2d::SetThickness(2) <<
        classifiedLines
        << vis::manip2d::Show(0);
}


TEST(Feature, LocalManhattanVanishingPointDetector) {

    using namespace core;

    // forged experiment for panorama
    core::Image im = cv::imread(ProjectDataDirStrings::PanoramaIndoor + "/14.jpg");
    core::ResizeToMakeWidthUnder(im, 2000);
    //vis::Visualizer2D(im) << vis::manip2d::Show();

    core::PanoramicCamera ocam(im.cols / M_PI / 2.0);
    core::PerspectiveCamera cam(800, 800, ocam.focal(), { 0, 0, 0 }, { -2, 0, -0.5 }, { 0, 0, -1 });

    auto pim = core::MakeCameraSampler(cam, ocam)(im);
    //vis::Visualizer2D(pim) << vis::manip2d::Show();

    core::LineSegmentExtractor lineseg;
    lineseg.params().minLength = 5;
    lineseg.params().algorithm = core::LineSegmentExtractor::LSD;
    std::vector<Line2> line2s = lineseg(pim);

    Vec3 vp1 = { 0, 0, 1 };
    std::vector<Classified<Line3>> line3s(line2s.size());
    std::vector<Vec3> line3norms(line2s.size());
    for (int i = 0; i < line2s.size(); i++) {
        line3s[i].component.first = normalize(cam.spatialDirection(line2s[i].first));
        line3s[i].component.second = normalize(cam.spatialDirection(line2s[i].second));
        line3norms[i] = line3s[i].component.first.cross(line3s[i].component.second);
        line3s[i].claz = abs(line3norms[i].dot(vp1)) < 0.006 ? 0 : -1;
    }

    std::vector<std::pair<int, int>> pairs;
    for (int i = 0; i < line2s.size(); i++){
        if (line3s[i].claz == 0)
            continue;
        if (abs(line3norms[i].dot(vp1)) < 0.01)
            continue;
        for (int j = i + 1; j < line2s.size(); j++){
            if (line3s[i].claz == 0)
                continue;
            if (abs(line3norms[j].dot(vp1)) < 0.01)
                continue;
            double dist = DistanceBetweenTwoLines(line2s[i], line2s[j]).first;
            auto & n1 = line3norms[i];
            auto & n2 = line3norms[j];
            auto inter = n1.cross(n2);
            auto interp = cam.screenProjection(inter);
            double dd = 40;
            if (dist < dd &&
                DistanceFromPointToLine(interp, line2s[i]).first < dd &&
                DistanceFromPointToLine(interp, line2s[j]).first < dd){
                pairs.emplace_back(i, j);
            }
        }
    }

    std::vector<std::pair<int, int>> orthoPairs;
    for (auto & p : pairs){
        auto & n1 = line3norms[p.first];
        auto & n2 = line3norms[p.second];
        auto p1 = normalize(n1.cross(vp1));
        auto p2 = normalize(n2.cross(vp1));
        if (abs(p1.dot(p2)) < 0.02)
            orthoPairs.push_back(p);
    }

    vis::Visualizer2D viz(pim);
    viz = viz << vis::manip2d::SetThickness(2);
    for (int i = 0; i < line2s.size(); i++){
        if (line3s[i].claz == 0){
            viz << vis::manip2d::SetColor(vis::ColorTag::Red) << line2s[i];
        }
        else {
            //viz << vis::manip2d::SetColor(vis::ColorTag::Black) << line2s[i];
        }
    }
    for (auto & op : orthoPairs){
        auto & n1 = line3norms[op.first];
        auto & n2 = line3norms[op.second];
        auto inter = n1.cross(n2);
        auto interp = cam.screenProjection(inter);
        viz << vis::manip2d::SetColor(vis::ColorTag::LightGray) << vis::manip2d::SetThickness(1)
            << Line2(line2s[op.first].center(), interp)
            << Line2(line2s[op.second].center(), interp);
        viz << vis::manip2d::SetColor(vis::ColorTag::White) << vis::manip2d::SetThickness(2) 
            << line2s[op.first] << line2s[op.second];
    }

    viz << vis::manip2d::Show();

    return;


    //core::LineSegmentExtractor lineseg;
    //lineseg.params().useLSD = true;
    //core::LocalManhattanVanishingPointsDetector::Params params;

    //std::vector<core::Image> ims = {
    //    cv::imread(ProjectDataDirStrings::LocalManhattan + "/cuboids.png")//,
    //    //cv::imread(ProjectDataDirStrings::LocalManhattan + "/buildings3.jpg"),
    //    //cv::imread(ProjectDataDirStrings::LocalManhattan + "/buildings2.jpg")
    //};

    //for (auto & im : ims) {
    //    core::ResizeToMakeWidthUnder(im, 1000);
    //    params.image = im;

    //    core::LocalManhattanVanishingPointsDetector vpdetector(params);

    //    std::vector<core::Line2> lines;
    //    core::LocalManhattanVanishingPointsDetector::Result result;

    //    lines = lineseg(im);
    //    //core::LoadFromDisk(ProjectDataDirStrings::Serialization + "/temp.state", lines);
    //    result = vpdetector(lines, core::Point2(im.cols / 2, im.rows / 2));

    //    for (auto & op : result.horizontalVanishingPointIds){
    //        std::cout << "[ " << op.first << ", " << op.second << "]" << "---"
    //            << result.vanishingPoints[op.first].value()
    //            << result.vanishingPoints[op.second].value() << std::endl;
    //    }

    //    std::vector<int> vpPairIds(result.vanishingPoints.size(), -1);
    //    for (int i = 0; i < result.horizontalVanishingPointIds.size(); i++){
    //        vpPairIds[result.horizontalVanishingPointIds[i].first] = i;
    //        vpPairIds[result.horizontalVanishingPointIds[i].second] = i;
    //    }

    //    auto ctable = vis::CreateRandomColorTableWithSize(result.vanishingPoints.size());
    //    ctable.exceptoinalColor() = vis::ColorFromTag(vis::ColorTag::White);

    //    auto viz = vis::Visualizer2D(im)
    //        << vis::manip2d::SetColorTable(ctable)
    //        << vis::manip2d::SetThickness(4);

    //    for (int i = 0; i < lines.size(); i++) {
    //        int claz = result.lineClasses[i];
    //        viz << vis::manip2d::SetThickness(claz == -1 ? 1 : vpPairIds[claz] + 3)
    //            << vis::manip2d::SetColor(ctable[claz])
    //            << core::NoteAs(lines[i], std::to_string(i) + "." + std::to_string(result.lineClasses[i]));
    //    }

    //    viz << vis::manip2d::SetColor(vis::ColorTag::Red)
    //        << result.vanishingPoints;

    //    viz /*<< vis::manip2d::SetColor(vis::ColorTag::Black)
    //        << result.hlineCands
    //        << vis::manip2d::SetColor(vis::ColorTag::Red)
    //        << result.horizon*/
    //        << vis::manip2d::Show(0);
    //}
}



TEST(Feature, FeatureExtractor) {
    core::SegmentationExtractor segmenter;
    core::LineSegmentExtractor::Params params;
    params.algorithm = core::LineSegmentExtractor::LSD;
    core::LineSegmentExtractor lineSegmentExtractor(params);
    
    for (int i = 0; i < 4; i++) {
        std::string name = ProjectDataDirStrings::Normal + "/" + "sampled_" + std::to_string(i) + ".png";
        cv::Mat im = cv::imread(name);
        auto segs = segmenter(im);
        vis::Visualizer2D(im) 
            << [](vis::Visualizer2D & viz) { viz.params.winName = "haha"; }
        << vis::manip2d::SetColorTable(vis::CreateRandomColorTableWithSize(segs.second))
            << segs.first
            << lineSegmentExtractor(im) 
            << vis::manip2d::Show();
    }
}

TEST(Feature, GeometricContextMatlab) {
    core::GeometricContextEstimator gcEstimator;
    gcEstimator.params().useMatlab = true;

    core::Image im = cv::imread(ProjectDataDirStrings::LocalManhattan + "/buildings2.jpg");
    auto gc = gcEstimator(im, core::SceneClass::Indoor);
}

TEST(Feature, GeometricContextPanorama) {
    //core::GeometricContextEstimator gcEstimator;
    //gcEstimator.params().useMatlab = true;

    //core::Image pan = cv::imread(ProjectDataDirStrings::PanoramaIndoor + "/13.jpg");
    //core::ResizeToMakeHeightUnder(pan, 500);
    //auto panview = core::CreatePanoramicView(pan);

    //auto pviews = core::SampleHorizontalPerspectiveViews(pan, 10);
    //std::vector<core::ViewX<core::PerspectiveCamera>> pviewGCs(pviews.size());
    ////for (int i = 0; i < pviews.size(); i++){
    ////    pviewGCs[i].camera = pviews[i].camera;
    ////    auto gc = gcEstimator(pviews[i].image);
    ////    pviewGCs[i].image.resize(gc.size());
    ////    for (int k = 0; k < gc.size(); k++)
    ////        pviewGCs[i].image[k] = gc[k];
    ////}
    //core::LoadFromDisk("./cache/Feature.GeometricContextPanorama.pviewGCs", pviewGCs);

    //auto gcs = core::MergePerspectiveViewsBack(pviewGCs, panview.camera).image;
    //core::Image wholeGC;
    //cv::merge(gcs, wholeGC);
    //vis::Visualizer2D(wholeGC) << vis::manip2d::Show();

}