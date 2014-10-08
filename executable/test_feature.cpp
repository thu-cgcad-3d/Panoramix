#include "../src/core/feature.hpp"
#include "../src/vis/visualize2d.hpp"

#include "test_config.hpp"

using namespace panoramix;
using namespace test;

DEBUG_TEST(Feature, SegmentationExtractor) {
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
        p.useSLIC = true;
        p.superpixelSizeSuggestion = 3000;
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
    params.useLSD = true;
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

    core::Image im = cv::imread(ProjectDataDirStrings::LocalManhattan + "/buildings2.jpg");
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
    core::LineSegmentExtractor::Params lsparams;
    lsparams.useLSD = true;
    core::LineSegmentExtractor lineseg(lsparams);    
    core::LocalManhattanVanishingPointsDetector::Params params;

    std::vector<core::Image> ims = {
        cv::imread(ProjectDataDirStrings::LocalManhattan + "/cuboids.png")//,
        //cv::imread(ProjectDataDirStrings::LocalManhattan + "/buildings3.jpg"),
        //cv::imread(ProjectDataDirStrings::LocalManhattan + "/buildings2.jpg")
    };

    for (auto & im : ims) {
        core::ResizeToMakeWidthUnder(im, 1000);
        params.image = im;

        core::LocalManhattanVanishingPointsDetector vpdetector(params);

        std::vector<core::Line2> lines;
        core::LocalManhattanVanishingPointsDetector::Result result;

        lines = lineseg(im);
        //core::LoadFromDisk(ProjectDataDirStrings::Serialization + "/temp.state", lines);
        result = vpdetector(lines, core::Point2(im.cols / 2, im.rows / 2));

        for (auto & op : result.horizontalVanishingPointIds){
            std::cout << "[ " << op.first << ", " << op.second << "]" << "---"
                << result.vanishingPoints[op.first].value()
                << result.vanishingPoints[op.second].value() << std::endl;
        }

        std::vector<int> vpPairIds(result.vanishingPoints.size(), -1);
        for (int i = 0; i < result.horizontalVanishingPointIds.size(); i++){
            vpPairIds[result.horizontalVanishingPointIds[i].first] = i;
            vpPairIds[result.horizontalVanishingPointIds[i].second] = i;
        }

        auto ctable = vis::CreateRandomColorTableWithSize(result.vanishingPoints.size());
        ctable.exceptoinalColor() = vis::ColorFromTag(vis::ColorTag::White);

        auto viz = vis::Visualizer2D(im)
            << vis::manip2d::SetColorTable(ctable)
            << vis::manip2d::SetThickness(4);

        for (int i = 0; i < lines.size(); i++) {
            int claz = result.lineClasses[i];
            viz << vis::manip2d::SetThickness(claz == -1 ? 1 : vpPairIds[claz] + 3)
                << vis::manip2d::SetColor(ctable[claz])
                << core::NoteAs(lines[i], std::to_string(i) + "." + std::to_string(result.lineClasses[i]));
        }

        viz << vis::manip2d::SetColor(vis::ColorTag::Red)
            << result.vanishingPoints;

        viz /*<< vis::manip2d::SetColor(vis::ColorTag::Black)
            << result.hlineCands
            << vis::manip2d::SetColor(vis::ColorTag::Red)
            << result.horizon*/
            << vis::manip2d::Show(0);
    }
}



TEST(Feature, FeatureExtractor) {
    core::SegmentationExtractor segmenter;
    core::LineSegmentExtractor::Params params;
    params.useLSD = true;
    core::LineSegmentExtractor lineSegmentExtractor(params);
    
    core::CVFeatureExtractor<cv::SIFT> sift;
    core::CVFeatureExtractor<cv::SURF> surf(300.0);
    for (int i = 0; i < 4; i++) {
        std::string name = ProjectDataDirStrings::Normal + "/" + "sampled_" + std::to_string(i) + ".png";
        cv::Mat im = cv::imread(name);
        auto segs = segmenter(im);
        vis::Visualizer2D(im) 
            << [](vis::Visualizer2D & viz) { viz.params.winName = "haha"; }
        << vis::manip2d::SetColorTable(vis::CreateRandomColorTableWithSize(segs.second))
            << segs.first
            << lineSegmentExtractor(im) 
            << sift(im) << surf(im)
            << vis::manip2d::Show();
    }
}


int main(int argc, char * argv[], char * envp[])
{
    testing::InitGoogleTest(&argc, argv);
    return DEBUG_RUN_ALL_TESTS();
}
