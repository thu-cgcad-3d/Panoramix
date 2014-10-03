#include "../src/core/feature.hpp"
#include "../src/vis/visualize2d.hpp"

#include "test_config.hpp"

using namespace panoramix;
using namespace test;

TEST(Feature, SegmentationExtractor) {
    core::SegmentationExtractor::Params p;
    p.c = 5;
    p.minSize = 400;
    p.sigma = 1;
    core::SegmentationExtractor seg(p);
    core::Image im = cv::imread(ProjectDataDirStrings::Normal + "/75.jpg");
    vis::Visualizer2D(im) << vis::manip2d::Show();
    core::Image segim = seg(im, true);
    vis::Visualizer2D(segim) << vis::manip2d::Show();
}

TEST(Feature, LineSegmentExtractor) {
    core::LineSegmentExtractor lineseg;
    core::Image im = cv::imread(ProjectDataDirStrings::Normal + "/2148.jpg");
    vis::Visualizer2D(im) 
        << vis::manip2d::SetColor(vis::ColorTag::Yellow) 
        << vis::manip2d::SetThickness(2) << 
        lineseg(im) 
        << vis::manip2d::Show();
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


DEBUG_TEST(Feature, LocalManhattanVanishingPointDetector) {
    core::LineSegmentExtractor lineseg;
    
    core::LocalManhattanVanishingPointsDetector::Params params;
	core::Image im = cv::imread(ProjectDataDirStrings::LocalManhattan + "/cuboids.png");
    core::ResizeToMakeWidthUnder(im, 1000);
    params.image = im;

    core::LocalManhattanVanishingPointsDetector vpdetector(params);	

	std::vector<core::Line2> lines;
    core::LocalManhattanVanishingPointsDetector::Result result;
    
    lines = lineseg(im);
    //core::LoadFromDisk(ProjectDataDirStrings::Serialization + "/temp.state", lines);
    result = vpdetector(lines, core::Point2(im.cols / 2, im.rows / 2));

    for (auto & op : result.horizontalVanishingPointIds){
        std::cout << "[ " << op.first << ", " << op.second << "]" << std::endl;
    }
    
    auto ctable = vis::CreateGreyColorTableWithSize(result.vanishingPoints.size());
    ctable.randomize();

    auto viz = vis::Visualizer2D(im)
        << vis::manip2d::SetColorTable(ctable)
        << vis::manip2d::SetThickness(4);
  
    for (int i = 0; i < lines.size(); i++) {
        viz //<< vis::manip2d::SetThickness(i + 1)
            << vis::manip2d::SetColor(ctable[result.lineClasses[i]])
            << core::NoteAs(lines[i], std::to_string(result.lineClasses[i]));
    }

    viz << vis::manip2d::SetColor(vis::ColorTag::Red)
        << result.vanishingPoints;

    viz /*<< vis::manip2d::SetColor(vis::ColorTag::Black)
        << result.hlineCands
        << vis::manip2d::SetColor(vis::ColorTag::Red)
        << result.horizon*/
        << vis::manip2d::Show(0);
    
    return;
    

   /* std::vector<core::Classified<core::Line2>> classifiedLines;
    for (int i = 0; i < lines.size(); i++) {
        int c = -1;
        if (result.lineClasses[i] == 0){
            c = 0;
        }
        else{
            for (int j = 0; j < result.horizontalVanishingPointIds.size(); j++){
                if (result.lineClasses[i] == result.horizontalVanishingPointIds[j].first ||
                    result.lineClasses[i] == result.horizontalVanishingPointIds[j].second){
                    c = j+1;
                    break;
                }
            }
        }
        
        classifiedLines.push_back({ c, lines[i] });
	}

	
    vis::Visualizer2D(im)
		<< vis::manip2d::SetColorTable(colorTable)
        << vis::manip2d::SetThickness(2) 
		<< classifiedLines
        << vis::manip2d::Show(0);*/
}


TEST(Feature, PerspectiveCamera){
    core::PerspectiveCamera cam(1000, 1000, 500, 
        core::Vec3(0, 0, 0), 
        core::Vec3(1, 0, 0));

    for (int i = 0; i < 100; i++){
        core::Vec2 v(abs(rand()), abs(rand()));
        auto p = cam.spatialDirection(v);
        auto v2 = cam.screenProjection(p);
        double dist = core::norm(v - v2);
        ASSERT_LT(dist, 0.01);
    }
    auto c = cam.screenProjection(cam.center());
    double dist = core::norm(c - core::Vec2(cam.screenSize().width / 2, cam.screenSize().height / 2));
    ASSERT_LT(dist, 2);
}

TEST(Feature, PerspectiveCameraRandom){
    for (int k = 0; k < 100; k++) {
        core::PerspectiveCamera cam(abs(rand()) % 500, abs(rand()) % 400, abs(rand()) % 600,
            core::Vec3(rand(), rand(), rand()),
            core::Vec3(rand(), rand(), rand()));
        for (int i = 0; i < 100; i++){
            core::Vec2 v(rand(), rand());
            auto p = cam.spatialDirection(v);
            auto v2 = cam.screenProjection(p);
            double dist = core::norm(v - v2);
            if (!std::isnan(dist) && !std::isinf(dist))
                ASSERT_LT(dist, 0.01);
        }
        auto c = cam.screenProjection(cam.center());
        double dist = core::norm(c - core::Vec2(cam.screenSize().width / 2, cam.screenSize().height / 2));
        if (!std::isnan(dist) && !std::isinf(dist))
            ASSERT_LT(dist, 2);
    }
}


TEST(Feature, CameraSampler) {
    cv::Mat im = cv::imread(ProjectDataDirStrings::PanoramaOutdoor + "/panofactory.jpg");

    EXPECT_EQ(4000, im.cols);
    EXPECT_EQ(2000, im.rows);
    cv::resize(im, im, cv::Size(1000, 500));
    vis::Visualizer2D viz(im);
    viz.params.alphaForNewImage = 0.3;

    core::PanoramicCamera originCam(im.cols / M_PI / 2.0);
    core::PanoramicCamera newCam(im.cols / M_PI / 2.0, 
        core::Vec3(0, 0, 0), 
        core::Vec3(0, 0, 1),
        core::Vec3(0, 1, 0));
    viz << core::CameraSampler<core::PanoramicCamera, core::PanoramicCamera>(newCam, originCam)(im)
        << core::SegmentationExtractor()(im, true)
        << vis::manip2d::Show();

    float camPositions[4][3] = {
        {1, 0, 0},
        {0, 1, 0},
        {-1, 0, 0},
        {0, -1, 0}
    };

    for (int i = 0; i < 4; i++){
        core::PerspectiveCamera cam(500, 600, 150,
            core::Vec3(0, 0, 0),
            core::Vec3(camPositions[i][0], camPositions[i][1], camPositions[i][2]),
            core::Vec3(0, 0, -1));
        core::CameraSampler<core::PerspectiveCamera, core::PanoramicCamera> sampler(cam,
            core::PanoramicCamera(im.cols / M_PI / 2.0));
        cv::Mat sampledIm = sampler(im);
        vis::Visualizer2D(sampledIm) << vis::manip2d::Show();
    }
}


TEST(Feature, FeatureExtractor) {
    core::SegmentationExtractor segmenter;
    core::LineSegmentExtractor::Params params;
    params.useExperimentalAlgorithm = false;
    core::LineSegmentExtractor lineSegmentExtractor(params);
    
    core::CVFeatureExtractor<cv::SIFT> sift;
    core::CVFeatureExtractor<cv::SURF> surf(300.0);
    for (int i = 0; i < 4; i++) {
        std::string name = ProjectDataDirStrings::Normal + "/" + "sampled_" + std::to_string(i) + ".png";
        cv::Mat im = cv::imread(name);
        vis::Visualizer2D(im) 
            << [](vis::Visualizer2D & viz) { viz.params.winName = "haha"; }
            << segmenter(im, true)    
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
