#include "../src/vis/renderable_object_tree.hpp"
#include "../src/vis/visualize3d.hpp"
#include "gtest/gtest.h"

#include <iostream>
#include <string>
#include <random>

#include "test_config.hpp"

using namespace panoramix;
using namespace test;

static_assert(vis::CanMakeRenderable<core::Line3>::value, "Line3 is not renderable!");
static_assert(vis::CanMakeRenderable<core::Point3>::value, "Point3 is not renderable!");
static_assert(vis::CanMakeRenderable<core::Classified<core::Line3>>::value, "Classified<Line3> is not renderable!");

TEST(Visualizer3D, RenderPoints) {
    std::vector<core::Point3> points = {
        {-1, -1, -1},
        {-1, 1, 1},
        {1, -1, 1}
    };
    vis::Visualizer3D()
        << vis::manip3d::SetDefaultForegroundColor(vis::ColorTag::DimGray)
        << vis::manip3d::SetDefaultPointSize(20.0)
        << points
        << vis::manip3d::SetBackgroundColor(vis::ColorTag::White)
        << vis::manip3d::Show();
}

TEST(Visualizer3D, RenderLines) { 
    std::vector<core::Line3> lines = {
        { {-1, 1}, {5, 5} },
        { { -2, -2 }, {8, 6} },
        { { 5, -3, -4 }, {-3, 2, 8} }
    };    
    vis::Visualizer3D()
        << vis::manip3d::SetDefaultForegroundColor(vis::ColorTag::Yellow)
        << vis::manip3d::SetDefaultLineWidth(3.5)
        << lines
        << vis::manip3d::SetBackgroundColor(vis::ColorTag::Red)
        << vis::manip3d::Show();
}

TEST(Visualizer3D, RenderClassifiedLines) {
    std::vector<core::Classified<core::Line3>> lines = {
        { -1, { { -1, 1 }, { 5, 5 } } },
        { 0, { { -2, -2 }, { 8, 6 } } },
        { 1, { { 5, -3, -4 }, { -3, 2, 8 } } }
    };
    vis::Visualizer3D()
        << vis::manip3d::SetDefaultColorTable(vis::ColorTable({vis::ColorTag::Blue, vis::ColorTag::Yellow}, vis::ColorTag::Black))
        << vis::manip3d::SetDefaultLineWidth(4)
        << lines
        << vis::manip3d::SetBackgroundColor(vis::ColorTag::White)
        << vis::manip3d::Show();
}

TEST(Visualizer3D, RenderSpatialProjectedPolygon) {
    vis::SpatialProjectedPolygon p[2];
    p[0].plane = core::Plane3({ 0, -1, 0 }, { 0, 1, 0 });
    p[0].projectionCenter = { 0, 0, 0 };
    for (int i = 0; i < 5; i++){
        p[0].corners.emplace_back(cos(i / 5.0*M_PI * 2), -1.0, sin(i / 5.0*M_PI * 2));
    }
    p[1].plane = core::Plane3({ 0, 1, 0 }, { 0, 1, 0 });
    p[1].projectionCenter = { 0, 0, 0 };
    for (int i = 0; i < 5; i++){
        p[1].corners.emplace_back(cos(i / 5.0*M_PI * 2), -1.0, sin(i / 5.0*M_PI * 2));
    }

    core::Image im = cv::imread(ProjectDataDirStrings::PanoramaIndoor + "/13.jpg");
    vis::Visualizer3D()
        << vis::manip3d::Begin(std::vector<vis::SpatialProjectedPolygon>(p, p + 2))
        << vis::manip3d::SetTexture(im)
        << vis::manip3d::End
        << vis::manip3d::SetBackgroundColor(vis::ColorTag::White)
        << vis::manip3d::Show();
}



int main(int argc, char * argv[], char * envp[])
{
    testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}
