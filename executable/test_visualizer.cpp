#include "../src/vis/visualizers.hpp"
#include "gtest/gtest.h"

#include <iostream>
#include <string>
#include <random>

#include "test_config.hpp"

using namespace panoramix;
using namespace test;

static_assert(vis::IsDiscretizable<core::Line3>::value, "Line3 is not renderable!");
static_assert(vis::IsDiscretizable<core::Point3>::value, "Point3 is not renderable!");
static_assert(vis::IsDiscretizable<core::Classified<core::Line3>>::value, "Classified<Line3> is not renderable!");
static_assert(vis::IsDiscretizable<std::vector<core::Classified<core::Line3>>>::value, 
    "std::vector<core::Classified<core::Line3>> is not renderable!");
static_assert(!vis::IsDiscretizable<int>::value, "");

TEST(Visualizer, _1){

    using namespace vis;
    std::vector<core::Point3> pts(10);
    for (int i = 0; i < pts.size(); i++)
        pts[i] = core::Point3(i / 10.0, i / 10.0, i / 10.0);
    core::Sphere3 s = { core::Point3(0, 0, 0), 0.3f };
    vis::ResourceStore::set("texture", cv::imread(ProjectDataDirStrings::PanoramaIndoor + "/13.jpg"));
    Visualizer("_1")
        .begin(core::Line3(core::Point3(1, 0, 1), core::Point3(0, 1, 0)))
            .shaderSource(vis::PredefinedShaderSource(vis::OpenGLShaderSourceDescriptor::XLines))
        .end()
        .beginWithBinding(s, [](vis::InteractionID iid, core::Sphere3 & s){
            if (iid == vis::ClickLeftButton)
                std::cout << "clicked on the sphere" << std::endl;
            })
            .resource("texture")
            .shaderSource(vis::PredefinedShaderSource(vis::OpenGLShaderSourceDescriptor::XPanorama))
        .end()
        .add(pts)
        .show();
}


int main(int argc, char * argv[], char * envp[]) {
    testing::InitGoogleTest(&argc, argv);
    testing::GTEST_FLAG(catch_exceptions) = false;
    testing::GTEST_FLAG(throw_on_failure) = true;
    return RUN_ALL_TESTS();
}
