#include "../src/vis/visualizers.hpp"
#include "gtest/gtest.h"

#include <iostream>
#include <string>
#include <random>

#include "config.hpp"

using namespace panoramix;
using namespace test;



TEST(Visualizer, Color) {

    int m = 600;
    core::ImageWithType<core::Vec<double, 4>> hs(m, m);
    core::ImageWithType<core::Vec<double, 4>> sv(m, m);
    for (int i = 0; i < m; i++){
        for (int j = 0; j < m; j++){
            hs(i, j) = vis::ColorFromHSV(i / (double)m, j / (double)m, 0.5);
            sv(i, j) = vis::ColorFromHSV(0.5, i / (double)m, j / (double)m);
        }
    }
    cv::imshow("[HS]V color", hs);
    cv::imshow("H[SV] color", sv);
    cv::waitKey();

}


static_assert(vis::IsDiscretizable<core::Line3>::value, "Line3 is not renderable!");
static_assert(vis::IsDiscretizable<core::Point3>::value, "Point3 is not renderable!");
static_assert(vis::IsDiscretizable<core::Classified<core::Line3>>::value, "Classified<Line3> is not renderable!");
static_assert(vis::IsDiscretizable<std::vector<core::Classified<core::Line3>>>::value, 
    "std::vector<core::Classified<core::Line3>> is not renderable!");
static_assert(!vis::IsDiscretizable<int>::value, "");

static auto a = [](vis::InteractionID iid, core::Sphere3 & s){
    if (iid == vis::ClickLeftButton)
        std::cout << "clicked on the sphere" << std::endl;
};

static auto b = [](vis::InteractionID iid, const vis::VisualObjectTree & tree, const vis::VisualObjectMeshTriangle & omt){
    if (iid == vis::ClickLeftButton)
        std::cout << "clicked on the sphere" << std::endl;
};

static_assert(vis::IsSimpleCallbackFunctionImp<std::pair<decltype(a), core::Sphere3>>::value, "");
static_assert(!vis::IsSimpleCallbackFunctionImp<std::pair<int, core::Sphere3>>::value, "");

static_assert(vis::CallbackFunctionTraits<decltype(a), core::Sphere3>::value == vis::CallbackFunctionType::Simple, "");
static_assert(vis::CallbackFunctionTraits<decltype(b), core::Sphere3>::value == vis::CallbackFunctionType::Complex, "");
static_assert(vis::CallbackFunctionTraits<int, core::Sphere3>::value == vis::CallbackFunctionType::None, "");

void Print(const core::PerspectiveCamera & cam) {
    std::cout << std::setprecision(6) << "======================================" << std::endl;
    std::cout << std::setprecision(6) << "eye:\t" << cam.eye() << std::endl;
    std::cout << std::setprecision(6) << "center:\t" << cam.center() << std::endl;
    std::cout << std::setprecision(6) << "up:\t" << cam.up() << std::endl;
    std::cout << std::setprecision(6) << "size:\t" << cam.screenSize() << std::endl;
    std::cout << std::setprecision(6) << "focal:\t" << cam.focal() << std::endl;
    std::cout << std::setprecision(6) << "near:\t" << cam.nearPlane() << std::endl;
    std::cout << std::setprecision(6) << "far:\t" << cam.farPlane() << std::endl;
    std::cout << std::setprecision(6) << "viewMat:\t" << std::endl << cam.viewMatrix() << std::endl;
    std::cout << std::setprecision(6) << "projectMat:\t" << std::endl << cam.projectionMatrix() << std::endl;
    std::cout << std::setprecision(6) << "viewProjMat:\t" << std::endl << cam.viewProjectionMatrix() << std::endl;
    std::cout << std::setprecision(6) << "======================================" << std::endl;
}

TEST(Visualizer, Interaction){

    using namespace vis;    
   
    std::vector<core::Sphere3> ss = { { core::Point3(1, 1, 1), 2.0 }, { core::Point3(-1, -1, -1), 2.0 } };

    auto t = core::MakeTriangle(core::Point3(0, 0, 1), core::Point3(1, 0, 0), core::Point3(0, 1, 0));
    
    vis::ResourceStore::set("texture", cv::imread(ProjectDataDirStrings::PanoramaIndoor + "/13.jpg"));

    int clickedCount = 0;
    Visualizer("_1")
        .begin(ss, [&clickedCount](vis::InteractionID iid, core::Sphere3 & s){
            if (iid == vis::ClickLeftButton)
                std::cout << "clicked on the spheres, its center is at " << s.center << std::endl;
            else
                std::cout << "pressed on the spheres, its center is at " << s.center << std::endl;
            })
            .resource("texture")
            .shaderSource(vis::PredefinedShaderSource(vis::OpenGLShaderSourceDescriptor::XPanorama))
        .end()
        .renderMode(vis::RenderModeFlag::Lines | vis::RenderModeFlag::Triangles)
        .show();
    
}


//int main(int argc, char * argv[], char * envp[]) {
//    testing::InitGoogleTest(&argc, argv);
//    testing::GTEST_FLAG(catch_exceptions) = false;
//    testing::GTEST_FLAG(throw_on_failure) = true;
//    return RUN_ALL_TESTS();
//}
