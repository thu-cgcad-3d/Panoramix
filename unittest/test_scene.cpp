#include "../src/gui/canvas.hpp"
#include "../src/gui/discretization.hpp"
#include "../src/gui/scene.hpp"
#include "gtest/gtest.h"

#include <iostream>
#include <string>
#include <random>

#include "config.hpp"

using namespace panoramix;
using namespace test;



TEST(Canvas, Color) {

    int m = 600;
    core::ImageOfType<core::Vec<double, 4>> hs(m, m);
    core::ImageOfType<core::Vec<double, 4>> sv(m, m);
    for (int i = 0; i < m; i++){
        for (int j = 0; j < m; j++){
            hs(i, j) = gui::ColorFromHSV(i / (double)m, j / (double)m, 0.5);
            sv(i, j) = gui::ColorFromHSV(0.5, i / (double)m, j / (double)m);
        }
    }
    cv::imshow("[HS]V color", hs);
    cv::imshow("H[SV] color", sv);
    cv::waitKey();

}


static_assert(gui::IsDiscretizable<core::Line3>::value, "Line3 is not renderable!");
static_assert(gui::IsDiscretizable<core::Point3>::value, "Point3 is not renderable!");
static_assert(gui::IsDiscretizable<core::Classified<core::Line3>>::value, "Classified<Line3> is not renderable!");
static_assert(gui::IsDiscretizable<core::LayeredShape3>::value, "");
static_assert(!gui::IsDiscretizable<int>::value, "");

static auto a = [](gui::InteractionID iid, core::Sphere3 & s){
    if (iid == gui::ClickLeftButton)
        std::cout << "clicked on the sphere" << std::endl;
};

static auto b = [](gui::InteractionID iid, const gui::SceneObjectTree & tree, const gui::SceneObjectMeshTriangle & omt){
    if (iid == gui::ClickLeftButton)
        std::cout << "clicked on the sphere" << std::endl;
};

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

TEST(Scene, LayeredShape3){

    using namespace core;
    using namespace gui;

    core::LayeredShape3 ls;
    ls.layers = {
        { Point3(0, 0, 0), Point3(1, 0, 0), Point3(1, 1, 0), Point3(0, 1, 0) },
        { Point3(0, 0, 0.5), Point3(1, 0, 0.5), Point3(1, 1, 0.5), Point3(0, 1, 0.5) },
        { Point3(0, 0, 1)/*, Point3(1, 0, 1)*/, Point3(1, 1, 1), Point3(0, 1, 1) }
    };
    ls.normal = Vec3(0, 0, 1);
    
    gui::ResourceStore::set("texture", cv::imread(ProjectDataDirStrings::PanoramaIndoor + "/13.jpg"));
    SceneBuilder().begin(ls)
        .resource("texture").shaderSource(gui::OpenGLShaderSourceDescriptor::XPanorama)
        .end()
        .show(true, true, gui::RenderOptions().renderMode(gui::RenderModeFlag::Lines | gui::RenderModeFlag::Triangles));

}

TEST(Scene, Interaction){

    using namespace gui;    
   
    std::list<core::Sphere3> ss = { { core::Point3(1, 1, 1), 2.0 }, { core::Point3(-1, -1, -1), 2.0 } };
    core::Line3 ll[] = { { core::Point3(-1, 1, 2), core::Point3(1, -1, 0) }, { core::Point3(0, 2, 3), core::Point3(-1, -2, -3) } };
    //auto t = core::MakeTriangle(core::Point3(0, 0, 1), core::Point3(1, 0, 0), core::Point3(0, 1, 0));
    
    gui::ResourceStore::set("texture", cv::imread(ProjectDataDirStrings::PanoramaIndoor + "/13.jpg"));

    int clickedCount = 0;
    SceneBuilder()
        .begin(ss, [&clickedCount](gui::InteractionID iid, core::Sphere3 & s){
            if (iid == gui::ClickLeftButton)
                std::cout << "clicked on the spheres, its center is at " << s.center << std::endl;
            else
                std::cout << "pressed on the spheres, its center is at " << s.center << std::endl;
            })
            .resource("texture")
            .shaderSource(gui::OpenGLShaderSourceDescriptor::XPanorama)
        .end()
        .begin(ll, [&clickedCount](gui::InteractionID iid, core::Line3 & l){
                if (iid == gui::ClickLeftButton)
                    std::cout << "clicked on the line (" << l.first << ", " << l.second << ")" << std::endl;
                else
                    std::cout << "pressed on the line (" << l.first << ", " << l.second << ")" << std::endl;
            })
            .shaderSource(gui::OpenGLShaderSourceDescriptor::XLines)
        .end()
        .show(true, true, gui::RenderOptions().renderMode(gui::RenderModeFlag::Lines | gui::RenderModeFlag::Triangles));
    
}