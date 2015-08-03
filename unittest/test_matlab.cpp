
#include <iostream>
#include <random>
#include <thread>

//#include "../src/misc/matlab_engine.hpp"

#include "config.hpp"

using namespace pano;
using namespace test;


//TEST(MatlabEngine, CVX) {
//    if (misc::MatlabEngine::IsBuilt()){
//        misc::MatlabEngine::Start();
//        misc::MatlabEngine::RunScript("clear");
//        misc::MatlabEngine::RunScript("cvx_setup");
//        int m = 16, n = 8;
//        std::vector<double> b(m);
//        std::generate(b.begin(), b.end(), std::rand);
//        misc::MatlabEngine::PutVariable("m", m);
//        misc::MatlabEngine::PutVariable("n", n);
//        misc::MatlabEngine::PutVariable("b", b);
//        misc::MatlabEngine::RunScript("b = b';");
//        misc::MatlabEngine::RunScript("[m n]");
//        misc::MatlabEngine matlab;
//        matlab
//            << "A = randn(m, n);"
//            //<< "b = randn(m, 1);"
//            << "save tempfile;"
//            << "cvx_begin"
//            << "   variable x(n)"
//            << "   minimize(norm(A*x - b))"
//            << "cvx_end";
//        std::vector<double> x;
//        misc::MatlabEngine::GetVariable("x", x, false);
//        ASSERT_EQ(x.size(), n);
//        misc::MatlabEngine::Close();
//    }
//
//    misc::MatlabEngine matlab;
//    matlab << "clear"
//        << "cvx_setup";
//    int m = 16, n = 8;
//    std::vector<double> b(m);
//    std::generate(b.begin(), b.end(), std::rand);
//    matlab.setVar("m", m);
//    matlab.setVar("n", n);
//    matlab.setVar("b", misc::MXArray(b));
//    matlab << "b = b'; [m, n]";
//    matlab
//        << "A = randn(m, n);"
//        //<< "b = randn(m, 1);"
//        << "save tempfile;"
//        << "cvx_begin"
//        << "   variable x(n)"
//        << "   minimize(norm(A*x - b))"
//        << "cvx_end";
//    std::vector<double> x;
//
//}
//
//
//TEST(MatlabEngine, CDAndAddPath){
//
//    if (misc::MatlabEngine::IsBuilt()){
//        misc::MatlabEngine::Start();
//        misc::MatlabEngine::RunScript("clear;");
//        misc::MatlabEngine::RunScript("clear functions;");
//        misc::MatlabEngine::RunScript("cd " + ProjectDataDirStrings::Scripts + "/GC");
//        misc::MatlabEngine::Close();
//    }
//}
//
//TEST(MatlabEngine, ImageConversion){
//
//    if (misc::MatlabEngine::IsBuilt()){
//        misc::MatlabEngine::Start();
//        ASSERT_TRUE(misc::MatlabEngine::RunScript("x = 1;"));
//
//        auto image = core::ImageRead(ProjectDataDirStrings::Normal + "/75.jpg");
//        ASSERT_TRUE(misc::MatlabEngine::PutVariable("im", image));
//        ASSERT_TRUE(misc::MatlabEngine::RunScript("imshow(im);"));
//        ASSERT_TRUE(misc::MatlabEngine::RunScript("im = horzcat(im * 0.5, im);"));
//        ASSERT_TRUE(misc::MatlabEngine::RunScript("imshow(im);"));
//        core::Image newImage;
//        ASSERT_TRUE(misc::MatlabEngine::GetVariable("im", newImage));
//
//        cv::imshow("new image", newImage);
//        cv::waitKey();
//
//        core::ImageOf<core::Vec<int, 3>> all123s(500, 500, core::Vec<int, 3>(1, 2, 3));
//        ASSERT_TRUE(misc::MatlabEngine::PutVariable("all123s", all123s));
//        ASSERT_TRUE(misc::MatlabEngine::RunScript("all246s = all123s * 2;"));
//        core::ImageOf<core::Vec<int, 3>> all246s;
//        ASSERT_TRUE(misc::MatlabEngine::GetVariable("all246s", all246s));
//
//        for (auto & i : all246s){
//            ASSERT_TRUE(i == (core::Vec<int, 3>(2, 4, 6)));
//        }
//        misc::MatlabEngine::Close();
//    }
//
//}
//
//TEST(MatlabEngine, GeneralConversion){
//
//    if (misc::MatlabEngine::IsBuilt()){
//        misc::MatlabEngine::Start();
//        auto points = std::vector<core::Point3>{core::Point3(1, 2, 3), core::Point3(4, 5, 6)};
//        ASSERT_TRUE(misc::MatlabEngine::PutVariable("x", points));
//        ASSERT_TRUE(misc::MatlabEngine::RunScript("x = x * 2;"));
//        std::vector<core::Point3> npoints;
//        ASSERT_TRUE(misc::MatlabEngine::GetVariable("x", npoints, true));
//
//        ASSERT_EQ(points.size(), npoints.size());
//        for (int i = 0; i < points.size(); i++){
//            ASSERT_TRUE(points[i] * 2 == npoints[i]);
//        }
//        misc::MatlabEngine::Close();
//    }
//
//}
//
//TEST(MatlabEngine, SparseMatrix){
//
//
//    if (misc::MatlabEngine::IsBuilt()){
//        misc::MatlabEngine::Start();
//        core::SparseMat<double> mat(2, (std::initializer_list<int>{ 15, 10 }).begin());
//        mat.ref(0, 0) = -10.0;
//        mat.ref(1, 2) = 5.0;
//        mat.ref(3, 5) = 7.0;
//        mat.ref(6, 5) = 10.0;
//        mat.ref(6, 6) = 11.0;
//        mat.ref(7, 6) = 12.0;
//        mat.ref(14, 9) = 20.0;
//
//        ASSERT_TRUE(misc::MatlabEngine::PutVariable("mat", mat));
//        ASSERT_TRUE(misc::MatlabEngine::RunScript("mat"));
//        misc::MatlabEngine::Close();
//    }
//}