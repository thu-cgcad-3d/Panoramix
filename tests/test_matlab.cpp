
#include <iostream>
#include <random>
#include <thread>

#include "../src/core/matlab.hpp"

#include "config.hpp"

using namespace panoramix;
using namespace test;


TEST(Matlab, CVX) {
#ifdef USE_MATLAB

    core::Matlab::RunScript("clear");
    core::Matlab::RunScript("cvx_setup");
    int m = 16, n = 8;
    core::Matlab::PutVariable("m", m);
    core::Matlab::PutVariable("n", n);
    core::Matlab::RunScript("[m n]");
    core::Matlab matlab;
    matlab
        << "A = randn(m, n);"
        << "b = randn(m, 1);"
        << "save tempfile;"
        << "cvx_begin"
        << "   variable x(n)"
        << "   minimize(norm(A*x - b))"
        << "cvx_end";
    std::vector<double> x;
    core::Matlab::GetVariable("x", x, false);
    ASSERT_EQ(x.size(), n);

#endif 
}


TEST(Matlab, CDAndAddPath){

#ifdef USE_MATLAB
    EXPECT_TRUE(core::Matlab::IsBuilt());
    ASSERT_TRUE(core::Matlab::IsUsable());
    core::Matlab::RunScript("clear;");
    core::Matlab::RunScript("clear functions;");
    core::Matlab::RunScript("cd " + ProjectDataDirStrings::Scripts + "/GC");
    

#else
    EXPECT_FALSE(core::Matlab::IsBuilt());
#endif   

}

TEST(Matlab, ImageConversion){

#ifdef USE_MATLAB
    EXPECT_TRUE(core::Matlab::IsBuilt());
    ASSERT_TRUE(core::Matlab::IsUsable());
    ASSERT_TRUE(core::Matlab::RunScript("x = 1;"));

    auto image = cv::imread(ProjectDataDirStrings::Normal + "/75.jpg");
    ASSERT_TRUE(core::Matlab::PutVariable("im", image));
    ASSERT_TRUE(core::Matlab::RunScript("imshow(im);"));
    ASSERT_TRUE(core::Matlab::RunScript("im = horzcat(im * 0.5, im);"));
    ASSERT_TRUE(core::Matlab::RunScript("imshow(im);"));
    core::Image newImage;
    ASSERT_TRUE(core::Matlab::GetVariable("im", newImage));

    cv::imshow("new image", newImage);
    cv::waitKey();

    core::ImageOfType<core::Vec<int, 3>> all123s(500, 500, core::Vec<int, 3>(1, 2, 3));
    ASSERT_TRUE(core::Matlab::PutVariable("all123s", all123s));
    ASSERT_TRUE(core::Matlab::RunScript("all246s = all123s * 2;"));
    core::ImageOfType<core::Vec<int, 3>> all246s;
    ASSERT_TRUE(core::Matlab::GetVariable("all246s", all246s));

    for (auto & i : all246s){
        ASSERT_TRUE(i == (core::Vec<int, 3>(2, 4, 6)));
    }
#else
    EXPECT_FALSE(core::Matlab::IsBuilt());
#endif   

}

TEST(Matlab, GeneralConversion){

#ifdef USE_MATLAB
    EXPECT_TRUE(core::Matlab::IsBuilt());
    ASSERT_TRUE(core::Matlab::IsUsable());
    auto points = std::vector<core::Point3>{core::Point3(1, 2, 3), core::Point3(4, 5, 6)};
    ASSERT_TRUE(core::Matlab::PutVariable("x", points));
    ASSERT_TRUE(core::Matlab::RunScript("x = x * 2;"));
    std::vector<core::Point3> npoints;
    ASSERT_TRUE(core::Matlab::GetVariable("x", npoints, true));

    ASSERT_EQ(points.size(), npoints.size());
    for (int i = 0; i < points.size(); i++){
        ASSERT_TRUE(points[i] * 2 == npoints[i]);
    }

#else
    EXPECT_FALSE(core::Matlab::IsBuilt());
#endif  

}

TEST(Matlab, SparseMatrix){

#ifdef USE_MATLAB
    EXPECT_TRUE(core::Matlab::IsBuilt());
    ASSERT_TRUE(core::Matlab::IsUsable());
    
    core::SparseMat<double> mat(2, (std::initializer_list<int>{ 15, 10 }).begin());
    mat.ref(0, 0) = -10.0;
    mat.ref(1, 2) = 5.0;
    mat.ref(3, 5) = 7.0;
    mat.ref(6, 5) = 10.0;
    mat.ref(6, 6) = 11.0;
    mat.ref(7, 6) = 12.0;
    mat.ref(14, 9) = 20.0;

    ASSERT_TRUE(core::Matlab::PutVariable("mat", mat));
    ASSERT_TRUE(core::Matlab::RunScript("mat"));

#else
    EXPECT_FALSE(core::Matlab::IsBuilt());
#endif  

}