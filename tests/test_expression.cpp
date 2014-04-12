#include "../src/core/basic_types.hpp"
#include "../src/core/expression.hpp"
#include "gtest/gtest.h"

#include <iostream>
#include <random>
#include <chrono>

using namespace panoramix::core;

TEST(Expression, Basic) {   
    
    Mat4 m1 = Mat4::eye();
    Mat4 m2 = Mat4::ones();
    auto e1 = MakeConstantValue(Mat4::eye());
    auto e2 = MakeConstantValue(Mat4::ones());
    auto fun = [](double d1, double d2){return d1 + d2 * 2; };
    auto e12 = PerformElementWiseOperation(fun, e1, e2);
    
    for (int i = 0; i < 4; i++){
        for (int j = 0; j < 4; j++){
            ASSERT_EQ(fun(m1(i, j), m2(i, j)), e12(i, j));
        }
    }

}

TEST(Expression, Speed1) {

    //Mat4 m1 = Mat4::randn(0, 1);
    //Mat4 m2 = Mat4::randn(0, 1);
    //auto e1 = MakeConstantValue(m1);
    //auto e2 = MakeConstantValue(m2);
    //auto fun = [](double d1, double d2){return d1 + d2 * 2 + sin(d1 - d2); };
    ////auto e12 = PerformElementwiseOperation(fun, e1, e2);

    //using namespace std::chrono;
    //auto start = high_resolution_clock::now();

    //Mat4 s1 = Mat4::zeros();
    //for (int k = 0; k < 10000; k++){
    //    Mat4 m;
    //    for (int i = 0; i < 4; i++){
    //        for (int j = 0; j < 4; j++){
    //            m(i, j) = fun(m1(i, j), m2(i, j));
    //        }
    //    }
    //    s1 += m;
    //}

    //std::cout << "time cost for direct computation: " << 
    //    duration_cast<milliseconds>(high_resolution_clock::now() - start).count() << std::endl;
}

TEST(Expression, Speed2) {

   /* Mat4 m1 = Mat4::randn(0, 1);
    Mat4 m2 = Mat4::randn(0, 1);
    auto e1 = MakeConstantValue(m1);
    auto e2 = MakeConstantValue(m2);
    auto fun = [](double d1, double d2){return d1 + d2 * 2 + sin(d1 - d2); };
    auto e12 = PerformElementwiseOperation(fun, e1, e2);

    using namespace std::chrono;
    auto start = high_resolution_clock::now();

    Mat4 s1 = Mat4::zeros();
    auto es = MakeConstantValue(s1);

    for (int k = 0; k < 10000; k++){
        es = PerformElementwiseOperation([](double d1, double d2){return d1 + d2; }, es, e12);
    }

    std::cout << "time cost for delayed computation: " <<
        duration_cast<milliseconds>(high_resolution_clock::now() - start).count() << std::endl;*/
}



TEST(Expression, Derivative) {

   

}

int main(int argc, char * argv[], char * envp[])
{
	testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}



