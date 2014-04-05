#include "../src/core/basic_types.hpp"
#include "../src/core/expression.hpp"
#include "gtest/gtest.h"

#include <iostream>
#include <random>

using namespace panoramix;


TEST(Expression, Basic) {   

    core::Expression<float> e = 1.0f;
    core::Expression<int> f = 1;
    auto fe = f + e;
    auto c = fe + core::MakeConstantValue(2);
    auto d = -c;
    auto g = core::sqrt(-d);

    ASSERT_FLOAT_EQ(-4, d.eval());
    ASSERT_FLOAT_EQ(2, g.eval());

    core::Expression<core::Vec3> v = core::Vec3(0, 0, 1);
    auto vvv = v * core::MakeConstantValue(3);
    auto vv = vvv - v;
    
    ASSERT_FLOAT_EQ(2, cv::norm(vv.eval()));

    int a = 0;
    auto s = core::MakeConstantValue(a);

    int v1, v2;
    auto var1 = core::MakeVariableAt(&v1);
    auto var2 = core::MakeVariableAt(&v2);
    auto var12 = var1 + var2;

    for (int i = 0; i < 1000; i++){
        v1 = rand();
        v2 = rand();
        ASSERT_EQ(v1 + v2, var12.eval());
    }
    
}

TEST(Expression, Derivative) {

    using namespace core;

    auto a = MakeConstantValue(Vec4(1.0, 1.0, 1.0));
    Vec4 xdata, ydata;
    auto x = MakeVariableAt(&xdata);
    auto y = MakeVariableAt(&ydata);
    
    auto axy = a * x * y;
    //auto axy_dx = axy.derivativeOn(x);

}

int main(int argc, char * argv[], char * envp[])
{
	testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}



