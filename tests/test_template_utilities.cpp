#include "../src/core/template_utilities.hpp"
#include "gtest/gtest.h"

#include <iostream>
#include <random>

using namespace panoramix;

TEST(Tuples, Invoke){
    
    auto args = std::make_tuple(1, 2, 3, 4, 5.0);
    auto fun = [](int a, int b, int c, int d, double e) {
        return (a + b + c + d + e) * 2.0;
    };
    
    auto sumDoubled = core::Invoke(fun, args);
    ASSERT_EQ(30, sumDoubled);

    auto fun2 = [](int a, int b, int c, int d, double e) {
        std::cout << a << b << c << d << e << std::endl;
    };

    core::InvokeWithoutReturn(fun, args);

}

int main(int argc, char * argv[], char * envp[])
{
	testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}

