#include "../src/core/template_utilities.hpp"
#include "gtest/gtest.h"

#include <iostream>
#include <random>

using namespace panoramix;

TEST(Tuples, Tuples){
    auto args = std::make_tuple(1, 2, 3, 4, 5.0);
    std::cout << core::Invoke([](int a, int b, int c, int d, double e) {
        return (a + b + c + d) * 2.0; 
    }, args) << std::endl;

}

int main(int argc, char * argv[], char * envp[])
{
	testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}

