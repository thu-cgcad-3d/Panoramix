#include "../src/core/util.hpp"

#include <random>
#include "gtest/gtest.h"

using namespace panoramix;

TEST(UtilTest, WrapBetween) {
	for (int i = 0; i < 10000; i++){
		double x = double(rand()) / rand() + rand();
		double a = double(rand()) / rand() + rand();
		double b = a + abs(double(rand())/rand());
		if (a == b || std::isnan(x) || std::isnan(a) || std::isnan(b) || 
			std::isinf(x) || std::isinf(a) || std::isinf(b))
			continue;
		double xx = core::WrapBetween(x, a, b);
		double rem = (xx - x) / (b - a) - std::round((xx - x) / (b - a));
		if (std::isnan(rem)){
			assert(0);
		}
		ASSERT_NEAR(0, rem, 1e-5);
		ASSERT_LE(a, xx);
		ASSERT_LT(xx, b);
	}
}


int main(int argc, char * argv[], char * envp[])
{
	testing::InitGoogleTest(&argc, argv);
	return RUN_ALL_TESTS();
	return 0;
}