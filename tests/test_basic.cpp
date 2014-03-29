#include "../src/core/version.hpp"
#include "../src/core/basic_types.hpp"
#include "gtest/gtest.h"

#include <iostream>
#include <random>

using namespace panoramix;

TEST(ConfigTest, Version) {
	EXPECT_EQ(PANORAMIX_VERSION_MAJOR, core::GetVersion().major);
	EXPECT_EQ(PANORAMIX_VERSION_MINOR, core::GetVersion().minor);
}

TEST(BasicTypeTest, HPoint) {
	for (int i = 0; i < 1000; i++){
		core::Vec4 v4;
		std::generate(v4.val, v4.val + 4, std::rand);
		auto hp = core::HPointFromVector(v4);
		auto p = hp.toVector();
		ASSERT_LT(core::norm(p - v4), 1e-5);
		core::HPoint<double, 4> hp5 = v4;
		auto p5 = hp5.toPoint();
		ASSERT_LT(core::norm(p5 - v4), 1e-5);
	}
}

int main(int argc, char * argv[], char * envp[])
{
	for (int i = 0; i < argc; i++) {
		std::cout << "[INPUT]:" << argv[i] << std::endl;
	}
	char** env;
	for (env = envp; *env != 0; env++) {
		char* thisEnv = *env;
		std::cout << "[ENV]:" << thisEnv << std::endl;
	}
	testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}



