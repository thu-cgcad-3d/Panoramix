#include "../src/config.hpp"
#include "../src/core/Dummy.hpp"
#include "gtest/gtest.h"

#include <iostream>

using namespace panoramix;

TEST(DummyTest, Zero) {
	EXPECT_EQ(PANORAMIX_VERSION_MAJOR, core::MakeAnIntDummy().a);
}

TEST(DummyTest, One) {
	EXPECT_EQ(PANORAMIX_VERSION_MINOR, core::MakeAnIntDummy().b);
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
	return 0;
}



