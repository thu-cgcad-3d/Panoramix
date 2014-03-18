#include "../src/core/version.hpp"
#include "gtest/gtest.h"

#include <iostream>

using namespace panoramix;

TEST(ConfigTest, Version) {
	EXPECT_EQ(PANORAMIX_VERSION_MAJOR, core::GetVersion().major);
	EXPECT_EQ(PANORAMIX_VERSION_MINOR, core::GetVersion().minor);
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



