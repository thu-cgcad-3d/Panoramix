#include "../src/core/version.hpp"
#include "../src/core/basic_types.hpp"
#include "../src/vis/qttest.hpp"
#include "gtest/gtest.h"

#include <iostream>
#include <random>
#include <QApplication>

using namespace panoramix;

TEST(ConfigTest, Version) {
	EXPECT_EQ(PANORAMIX_VERSION_MAJOR, core::GetVersion().major);
	EXPECT_EQ(PANORAMIX_VERSION_MINOR, core::GetVersion().minor);
}

TEST(Qt, Basic) {
	vis::TestWidget * w = new vis::TestWidget;
	w->resize(400, 400);
	w->show();
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
	QApplication app(argc, argv);
    RUN_ALL_TESTS();
    return app.exec();
}



