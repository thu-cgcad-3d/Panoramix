#include "../src/core/feature.hpp"
#include "gtest/gtest.h"

#include <iostream>
#include <string>

using namespace panoramix;

static const std::string ProjectTestDataDirStr = PROJECT_TEST_DATA_DIR_STR;

TEST(Feature, ExtractLines) {
	cv::Mat im = cv::imread(ProjectTestDataDirStr + "/" + "panofactory.jpg");
	
    EXPECT_EQ(4000, im.cols);
	EXPECT_EQ(2000, im.rows);
	
    std::list<core::LineData<double, 2>> lines;
    core::ExtractLines(im, lines, 150);
    
    EXPECT_GT(lines.size(), 0);
    
    
    
}


int main(int argc, char * argv[], char * envp[])
{
	testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
	return 0;
}
