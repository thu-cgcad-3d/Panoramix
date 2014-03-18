#include "../src/core/feature.hpp"
#include "gtest/gtest.h"

#include <iostream>
#include <string>
#include <random>

using namespace panoramix;

// PROJECT_TEST_DATA_DIR_STR is predefined using CMake
static const std::string ProjectTestDataDirStr = PROJECT_TEST_DATA_DIR_STR;

TEST(Feature, ExtractLines) {
	cv::Mat im = cv::imread(ProjectTestDataDirStr + "/" + "panofactory.jpg");
	
    EXPECT_EQ(4000, im.cols);
	EXPECT_EQ(2000, im.rows);
	
    //std::list<core::LineData<double, 2>> lines;
    //core::ExtractLines(im, lines, 150, 10, 20, 8);
    
    //EXPECT_EQ(254, lines.size());
}

TEST(Feature, Camera){
	for (int k = 0; k < 100; k++) {
		core::Camera cam(rand() % 500, rand() % 400, rand() % 600, 
			core::Camera::Vec3(rand(), rand(), rand()), 
			core::Camera::Vec3(rand(), rand(), rand()));
		for (int i = 0; i < 100; i++){
			core::Camera::Vec2 v(rand(), rand());
			auto p = cam.spatialDirection(v);
			auto v2 = cam.screenProjection(p);
			double dist = (v - v2).norm();
			if (!std::isnan(dist))
				EXPECT_LT(dist, 0.01);
		}
	}
}


int main(int argc, char * argv[], char * envp[])
{
	testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
	return 0;
}
