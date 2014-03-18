#include "../src/core/feature.hpp"
#include "gtest/gtest.h"

#include <iostream>
#include <string>
#include <random>

using namespace panoramix;

// PROJECT_TEST_DATA_DIR_STR is predefined using CMake
static const std::string ProjectTestDataDirStr = PROJECT_TEST_DATA_DIR_STR;

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
	core::Camera cam(1000, 1000, 500, core::Camera::Vec3(1, 0, 0));
	for (int i = 0; i < 100; i++){
		core::Camera::Vec2 v(rand(), rand());
		auto p = cam.spatialDirection(v);
		auto v2 = cam.screenProjection(p);
		double dist = (v - v2).norm();
		EXPECT_LT(dist, 0.01);
	}
}


TEST(Feature, PerspectiveSampler) {
	cv::Mat im = cv::imread(ProjectTestDataDirStr + "/" + "panofactory.jpg");

	EXPECT_EQ(4000, im.cols);
	EXPECT_EQ(2000, im.rows);
	cv::resize(im, im, cv::Size(1000, 500));
	cv::imshow("panorama", im);

	float camPositions[4][3] = {
		{1, 0, 0},
		{0, 1, 0},
		{-1, 0, 0},
		{0, -1, 0}
	};
	for (int i = 0; i < 4; i++){
		core::Camera cam(200, 800, 30, 
			core::Camera::Vec3(0, 0, 0),
			core::Camera::Vec3(camPositions[i][0], camPositions[i][1], camPositions[i][2]),
			core::Camera::Vec3(0, 0, 1));
		core::PerspectiveSampler sampler(cam);
		cv::Mat sampledIm = sampler(im);
		cv::imshow("sampled " + std::to_string(i), sampledIm);
	}

	cv::waitKey();
}


int main(int argc, char * argv[], char * envp[])
{
	testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}
