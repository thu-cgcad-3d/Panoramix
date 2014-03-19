#include "../src/core/feature.hpp"
#include "gtest/gtest.h"

#include <iostream>
#include <string>
#include <random>

using namespace panoramix;

// PROJECT_TEST_DATA_DIR_STR is predefined using CMake
static const std::string ProjectTestDataDirStr = PROJECT_TEST_DATA_DIR_STR;

TEST(Feature, Camera){
	core::Camera cam(1000, 1000, 500, core::Camera::Vec3(0, 0, 0), core::Camera::Vec3(1, 0, 0));

	for (int i = 0; i < 100; i++){
		core::Camera::Vec2 v(abs(rand()), abs(rand()));
		auto p = cam.spatialDirection(v);
		auto v2 = cam.screenProjection(p);
		double dist = (v - v2).norm();
		ASSERT_LT(dist, 0.01);
	}
	auto c = cam.screenProjection(cam.center());
	double dist = (c - core::Camera::Vec2(cam.screenSize().width / 2, cam.screenSize().height / 2)).norm();
	ASSERT_LT(dist, 2);
}

TEST(Feature, CameraRandom){
	for (int k = 0; k < 100; k++) {
		core::Camera cam(abs(rand()) % 500, abs(rand()) % 400, abs(rand()) % 600,
			core::Camera::Vec3(rand(), rand(), rand()),
			core::Camera::Vec3(rand(), rand(), rand()));
		for (int i = 0; i < 100; i++){
			core::Camera::Vec2 v(rand(), rand());
			auto p = cam.spatialDirection(v);
			auto v2 = cam.screenProjection(p);
			double dist = (v - v2).norm();
			if (!std::isnan(dist) && !std::isinf(dist))
				ASSERT_LT(dist, 0.01);
		}
		auto c = cam.screenProjection(cam.center());
		double dist = (c - core::Camera::Vec2(cam.screenSize().width / 2, cam.screenSize().height / 2)).norm();
		if (!std::isnan(dist) && !std::isinf(dist))
			ASSERT_LT(dist, 2);
	}
}


TEST(Feature, CameraSampler) {
	cv::Mat im = cv::imread(ProjectTestDataDirStr + "/" + "panofactory.jpg");

	EXPECT_EQ(4000, im.cols);
	EXPECT_EQ(2000, im.rows);
	cv::resize(im, im, cv::Size(1000, 500));
	cv::imshow("panorama", im);

	core::PanoramicCamera originCam(im.cols / M_PI / 2.0);
	core::PanoramicCamera newCam(im.cols / M_PI / 2.0, 
		core::PanoramicCamera::Vec3(0, 0, 0), 
		core::PanoramicCamera::Vec3(0, 0, 1),
		core::PanoramicCamera::Vec3(0, 1, 0));
	cv::imshow("panorama ii", 
		core::CameraSampler<core::PanoramicCamera, core::PanoramicCamera>(newCam, originCam)(im));
	cv::waitKey();

	float camPositions[4][3] = {
		{1, 0, 0},
		{0, 1, 0},
		{-1, 0, 0},
		{0, -1, 0}
	};
	for (int i = 0; i < 4; i++){
		core::Camera cam(300, 800, 150, 
			core::Camera::Vec3(0, 0, 0),
			core::Camera::Vec3(camPositions[i][0], camPositions[i][1], camPositions[i][2]),
			core::Camera::Vec3(0, 0, -1));
		core::CameraSampler<core::Camera, core::PanoramicCamera> sampler(cam, 
			core::PanoramicCamera(im.cols / M_PI / 2.0));
		cv::Mat sampledIm = sampler(im);
		cv::imshow("sampled " + std::to_string(i), sampledIm);
		cv::waitKey();
	}
}


int main(int argc, char * argv[], char * envp[])
{
	testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}
