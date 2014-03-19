#include "../src/core/feature.hpp"
#include "../src/core/feature_visualize.hpp"
#include "gtest/gtest.h"

#include <iostream>
#include <string>
#include <random>

using namespace panoramix;

// PROJECT_TEST_DATA_DIR_STR is predefined using CMake
static const std::string ProjectTestDataDirStr = PROJECT_TEST_DATA_DIR_STR;

TEST(Feature, PerspectiveCamera){
	core::PerspectiveCamera cam(1000, 1000, 500, 
		core::PerspectiveCamera::Vec3(0, 0, 0), 
		core::PerspectiveCamera::Vec3(1, 0, 0));

	for (int i = 0; i < 100; i++){
		core::PerspectiveCamera::Vec2 v(abs(rand()), abs(rand()));
		auto p = cam.spatialDirection(v);
		auto v2 = cam.screenProjection(p);
		double dist = (v - v2).norm();
		ASSERT_LT(dist, 0.01);
	}
	auto c = cam.screenProjection(cam.center());
	double dist = (c - core::PerspectiveCamera::Vec2(cam.screenSize().width / 2, cam.screenSize().height / 2)).norm();
	ASSERT_LT(dist, 2);
}

TEST(Feature, PerspectiveCameraRandom){
	for (int k = 0; k < 100; k++) {
		core::PerspectiveCamera cam(abs(rand()) % 500, abs(rand()) % 400, abs(rand()) % 600,
			core::PerspectiveCamera::Vec3(rand(), rand(), rand()),
			core::PerspectiveCamera::Vec3(rand(), rand(), rand()));
		for (int i = 0; i < 100; i++){
			core::PerspectiveCamera::Vec2 v(rand(), rand());
			auto p = cam.spatialDirection(v);
			auto v2 = cam.screenProjection(p);
			double dist = (v - v2).norm();
			if (!std::isnan(dist) && !std::isinf(dist))
				ASSERT_LT(dist, 0.01);
		}
		auto c = cam.screenProjection(cam.center());
		double dist = (c - core::PerspectiveCamera::Vec2(cam.screenSize().width / 2, cam.screenSize().height / 2)).norm();
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
		core::PerspectiveCamera cam(500, 600, 150,
			core::PerspectiveCamera::Vec3(0, 0, 0),
			core::PerspectiveCamera::Vec3(camPositions[i][0], camPositions[i][1], camPositions[i][2]),
			core::PerspectiveCamera::Vec3(0, 0, -1));
		core::CameraSampler<core::PerspectiveCamera, core::PanoramicCamera> sampler(cam,
			core::PanoramicCamera(im.cols / M_PI / 2.0));
		cv::Mat sampledIm = sampler(im);
		cv::imshow("sampled_" + std::to_string(i), sampledIm);
		cv::waitKey();
	}
}


TEST(Feature, FeatureExtractor) {
	core::LineSegmentExtractor lineSegmentExtractor;
	core::CVFeatureExtractor<cv::SIFT> sift;
	core::CVFeatureExtractor<cv::SURF> surf(300.0);
	for (int i = 0; i < 4; i++) {
		std::string name = ProjectTestDataDirStr + "/" + "sampled_" + std::to_string(i) + ".png";
		cv::Mat im = cv::imread(name);
		core::ImageFeatureVisualizer viz(im, name);
		viz << lineSegmentExtractor(im);
		viz << sift(im);
		viz << surf(im);
	}
}


int main(int argc, char * argv[], char * envp[])
{
	testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}
