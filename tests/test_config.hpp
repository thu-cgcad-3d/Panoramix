#ifndef TEST_CONFIG_HPP
#define TEST_CONFIG_HPP
 
#include <iostream>
#include <chrono>
#include <random>

#include "gtest/gtest.h"

namespace test {

	static const std::string ProjectTestDataDirStr = PROJECT_TEST_DATA_DIR_STR;
	static const std::string ProjectTestDataDirStr_Normal = ProjectTestDataDirStr + "/normal";
	static const std::string ProjectTestDataDirStr_PanoramaIndoor = ProjectTestDataDirStr + "/panorama/indoor";
	static const std::string ProjectTestDataDirStr_PanoramaOutdoor = ProjectTestDataDirStr + "/panorama/outdoor";
    static const std::string ProjectTestDataDirStr_Serialization = ProjectTestDataDirStr + "/serialization";
 
}

#endif