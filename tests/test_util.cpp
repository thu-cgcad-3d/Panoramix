#include "../src/core/util.hpp"

#include <random>

#include <Eigen/Core>

#include "gtest/gtest.h"

using namespace panoramix;
using Mat4 = Eigen::Matrix4d;
using Vec3 = Eigen::Vector3d;
using Vec4 = Eigen::Vector4d;

TEST(UtilTest, WrapBetween) {
	for (int i = 0; i < 10000; i++){
		double x = double(rand()) / rand() + rand();
		double a = double(rand()) / rand() + rand();
		double b = a + abs(double(rand())/rand());
		if (a == b || std::isnan(x) || std::isnan(a) || std::isnan(b) || 
			std::isinf(x) || std::isinf(a) || std::isinf(b))
			continue;
		double xx = core::WrapBetween(x, a, b);
		double rem = (xx - x) / (b - a) - std::round((xx - x) / (b - a));
		if (std::isnan(rem)){
			assert(0);
		}
		ASSERT_NEAR(0, rem, 1e-5);
		ASSERT_LE(a, xx);
		ASSERT_LT(xx, b);
	}
}

TEST(UtilTest, MatrixLookAt) {
	Mat4 m;
	m.setIdentity();
	m = core::Matrix4MakeLookAt(Vec3(0, 0, 0), Vec3(1, 0, 0), Vec3(0, 0, 1), m);
	Vec4 p = m * Vec4(1, 0, 0, 1);
	Vec3 pj = Vec3(p(0), p(1), p(2)) / p(3);
	ASSERT_LT((pj - Vec3(0, 0, 1)).norm(), 2);
}


int main(int argc, char * argv[], char * envp[])
{
	testing::InitGoogleTest(&argc, argv);
	return RUN_ALL_TESTS();
}