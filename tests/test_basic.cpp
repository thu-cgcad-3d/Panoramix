#include "../src/core/version.hpp"
#include "../src/core/basic_types.hpp"
#include "gtest/gtest.h"

#include <iostream>
#include <random>

using namespace panoramix;

TEST(ConfigTest, Version) {
    EXPECT_EQ(PANORAMIX_VERSION_MAJOR, core::GetVersion().major);
    EXPECT_EQ(PANORAMIX_VERSION_MINOR, core::GetVersion().minor);
}

TEST(BasicType, Vec) {
    {
        core::Vec3 v1(0, 0, 0), v2(1, 1, 1);
        core::Vec3 v12 = v1 + v2 * 3.0;
        ASSERT_TRUE(v12 == core::Vec3(3, 3, 3));
    }
    {
        core::Vec4 v1(0, 0, 0, 0), v2(1, 1, 1, 1);
        core::Vec4 v12 = v1 + v2 * 3.0;
        ASSERT_TRUE(v12 == core::Vec4(3, 3, 3, 3));
    }
}

TEST(BasicType, HPoint) {
    for (int i = 0; i < 1000; i++){
        core::Vec4 v4;
        std::generate(v4.val, v4.val + 4, std::rand);
        auto hp = core::HPointFromVector(v4);
        auto p = VectorFromHPoint(hp);
        ASSERT_LT(core::norm(p - v4), 1e-5);
        core::HPoint<double, 4> hp5 = v4;
        auto p5 = hp5.value();
        ASSERT_LT(core::norm(p5 - v4), 1e-5);
    }
}


