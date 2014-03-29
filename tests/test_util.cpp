#include "../src/core/utilities.hpp"

#include <random>
#include <list>

#include <Eigen/Core>

#include "gtest/gtest.h"

using namespace panoramix;

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

TEST(UtilTest, AngleBetweenDirections) {
    core::Vec2 v1(1, 0), v2(1, 1);
    ASSERT_DOUBLE_EQ(M_PI_4, core::AngleBetweenDirections(v1, v2));
    ASSERT_DOUBLE_EQ(M_PI_4, core::SignedAngleBetweenDirections(v1, v2, false));
    ASSERT_DOUBLE_EQ(-M_PI_4, core::SignedAngleBetweenDirections(v1, v2, true));
    ASSERT_DOUBLE_EQ(-M_PI_4, core::SignedAngleBetweenDirections(v1, v2));
    core::Vec2 v3(-1, -1);
    ASSERT_DOUBLE_EQ(-M_PI_4 * 3, core::SignedAngleBetweenDirections(v1, v3, false));
    ASSERT_DOUBLE_EQ(M_PI_4 * 3, core::SignedAngleBetweenDirections(v1, v3, true));
}

TEST(UtilTest, DistanceFromPointToLine) {
    core::Line3 l;
    l.first = { 1, 0, 0 };
    l.second = { -1, 0, 0 };
    for (double x = -3; x <= 3; x += 0.5) {
        core::Point3 p(x, 1, 0);
        if (x < -1){
            ASSERT_DOUBLE_EQ(core::norm(l.second - p), core::DistanceFromPointToLine(p, l).first);
        }
        else if (x > 1){
            ASSERT_DOUBLE_EQ(core::norm(l.first - p), core::DistanceFromPointToLine(p, l).first);
        }
        else{
            ASSERT_DOUBLE_EQ(1, core::DistanceFromPointToLine(p, l).first);
        }
    }
}


TEST(UtilTest, MergeNear) {
    std::list<double> arr1;
    arr1.resize(1000);
    std::generate(arr1.begin(), arr1.end(), std::rand);
    std::vector<double> arr2(arr1.begin(), arr1.end());

    double thres = 10;
    auto gBegins1 = core::MergeNear(std::begin(arr1), std::end(arr1), std::false_type(), thres);
    auto gBegins2 = core::MergeNear(std::begin(arr2), std::end(arr2), std::true_type(), thres);
    ASSERT_EQ(gBegins1.size(), gBegins2.size());
    auto i = gBegins1.begin();
    auto j = gBegins2.begin();
    for (; i != gBegins1.end(); ++i, ++j){
        EXPECT_EQ(**i, **j);
    }
    for (auto i = gBegins2.begin(); i != gBegins2.end(); ++i){
        auto inext = std::next(i);
        auto begin = *i;
        auto end = inext == gBegins2.end() ? std::end(arr2) : *inext;
        auto beginVal = *begin;
        for (auto j = begin; j != end; ++j){
            EXPECT_NEAR(*j, beginVal, thres);
        }
    }
}


int main(int argc, char * argv[], char * envp[])
{
    testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}