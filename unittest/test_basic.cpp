#include "../src/core/version.hpp"
#include "../src/core/basic_types.hpp"
#include "../src/core/ring.hpp"
#include "../src/core/utility.hpp"
#include "config.hpp"

#include <iostream>
#include <random>
#include <Eigen/Dense>

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

TEST(BasicType, Line) {

    std::vector<core::Line2> lines = {
        { core::Point2(1, 2), core::Point2(3, 4) },
        { core::Point2(5, 6), core::Point2(7, 8) },
        { core::Point2(9, 10), core::Point2(11, 12) },
        { core::Point2(13, 14), core::Point2(15, 16) },
        { core::Point2(17, 18), core::Point2(19, 20) },
        { core::Point2(21, 22), core::Point2(23, 24) }
    };

    Eigen::Map<Eigen::Matrix<double, 4, Eigen::Dynamic>> linesData((double*)(lines.data()), 4, lines.size());
    std::cout << linesData << std::endl;

}

TEST(BasicType, VecCast) {

    core::Imaged im(100, 100, 0.0);
    auto imi = core::vec_cast<int>(im);

    for (auto & i : imi){
        ASSERT_EQ(i, 0);
    }

    core::Image3d im3(100, 100, core::Vec3(1, 2, 3));
    auto imi3 = core::vec_cast<int>(im3);
    auto imd3 = core::vec_cast<double>(im3);

    for (auto & i : imi3){
        ASSERT_TRUE(i == core::Vec3i(1, 2, 3));
    }
    for (auto & i : imd3){
        ASSERT_TRUE(i == core::Vec3(1, 2, 3));
    }


}

TEST(BasicType, Ring) {

    //core::Radian r = M_PI;
    //double rr = r;
    //ASSERT_FLOAT_EQ(rr, M_PI);

    //r += M_PI;
    //double rr2 = r;
    //ASSERT_FLOAT_EQ(rr2, 0);

    core::Radian rrr = M_PI * 2;
    double rrr2 = rrr;
    ASSERT_FLOAT_EQ(rrr2, 0);

}

TEST(BasicType, RingPerformance){
    auto r = core::Radian::toRep(M_PI);
    for (uint64_t i = 0; i < 1e8; i++){
        r = r + core::Radian::toRep(M_PI) * 2.5 + core::Radian::toRep(1.0);
    }
    std::cout << core::Radian::toValue(r) << std::endl;
}

TEST(BasicType, RingPerformance2){
    core::Radian r = M_PI;
    for (uint64_t i = 0; i < 1e8; i++){
        r = r + M_PI * 2.5 + 1.0;
    }
    std::cout << double(r) << std::endl;
}

TEST(BasicType, RingPerformanceBaseline){
    double r = M_PI;
    for (uint64_t i = 0; i < 1e8; i++){
        r = r + M_PI * 2.5 + 1.0;
        r = core::WrapBetween(r, 0.0, M_PI * 2);
    }
    std::cout << double(r) << std::endl;
}


