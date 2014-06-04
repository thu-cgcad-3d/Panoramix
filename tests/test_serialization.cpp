#include "../src/core/basic_types_serialization.hpp"

#include <random>

#include "gtest/gtest.h"

#include "test_config.hpp"

using namespace panoramix;
using namespace test;

double randf(){
    return (std::rand() % 100000) / 100000.0;
}

TEST(Serialization, BasicTypesOutput) {

    std::vector<core::Classified<std::pair<core::Line3, core::HPoint2>>> data1(1000);
    std::vector<core::Classified<std::pair<core::Line3, core::HPoint2>>> data1c;
    for (auto & d : data1) {
        d.claz = std::rand();
        d.component.first = core::Line3(core::Point3(randf(), randf(), randf()), core::Point3(randf(), randf(), randf()));
        d.component.second = core::HPoint2(core::Point2(randf(), randf()), randf());
    }
    std::map<std::string, std::pair<core::Box2, core::GeoCoord>> data2;
    std::map<std::string, std::pair<core::Box2, core::GeoCoord>> data2c;
    for (int i = 0; i < 1000; i++) {
        data2[std::to_string(i)] =
            std::make_pair(core::Box2(core::Point2(randf(), randf()), core::Point2(randf(), randf())),
            core::GeoCoord(randf(), randf()));
    }

    {
        std::ofstream os(ProjectTestDataDirStr_Serialization + "/data1data2.cereal", std::ios::binary);
        cereal::BinaryOutputArchive archive(os);
        archive(data1, data2);
    }
    {
        std::ifstream is(ProjectTestDataDirStr_Serialization + "/data1data2.cereal", std::ios::binary);
        cereal::BinaryInputArchive archive(is);
        archive(data1c, data2c);
    } 

    EXPECT_EQ(data1.size(), data1c.size());
    EXPECT_EQ(data2.size(), data2c.size());

    for (int i = 0; i < data1.size(); i++) {
        EXPECT_TRUE(data1[i] == data1c[i]);
    }

    for (int i = 0; i < 1000; i++) {
        EXPECT_TRUE(data2[std::to_string(i)] == data2c[std::to_string(i)]);
    }

}




int main(int argc, char * argv[], char * envp[])
{
    testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();

}