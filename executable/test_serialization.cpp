#include "../src/core/basic_types.hpp"
#include "../src/core/graphical_model.hpp"

#include <random>

#include "gtest/gtest.h"

#include "test_config.hpp"

using namespace panoramix;
using namespace test;

double randf(){
    return (std::rand() % 100000) / 100000.0;
}

TEST(Serialization, BasicTypes) {

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

TEST(Serialization, GraphicalModel) {

    using CGraph = core::GraphicalModel<core::Dummy, core::LayerConfig<core::Dummy, core::Dynamic>, core::LayerConfig<core::Point3, 4>>;

    CGraph cgraph, cgraph2;
    auto c0 = cgraph.add();
    auto c1 = cgraph.add();
    auto c2 = cgraph.add();
    auto c3 = cgraph.add();

    auto cc012 = cgraph.add<1>({ c0, c1, c2 });
    auto cc123 = cgraph.add<1>({ c1, c2, c3 });
    auto cc230 = cgraph.add<1>({ c2, c3, c0 });
    auto cc301 = cgraph.add<1>({ c3, c0, c1 });

    auto ccc = cgraph.add<2>({ cc012, cc123, cc230, cc301 }, core::Point3(1, 2, 3));

    {
        std::ofstream os(ProjectTestDataDirStr_Serialization + "/gm.cereal", std::ios::binary);
        cereal::BinaryOutputArchive archive(os);
        archive(cgraph);
    }
    {
        std::ifstream is(ProjectTestDataDirStr_Serialization + "/gm.cereal", std::ios::binary);
        cereal::BinaryInputArchive archive(is);
        archive(cgraph2);
    }

    EXPECT_EQ(cgraph.internalElements<0>().size(), cgraph2.internalElements<0>().size());
    EXPECT_EQ(cgraph.internalElements<1>().size(), cgraph2.internalElements<1>().size());
    EXPECT_EQ(cgraph.internalElements<2>().size(), cgraph2.internalElements<2>().size());

    EXPECT_TRUE(cgraph.data(ccc) == cgraph2.data(ccc));

}

TEST(Serialization, FileTime) {
    
    core::Line3 testEntity;
    {
        std::ofstream os(ProjectTestDataDirStr_Serialization + "/line.cereal", std::ios::binary);
        cereal::BinaryOutputArchive archive(os);
        archive(testEntity);
    }

    auto t1 = core::LastModifiedTimeOfFile(ProjectTestDataDirStr_Serialization + "/gm.cereal");
    auto t2 = core::LastModifiedTimeOfFile(ProjectTestDataDirStr_Serialization + "/line.cereal");
    auto t3 = core::CurrentTime();

    ASSERT_LE(t1, t2);
    ASSERT_LE(t2, t3);

}





int main(int argc, char * argv[], char * envp[])
{
    testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();

}