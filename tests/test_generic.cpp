#include "../src/core/generic.hpp"

#include <random>
#include <list>

#include "gtest/gtest.h"

using namespace panoramix;


TEST(BoundingBoxTest, Box) {
    
    using namespace core;
    Line3 l1(Point3(0.5, 0.1, 1), Point3(1, 0.4, 0.7));
    Line3 l2(Point3(0.6, 1, 0.9), Point3(0.2, -1, 0.5));

    Line3 lines[] = { l1, l2 };
    auto box = BoundingBoxOfContainer(lines);

    ASSERT_EQ(0.2, box.minCorner[0]);
    ASSERT_EQ(1, box.maxCorner[0]);

    ASSERT_EQ(-1, box.minCorner[1]);
    ASSERT_EQ(1, box.maxCorner[1]);

    ASSERT_EQ(0.5, box.minCorner[2]);
    ASSERT_EQ(1, box.maxCorner[2]);

}



int main(int argc, char * argv[], char * envp[])
{
    testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}