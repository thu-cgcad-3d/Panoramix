#include "gtest/gtest.h"

int main(int argc, char * argv[], char * envp[]) {
   testing::InitGoogleTest(&argc, argv);
   testing::GTEST_FLAG(catch_exceptions) = false;
   testing::GTEST_FLAG(throw_on_failure) = true;
   testing::GTEST_FLAG(filter) = "MixedGraph.Test";
   return RUN_ALL_TESTS();
}
