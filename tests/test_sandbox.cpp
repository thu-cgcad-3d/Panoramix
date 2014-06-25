#include <type_traits>

#include "../src/sandbox/for_expression.hpp"
#include "../src/sandbox/for_expression2.hpp"

#include "../src/deriv/data_traits_definitions.hpp"

#include "gtest/gtest.h"

struct B {
    int v;
};

struct A {
//private:
    inline A & operator = (const A & a) {
        v = a.v;
        return *this;
    }
    int v;
};

TEST(Concept, Concepts) {
    A a;
    A aa;
    //aa = a;
    ASSERT_TRUE((std::is_trivially_assignable<A&, A&&>::value));
    
}


int main(int argc, char * argv[], char * envp[])
{
    testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}