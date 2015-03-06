#include "../src/core/cons_graph.hpp"
#include "../src/core/optimization.hpp"
#include "../src/vis/visualize2d.hpp"

#include "config.hpp"

using namespace panoramix;
using namespace test;

using CG = core::ConstraintGraph<
        std::tuple<bool, int, double, float, core::Vec3>, 
        std::tuple<
            core::ConstraintConfig<std::string, int, double>,
            core::ConstraintConfig<long long, bool, float, core::Vec3>
        >
    >;

struct Print {
    template <class T>
    inline void operator()(const T & t) const {
        if (t.exists)
            std::cout << t.data << std::endl;
    }
};


TEST(ConstraintGraph, Basic){

    CG cg;

    auto h1 = cg.addComponent(1);
    auto h2 = cg.addComponent(2.0);
    auto h3 = cg.addComponent(true);
    auto h4 = cg.addComponent(4.0f);
    auto h5 = cg.addComponent(5);
    auto h6 = cg.addComponent(core::Vec3(1, 2, 3));
    auto h7 = cg.addComponent(core::Vec3(2, 3, 4));
    auto h8 = cg.addComponent(core::Vec3(3, 4, 5));
    auto h9 = cg.addComponent(9.0);

    auto c1 = cg.addConstraint(std::string("a constraint"), h1, h2);
    auto c2 = cg.addConstraint(64ll, h3, h4, h6);
    auto c3 = cg.addConstraint(std::string("another constraint"), h1, h9);

    core::IterateOver(std::make_pair(cg.allComponents(), cg.allConstraints()), Print());

    cg.remove(h1);
    ASSERT_TRUE(cg.removed(c1));
    ASSERT_FALSE(cg.removed(c2));
    ASSERT_TRUE(cg.removed(c3));

    auto cg2 = cg;
    cg2.merge(cg);

}




