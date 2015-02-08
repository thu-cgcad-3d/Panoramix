#include "../src/core/cons_graph.hpp"
#include "../src/core/optimization.hpp"
#include "../src/vis/visualize2d.hpp"

#include "config.hpp"

using namespace panoramix;
using namespace test;

using CG = core::ConstraintGraph<
        std::tuple<int, double, std::string, core::Vec3>, 
        std::tuple<
            core::ConstraintConfig<std::string, int, double, std::string, core::Vec3>
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
    auto h2 = cg.addComponent(std::string("2_string"));
    auto h3 = cg.addComponent(3.0);
    auto h4 = cg.addComponent(4.0);
    auto h5 = cg.addComponent(core::Vec3(1, 2, 3));
    auto h6 = cg.addComponent(core::Vec3(2, 3, 4));
    auto h7 = cg.addComponent(core::Vec3(3, 4, 5));

    auto c = cg.addConstraint(std::string("a constraint"), 
        h1, h3, h2, h5);

    core::IterateOver(std::make_pair(cg.allComponents(), cg.allConstraints()), Print());

    ASSERT_TRUE(!cg.removed(c));
    cg.remove(h1);    
    ASSERT_TRUE(cg.removed(c));

    auto cg2 = cg;
    cg2.merge(cg);

}




