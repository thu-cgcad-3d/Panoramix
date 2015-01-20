#include "../src/core/cons_graph.hpp"
#include "../src/core/optimization.hpp"
#include "../src/core/p3graph.hpp"
#include "../src/vis/visualize2d.hpp"

#include "config.hpp"

using namespace panoramix;
using namespace test;

using CG = core::ConstraintGraph<
        std::tuple<int, double, std::string, core::Vec3>, 
        std::tuple<
            core::ConstraintConfig<std::string, 
                core::ComponentOccupation<int, 2>, 
                core::ComponentOccupation<double, core::Dynamic>,
                core::ComponentOccupation<std::string, 1>,
                core::ComponentOccupation<core::Vec3, 3>
            >
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
        core::Depends(h1),
        core::Depends(h3, h4),
        core::Depends(h2),
        core::Depends(h5, h6, h7));

    core::IterateOver(std::make_pair(cg.allComponents(), cg.allConstraints()), Print());

    ASSERT_TRUE(!cg.removed(c));
    cg.remove(h1);    
    ASSERT_TRUE(cg.removed(c));
}


//TEST(ConstraintGraph, P3Graph) {
//
//    core::P3Graph p3g;
//    auto h1 = p3g.addComponent(core::P3Line());
//    auto h2 = p3g.addComponent(core::P3Region());
//
//    core::P3Environment env;
//
//    auto varReg = core::RegisterComponents<double>(p3g, env);
//    for (int i = 0; i < 3; i++)
//        varReg(h2, i) = i+1;
//
//    core::UpdateComponents(p3g, varReg, env);
//    auto & reg = p3g.data(h2);
//    ASSERT_TRUE(reg.plane.component == core::Plane3FromEquation(1.0, 2.0, 3.0));
//
//    auto consReg = core::RegisterConstraints(p3g, varReg, env);
//
//}




