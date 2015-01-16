#include "../src/core/cons_graph.hpp"
#include "../src/vis/visualize2d.hpp"

#include "config.hpp"

using namespace panoramix;
using namespace test;





TEST(ConstraintGraph, Basic){

    core::ConstraintGraph<
        std::tuple<int, double, std::string>, 
        std::tuple<
            core::ConstraintConfig<std::string, 
                core::ComponentOccupation<int, 2>, 
                core::ComponentOccupation<double, core::Dynamic>
            >
        >
    > cg;

    const double a = 0.0;
    auto h1 = cg.addComponent(12);
    auto h2 = cg.addComponent(std::string("hahaha"));
    auto h3 = cg.addComponent(10.0);
    auto h4 = cg.addComponent(a);

    auto aa = core::Depends<2>(h1, h1);

    auto c = cg.addConstraint(std::string("a constraint"), core::Depends<2>(h1, h1), core::Depends<-1>(h3, h4));
    
    cg.topo(h1);
    cg.topo(c);
    //ASSERT_EQ(cg.topo(c).components<int>().size(), 1);

    std::cout << cg.data(h1) << ";" << cg.data(h2) << std::endl;

}