#include "../src/core/optimization.hpp"

#include "config.hpp"

#include <iostream>
#include <random>

using namespace panoramix;

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

//using VarReg = core::VarRegistrationFromConstraintGraph<double, CG>::type;

TEST(Optimization, Basic){

    

}