#include "../src/core/version.hpp"
#include "../src/core/basic_types.hpp"
#include "gtest/gtest.h"

//#include <armadillo>
//#include <armadillo_bits/config.hpp>

#include <iostream>
#include <random>

using namespace panoramix;

TEST(ConfigTest, Version) {
	EXPECT_EQ(PANORAMIX_VERSION_MAJOR, core::GetVersion().major);
	EXPECT_EQ(PANORAMIX_VERSION_MINOR, core::GetVersion().minor);
}

TEST(BasicTypeTest, HPoint) {
	for (int i = 0; i < 1000; i++){
		core::Vec4 v4;
		std::generate(v4.val, v4.val + 4, std::rand);
		auto hp = core::HPointFromVector(v4);
		auto p = hp.toVector();
		ASSERT_LT(core::norm(p - v4), 1e-5);
		core::HPoint<double, 4> hp5 = v4;
		auto p5 = hp5.toPoint();
		ASSERT_LT(core::norm(p5 - v4), 1e-5);
	}
}

//TEST(Arma, A) {
//
    //using namespace arma;
    //using namespace std;
    //// points to which we will fit the line
    //mat data = "1 6; 2 5; 3 7; 4 10";
//
    //mat aa = "1 2 3; 3 4 5";
    //mat bb = "2 3 4; 5 6 7";
    //auto ab = aa % bb;
    //cout << ab << endl;
//
    //cout << "Points used for the estimation:" << endl;
    //cout << data << endl;
//
    //// Build matrices to solve Ax = b problem:
    //vec b(data.n_rows);
    //mat C(data.n_rows, 2);
//
    //for (u32 i = 0; i<data.n_rows; ++i)
    //{
        //b(i) = data(i, 1);
//
        //C(i, 0) = 1;
        //C(i, 1) = data(i, 0);
    //}
//
    //cout << "b:" << endl;
    //cout << b << endl;
//
    //cout << "Constraint matrix:" << endl;
    //cout << C << endl;
//
    //// Compute least-squares solution:
    //vec solution = solve(C, b);
//
    //// solution should be "3.5; 1.4"
    //cout << "solution:" << endl;
    //cout << solution << endl;
//
//
    //cout << "Reprojection error:" << endl;
//
    //for (u32 i = 0; i<data.n_rows; ++i)
    //{
        //cout << "  residual: " << (data(i, 1) - (solution(0) + solution(1) * data(i, 0))) << endl;
    //}
//
//}



int main(int argc, char * argv[], char * envp[])
{
	for (int i = 0; i < argc; i++) {
		std::cout << "[INPUT]:" << argv[i] << std::endl;
	}
	char** env;
	for (env = envp; *env != 0; env++) {
		char* thisEnv = *env;
		std::cout << "[ENV]:" << thisEnv << std::endl;
	}
	testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}



