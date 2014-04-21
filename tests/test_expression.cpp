#include "../src/core/basic_types.hpp"
#include "../src/core/utilities.hpp"
#include "../src/core/expression.hpp"
#include "gtest/gtest.h"

#include <iostream>
#include <random>
#include <chrono>

using namespace panoramix::core;

TEST(Expression, Shape) {

    Shape<1, 2, Dynamic, 3, Dynamic, Dynamic> s(5, 6, 7);

    EXPECT_EQ(6, s.Rank);
    EXPECT_EQ(3, s.DynamicNum);
    
    // sizes
    EXPECT_EQ(1, s.size<0>());
    EXPECT_EQ(2, s.size<1>());
    EXPECT_EQ(5, s.size<2>());
    EXPECT_EQ(3, s.size<3>());
    EXPECT_EQ(6, s.size<4>());
    EXPECT_EQ(7, s.size<5>());

    // resize
    s.resize<2>(500);
    s.resize<4>(600);
    s.resize<5>(700);

    EXPECT_EQ(500, s.size<2>());
    EXPECT_EQ(600, s.size<4>());
    EXPECT_EQ(700, s.size<5>());

    s.resize<2>(5);
    s.resize<4>(6);
    s.resize<5>(7);

    // dynamic sizes
    EXPECT_EQ(5, s.dsize<0>());
    EXPECT_EQ(6, s.dsize<1>());
    EXPECT_EQ(7, s.dsize<2>());

    // dynamic resize
    s.dresize<0>(50);
    s.dresize<1>(60);
    s.dresize<2>(70);

    EXPECT_EQ(50, s.dsize<0>());
    EXPECT_EQ(60, s.dsize<1>());
    EXPECT_EQ(70, s.dsize<2>());

    s.dresize<0>(5);
    s.dresize<1>(6);
    s.dresize<2>(7);

    // static sizes
    EXPECT_EQ(1, s.ssize<0>());
    EXPECT_EQ(2, s.ssize<1>());
    EXPECT_EQ(3, s.ssize<2>());

    // volume
    EXPECT_EQ(1 * 2 * 5 * 3 * 6 * 7, s.volume());

    // shape to array
    auto arr = MakeArrayFromShape(s);
    EXPECT_EQ(1, arr[0]);
    EXPECT_EQ(2, arr[1]);
    EXPECT_EQ(5, arr[2]);
    EXPECT_EQ(3, arr[3]);
    EXPECT_EQ(6, arr[4]);
    EXPECT_EQ(7, arr[5]);

    auto tp = MakeTupleFromShape(s);
    EXPECT_TRUE(std::make_tuple(1, 2, 5, 3, 6, 7) == tp);

    // sub2ind & ind2sub
    EXPECT_EQ(0, IndexFromSubs(s, 0, 0, 0, 0, 0, 0));
    EXPECT_EQ(1, IndexFromSubs(s, 0, 1, 0, 0, 0, 0));

    int correctInd = 0;
    for (int f = 0; f < s.size<5>(); f++){
        for (int e = 0; e < s.size<4>(); e++){
            for (int d = 0; d < s.size<3>(); d++){
                for (int c = 0; c < s.size<2>(); c++){
                    for (int b = 0; b < s.size<1>(); b++){
                        for (int a = 0; a < s.size<0>(); a++){
                            int ind = IndexFromSubs(s, a, b, c, d, e, f);
                            int ind2 = IndexFromSubs(s, std::array<int, 6>{{ a, b, c, d, e, f }});
                            ASSERT_EQ(correctInd, ind);
                            ASSERT_EQ(correctInd, ind2);
                            correctInd++;

                            int aa, bb, cc, dd, ee, ff;
                            SubsFromIndex(s, ind, aa, bb, cc, dd, ee, ff);
                            std::array<int, 6> subs;
                            SubsFromIndex(s, ind, subs);
                            ASSERT_TRUE(std::tie(a, b, c, d, e, f) == std::tie(aa, bb, cc, dd, ee, ff));
                            ASSERT_TRUE(std::tie(a, b, c, d, e, f) == std::tie(subs[0], subs[1], subs[2], subs[3], subs[4], subs[5]));
                        }
                    }
                }
            }
        }
    }
   


    // make shape from dynamic sizes of a existing shape
    auto s2 = MakeShapeFromDynamicSizes<Dynamic, Dynamic, Dynamic>(s);

    EXPECT_EQ(3, s2.Rank);
    EXPECT_EQ(3, s2.DynamicNum);
    EXPECT_EQ(5, s2.size<0>());
    EXPECT_EQ(6, s2.size<1>());
    EXPECT_EQ(7, s2.size<2>());
    EXPECT_EQ(5 * 6 * 7, s2.volume());

    // make descartes product shape
    auto s_s2 = MakeDescartesProductShape(s, s2);

    EXPECT_EQ(s.Rank + s2.Rank, s_s2.Rank);
    EXPECT_EQ(s.DynamicNum + s2.DynamicNum, s_s2.DynamicNum);

    EXPECT_EQ(1, s_s2.size<0>());
    EXPECT_EQ(2, s_s2.size<1>());
    EXPECT_EQ(5, s_s2.size<2>());
    EXPECT_EQ(3, s_s2.size<3>());
    EXPECT_EQ(6, s_s2.size<4>());
    EXPECT_EQ(7, s_s2.size<5>());
    EXPECT_EQ(5, s_s2.size<6>());
    EXPECT_EQ(6, s_s2.size<7>());
    EXPECT_EQ(7, s_s2.size<8>());

    // make matrix product shape
    Shape<5, 6, Dynamic> aa(10);
    Shape<6, 2, 4> bb;
    auto aabb = MakeMatrixProductShape(aa, bb);

    EXPECT_EQ(5, aabb.size<0>());
    EXPECT_EQ(2, aabb.size<1>());
    EXPECT_EQ(10, aabb.size<2>());
    EXPECT_EQ(4, aabb.size<3>());


    // directly instantiate a static shape
    Shape<5, 6, 7> s3;

    EXPECT_EQ(3, s3.Rank);
    EXPECT_EQ(0, s3.DynamicNum);
    EXPECT_EQ(5, s3.size<0>());
    EXPECT_EQ(6, s3.size<1>());
    EXPECT_EQ(7, s3.size<2>());
    EXPECT_EQ(5 * 6 * 7, s3.volume());

    // copy
    auto ss = s;

    EXPECT_EQ(6, ss.Rank);
    EXPECT_EQ(3, ss.DynamicNum);
    EXPECT_EQ(1, ss.size<0>());
    EXPECT_EQ(2, ss.size<1>());
    EXPECT_EQ(5, ss.size<2>());
    EXPECT_EQ(3, ss.size<3>());
    EXPECT_EQ(6, ss.size<4>());
    EXPECT_EQ(7, ss.size<5>());
    EXPECT_EQ(1 * 2 * 5 * 3 * 6 * 7, ss.volume());    

}

TEST(Expression, Basic) {

    ExpressionGraph<double> graph;
    auto e1 = graph.addExpression<Const<double>>(Shape<2, 3>(), {}, 5);
    auto e2 = graph.addExpression<Const<double>>(Shape<2, 3>(), {}, 10);
    auto e3 = graph.addExpression<Const<double>>(Shape<2, 3>(), {}, 30);
    auto e123 = graph.addExpression<Plus<double>>(Shape<2, Dynamic>(3), {e1.handle, e2.handle, e3.handle});
    
    EXPECT_EQ(45, e123.eval(0, 1));

}


TEST(Expression, NN) {   
    
    //Expr<double, Dynamic> in(10);
    //Expr<double, Dynamic, Dynamic> w(10, 10);
    
    //Expr<double, Dynamic, 1> p = MatrixMult(w, in.promote());
    //Expr<double, Dynamic, Dynamic> j = JacobianMatrix(p, w);



    /*
    Expr<double, Dynamic> nnOutputs = 1.0 / (1.0 + exp(-nnDirectOutputs));
    
    Expr<double, Dynamic> expectedOutputs;

    Expr<double> squaredSumError = SumAll(Square(expectedOutputs - nnOutputs));
    Expr<double, Dynamic, Dynamic> nnDeltaWeights = JacobianMatrix(squaredSumError, nnWeights).reshape(nnWeights.shape()) * 0.3;

    Expr<double, Dynamic, Dynamic> newNNWeights = nnWeights - nnDeltaWeights;*/

    

}



int main(int argc, char * argv[], char * envp[])
{
	testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}



