#include "../src/core/expression.hpp"
#include "../src/core/expression_ops.hpp"

#include "gtest/gtest.h"


#include <iostream>
#include <random>
#include <chrono>

using namespace panoramix::core;
using namespace Eigen;



TEST(Expression, addConst) {

    ExpressionGraph graph;
    
    bool b = traits<decltype(Matrix<double, 4, 5>::Ones())>::should_be_cached;

    auto m1 = graph.addConst(Matrix<double, 4, 5>::Ones());
    EXPECT_TRUE(m1.result() == (Matrix<double, 4, 5>::Ones()));

    auto m2 = graph.addConst(Matrix<float, 5, 10>::Identity() * Matrix<float, 10, 6>::Ones());
    EXPECT_TRUE(m2.result() == (Matrix<float, 5, 10>::Identity() * Matrix<float, 10, 6>::Ones()));
    
    
    MatrixXd m3d(5, 2);
    auto m3f = m3d.cast<float>().eval();

    m3d << 1, 2, 3, 4, 5, 6, 7, 8, 9, 0;
    auto m3 = graph.addConst(m3d * Matrix<double, 2, 6>::Ones());
    EXPECT_TRUE(m3.result() == (m3d * Matrix<double, 2, 6>::Ones()));
}


TEST(Expression, addRef) {

    ExpressionGraph graph;
    MatrixXd a(4, 5);
    a.fill(5.5);   

    auto aa = graph.addRef(a);
    EXPECT_TRUE(aa.result() == a);

    a.fill(2.2);
    std::cout << aa.result() << std::endl;
    EXPECT_TRUE(aa.result() == a);

    a.resize(7, 4);
    EXPECT_TRUE(aa.result() == a);

    a.transposeInPlace();
    EXPECT_TRUE(aa.result() == a);

    a *= Matrix<double, 7, 8>::Identity();
    EXPECT_TRUE(aa.result() == a);
    
}


TEST(Expression, Scalars) {

    ExpressionGraph graph;

    double x = 1.0;
    auto xx = graph.addRef(x);

    EXPECT_EQ(x, xx.result());
    x = 3.0;
    EXPECT_EQ(x, xx.result());

    float y = 2.0;
    auto yy = graph.addRef(y);

    EXPECT_EQ(y, yy.result());
    y = 10;
    EXPECT_EQ(y, yy.result());

    auto zz = xx + yy;
    zz.execute();
    EXPECT_EQ(zz.result(), x + y);

    auto dzz = zz.derivatives(xx);
    auto dzzxx = std::get<0>(dzz);
    dzzxx.execute();
    EXPECT_EQ(dzzxx.result(), 1);

}


TEST(Expression, Op) {
//void run(){

    ExpressionGraph graph;

    double xv = 9.0;
    double yv = 5.0;
    
    auto x = graph.addRef(xv, "x");
    auto y = graph.addRef(yv, "y");

    {
        // plus
        auto f = x + y;

        std::cout << f << std::endl;

        // plus deriv
        auto df = f.derivatives(x, y);
        auto dfdx = std::get<0>(df);
        auto dfdy = std::get<1>(df);

        for (int i = 0; i < 10; i++) {
            xv = (std::rand() % 100000) / 100000.0;
            yv = (std::rand() % 100000) / 100000.0;
            EXPECT_DOUBLE_EQ(xv + yv, f.execute());
            EXPECT_DOUBLE_EQ(1.0, dfdx.execute());
            EXPECT_DOUBLE_EQ(1.0, dfdy.execute());
        }
    }

    {
        // multiplication
        auto f = x * y;

        std::cout << f << std::endl;

        auto df = f.derivatives(x, y);
        auto dfdx = std::get<0>(df);
        auto dfdy = std::get<1>(df);
        for (int i = 0; i < 10; i++) {
            xv = (std::rand() % 100000) / 100000.0;
            yv = (std::rand() % 100000) / 100000.0;
            EXPECT_NEAR(xv * yv, f.execute(), 0.01);
            EXPECT_NEAR(y.execute(), dfdx.execute(), 0.01);
            EXPECT_NEAR(x.execute(), dfdy.execute(), 0.01);
        }
    }

    {
        // exp
        auto f = exp(x);

        std::cout << f << std::endl;

        auto df = f.derivatives(x);
        auto dfdx = std::get<0>(df);
        for (int i = 0; i < 10; i++) {
            xv = (std::rand() % 100000) / 100000.0;
            EXPECT_NEAR(exp(xv), f.execute(), 0.01);
            EXPECT_NEAR(f.execute(), dfdx.execute(), 0.01);
        }
    }

    {
        // exp(x+y)
        auto f = exp(x + y);

        std::cout << f << std::endl;

        auto df = f.derivatives(x, y);
        auto dfdx = std::get<0>(df);
        auto dfdy = std::get<1>(df);

        for (int i = 0; i < 10; i++) {
            xv = (std::rand() % 100000) / 100000.0;
            yv = (std::rand() % 100000) / 100000.0;
            EXPECT_NEAR(exp(xv + yv), f.execute(), 0.01);
            EXPECT_NEAR(f.execute(x, y), dfdx.execute(x, y), 0.01);
            EXPECT_NEAR(f.execute(x, y), dfdy.execute(x, y), 0.01);
        }
    }

    {
        auto f = x + x;

        std::cout << f << std::endl;
        auto df = f.derivatives(x);
        auto dfdx = std::get<0>(df);

        std::cout << dfdx << std::endl;

        for (int i = 0; i < 10; i++) {
            xv = (std::rand() % 100000) / 100000.0;
            EXPECT_NEAR(xv + xv, f.execute(), 0.01);
            EXPECT_NEAR(2.0, dfdx.execute(), 0.01);
        }
    }

    {
        // exp(x*y)
        auto f = exp(x * y);

        std::cout << f << std::endl;

        auto df = f.derivatives(x, y);

        auto dfdx = std::get<0>(df);
        auto dfdy = std::get<1>(df);

        std::cout << dfdx << std::endl;
        std::cout << dfdy << std::endl;

        for (int i = 0; i < 10; i++) {
            xv = (std::rand() % 100000) / 100000.0;
            yv = (std::rand() % 100000) / 100000.0;
            EXPECT_NEAR(exp(xv * yv), f.execute(), 0.01);
            EXPECT_NEAR(exp(xv * yv) * yv, dfdx.execute(), 0.01);
            EXPECT_NEAR(exp(xv * yv) * xv, dfdy.execute(), 0.01);
        }
    }

    {
        // exp(x)+x
        auto f = exp(x) + x;

        std::cout << "=============================" << std::endl;
        std::cout << f << std::endl;

        auto df = f.derivatives(x);
        auto dfdx = std::get<0>(df);
        std::cout << dfdx << std::endl;

        std::cout << "=============================" << std::endl;

        for (int i = 0; i < 10; i++) {
            xv = (std::rand() % 100000) / 100000.0;
            EXPECT_NEAR(exp(xv) + xv, f.execute(), 0.01);
            EXPECT_NEAR(exp(xv) + 1, dfdx.execute(), 0.01);
        }
    }

    {
        // exp(x)*x
        auto f = cwise_product(exp(x), x);

        std::cout << "=============================" << std::endl;
        std::cout << f << std::endl;

        auto df = f.derivatives(x);

        auto dfdx = std::get<0>(df);

        std::cout << dfdx << std::endl;

        std::cout << "=============================" << std::endl;

        for (int i = 0; i < 10; i++) {
            xv = (std::rand() % 100000) / 100000.0;
            ASSERT_NEAR(exp(xv) * xv, f.execute(), 0.01);
            ASSERT_NEAR(exp(xv)*xv + exp(xv), dfdx.execute(), 0.01);
        }
    }




    {
        // exp(x)+x+y
        auto f = exp(x) + x + y;

        std::cout << "=============================" << std::endl;
        std::cout << f << std::endl;

        auto df = f.derivatives(x, y);

        auto dfdx = std::get<0>(df);
        auto dfdy = std::get<1>(df);

        std::cout << dfdx << std::endl;
        std::cout << dfdy << std::endl;

        std::cout << "=============================" << std::endl;

        for (int i = 0; i < 10; i++) {
            xv = (std::rand() % 100000) / 100000.0;
            yv = (std::rand() % 100000) / 100000.0;
            EXPECT_NEAR(exp(xv) + xv + yv, f.execute(), 0.01);
            EXPECT_NEAR(exp(xv) + 1, dfdx.execute(), 0.01);
            EXPECT_NEAR(1, dfdy.execute(), 0.01);
        }
    }

    {
        // exp(x*y)+x+y
        auto f = exp(x * y) + x + y;

        std::cout << "=============================" << std::endl;
        std::cout << f << std::endl;

        auto df = f.derivatives(x, y);

        auto dfdx = std::get<0>(df);
        auto dfdy = std::get<1>(df);

        std::cout << dfdx << std::endl;
        std::cout << dfdy << std::endl;

        std::cout << "=============================" << std::endl;

        for (int i = 0; i < 10; i++) {
            xv = (std::rand() % 100000) / 100000.0;
            yv = (std::rand() % 100000) / 100000.0;
            EXPECT_NEAR(exp(xv * yv) + xv + yv, f.execute(), 0.01);
            EXPECT_NEAR(exp(xv * yv) * yv + 1, dfdx.execute(), 0.01);
            EXPECT_NEAR(exp(xv * yv) * xv + 1, dfdy.execute(), 0.01);
        }
    }

    {
        // exp(x*exp(y))
        auto f = exp(x * exp(y));
        auto df = f.derivatives(x, y);
        auto dfdx = std::get<0>(df);
        auto dfdy = std::get<1>(df);
        for (int i = 0; i < 10; i++) {
            xv = (std::rand() % 100000) / 100000.0;
            yv = (std::rand() % 100000) / 100000.0;
            EXPECT_NEAR(exp(xv * exp(yv)), f.execute(), 0.01);
            EXPECT_NEAR(exp(xv * exp(yv)) * exp(yv), dfdx.execute(), 0.01);
            EXPECT_NEAR(exp(xv * exp(yv)) * xv * exp(yv), dfdy.execute(), 0.01);
        }
    }

}

namespace Helper {

    template<typename DerivedA, typename DerivedB>
    bool allclose(const Eigen::DenseBase<DerivedA>& a,
        const Eigen::DenseBase<DerivedB>& b,
        const typename DerivedA::RealScalar& rtol
        = Eigen::NumTraits<typename DerivedA::RealScalar>::dummy_precision(),
        const typename DerivedA::RealScalar& atol
        = Eigen::NumTraits<typename DerivedA::RealScalar>::epsilon())
    {
        return ((a.derived() - b.derived()).array().abs()
            <= (atol + rtol * b.derived().array().abs())).all();
    }

#define EXPECT_MATRIX_EQ(a, b) EXPECT_TRUE(Helper::allclose(a, b))

}


TEST(Expression, MatrixPlusSum) {
    //void run(){
    ExpressionGraph graph;

    Matrix<double, Dynamic, Dynamic> xv;
    Matrix<double, Dynamic, Dynamic> yv;
    decltype(std::declval<TransposeResultType<MatrixXd>>().eval()) a;

    auto x = graph.addRef(xv);
    auto y = graph.addRef(yv);

    {
        // +
        // plus
        auto f = (x + y).sum();

        // plus deriv
        auto df = f.derivatives(x, y);
        auto dfdx = std::get<0>(df);
        auto dfdy = std::get<1>(df);

        for (int i = 0; i < 10; i++) {
            int r = std::rand() % 100 + 1;
            int c = std::rand() % 100 + 1;
            xv.setRandom(r, c);
            yv.setRandom(r, c);
            EXPECT_DOUBLE_EQ((xv + yv).sum(), f.execute());
            EXPECT_MATRIX_EQ(MatrixXd::Ones(r, c), dfdx.execute());
            EXPECT_MATRIX_EQ(MatrixXd::Ones(r, c), dfdy.execute());
        }
    }

}


TEST(Expression, ArrayOp){
//void run(){
    ExpressionGraph graph;
    ArrayXXd xv;
    ArrayXXd yv;
    auto x = graph.addRef(xv);
    auto y = graph.addRef(yv);

    {       
        auto xy = cwise_product(x, y);
        auto f = xy.sum();

        // mult deriv
        auto df = f.derivatives(x, y);
        auto dfdx = std::get<0>(df);
        auto dfdy = std::get<1>(df);

        for (int i = 0; i < 10; i++) {
            int a = std::rand() % 100 + 1;
            int b = std::rand() % 100 + 1;
            xv.setRandom(a, b);
            yv.setRandom(a, b);
            EXPECT_DOUBLE_EQ((xv.cwiseProduct(yv)).sum(), f.execute());
            EXPECT_MATRIX_EQ(yv, dfdx.execute());
            EXPECT_MATRIX_EQ(xv, dfdy.execute());
        }
    }

    {
        auto x2 = pow(x, 2.0);
        auto f = x2.sum();
        auto df = std::get<0>(f.derivatives(x));

        for (int i = 0; i < 10; i++) {
            int a = std::rand() % 100 + 1;
            int b = std::rand() % 100 + 1;
            xv.setRandom(a, b);

            EXPECT_DOUBLE_EQ(xv.pow(2.0).sum(), f.execute());
            EXPECT_MATRIX_EQ(xv * 2.0, df.execute());
        }

    }

    {
        // exp
        auto f = exp(x).sum();
        auto df = f.derivatives(x);
        auto dfdx = std::get<0>(df);
        for (int i = 0; i < 10; i++) {
            int a = std::rand() % 100 + 1;
            int b = std::rand() % 100 + 1;
            xv.setRandom(a, b);

            EXPECT_DOUBLE_EQ(exp(xv).sum(), f.execute());
            EXPECT_MATRIX_EQ(exp(xv), dfdx.execute());
        }
    }

    {
        // exp(x*y)+x+y
        auto f = (exp(x * y) + x + y).sum();
        auto df = f.derivatives(x, y);

        auto dfdx = std::get<0>(df);
        auto dfdy = std::get<1>(df);

        for (int i = 0; i < 10; i++) {
            int a = std::rand() % 100 + 1;
            int b = std::rand() % 100 + 1;
            xv.setRandom(a, b);
            yv.setRandom(a, b);
            EXPECT_DOUBLE_EQ(((xv * yv).exp() + xv + yv).sum(), f.execute());
            EXPECT_MATRIX_EQ((xv * yv).exp() * yv + 1.0, dfdx.execute());
            EXPECT_MATRIX_EQ((xv * yv).exp() * xv + 1.0, dfdy.execute());
        }
    }

}


TEST(Expression, MatrixTranspose) {
//void run(){

    ExpressionGraph graph;
    MatrixXd xv;
    auto x = graph.addRef(xv);
    auto tx = transpose(x);
    auto ttx = transpose(tx);
    auto tttx = transpose(ttx);

    auto f = tttx.sum();
    auto df = std::get<0>(f.derivatives(x));

    for (int i = 0; i < 10; i++){
        int a = std::rand() % 100 + 1;
        int b = std::rand() % 100 + 1;
        xv.setRandom(a, b);
        
        EXPECT_DOUBLE_EQ(xv.sum(), f.execute());
        EXPECT_MATRIX_EQ(MatrixXd::Ones(a, b), df.execute());
    }

}


TEST(Expression, MatrixMultSum){
//void run(){
    ExpressionGraph graph;
    MatrixXd xv;
    MatrixXd yv;
    auto x = graph.addRef(xv);
    auto y = graph.addRef(yv);    

    {
        auto xy = x * y;
        auto f = xy.sum();
        auto df = f.derivatives(x, y);
        auto dfdx = std::get<0>(df);
        auto dfdy = std::get<1>(df);

        for (int i = 0; i < 10; i++) {
            int a = std::rand() % 100 + 1;
            int b = std::rand() % 100 + 1;
            int c = std::rand() % 100 + 1;
            xv.setRandom(a, b);
            yv.setRandom(b, c);

            EXPECT_DOUBLE_EQ((xv * yv).sum(), f.execute());
            EXPECT_MATRIX_EQ(MatrixXd::Ones(a, c) * yv.transpose(), dfdx.execute());
            EXPECT_MATRIX_EQ(xv.transpose() * MatrixXd::Ones(a, c), dfdy.execute());
        }
    }

    {
        
        
    }
}





int main(int argc, char * argv[], char * envp[])
{
	testing::InitGoogleTest(&argc, argv);
    testing::GTEST_FLAG(catch_exceptions) = true;
    testing::GTEST_FLAG(show_internal_stack_frames) = true;
    return RUN_ALL_TESTS();
    //run();
    return 1;
}



