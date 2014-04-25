#include "../src/core/expression.hpp"
#include "../src/core/expression_ops.hpp"

#include "gtest/gtest.h"

#include <iostream>
#include <random>
#include <chrono>

using namespace panoramix::core;


TEST(Expression, addConst) {

    ExpressionGraph graph;
    
    auto m1 = graph.addConst(Matrix<double, 4, 5>::Ones());
    EXPECT_TRUE(m1.eval() == (Matrix<double, 4, 5>::Ones()));

    auto m2 = graph.addConst(Matrix<float, 5, 10>::Identity() * Matrix<float, 10, 6>::Ones());
    EXPECT_TRUE(m2.eval() == (Matrix<float, 5, 10>::Identity() * Matrix<float, 10, 6>::Ones()));

    MatrixXd m3d(5, 2);
    m3d << 1, 2, 3, 4, 5, 6, 7, 8, 9, 0;
    auto m3 = graph.addConst(m3d * Matrix<double, 2, 6>::Ones());
    EXPECT_TRUE(m3.eval() == (m3d * Matrix<double, 2, 6>::Ones()));

    SparseMatrix<double> m4d;
    m4d.resize(5, 5);
    m4d.setIdentity();
    Fill(m4d, 5);
    auto m4 = graph.addConst(m4d);
    EXPECT_TRUE(MatrixXd(m4d) == MatrixXd(m4.eval()));

}

TEST(Expression, addRef) {

    ExpressionGraph graph;
    MatrixXd a(4, 5);
    a.fill(5.5);   

    auto aa = graph.addRef(a);
    EXPECT_TRUE(aa.eval() == a);

    a.fill(2.2);
    EXPECT_TRUE(aa.eval() == a);

    a.resize(7, 4);
    EXPECT_TRUE(aa.eval() == a);

    a.transposeInPlace();
    EXPECT_TRUE(aa.eval() == a);

    a *= Matrix<double, 7, 8>::Identity();
    EXPECT_TRUE(aa.eval() == a);
    
}


TEST(Expression, Scalars) {

    ExpressionGraph graph;

    double x = 1.0;
    auto xx = graph.addRef(x);

    EXPECT_EQ(x, xx.eval());
    x = 3.0;
    EXPECT_EQ(x, xx.eval());

    float y = 2.0;
    auto yy = graph.addRef(y);

    EXPECT_EQ(y, yy.eval());
    y = 10;
    EXPECT_EQ(y, yy.eval());

    std::complex<double> c(4, 3);
    auto cc = graph.addRef(c);

    EXPECT_TRUE(cc.eval() == c);
    c.imag(5);
    EXPECT_TRUE(cc.eval() == c);

}


TEST(Expression, ScalarOp) {
    ExpressionGraph graph;

    double xv = 9.0;
    double yv = 5.0;
    
    auto x = graph.addRef(xv, "x");
    auto y = graph.addRef(yv, "y");

    {
        // plus
        auto f = x + y;
        // plus deriv
        auto df = f.derivatives(x, y);
        auto dfdx = std::get<0>(df);
        auto dfdy = std::get<1>(df);

        for (int i = 0; i < 10; i++) {
            xv = (std::rand() % 100000) / 100000.0;
            yv = (std::rand() % 100000) / 100000.0;
            EXPECT_DOUBLE_EQ(xv + yv, f.eval());
            EXPECT_DOUBLE_EQ(1.0, dfdx.eval());
            EXPECT_DOUBLE_EQ(1.0, dfdy.eval());
        }
    }

    {
        // multiplication
        auto f = x * y;
        auto df = f.derivatives(x, y);
        auto dfdx = std::get<0>(df);
        auto dfdy = std::get<1>(df);
        for (int i = 0; i < 10; i++) {
            xv = (std::rand() % 100000) / 100000.0;
            yv = (std::rand() % 100000) / 100000.0;
            EXPECT_NEAR(xv * yv, f.eval(), 0.01);
            EXPECT_NEAR(y.eval(), dfdx.eval(), 0.01);
            EXPECT_NEAR(x.eval(), dfdy.eval(), 0.01);
        }
    }

    {
        // exp
        auto f = exp(x);
        auto df = f.derivatives(x, y);
        auto dfdx = std::get<0>(df);
        for (int i = 0; i < 10; i++) {
            xv = (std::rand() % 100000) / 100000.0;
            EXPECT_NEAR(exp(xv), f.eval(), 0.01);
            EXPECT_NEAR(f.eval(), dfdx.eval(), 0.01);
        }
    }

    {
        // exp(x+y)
        auto f = exp(x + y);
        auto df = f.derivatives(x, y);
        auto dfdx = std::get<0>(df);
        auto dfdy = std::get<1>(df);
        for (int i = 0; i < 10; i++) {
            xv = (std::rand() % 100000) / 100000.0;
            yv = (std::rand() % 100000) / 100000.0;
            EXPECT_NEAR(exp(xv + yv), f.eval(), 0.01);
            EXPECT_NEAR(f.eval(), dfdx.eval(), 0.01);
            EXPECT_NEAR(f.eval(), dfdy.eval(), 0.01);
        }
    }

    {
        // exp(x*y)
        auto f = exp(x * y);
        auto df = f.derivatives(x, y);
        auto dfdx = std::get<0>(df);
        auto dfdy = std::get<1>(df);
        for (int i = 0; i < 10; i++) {
            xv = (std::rand() % 100000) / 100000.0;
            yv = (std::rand() % 100000) / 100000.0;
            EXPECT_NEAR(exp(xv * yv), f.eval(), 0.01);
            EXPECT_NEAR(exp(xv * yv) * yv, dfdx.eval(), 0.01);
            EXPECT_NEAR(exp(xv * yv) * xv, dfdy.eval(), 0.01);
        }
    }

    {
        // exp(x*y)+x+y
        auto f = exp(x * y) + x + y;
        auto df = f.derivatives(x, y);
        auto dfdx = std::get<0>(df);
        auto dfdy = std::get<1>(df);
        for (int i = 0; i < 10; i++) {
            xv = (std::rand() % 100000) / 100000.0;
            yv = (std::rand() % 100000) / 100000.0;
            EXPECT_NEAR(exp(xv * yv) + xv + yv, f.eval(), 0.01);
            EXPECT_NEAR(exp(xv * yv) * yv + 1, dfdx.eval(), 0.01);
            EXPECT_NEAR(exp(xv * yv) * xv + 1, dfdy.eval(), 0.01);
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
            EXPECT_NEAR(exp(xv * exp(yv)), f.eval(), 0.01);
            EXPECT_NEAR(exp(xv * exp(yv)) * exp(yv), dfdx.eval(), 0.01);
            EXPECT_NEAR(exp(xv * exp(yv)) * xv * exp(yv), dfdy.eval(), 0.01);
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

template <class T>
void foo(const T & t){
    std::cout << "default" << std::endl;
}

template <class Derived>
void foo(const Eigen::MatrixBase<Derived> & t) {
    std::cout << "dense version" << std::endl;
}

TEST(Expression, TemplateMatching) {
    MatrixXd xv;
    foo(xv);
}

template <class T, class T1, class T2>
T passValue(T1 && t1, T2 && t2) {
    return t1 + t2 * 2;
}

TEST(Expression, EigenExprTypePassing) {
    MatrixXd x;
    MatrixXd y;
    x.setRandom(2, 2);
    y.setRandom(2, 2);
    using T = decltype(x + y * 2);
    auto m = passValue<T>(x, y);

    MatrixXd sum = m;
    EXPECT_MATRIX_EQ(sum, x + y * 2);
}


TEST(Expression, MatrixOp) {

    ExpressionGraph graph;

    Matrix<double, Dynamic, Dynamic> xv;
    Matrix<double, Dynamic, Dynamic> yv;

    auto x = graph.addRef(xv);
    auto y = graph.addRef(yv);

    {
        // +
        // plus
        auto f = SumElements(x + y);
        // plus deriv
        auto df = f.derivatives(x, y);
        auto dfdx = std::get<0>(df);
        auto dfdy = std::get<1>(df);

        for (int i = 0; i < 10; i++) {
            int r = std::rand() % 100 + 1;
            int c = std::rand() % 100 + 1;
            xv.setRandom(r, c);
            yv.setRandom(r, c);
            EXPECT_DOUBLE_EQ((xv + yv).sum(), f.eval());
            EXPECT_MATRIX_EQ(MatrixXd::Ones(r, c), dfdx.eval());
            EXPECT_MATRIX_EQ(MatrixXd::Ones(r, c), dfdy.eval());
        }
    }



    //{
    //    // *
    //    // mult
    //    auto f = SumElements(x * y);
    //    // plus deriv
    //    auto df = f.derivatives(x, y);
    //    auto dfdx = std::get<0>(df);
    //    auto dfdy = std::get<1>(df);

    //    for (int i = 0; i < 10; i++) {
    //        int a = std::rand() % 100 + 1;
    //        int b = std::rand() % 100 + 1;
    //        int c = std::rand() % 100 + 1;
    //        xv.setRandom(a, b);
    //        yv.setRandom(b, c);
    //        EXPECT_DOUBLE_EQ((xv * yv).sum(), f.eval());
    //        EXPECT_MATRIX_EQ(yv.transpose(), dfdx.eval());
    //        EXPECT_MATRIX_EQ(xv.transpose(), dfdy.eval());
    //    }
    //}

}








int main(int argc, char * argv[], char * envp[])
{
	testing::InitGoogleTest(&argc, argv);
    testing::GTEST_FLAG(catch_exceptions) = true;
    testing::GTEST_FLAG(show_internal_stack_frames) = true;
    return RUN_ALL_TESTS();
}



