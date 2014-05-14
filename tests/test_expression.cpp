#include <iostream>
#include <random>
#include <chrono>

#include "gtest/gtest.h"

#include "../src/deriv/derivative.hpp"

using namespace panoramix::deriv;
using namespace Eigen;



TEST(Expression, addConst) {

    using namespace panoramix::deriv;

    ExpressionGraph graph;

    auto m1 = graph.addConst(Matrix<double, 4, 5>::Ones());

    EXPECT_TRUE(m1.result() == (Matrix<double, 4, 5>::Ones()));

    auto m2 = graph.addConst(Matrix<float, 5, 10>::Identity() * Matrix<float, 10, 6>::Ones());
    EXPECT_TRUE(m2.result() == (Matrix<float, 5, 10>::Identity() * Matrix<float, 10, 6>::Ones()));

    auto aaa = graph.addConst(1.0);
    auto bbb = graph.addConst(2.0);
    auto ccc = generalSum(aaa, bbb);
    EXPECT_EQ(3.0, ccc.execute());
    auto ddd = generalProd(aaa, bbb);
    EXPECT_EQ(2.0, ddd.execute());
    auto n1 = -m1;
    auto dn1 = n1.sum().derivatives(m1);

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


TEST(Expression, ScalarOp) {
//void run(){

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
            EXPECT_DOUBLE_EQ(xv + yv, f.execute());
            EXPECT_DOUBLE_EQ(1.0, dfdx.execute());
            EXPECT_DOUBLE_EQ(1.0, dfdy.execute());
        }
    }

    {
        // multiplication
        auto f = cwiseProd(x, y);

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
        // pow

        auto f = pow(x, 2.0);
        auto df = f.derivatives(x);
        auto dfdx = std::get<0>(df);
        for (int i = 0; i < 10; i++) {
            xv = (std::rand() % 100000) / 100000.0;
            EXPECT_NEAR(pow(xv,2), f.execute(), 0.01);
            EXPECT_NEAR(x.execute()*2, dfdx.execute(), 0.01);
        }

    }



    {
        // exp
        auto f = exp(x);

        auto df = f.derivatives(x);
        auto dfdx = std::get<0>(df);
        for (int i = 0; i < 10; i++) {
            xv = (std::rand() % 100000) / 100000.0;
            EXPECT_NEAR(exp(xv), f.execute(), 0.01);
            EXPECT_NEAR(f.execute(), dfdx.execute(), 0.01);
        }
    }

    {
        // log
        auto f = log(x);
        auto dfdx = std::get<0>(f.derivatives(x));
        for (int i = 0; i < 10; i++) {
            xv = (std::rand() % 100000) / 100000.0;
            EXPECT_NEAR(log(xv), f.execute(), 0.01);
            EXPECT_NEAR(1.0 / xv, dfdx.execute(), 0.01);
        }
    }

    {
        // select
        auto fv = common::CWiseSelect(xv, xv, -xv);
        using TT = CWiseSelectWhenElseRetIsConstTraits<const double &, const double &, const double &>::OutputType;
        using TTT = decltype(TAG<TT>()); 

        using KK = decltype(-x);
        
        auto f = cwiseSelect(x, x, -x);
        auto dfdx = std::get<0>(f.derivatives(x));
        for (int i = 0; i < 10; i++) {
            xv = (std::rand() % 100000) / 100000.0;
            auto v = f.execute();
            EXPECT_NEAR(xv > 0 ? xv : -xv, v, 0.01);
            EXPECT_NEAR(xv > 0 ? 1 : -1, dfdx.execute(), 0.01);
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
            EXPECT_NEAR(exp(xv + yv), f.execute(), 0.01);
            EXPECT_NEAR(f.execute(x, y), dfdx.execute(x, y), 0.01);
            EXPECT_NEAR(f.execute(x, y), dfdy.execute(x, y), 0.01);
        }
    }

    {
        auto f = x + x;
        auto df = f.derivatives(x);
        auto dfdx = std::get<0>(df);

        for (int i = 0; i < 10; i++) {
            xv = (std::rand() % 100000) / 100000.0;
            EXPECT_NEAR(xv + xv, f.execute(), 0.01);
            EXPECT_NEAR(2.0, dfdx.execute(), 0.01);
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
            EXPECT_NEAR(exp(xv * yv), f.execute(), 0.01);
            EXPECT_NEAR(exp(xv * yv) * yv, dfdx.execute(), 0.01);
            EXPECT_NEAR(exp(xv * yv) * xv, dfdy.execute(), 0.01);
        }
    }

    {
        // exp(x)+x
        auto f = exp(x) + x;
        auto df = f.derivatives(x);
        auto dfdx = std::get<0>(df);

        for (int i = 0; i < 10; i++) {
            xv = (std::rand() % 100000) / 100000.0;
            EXPECT_NEAR(exp(xv) + xv, f.execute(), 0.01);
            EXPECT_NEAR(exp(xv) + 1, dfdx.execute(), 0.01);
        }
    }

    {
        // exp(x)*x
        auto f = exp(x) * x;
        auto df = f.derivatives(x);
        auto dfdx = std::get<0>(df);

        for (int i = 0; i < 10; i++) {
            xv = (std::rand() % 100000) / 100000.0;
            ASSERT_NEAR(exp(xv) * xv, f.execute(), 0.01);
            ASSERT_NEAR(exp(xv)*xv + exp(xv), dfdx.execute(), 0.01);
        }
    }




    {
        // exp(x)+x+y
        auto f = exp(x) + x + y;
        auto df = f.derivatives(x, y);

        auto dfdx = std::get<0>(df);
        auto dfdy = std::get<1>(df);

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

        auto df = f.derivatives(x, y);

        auto dfdx = std::get<0>(df);
        auto dfdy = std::get<1>(df);

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


TEST(Expression, ArrayOp){
//void run(){
    ExpressionGraph graph;
    ArrayXXd xv;
    ArrayXXd yv;
    auto x = graph.addRef(xv, "x");
    auto y = graph.addRef(yv, "y");

    {       
        auto xy = x * y;
        auto f = sumElements(xy);

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
        auto xy = cwiseProd(x, y);
        auto f = sumElements(xy);

        // mult deriv
        auto df = f.derivatives(x, y);
        auto dfdx = std::get<0>(df);
        auto dfdy = std::get<1>(df);

        for (int i = 0; i < 10; i++) {
            int a = std::rand() % 100 + 1;
            int b = std::rand() % 100 + 1;
            xv.setRandom(a, b);
            yv.setRandom(a, b);

            EXPECT_MATRIX_EQ(xv.cwiseProduct(yv).matrix(), xy.execute().matrix());
            EXPECT_DOUBLE_EQ((xv.cwiseProduct(yv)).sum(), f.execute());
            EXPECT_MATRIX_EQ(yv, dfdx.execute());
            EXPECT_MATRIX_EQ(xv, dfdy.execute());
        }
    }

    {
        // select
        auto f = cwiseSelect(x, x, -x).sum();
        auto dfdx = std::get<0>(f.derivatives(x));
        for (int i = 0; i < 10; i++) {
            int a = std::rand() % 100 + 1;
            int b = std::rand() % 100 + 1;
            xv.setRandom(a, b);
            auto v = f.execute();
            EXPECT_NEAR(xv.abs().sum(), v, 0.01);
            EXPECT_MATRIX_EQ((xv > 0).select(ArrayXXd::Ones(xv.rows(), xv.cols()), -1), dfdx.execute());
        }
    }

    {
        auto x2 = pow(x, 2.0);
        auto f = sumElements(x2);
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
        auto f = sumElements(exp(x));
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
        // exp(x*y)+x*2+y+x*x
        auto f = sumElements(exp(x * y) + x * 2 + y + x * x);
        auto df = f.derivatives(x, y);

        auto dfdx = std::get<0>(df);
        auto dfdy = std::get<1>(df);

        for (int i = 0; i < 10; i++) {
            int a = std::rand() % 100 + 1;
            int b = std::rand() % 100 + 1;
            xv.setRandom(a, b);
            yv.setRandom(a, b);
            EXPECT_DOUBLE_EQ(((xv * yv).exp() + xv * 2 + yv + xv * xv).sum(), f.execute());
            EXPECT_MATRIX_EQ((xv * yv).exp() * yv + 2 + 2 * xv, dfdx.execute());
            EXPECT_MATRIX_EQ((xv * yv).exp() * xv + 1, dfdy.execute());
        }
    }

    {
        // log(exp(x*y)+x*2+y+x*x)
        auto f = sumElements(log(exp(x * y) + x * 2 + y + x * x));
        auto df = f.derivatives(x, y);

        auto dfdx = std::get<0>(df);
        auto dfdy = std::get<1>(df);

        for (int i = 0; i < 10; i++) {
            int a = std::rand() % 100 + 1;
            int b = std::rand() % 100 + 1;
            xv.setRandom(a, b);
            xv = xv.abs() + 1;
            yv.setRandom(a, b);
            yv = yv.abs() + 1;
            EXPECT_DOUBLE_EQ(log((xv * yv).exp() + xv * 2 + yv + xv * xv).sum(), f.execute());
            EXPECT_MATRIX_EQ(((xv * yv).exp() * yv + 2 + 2 * xv).cwiseQuotient((xv * yv).exp() + xv * 2 + yv + xv * xv), dfdx.execute());
            EXPECT_MATRIX_EQ(((xv * yv).exp() * xv + 1).cwiseQuotient((xv * yv).exp() + xv * 2 + yv + xv * xv), dfdy.execute());
        }
    }

    {
        // tanh
        auto f = tanh(x).sum();
        auto df = std::get<0>(f.derivatives(x));
        for (int i = 0; i < 10; i++) {
            int a = std::rand() % 100 + 1;
            int b = std::rand() % 100 + 1;
            xv.setRandom(a, b);
            EXPECT_DOUBLE_EQ(tanh(xv).sum(), f.execute());
            EXPECT_MATRIX_EQ(((1 - tanh(xv))*(1 - tanh(xv))).matrix(), df.execute().matrix());
        }
    }

}


TEST(Expression, MatrixOp1) {
//void run(){

    ExpressionGraph graph;
    MatrixXd xv;
    auto x = graph.addRef(xv);
    auto tx = transpose(x);
    auto ttx = transpose(tx);
    auto tttx = transpose(ttx);

    auto f = sumElements(tttx);
    auto df = std::get<0>(f.derivatives(x));

    for (int i = 0; i < 10; i++){
        int a = std::rand() % 100 + 1;
        int b = std::rand() % 100 + 1;
        xv.setRandom(a, b);
        
        EXPECT_DOUBLE_EQ(xv.sum(), f.execute());
        EXPECT_MATRIX_EQ(MatrixXd::Ones(a, b), df.execute());
    }

}


TEST(Expression, MatrixOp2) {

    // plus

    //void run(){
    ExpressionGraph graph;

    Matrix<float, Dynamic, Dynamic> xv;
    Matrix<float, Dynamic, Dynamic> yv;

    auto x = graph.addRef(xv);
    auto y = graph.addRef(yv);

    {
        // +
        // plus
        auto f = sumElements(x + y * 2);

        // plus deriv
        auto df = f.derivatives(x, y);
        auto dfdx = std::get<0>(df);
        auto dfdy = std::get<1>(df);

        for (int i = 0; i < 10; i++) {
            int r = std::rand() % 100 + 1;
            int c = std::rand() % 100 + 1;
            xv.setRandom(r, c);
            yv.setRandom(r, c);
            EXPECT_FLOAT_EQ((xv + yv * 2).sum(), f.execute());
            EXPECT_MATRIX_EQ(MatrixXf::Ones(r, c), dfdx.execute());
            EXPECT_MATRIX_EQ(MatrixXf::Ones(r, c) * 2, dfdy.execute());
        }
    }

}

TEST(Expression, MatrixOp3){
//void run(){

    ExpressionGraph graph;
    MatrixXd xv;
    MatrixXd yv;
    auto x = graph.addRef(xv, "x");
    auto y = graph.addRef(yv, "y");

    {
        // mult
        auto xy = x * y;
        auto f = sumElements(xy);
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
       /* xv.setRandom(3, 4);
        yv.setRandom(4, 4);
        auto z = xv * (yv * yv);
        auto zz = yv * yv;
        auto zzz = xv * zz;
        std::cout << z << std::endl;
        std::cout << zz << std::endl;
        std::cout << zzz << std::endl;*/

    }

    {
        //// general product -> eval
        //xv.setRandom(3, 4);
        //yv.setRandom(4, 4);
        //auto xy = generalProd(x, y).eval();
        //std::cout << xy << std::endl;
        //std::cout << xy.execute() << std::endl;

        ////auto a = xv * yv * yv * yv;       
        /*xv.setRandom(3, 4);
        yv.setRandom(4, 4);
        auto z = ProdAll(xv, yv, yv, yv);
        std::cout << z << std::endl;*/
    }

    { // matrix array conversion
        auto xx = arrayToMatrix(matrixToArray(x));
        for (int i = 0; i < 10; i++){
            int a = std::rand() % 100 + 1;
            int b = std::rand() % 100 + 1;
            int c = std::rand() % 100 + 1;
            xv.setRandom(a, b);
            yv.setRandom(b, c);

            EXPECT_MATRIX_EQ(xv, xv.array().matrix());
            EXPECT_MATRIX_EQ(xv, xx.execute());
        }

        auto xxp = pow(matrixToArray(x).eval(), 2);
        for (int i = 0; i < 10; i++) {
            int a = std::rand() % 100 + 1;
            int b = std::rand() % 100 + 1;
            xv.setRandom(a, b);

            auto dxxp = std::get<0>(sumElements(xxp).derivatives(x));

            EXPECT_MATRIX_EQ(xv.array().pow(2).matrix(), xxp.execute().matrix());
            EXPECT_MATRIX_EQ(xv * 2, dxxp.execute().matrix());
        }
    }

    {
        // matrix inverse
        auto xinv = inverseMatrix(x);
        for (int i = 0; i < 10; i++) {
            int a = std::rand() % 100 + 1;
            xv.setRandom(a, a);
            auto dxinv = std::get<0>(sumElements(xinv).derivatives(x));
            EXPECT_MATRIX_EQ(xv.inverse(), xinv.execute());
            auto grad = dxinv.execute();
            EXPECT_MATRIX_EQ((xv.inverse() * -MatrixXd::Ones(a, a) * xv.inverse()).transpose(), dxinv.execute());
        }
    }

}


TEST(Expression, VectorOp) {

    ExpressionGraph graph;
    VectorXd xv, yv;
    auto x = graph.addRef(xv, "x");
    auto y = graph.addRef(yv, "y");

    {
        // dot product
        auto d = dotProd(x, y);
        auto dd = d.derivatives(x, y);
        auto dddx = std::get<0>(dd);
        auto dddy = std::get<1>(dd);
        for (int i = 0; i < 10; i++){
            int a = std::rand() % 100 + 1;
            xv.setRandom(a);
            yv.setRandom(a);
            EXPECT_NEAR(xv.dot(yv), d.execute(), 0.01);
            EXPECT_MATRIX_EQ(xv, dddy.execute());
            EXPECT_MATRIX_EQ(yv, dddx.execute());
        }
    }

    {
        // cross product
        auto c = makeCross3ProductMatrix<Matrix3d>(x.assign<Vector3d>()) * y.assign<Vector3d>();
        std::cout << c << std::endl;
        auto dc = c.sum().derivatives(x, y);
        auto dcdx = std::get<0>(dc);
        auto dcdy = std::get<1>(dc);
        for (int i = 0; i < 10; i++){
            xv.setRandom(3);
            yv.setRandom(3);
            EXPECT_MATRIX_EQ(x.assign<Vector3d>().execute().cross(y.assign<Vector3d>().execute()), 
                c.execute());
            dcdx.execute();
            dcdy.execute();
        }
    }

    {
        // norm
        auto n = norm(x);
        auto nx = x / n;
        auto dn = std::get<0>(n.derivatives(x));
        auto dnx = std::get<0>(nx.sum().derivatives(x));
        for (int i = 0; i < 10; i++){
            int a = std::rand() % 100 + 1;
            xv.setRandom(a);
            EXPECT_NEAR(xv.norm(), n.execute(), 0.001);
            EXPECT_MATRIX_EQ(xv.normalized(), nx.execute());
            EXPECT_MATRIX_EQ(xv.normalized(), dn.execute());
            dnx.execute();
        }
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



