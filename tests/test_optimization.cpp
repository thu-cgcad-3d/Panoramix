#include <iostream>
#include <random>

#include "../src/deriv/derivative.hpp"
#include "../src/core/feature.hpp"

#include "gtest/gtest.h"

using namespace Eigen;

using namespace panoramix::core;
using namespace panoramix::deriv;

TEST(Optimization, Experiment1){

    ExpressionGraph graph;
    MatrixXd xv;
    auto x = graph.addRef(xv, "x");

    // simple optimization
    xv.setRandom(100, 100);
    auto energy = inverseMatrix(x * x * x + arrayToMatrix(exp(matrixToArray(x))).eval()).sum();
    auto grad = std::get<0>(energy.derivatives(x));

    std::cout << "energy: " << energy << std::endl;
    //std::cout << "grad:   " << grad << std::endl;

    for (int i = 0; i < 50; i++) {
        std::cout << "energy: " << energy.execute(x) << std::endl;
        auto g = grad.execute(x, energy);
        xv -= g;
    }

    ArrayXd a, b;
    a.select(a, b);
    auto c = (a > b).eval();
    std::cout << typeid(c).name() << std::endl;
}


TEST(Optimization, ANN) {
//void run(){

    int inputn = 100, hidn = 80, outputn = 10;

    ExpressionGraph graph;
    MatrixXd weightVals1, weightVals2;
    MatrixXd biasVals1, biasVals2;
    
    auto inputs = graph.addConst(MatrixXd::Random(1, inputn));
    
    ArrayXXd arr = ArrayXXd::Random(1, outputn).cwiseAbs();
    arr = arr / arr.maxCoeff() / 2.0;
    auto labels = graph.addConst(arr);
    
    weightVals1.setRandom(inputn, hidn);
    weightVals2.setRandom(hidn, outputn);
    biasVals1.setRandom(1, hidn);
    biasVals2.setRandom(1, outputn);
    auto weights1 = graph.addRef(weightVals1);
    auto weights2 = graph.addRef(weightVals2);
    auto bias1 = graph.addRef(biasVals1);
    auto bias2 = graph.addRef(biasVals2);

    auto hidunits = sigmoid(matrixToArray(inputs * weights1 + bias1));
    auto diff = sigmoid(matrixToArray(arrayToMatrix(hidunits) * weights2 + bias2)) - labels;
    auto energy = (diff * diff).sum();

    Expression<MatrixXd> grad1, grad2, gradb1, gradb2;
    std::tie(grad1, grad2, gradb1, gradb2) = energy.derivatives(weights1, weights2, bias1, bias2);

    MatrixXd lastchange1 = MatrixXd::Zero(weightVals1.rows(), weightVals1.cols());
    MatrixXd lastchange2 = MatrixXd::Zero(weightVals2.rows(), weightVals2.cols());
    MatrixXd lastchangeb1 = MatrixXd::Zero(biasVals1.rows(), biasVals1.cols());
    MatrixXd lastchangeb2 = MatrixXd::Zero(biasVals2.rows(), biasVals2.cols());

    double delta = 5;

    for (int i = 0; i < 200; i++) {
        std::cout << "energy: " << energy.execute(weights1, weights2) << std::endl;
        MatrixXd g1 = grad1.execute(weights1, weights2, bias1, bias2, energy);
        MatrixXd g2 = grad2.execute(weights1, weights2, bias1, bias2, energy);
        MatrixXd gb1 = gradb1.execute(weights1, weights2, bias1, bias2, energy);
        MatrixXd gb2 = gradb2.execute(weights1, weights2, bias1, bias2, energy);
        MatrixXd change1 = (-g1*0.7 + lastchange1*0.3);
        MatrixXd change2 = (-g2*0.7 + lastchange2*0.3);
        MatrixXd changeb1 = (-gb1*0.7 + lastchangeb1*0.3);
        MatrixXd changeb2 = (-gb2*0.7 + lastchangeb2*0.3);
        weightVals1 += change1 * delta;
        weightVals2 += change2 * delta;
        biasVals1 += changeb1 * delta;
        biasVals2 += changeb2 * delta;
        lastchange1 = change1;
        lastchange2 = change2;
        lastchangeb1 = changeb1;
        lastchangeb2 = changeb2;
    }

}


TEST(Optimization, Camera) {

    PerspectiveCamera cam1, cam2;
    ExpressionGraph graph;



}





int main(int argc, char * argv[], char * envp[])
{
	testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
    //run();
    return 1;
}



