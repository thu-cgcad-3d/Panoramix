#include <iostream>
#include <random>

#include "../src/deriv/derivative.hpp"
#include "../src/core/feature.hpp"

#include "gtest/gtest.h"

using namespace Eigen;

using namespace panoramix::core;
using namespace panoramix::deriv;


TEST(Optimization, Quadratic) {

    ExpressionGraph graph;
    Matrix<double, 2, 1> xv;
    xv.setRandom();
    
    auto x = graph.addRef(xv);
    Matrix2d Av;
    Av << 1, 0, 0, 1;
    auto A = graph.addConst(Av);

    auto f = (transpose(x) * A * x).sum();
    auto df = std::get<0>(f.derivatives(x));

    for (int i = 0; i < 10; i++){
        std::cout << "f: " << f.execute() << std::endl;
        auto g = df.execute();
        xv -= g * 0.48;
    }

}


TEST(Optimization, InverseMatrix){

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
        xv -= g * 0.5;
    }
}

TEST(Optimization, DISABLED_ANN) {
//void run(){

    int inputn = 100, hidn = 80, outputn = 10;
    int samplesn = 100;

    ExpressionGraph graph;
    MatrixXd weightVals1, weightVals2;
    MatrixXd biasVals1, biasVals2;
    
    MatrixXd inputVals = MatrixXd::Random(samplesn, inputn);
    auto inputs = graph.addConst(inputVals);
    ArrayXXd labelVals = ArrayXXd::Zero(samplesn, outputn);
    ArrayXXd temp = inputVals.array() * (1.0 - inputVals.array());
    for (int i = 0; i < temp.rows(); i++){
        labelVals.col(i * outputn / inputn) = temp.col(i);
    }
    auto labels = graph.addConst(labelVals);
    
    weightVals1.setRandom(inputn, hidn);
    weightVals2.setRandom(hidn, outputn);
    biasVals1.setZero(1, hidn);
    biasVals2.setZero(1, outputn);
    auto weights1 = graph.addRef(weightVals1);
    auto weights2 = graph.addRef(weightVals2);
    //auto bias1 = graph.addRef(biasVals1);
    //auto bias2 = graph.addRef(biasVals2);


    auto hidunits = tanh(matrixToArray(inputs * weights1)); // samplesn x hidn
    auto output = tanh(matrixToArray(arrayToMatrix(hidunits) * weights2)); // samplesn x outputn
    auto diff = cwiseSelect(cwiseProd(output, labels), 0.0, (output - labels));
    auto energy = (diff * diff).sum();

    Expression<MatrixXd> grad1, grad2;
    std::tie(grad1, grad2) = energy.derivatives(weights1, weights2);

    MatrixXd lastchange1 = MatrixXd::Zero(weightVals1.rows(), weightVals1.cols());
    MatrixXd lastchange2 = MatrixXd::Zero(weightVals2.rows(), weightVals2.cols());
    //MatrixXd lastchangeb1 = MatrixXd::Zero(biasVals1.rows(), biasVals1.cols());
    //MatrixXd lastchangeb2 = MatrixXd::Zero(biasVals2.rows(), biasVals2.cols());

    double delta = 1;

    for (int i = 0; i < 200; i++) {
        std::cout << "energy: " << energy.execute(weights1, weights2) << std::endl;
        MatrixXd g1 = grad1.execute(weights1, weights2, energy);
        MatrixXd g2 = grad2.execute(weights1, weights2, energy);
        //MatrixXd gb1 = gradb1.execute(weights1, weights2, bias1, bias2, energy);
        //MatrixXd gb2 = gradb2.execute(weights1, weights2, bias1, bias2, energy);
        double a = 0.6, b = 0.3;
        MatrixXd change1 = (-g1 * a +  lastchange1*b);
        MatrixXd change2 = (-g2 * a +  lastchange2*b);
        //MatrixXd changeb1 = (-gb1*a + lastchangeb1*b);
        //MatrixXd changeb2 = (-gb2*a + lastchangeb2*b);
        weightVals1 += change1 * delta;
        weightVals2 += change2 * delta;
        //biasVals1 += changeb1 * delta;
        //biasVals2 += changeb2 * delta;
        lastchange1 = change1;
        lastchange2 = change2;
        //lastchangeb1 = changeb1;
        //lastchangeb2 = changeb2;
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



