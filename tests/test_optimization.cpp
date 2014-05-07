#include <iostream>
#include <random>

#include "../src/deriv/derivative.hpp"
#include "gtest/gtest.h"

using namespace Eigen;
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

    for (int i = 0; i < 20; i++) {
        std::cout << "energy: " << energy.execute(x) << std::endl;
        auto g = grad.execute(x, energy);
        xv -= g;
    }

}


TEST(Optimization, DISABLED_ANN) {
//void run(){

    int inputn = 100, hidn = 100, outputn = 40;

    ExpressionGraph graph;
    MatrixXd weightVals1, weightVals2;
    
    auto inputs = graph.addConst(MatrixXd::Random(1, inputn));
    
    ArrayXXd arr = ArrayXXd::Random(1, outputn);
    arr = arr / arr.cwiseAbs().maxCoeff();
    auto labels = graph.addConst(arr);
    
    weightVals1.setRandom(inputn, hidn);
    weightVals2.setRandom(hidn, outputn);
    auto weights1 = graph.addRef(weightVals1);
    auto weights2 = graph.addRef(weightVals2);

    auto hidunits = sigmoid(matrixToArray(inputs * weights1));
    auto diff = sigmoid(matrixToArray(arrayToMatrix(hidunits) * weights2)) - labels;
    auto energy = (diff * diff).sum();
    Expression<MatrixXd> grad1, grad2;
    std::tie(grad1, grad2) = energy.derivatives(weights1, weights2);

    MatrixXd lastchange1 = MatrixXd::Zero(weightVals1.rows(), weightVals1.cols());
    MatrixXd lastchange2 = MatrixXd::Zero(weightVals2.rows(), weightVals2.cols());
    for (int i = 0; i < 20; i++) {
        std::cout << "energy: " << energy.execute(weights1, weights2) << std::endl;
        MatrixXd g1 = grad1.execute(weights1, weights2, energy);
        MatrixXd g2 = grad2.execute(weights1, weights2, energy);
        MatrixXd change1 = (-g1*0.3 + lastchange1*0.7);
        MatrixXd change2 = (-g2*0.3 + lastchange2*0.7);
        weightVals1 += change1;
        weightVals2 += change2;
        lastchange1 = change1;
        lastchange2 = change2;
    }

}





int main(int argc, char * argv[], char * envp[])
{
	testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
    //run();
    return 1;
}



