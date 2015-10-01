#include "eval.hpp"

namespace panolyz {

    std::vector<Vec3> FibonacciDirections(int N) {
        std::vector<Vec3> directions;
        directions.reserve(N);
        GenerateFibonacciDirections(N, [&directions](double x, double y, double z) {
            directions.push_back(normalize(Vec3(x, y, z)));
        });
        return directions;
    }

    double AverageVariation(const ReconstructedModel & m1, const ReconstructedModel & m2, const std::vector<Vec3> & directions) {
        DenseMatd depths1(directions.size(), 1, 0.0), depths2(directions.size(), 1, 0.0);
        int invalidDirectionsNum = 0;
        for (int i = 0; i < directions.size(); i++) {
            if (m1.isValid(directions[i]) && m2.isValid(directions[i])) {
                depths1(i) = m1.depthAt(directions[i]);
                depths2(i) = m2.depthAt(directions[i]);
            } else {
                invalidDirectionsNum++;
            }
        }

        double dist = norm(normalize(depths1) - normalize(depths2));
        std::cout << "invalid directions num: " << invalidDirectionsNum << std::endl;
        return dist;
    }

}