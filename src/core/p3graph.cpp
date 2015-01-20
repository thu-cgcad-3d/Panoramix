
#include "iterators.hpp"
#include "p3graph.hpp"

namespace panoramix {
    namespace core {

        double InverseDepthAtDirection(const Vec3 & direction, double v1, double v2,
            const P3FreeLine & u){
            const auto & line = u.line;
            double theta = AngleBetweenDirections(line.first, line.center());
            double phi = AngleBetweenDirections(line.second, direction);
            /*           | sin(theta) | | p | | q |
            len:---------------------------------   => 1/len: [(1/q)sin(phi) - (1/p)sin(phi-theta)] / sin(theta)
            | p sin(phi) - q sin(phi - theta) |
            */
            // variables[0] -> 1/p
            // variables[1] -> 1/q
            double coeffFor1_p = -sin(phi - theta) / sin(theta);
            double coeffFor1_q = sin(phi) / sin(theta);
            assert(!IsInfOrNaN(coeffFor1_p) && !IsInfOrNaN(coeffFor1_q));
            assert(!IsInfOrNaN(coeffFor1_p) && !IsInfOrNaN(coeffFor1_q));
            return v1 * coeffFor1_p + v2 * coeffFor1_q;
        }

        //double InverseDepthAtDirection(const Vec3 & direction, double v1, double v2,
        //    const P3OrientedLine & u, const std::vector<Vec3> & vps){

        //}

    }
}