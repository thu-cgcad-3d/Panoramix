#ifndef PANORAMIX_DERIV_EXPRESSION_OPS_SIGNAL_HPP
#define PANORAMIX_DERIV_EXPRESSION_OPS_SIGNAL_HPP

#include <Eigen/Dense>
#include <Eigen/StdVector>

#include "expression_ops.hpp"

namespace panoramix {
    namespace deriv {

        template <class S, int R, int C>
        using CovolutionalFilterWeights = std::vector<Eigen::Matrix<S, R, C>>;

        

    }
}
 
#endif