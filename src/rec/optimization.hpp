#ifndef PANORAMIX_REC_OPTIMIZATION_HPP
#define PANORAMIX_REC_OPTIMIZATION_HPP
 
#include <Eigen/Core>
#include <Eigen/Dense>
#include <Eigen/Sparse>

#include <unsupported/Eigen/NonLinearOptimization>

#include <glpk.h>
#include <setjmp.h>

#include "../deriv/derivative.hpp"

#include "../core/basic_types.hpp"
#include "../core/utilities.hpp"


namespace panoramix {
    namespace rec {

        struct sinfo {
            char * text;
            jmp_buf * env;
        };

        void glpErrorHook(void * in){
            sinfo * info = (sinfo*)in;
            glp_free_env();
            longjmp(*(info->env), 1);
        }
    }
}
 
#endif