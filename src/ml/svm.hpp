#ifndef PANORAMIX_ML_SVM_HPP
#define PANORAMIX_ML_SVM_HPP

#include <opencv2/ml/ml.hpp>
#include "../core/basic_types.hpp"

namespace panoramix {
    namespace ml {

        class SVM {
        public:
            struct Params {

            };
            
        public:

            void train();
            void predict() const ;

        private:

        };

    }
}


#endif