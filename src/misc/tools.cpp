#include "matlab_engine.hpp"
#include "tools.hpp"

namespace panoramix {
    namespace misc {

        namespace dataset {

            namespace YorkUrbanDB {

                GroundTruth LoadGroundTruth(const std::string & camParams){
                    GroundTruth gt;
                    MatlabEngine matlab;
                    matlab << ("temp = load('" + camParams + "');");
                    matlab << "temp = temp.vp_orthogonal;";
                    matlab.GetVariable("temp", gt.vp_orthogonal);
                    return gt;
                }


            }



        }

    }
}