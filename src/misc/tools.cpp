#include "matlab.hpp"
#include "tools.hpp"

namespace panoramix {
    namespace misc {

        namespace dataset {

            namespace YorkUrbanDB {

                GroundTruth LoadGroundTruth(const std::string & camParams){
                    GroundTruth gt;
                    Matlab matlab;
                    matlab << ("temp = load('" + camParams + "');");
                    matlab << "temp = temp.vp_orthogonal;";
                    matlab.GetVariable("temp", gt.vp_orthogonal);
                    return gt;
                }


            }



        }

    }
}