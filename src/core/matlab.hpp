#ifndef PANORAMIX_CORE_MATLAB_HPP
#define PANORAMIX_CORE_MATLAB_HPP

#include "basic_types.hpp"

namespace panoramix {
    namespace core {

        using CVInputArray = cv::InputArray;
        using CVOutputArray = cv::OutputArray;

        class Matlab {
        public:
            static bool IsBuilt();

            static bool IsUsable();
            static bool RunScript(const char * cmd);
            static bool RunScript(const std::string & cmd) { return RunScript(cmd.data()); }
            static const char * LastMessage();

            static inline bool CDAndAddAllSubfolders(const std::string & dir){
                return RunScript("cd " + dir) && RunScript("addpath(genpath('.'));");
            }

            static bool PutVariable(const char * name, CVInputArray a);
            static bool GetVariable(const char * name, CVOutputArray a, bool lastDimIsChannel = true);
            
            static inline bool PutVariable(const std::string & name, CVInputArray a){
                PutVariable(name.data(), a); 
            }
            static inline bool GetVariable(const std::string & name, CVOutputArray a, bool lastDimIsChannel = true) {
                GetVariable(name.data(), a, lastDimIsChannel); 
            }           

        public:
            inline Matlab & operator << (const std::string & cmd) { 
                RunScript(cmd + ";");
                return *this;
            }
        };

    }
}
 
#endif