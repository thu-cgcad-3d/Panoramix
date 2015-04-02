#ifndef PANORAMIX_MISC_MATLAB_HPP
#define PANORAMIX_MISC_MATLAB_HPP

#include "../core/basic_types.hpp"


namespace panoramix {
    namespace misc {

        using CVInputArray = cv::InputArray;
        using CVOutputArray = cv::OutputArray;

        class Matlab {
        public:
            static bool IsBuilt();
            static std::string DefaultCodeDir();

            static bool IsUsable();
            static bool RunScript(const char * cmd);
            static bool RunScript(const std::string & cmd) { return RunScript(cmd.data()); }
            static const char * LastMessage();
            static inline void PrintLastMessage() { std::cout << LastMessage(); }

            static inline bool CDAndAddAllSubfolders(const std::string & dir){
                return RunScript("cd " + dir) && RunScript("addpath(genpath('.'));");
            }

            static bool PutVariable(const char * name, CVInputArray a);
            static void * PutVariable(CVInputArray a);
            static bool GetVariable(const char * name, CVOutputArray a, bool lastDimIsChannel = true);
            static bool GetVariable(const void * mxa, CVOutputArray a, bool lastDimIsChannel = true);
            static bool GetVariable(const char * name, double & a);
            static bool GetVariable(const void * mxa, double & a);

            static bool PutVariable(const char * name, const cv::SparseMat & mat);
            static void * PutVariable(const cv::SparseMat & mat);
            
            static inline bool PutVariable(const std::string & name, CVInputArray a){
                return PutVariable(name.data(), a); 
            }
            static inline bool GetVariable(const std::string & name, CVOutputArray a, bool lastDimIsChannel = true) {
                return GetVariable(name.data(), a, lastDimIsChannel); 
            }           

        public:
            inline Matlab & operator << (const std::string & cmd) { 
                RunScript(cmd);
                return *this;
            }
        };

    }
}
 
#endif