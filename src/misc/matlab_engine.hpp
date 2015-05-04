#ifndef PANORAMIX_MISC_MATLAB_ENGINE_HPP
#define PANORAMIX_MISC_MATLAB_ENGINE_HPP

#include "../core/basic_types.hpp"

namespace panoramix {
    namespace misc {


        using CVInputArray = cv::InputArray;
        using CVOutputArray = cv::OutputArray;


        class MatlabEngine {
        public:
            static bool IsBuilt();

            static bool Start();
            static void Close();
            static bool Started();

            static std::string DefaultCodeDir();

            static bool RunScript(const char * cmd);
            static bool RunScript(const std::string & cmd) { return RunScript(cmd.data()); }
            static const char * LastMessage();
            static inline void PrintLastMessage() { std::cout << LastMessage(); }

            static inline bool CDAndAddAllSubfolders(const std::string & dir){
                return RunScript("cd " + dir) && RunScript("addpath(genpath('.'));");
            }
            static inline bool Load(const std::string & matfile){
                return RunScript("load('" + matfile + "');");
            }

            static bool PutVariable(const char * name, CVInputArray a);
            static bool GetVariable(const char * name, CVOutputArray a, bool lastDimIsChannel = true);
            static bool GetVariable(const char * name, double & a);
            static bool GetVariable(const char * name, std::string & t);
            static bool PutVariable(const char * name, const cv::SparseMat & mat);
            
            static inline bool PutVariable(const std::string & name, CVInputArray a){
                return PutVariable(name.data(), a); 
            }
            static inline bool GetVariable(const std::string & name, CVOutputArray a, bool lastDimIsChannel = true) {
                return GetVariable(name.data(), a, lastDimIsChannel); 
            }           

        public:
            MatlabEngine() { Start(); }
            ~MatlabEngine() { Close(); }
            inline MatlabEngine & operator << (const std::string & cmd) { 
                RunScript(cmd);
                return *this;
            }
        };

    }
}
 
#endif