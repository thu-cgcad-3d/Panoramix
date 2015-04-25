#ifndef PANORAMIX_MISC_MATLAB_HPP
#define PANORAMIX_MISC_MATLAB_HPP

#include "../core/basic_types.hpp"

namespace panoramix {
    namespace misc {

        struct AnyPtr {
            uintptr_t data;
            AnyPtr() : data(0) {}
            template <class T> AnyPtr(T * d) : data(static_cast<uintptr_t>(d)) {}
            template <class T> operator T* () const { return static_cast<T*>(data); }
            bool operator == (nullptr_t) const { return data == 0; }
            bool operator != (nullptr_t) const { return data != 0; }
        };


        using MXAPtr = AnyPtr;


        using CVInputArray = cv::InputArray;
        using CVOutputArray = cv::OutputArray;


        class Matlab {
        public:
            static bool IsBuilt();

            static bool StartEngine();
            static void CloseEngine();
            static bool EngineStarted();

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
            static void * PutVariable(CVInputArray a);
            static bool GetVariable(const char * name, CVOutputArray a, bool lastDimIsChannel = true);
            static bool GetVariable(const void * mxa, CVOutputArray a, bool lastDimIsChannel = true);
            static bool GetVariable(const char * name, double & a);
            static bool GetVariable(const void * mxa, double & a);
            static bool GetVariable(const char * name, std::string & t);
            static bool GetVariable(const void * mxa, std::string & t);

            static bool PutVariable(const char * name, const cv::SparseMat & mat);
            static void * PutVariable(const cv::SparseMat & mat);
            
            static inline bool PutVariable(const std::string & name, CVInputArray a){
                return PutVariable(name.data(), a); 
            }
            static inline bool GetVariable(const std::string & name, CVOutputArray a, bool lastDimIsChannel = true) {
                return GetVariable(name.data(), a, lastDimIsChannel); 
            }           

        public:
            Matlab() { StartEngine(); }
            ~Matlab() { CloseEngine(); }
            inline Matlab & operator << (const std::string & cmd) { 
                RunScript(cmd);
                return *this;
            }
        };

    }
}
 
#endif