#ifndef PANORAMIX_CORE_MATLAB_HPP
#define PANORAMIX_CORE_MATLAB_HPP

#include "basic_types.hpp"

namespace panoramix {
    namespace core {

        class Matlab {
        public:
            static bool IsBuilt();

            static bool IsUsable();
            static bool RunScript(const char * cmd);

            static bool PutVariable(const char * name, const Image & im);
            static bool GetVariable(const char * name, Image & im);

        public:
            inline Matlab & operator << (const std::string & cmd) { 
                RunScript((cmd + ";").data());
                return *this; 
            }
        };

    }
}
 
#endif