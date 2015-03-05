#ifndef PANORAMIX_MISC_CMD_TOOLS_HPP
#define PANORAMIX_MISC_CMD_TOOLS_HPP

#include "../core/any.hpp"
#include "../core/basic_types.hpp"

namespace panoramix {
	namespace misc{

        struct CmdOption {
            std::string name;
            core::Any value;
            std::string description;
        };

        class CmdOptions {
        public:
            CmdOptions() {}
            CmdOptions(std::initializer_list<CmdOption> ilist);
            
            CmdOption & operator[](const std::string & name) { return _options[name]; }
            const CmdOption & operator[](const std::string & name) const { return _options.at(name); }

            template <class T>
            T value(const std::string & name) const { 
                return _options.at(name).value.cast<T>(); 
            }

            bool parseArguments(int argc, const char * const * argv);

        private:
            std::string _programName;
            std::map<std::string, CmdOption> _options;
        };


	}
}
 
#endif