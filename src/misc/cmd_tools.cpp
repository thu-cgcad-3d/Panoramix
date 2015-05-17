#include "../core/utility.hpp"
#include "cmd_tools.hpp"

namespace panoramix {
    namespace misc{

        CmdOptions::CmdOptions(std::initializer_list<CmdOption> ilist){
            for (auto && opt : ilist){
                _options[opt.name] = std::move(opt);
            }
        }

        bool CmdOptions::parseArguments(int argc, const char * const * argv) {
            if (argc <= 0){
                std::cout << "no arguments received!" << std::endl;
                return false;
            }
            _programName = argv[0];
            for (int i = 1; i < argc; i++){
                std::string arg = argv[i];
                if (arg.empty())
                    continue;
                if (arg[0] != '-' || arg.size() <= 1){
                    std::cout << "\"" << arg << "\" does not satisfy '-%name%=%value%' or -%name% (for boolean option)" << std::endl;
                    continue;
                }
                int eqPos = arg.find_first_of('=');
                std::string key = arg.substr(1, eqPos);
                if (!core::Contains(_options, key)){
                    std::cout << "\"" << arg << "\" invalid option name : " << key << std::endl;
                    continue;
                }
                CmdOption & opt = _options[key];
                if (eqPos == std::string::npos){
                    if (opt.value.is<bool>()){
                        opt.value = true;
                    }
                    else{
                        std::cout << "\"" << arg << "\" invalid option value" << std::endl;
                        continue;
                    }
                }
                std::string value = arg.substr(eqPos + 1);
                if (opt.value.is<std::string>()){
                    opt.value = value;
                }
                else if (opt.value.is<int>()){
                    opt.value = std::stoi(value);
                }
                else if (opt.value.is<double>()){
                    opt.value = std::stod(value);
                }
                else if (opt.value.is<bool>()){
                    if (value == "on"){
                        opt.value = true;
                    }
                    else if (value == "off"){
                        opt.value = false;
                    }
                    else{
                        std::cout << "\"" << arg << "\" invalid option value" << std::endl;
                        continue;
                    }
                }
                else{
                    std::cout << "\"" << arg << "\" option type is not supported" << std::endl;
                    continue;
                }
            }

            return true;
        }


	}
}