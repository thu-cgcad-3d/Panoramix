#ifndef PANORAMIX_CORE_CLOCK_HPP
#define PANORAMIX_CORE_CLOCK_HPP

#include <string>
#include <chrono>
#include <iostream>

namespace panoramix {
    namespace core {

        class Clock {
        public:
            Clock(const std::string & msg) : _message(msg) {
                _startTime = std::chrono::system_clock::now();
            }
            ~Clock(){
                auto duration = std::chrono::system_clock::now() - _startTime;
                std::cout << "[" << _message << "] Time Elapsed: "
                    << std::chrono::duration_cast<std::chrono::milliseconds>(duration).count()
                    << " ms" << std::endl;
            }
        private:
            std::chrono::system_clock::time_point _startTime;
            std::string _message;
        };

   	}
}

#endif