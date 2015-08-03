#pragma once

#include <string>
#include <chrono>
#include <iostream>

namespace pano {
    namespace core {

        class Clock {
        public:
            inline Clock(const std::string & msg) : _message(msg) {
                std::cout << "[" << _message << "] Started." << std::endl;
                _startTime = std::chrono::system_clock::now();
            }
            ~Clock(){
                auto duration = std::chrono::system_clock::now() - _startTime;
                std::cout << "[" << _message << "] Stopped. Time Elapsed: "
                    << std::chrono::duration_cast<std::chrono::milliseconds>(duration).count()
                    << " ms" << std::endl;
            }
        private:
            std::chrono::system_clock::time_point _startTime;
            std::string _message;
        };

#define SetClock() Clock clock##__COUNTER__(__FUNCTION__)

   	}
}

