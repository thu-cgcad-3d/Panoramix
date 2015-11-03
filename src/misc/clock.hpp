#pragma once

#include <string>
#include <chrono>
#include <ctime>
#include <iostream>

namespace pano {
    namespace misc {

        class Clock {
        public:
            Clock(const std::string & msg);
            ~Clock();
        private:
            std::chrono::system_clock::time_point _startTime;
            std::string _message;
        };

#define SetClock() pano::misc::Clock clock##__COUNTER__(__FUNCTION__)
        
        std::string CurrentTimeString(bool tagified = false);

   	}
}

