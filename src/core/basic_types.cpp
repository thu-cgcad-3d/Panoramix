#include "basic_types.hpp"

namespace panoramix {
    namespace core {

        std::pair<PixelLoc, PixelLoc> MinMaxLocOfImage(const Image & im) {
            PixelLoc minLoc, maxLoc;
            double minVal, maxVal;
            cv::minMaxLoc(im, &minVal, &maxVal, &minLoc, &maxLoc);
            return std::make_pair(minLoc, maxLoc);
        }
        
        std::pair<double, double> MinMaxValOfImage(const Image & im) {
            PixelLoc minLoc, maxLoc;
            double minVal, maxVal;
            cv::minMaxLoc(im, &minVal, &maxVal, &minLoc, &maxLoc);
            return std::make_pair(minVal, maxVal);
        }

    }
}