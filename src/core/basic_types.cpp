#include "basic_types.hpp"

namespace panoramix {
    namespace core {

		void ResizeToMakeWidthUnder(Image & im, int widthUpperBound) {
			if (im.cols <= widthUpperBound)
				return;
			cv::resize(im, im, cv::Size(widthUpperBound, static_cast<int>(im.rows * widthUpperBound / im.cols)));
		}

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