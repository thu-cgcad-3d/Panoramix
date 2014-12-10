#include "basic_types.hpp"

namespace panoramix {
    namespace core {

		void ResizeToMakeWidthUnder(Image & im, int widthUpperBound) {
			if (im.cols <= widthUpperBound)
				return;
			cv::resize(im, im, cv::Size(widthUpperBound, static_cast<int>(im.rows * widthUpperBound / im.cols)));
		}

        void ResizeToMakeHeightUnder(Image & im, int heightUpperBound) {
            if (im.rows <= heightUpperBound)
                return;
            cv::resize(im, im, cv::Size(static_cast<int>(im.cols * heightUpperBound / im.rows), heightUpperBound));
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



        PixelLoc PixelLocFromGeoCoord(const GeoCoord & p, int longidiv, int latidiv) {
            int longtid = static_cast<int>((p.longitude + M_PI) * longidiv / M_PI / 2);
            int latid = static_cast<int>((p.latitude + M_PI_2) * latidiv / M_PI);
            longtid = (longtid % longidiv + longidiv) % longidiv;
            latid = (latid % latidiv + latidiv) % latidiv;
            return PixelLoc(longtid, latid);
        }

        GeoCoord GeoCoordFromPixelLoc(const cv::Point & pixel, int longidiv, int latidiv) {
            return GeoCoord{ pixel.x * M_PI * 2 / longidiv - M_PI, pixel.y * M_PI / latidiv - M_PI_2 };
        }


    }
}