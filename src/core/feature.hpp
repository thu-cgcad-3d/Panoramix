#ifndef PANORAMIX_CORE_FEATURE_HPP
#define PANORAMIX_CORE_FEATURE_HPP

#include <opencv2/opencv.hpp>
 
namespace panoramix {
	namespace core {
        
        using cv::Mat;
        using cv::Vec2d;
        using cv::Vec3d;

        template <class ValueT, int dim>
        struct LineData {
        	cv::Vec<ValueT, dim> p1, p2;
        };
        
		void ExtractLines(const Mat& im,
                          std::list<LineData<double, 2>> & lines,
                          int minlen = 20, int xborderw = 10, int yborderw = 20,
                          int numdir = 8);
        
 
        
        
        
	}
}
 
#endif