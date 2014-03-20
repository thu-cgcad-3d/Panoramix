#ifndef PANORAMIX_CORE_BASIC_TYPES_HPP
#define PANORAMIX_CORE_BASIC_TYPES_HPP

#include <vector>
#include <list>
#include <deque>
#include <set>
#include <unordered_set>
#include <forward_list>
#include <array>

#include <opencv2/opencv.hpp>
#include <Eigen/Core>

namespace panoramix {
	namespace core {
 
		template <class T, int N> using Vec = cv::Vec<T, N>;
		using Vec2 = Vec<double, 2>;
		using Vec3 = Vec<double, 3>;
		using Vec4 = Vec<double, 4>;

		template <class T, int N> 
		struct Line {
			Vec<T, N> first, second;
		};
		using Line2 = Line<double, 2>;
		using Line3 = Line<double, 3>;

		using KeyPoint = cv::KeyPoint;
		
		template <class T, int N>
		struct Circle {
			Vec<T, N> center;
			T radius;
		};
		using Circle2 = Circle<double, 2>;
		using Circle3 = Circle<double, 3>;

		using Image = cv::Mat;

		using Color = cv::Scalar;
 
	}
}
 
#endif