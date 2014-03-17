#include "dummy.hpp"
#include "../config.hpp"

#include <iostream>

#include <Eigen/Core>
#include <opencv2/opencv.hpp>

namespace panoramix {
	namespace core {
 
		Dummy<int> MakeAnIntDummy() {
			using namespace std;

			cout << "Panoramix version: " << PANORAMIX_VERSION_MAJOR << "." << PANORAMIX_VERSION_MINOR << endl;

		    Eigen::Array3cf arr;
		    arr << 1.0f, 2.0f, 3.0f;
		    cout << "Eigen test" << endl;
		    cout << arr << endl;

		    cv::Mat im = cv::Mat::eye(500, 500, CV_8UC1);
		    im =  255 - im * 255;

			return Dummy<int>{PANORAMIX_VERSION_MAJOR, PANORAMIX_VERSION_MINOR};
		}
 
	}
}

