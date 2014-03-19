#ifndef PANORAMIX_CORE_FEATURE_HPP
#define PANORAMIX_CORE_FEATURE_HPP

#include <Eigen/Core>
#include <opencv2/opencv.hpp>
 
namespace panoramix {
	namespace core {
        
		// perspective camera
		class Camera {
		public:
			using Vec2 = Eigen::Vector2d;
			using Vec3 = Eigen::Vector3d;
			using Vec4 = Eigen::Vector4d;
			using Mat4 = Eigen::Matrix4d;
			using Size2 = cv::Size2i;

		public:
			explicit Camera(int w, int h, double focal, const Vec3 & eye = Vec3(0, 0, 0),
				const Vec3 & center = Vec3(1, 0, 0), const Vec3 & up = Vec3(0, 0, 1));

			inline Size2 screenSize() const { return Size2(_screenW, _screenH); }
			inline const Vec3 & eye() const { return _eye; }
			inline const Vec3 & center() const { return _center; }
			inline const Vec3 & up() const { return _up; }
			Vec2 screenProjection(const Vec3 & p3d) const;
			Vec3 spatialDirection(const Vec2 & p2d) const;

			inline const Mat4 & viewMatrix() const { return _viewMatrix; }
			inline const Mat4 & projectionMatrix() const { return _viewProjectionMatrix; }
			inline const Mat4 & viewProjectionMatrix() const { return _viewProjectionMatrix; }
			inline const Mat4 & viewProjectionMatrixInv() const { return _viewProjectionMatrixInv; }

		private:
			int _screenW, _screenH;
			double _focal;
			Vec3 _eye, _center, _up;
			Mat4 _viewMatrix, _projectionMatrix, _viewProjectionMatrix, _viewProjectionMatrixInv;
		};

		// panoramic camera
		class PanoramicCamera {
		public:
			using Vec2 = Eigen::Vector2d;
			using Vec3 = Eigen::Vector3d;
			using Vec4 = Eigen::Vector4d;
			using Mat4 = Eigen::Matrix4d;
			using Size2 = cv::Size2f;

		public:
			explicit PanoramicCamera(double focal, const Vec3 & eye = Vec3(0, 0, 0),
				const Vec3 & center = Vec3(1, 0, 0), const Vec3 & up = Vec3(0, 0, 1));

			inline Size2 screenSize() const { return Size2(_focal * 2 * M_PI, _focal * M_PI); }
			inline const Vec3 & eye() const { return _eye; }
			inline const Vec3 & center() const { return _center; }
			inline const Vec3 & up() const { return _up; }
			Vec2 screenProjection(const Vec3 & p3d) const;
			Vec3 spatialDirection(const Vec2 & p2d) const;

		private:
			double _focal;
			Vec3 _eye, _center, _up;
		};


		// sample image from image
		template <class OutCameraT, class InCameraT>
		class CameraSampler {
		public:
			using Image = cv::Mat;
		public:
			explicit inline CameraSampler(const OutCameraT & outCam, const InCameraT & inCam)
				: _outCam(outCam), _inCam(inCam) {
				assert(outCam.eye() == inCam.eye());
			}
			Image operator() (const Image & inputIm) const;
		private:
			OutCameraT _outCam;
			InCameraT _inCam;
		};

		template <class OutCameraT, class InCameraT>
		typename CameraSampler<OutCameraT, InCameraT>::Image CameraSampler<OutCameraT, InCameraT>::operator() (const Image & inputIm) const {
			Image undistorted_im;
			cv::Mat mapx = cv::Mat::zeros(_outCam.screenSize(), CV_32FC1);
			cv::Mat mapy = cv::Mat::zeros(_outCam.screenSize(), CV_32FC1);
			for (int i = 0; i < _outCam.screenSize().width; i++){
				for (int j = 0; j < _outCam.screenSize().height; j++){
					Camera::Vec2 screenp(i, j);
					Camera::Vec3 p3 = _outCam.spatialDirection(screenp);
					Camera::Vec2 screenpOnInCam = _inCam.screenProjection(p3);
					mapx.at<float>(j, i) = screenpOnInCam(0);
					mapy.at<float>(j, i) = screenpOnInCam(1);
				}
			}
			cv::remap(inputIm, undistorted_im, mapx, mapy, cv::INTER_LINEAR, cv::BORDER_REPLICATE);
			return undistorted_im;
		}


		





		

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

		void EstimateOrientationMap(const Mat& im, Mat & omap);

		void EstimateGeometricContext(const Mat& im, Mat& gcim);

		void EstimateManhattanJunctionDistribution(const Mat & im, Mat & mjim);

	}
}
 
#endif