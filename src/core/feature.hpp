#ifndef PANORAMIX_CORE_FEATURE_HPP
#define PANORAMIX_CORE_FEATURE_HPP

#include "basic_types.hpp"

#include <opencv2/features2d/features2d.hpp>
#include <opencv2/nonfree/features2d.hpp>
 
namespace panoramix {
	namespace core {
        
		// perspective camera
		class PerspectiveCamera {
		public:
			EIGEN_MAKE_ALIGNED_OPERATOR_NEW
			using Vec2 = Eigen::Vector2d;
			using Vec3 = Eigen::Vector3d;
			using Vec4 = Eigen::Vector4d;
			using Mat4 = Eigen::Matrix4d;
			using Size2 = cv::Size2f;

		public:
			explicit PerspectiveCamera(int w, int h, double focal, const Vec3 & eye = Vec3(0, 0, 0),
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
			EIGEN_MAKE_ALIGNED_OPERATOR_NEW
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
			Vec3 _xaxis, _yaxis, _zaxis;
		};


		// sample image from image using camera conversion
		template <class OutCameraT, class InCameraT>
		class CameraSampler {
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
		Image CameraSampler<OutCameraT, InCameraT>::operator() (const Image & inputIm) const {
			Image result;
			cv::Mat mapx = cv::Mat::zeros(_outCam.screenSize(), CV_32FC1);
			cv::Mat mapy = cv::Mat::zeros(_outCam.screenSize(), CV_32FC1);
			for (int i = 0; i < _outCam.screenSize().width; i++){
				for (int j = 0; j < _outCam.screenSize().height; j++){
					typename OutCameraT::Vec2 screenp(i, j);
					typename OutCameraT::Vec3 p3 = _outCam.spatialDirection(screenp);
					typename InCameraT::Vec2 screenpOnInCam = _inCam.screenProjection(p3);
					mapx.at<float>(j, i) = screenpOnInCam(0);
					mapy.at<float>(j, i) = screenpOnInCam(1);
				}
			}
			cv::remap(inputIm, result, mapx, mapy, cv::INTER_LINEAR, cv::BORDER_REPLICATE);
			return result;
		}




		

		/// features


		// line extractor
		class LineSegmentExtractor {
		public:
			using Feature = std::list<Line2>;
			struct Params {
				inline Params() : minLength(20), xBorderWidth(10), yBorderWidth(10), numDirs(8) {}
				int minLength;
				int xBorderWidth, yBorderWidth;
				int numDirs;
			};
		public:
			inline explicit LineSegmentExtractor(const Params & params = Params()) : _params(params){}
			Feature operator() (const Image & im) const;
		private:
			Params _params;
		};

		// point feature extractor
		template <class CVFeatureT>
		class CVFeatureExtractor {
		public:
			using Feature = std::vector<cv::KeyPoint>;
		public:
			template <class ... ArgT>
			explicit inline CVFeatureExtractor(ArgT ... args) : _feature2D(args ...) {}
			inline Feature operator() (const Image & im, const Image & mask = cv::Mat(), 
				cv::OutputArray descriptors = cv::noArray()) const {
				Feature fea;
				_feature2D(im, mask, fea, descriptors);
				return fea;
			}
		private:
			CVFeatureT _feature2D;
		};

		// region feature extractor
		// manhattan junction
		class ManhattanJunctionExtractor {
		public:
			using Feature = Image;
			struct Params {
				inline Params() {}
			};
		public:
			inline explicit ManhattanJunctionExtractor(const Params & params = Params()) : _params(params) {}
			Feature operator () (const Image & im) const { return Feature(); }
		private:
			Params _params;
		};

		// geometric context
		class GeometricContextExtractor {
		public:
			using Feature = Image;
			struct Params {
				inline Params() {}
			};
		public:
			inline explicit GeometricContextExtractor(const Params & params = Params()) : _params(params) {}
			Feature operator () (const Image & im) const { return Feature(); }
		private:
			Params _params;
		};

	}
}
 
#endif