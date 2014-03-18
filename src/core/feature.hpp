#ifndef PANORAMIX_CORE_FEATURE_HPP
#define PANORAMIX_CORE_FEATURE_HPP

#include <Eigen/Core>
#include <opencv2/opencv.hpp>
 
namespace panoramix {
	namespace core {
        
		class Camera {
		public:
			using Vec2 = Eigen::Vector2d;
			using Vec3 = Eigen::Vector3d;
			using Vec4 = Eigen::Vector4d;
			using Mat4 = Eigen::Matrix4d;

		public:
			Camera(int w, int h, double focal, const Vec3 & eye, 
				const Vec3 & center = Vec3(0, 0, 0), const Vec3 & up = Vec3(0, 1, 0));
			
			Vec2 screenProjection(const Vec3 & p3d) const;
			Vec3 spatialDirection(const Vec2 & p2d) const;

		private:
			int _screenW, _screenH;
			double _focal;
			Vec3 _eye, _center, _up;
			Mat4 _viewMatrix, _projectionMatrix, _viewProjectionMatrix, _viewProjectionMatrixInv;
		};


		class PerspectiveSampler {
		public:
			using Image = cv::Mat;
		public:
			inline PerspectiveSampler(const Camera & cam) : _cam(cam) {}
			Image operator() (const Image & panorama) const;
		private:
			Camera _cam;
		};

		

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