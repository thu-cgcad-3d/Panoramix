#ifndef PANORAMIX_CORE_FEATURE_HPP
#define PANORAMIX_CORE_FEATURE_HPP

#include <opencv2/features2d/features2d.hpp>
#include <opencv2/nonfree/features2d.hpp>

#include "basic_types.hpp"
 
namespace panoramix {
    namespace core {
        
        // perspective camera
        class PerspectiveCamera {
        public:
            explicit PerspectiveCamera(int w = 500, int h = 500, double focal = 250, 
                const Vec3 & eye = Vec3(0, 0, 0),
                const Vec3 & center = Vec3(1, 0, 0), 
                const Vec3 & up = Vec3(0, 0, 1),
                double nearPlane = 0.01, double farPlane = 1e4);

            inline Size screenSize() const { return Size(static_cast<float>(_screenW), static_cast<float>(_screenH)); }
            inline double fovRadians() const { return atan(_screenH / 2.0 / _focal) * 2; }
            inline double fovAngles() const { return fovRadians() * 180.0 / M_PI; }
            inline double aspect() const { return double(_screenW) / double(_screenH); }
            inline double focal() const { return _focal; }
            inline const Vec3 & eye() const { return _eye; }
            inline const Vec3 & center() const { return _center; }
            inline const Vec3 & up() const { return _up; }
            inline double nearPlane() const { return _near; }
            inline double farPlane() const { return _far; }
            Vec2 screenProjection(const Vec3 & p3d) const;
            Vec3 spatialDirection(const Vec2 & p2d) const;

            inline const Mat4 & viewMatrix() const { return _viewMatrix; }
            inline const Mat4 & projectionMatrix() const { return _projectionMatrix; }
            inline const Mat4 & viewProjectionMatrix() const { return _viewProjectionMatrix; }
            inline const Mat4 & viewProjectionMatrixInv() const { return _viewProjectionMatrixInv; }

            // operations
            void resizeScreen(const Size & sz, bool updateMat = true);
            void setFocal(double f, bool updateMat = true);
            void setEye(const Vec3 & e, bool updateMat = true);
            void setCenter(const Vec3 & c, bool updateMat = true);
            void setUp(const Vec3 & up, bool updateMat = true);
            void setNearAndFarPlanes(double nearPlane, double farPlane, bool updateMat = true);

        private:
            void updateMatrices();

        private:
            int _screenW, _screenH;
            double _focal;
            double _near, _far;
            Vec3 _eye, _center, _up;
            Mat4 _viewMatrix, _projectionMatrix, _viewProjectionMatrix, _viewProjectionMatrixInv;
        };

        // panoramic camera
        class PanoramicCamera {
        public:
            explicit PanoramicCamera(double focal = 250, 
                const Vec3 & eye = Vec3(0, 0, 0),
                const Vec3 & center = Vec3(1, 0, 0), 
                const Vec3 & up = Vec3(0, 0, 1));

            inline Size screenSize() const { return Size(static_cast<float>(_focal * 2 * M_PI), static_cast<float>(_focal * M_PI)); }
            inline double focal() const { return _focal; }
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
            
            Image operator() (const Image & inputIm) const {
                Image outputIm;
                cv::Mat mapx = cv::Mat::zeros(_outCam.screenSize(), CV_32FC1);
                cv::Mat mapy = cv::Mat::zeros(_outCam.screenSize(), CV_32FC1);
                for (int j = 0; j < _outCam.screenSize().height; j++){
                    for (int i = 0; i < _outCam.screenSize().width; i++){
                        Vec2 screenp(i, j);
                        Vec3 p3 = _outCam.spatialDirection(screenp);
                        Vec2 screenpOnInCam = _inCam.screenProjection(p3);
                        mapx.at<float>(j, i) = screenpOnInCam(0);
                        mapy.at<float>(j, i) = screenpOnInCam(1);
                    }
                }
                cv::remap(inputIm, outputIm, mapx, mapy, cv::INTER_LINEAR, cv::BORDER_REPLICATE);
                return outputIm;
            }

        private:
            OutCameraT _outCam;
            InCameraT _inCam;
        };






        

        /// features

        // line extractor
        class LineSegmentExtractor {
        public:
            using Feature = std::vector<Line2>;
            struct Params {
                inline Params() : minLength(15), xBorderWidth(3), yBorderWidth(3), numDirs(8), useExperimentalAlgorithm(false) {}
                int minLength;
                int xBorderWidth, yBorderWidth;
                int numDirs;
                bool useExperimentalAlgorithm;
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


        // segmentation
        class SegmentationExtractor {
        public:
            using Feature = Image; // CV_32SC1, from 0 to numRegion, use at<int32_t> to extract
            struct Params {
                inline Params() : sigma(0.8f), c(100.0f), minSize(100) {}
                float sigma; // for smoothing
                float c; // threshold function
                int minSize; // min component size
            };
        public:
            inline explicit SegmentationExtractor(const Params & params = Params()) : _params(params){}
            Feature operator() (const Image & im, bool forVisualization = false) const;
        private:
            Params _params;
        };

    }
}
 
#endif