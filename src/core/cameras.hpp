#ifndef PANORAMIX_CORE_CAMERAS_HPP
#define PANORAMIX_CORE_CAMERAS_HPP

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
            inline double screenWidth() const { return _screenW; }
            inline double screenHeight() const { return _screenH; }
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
            bool isVisibleOnScreen(const Vec3 & p3d) const;
            HPoint2 screenProjectionInHPoint(const Vec3 & p3d) const;
            Vec3 spatialDirection(const Vec2 & p2d) const;
            inline Vec3 spatialDirection(const PixelLoc & p) const { return spatialDirection(Vec2(p.x, p.y)); }

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

            template <class Archive> inline void serialize(Archive & ar) {
                ar(_screenW, _screenH);
                ar(_focal, _near, _far);
                ar(_eye, _center, _up);
                ar(_viewMatrix, _projectionMatrix, _viewProjectionMatrix, _viewProjectionMatrixInv);
            }

        private:
            void updateMatrices();

        private:
            double _screenW, _screenH;
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

            template <class Archive> inline void serialize(Archive & ar) {
                ar(_focal);
                ar(_eye, _center, _up);
                ar(_xaxis, _yaxis, _zaxis);
            }

        private:
            double _focal;
            Vec3 _eye, _center, _up;
            Vec3 _xaxis, _yaxis, _zaxis;
        };


        // sample image from image using camera conversion
        template <class OutCameraT, class InCameraT>
        class CameraSampler {
        public:
            CameraSampler(const OutCameraT & outCam, const InCameraT & inCam)
                : _outCam(outCam), _inCam(inCam) {
                assert(outCam.eye() == inCam.eye());
                _mapx = cv::Mat::zeros(_outCam.screenSize(), CV_32FC1);
                _mapy = cv::Mat::zeros(_outCam.screenSize(), CV_32FC1);
                for (int j = 0; j < _outCam.screenSize().height; j++) {
                    for (int i = 0; i < _outCam.screenSize().width; i++) {
                        Vec2 screenp(i, j);
                        Vec3 p3 = _outCam.spatialDirection(screenp);
                        Vec2 screenpOnInCam = _inCam.screenProjection(p3);
                        _mapx.at<float>(j, i) = static_cast<float>(screenpOnInCam(0));
                        _mapy.at<float>(j, i) = static_cast<float>(screenpOnInCam(1));
                    }
                }
            }

            Image operator() (const Image & inputIm,
                int borderMode = cv::BORDER_REPLICATE,
                const cv::Scalar & borderValue = cv::Scalar(0, 0, 0, 1)) const {
                Image outputIm;
                cv::remap(inputIm, outputIm, _mapx, _mapy, cv::INTER_LINEAR, borderMode, borderValue);
                return outputIm;
            }

            template <class Archive>
            inline void serialize(Archive & ar) {
                ar(_outCam, _inCam, _mapx, _mapy);
            }

        private:
            OutCameraT _outCam;
            InCameraT _inCam;
            cv::Mat _mapx, _mapy;
        };





    }
}
 
#endif