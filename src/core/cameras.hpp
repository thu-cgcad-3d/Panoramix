#ifndef PANORAMIX_CORE_CAMERAS_HPP
#define PANORAMIX_CORE_CAMERAS_HPP

#include "basic_types.hpp"
#include "feature.hpp"

namespace panoramix {
    namespace core {

        
        // perspective camera
        class PerspectiveCamera {
        public:
            explicit PerspectiveCamera(int w = 500, int h = 500, double focal = 250,
                const Vec3 & eye = Vec3(0, 0, 0),
                const Vec3 & center = Vec3(1, 0, 0),
                const Vec3 & up = Vec3(0, 0, -1),
                double nearPlane = 0.01, double farPlane = 1e4);

            inline SizeI screenSize() const { return Size(static_cast<int>(_screenW), static_cast<int>(_screenH)); }
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
            //inline const Mat4 & viewProjectionMatrixInv() const { return _viewProjectionMatrixInv; }

            // operations
            void resizeScreen(const Size & sz, bool updateMat = true);
            void setFocal(double f, bool updateMat = true);
            void setEye(const Vec3 & e, bool updateMat = true);
            void setCenter(const Vec3 & c, bool updateMat = true);
            void setUp(const Vec3 & up, bool updateMat = true);
            void setNearAndFarPlanes(double nearPlane, double farPlane, bool updateMat = true);

            // advanced functions
            inline Vec3 forward() const { return normalize(_center - _eye); }
            inline Vec3 backward() const { return normalize(_eye - _center); }
            inline Vec3 upward() const { return normalize(_up); }
            inline Vec3 downward() const { return -upward(); }
            inline Vec3 leftward() const { return normalize(_up.cross(forward())); }
            inline Vec3 rightward() const { return -leftward(); }

            void focusOn(const Sphere3 & target, bool updateMat = true);
            void translate(const Vec3 & t, const Sphere3 & target, bool updateMat = true);
            void moveEyeWithCenterFixed(const Vec3 & t, const Sphere3 & target, bool distanceFixed = false, bool updateMat = true);

        private:
            void updateMatrices();

        private:
            double _screenW, _screenH;
            double _focal;
            double _near, _far;
            Vec3 _eye, _center, _up;
            Mat4 _viewMatrix, _projectionMatrix, _viewProjectionMatrix/*, _viewProjectionMatrixInv*/;

            template <class Archive> inline void serialize(Archive & ar) {
                ar(_screenW, _screenH);
                ar(_focal, _near, _far);
                ar(_eye, _center, _up);
                ar(_viewMatrix, _projectionMatrix, _viewProjectionMatrix/*, _viewProjectionMatrixInv*/);
            }
            friend class cereal::access;
        };

        // panoramic camera
        class PanoramicCamera {
        public:
            explicit PanoramicCamera(double focal = 250,
                const Vec3 & eye = Vec3(0, 0, 0),
                const Vec3 & center = Vec3(1, 0, 0),
                const Vec3 & up = Vec3(0, 0, 1));

            inline SizeI screenSize() const { return SizeI(static_cast<int>(_focal * 2 * M_PI), static_cast<int>(_focal * M_PI)); }
            inline double focal() const { return _focal; }
            inline const Vec3 & eye() const { return _eye; }
            inline const Vec3 & center() const { return _center; }
            inline const Vec3 & up() const { return _up; }
            Vec2 screenProjection(const Vec3 & p3d) const;
            bool isVisibleOnScreen(const Vec3 & p3d) const { return true; }
            inline HPoint2 screenProjectionInHPoint(const Vec3 & p3d) const { return HPoint2(screenProjection(p3d), 1.0); }
            Vec3 spatialDirection(const Vec2 & p2d) const;
            inline Vec3 spatialDirection(const PixelLoc & p) const { return spatialDirection(Vec2(p.x, p.y)); }
            inline Vec3 spatialDirection(const HPoint2 & p) const { return spatialDirection(p.value()); }

        private:
            double _focal;
            Vec3 _eye, _center, _up;
            Vec3 _xaxis, _yaxis, _zaxis;

            template <class Archive> inline void serialize(Archive & ar) {
                ar(_focal);
                ar(_eye, _center, _up);
                ar(_xaxis, _yaxis, _zaxis);
            }
            friend class cereal::access;
        };


        // partial panoramic camera
        class PartialPanoramicCamera {
        public:
            explicit PartialPanoramicCamera(int w = 500, int h = 500, double focal = 250,
                const Vec3 & eye = Vec3(0, 0, 0),
                const Vec3 & center = Vec3(1, 0, 0),
                const Vec3 & up = Vec3(0, 0, -1));
            explicit PartialPanoramicCamera(const PanoramicCamera & panoCam, int w = 500, int h = 500);

            inline SizeI screenSize() const { return Size(static_cast<int>(_screenW), static_cast<int>(_screenH)); }
            inline double focal() const { return _focal; }
            inline const Vec3 & eye() const { return _eye; }
            inline const Vec3 & center() const { return _center; }
            inline const Vec3 & up() const { return _up; }
            Vec2 screenProjection(const Vec3 & p3d) const;
            bool isVisibleOnScreen(const Vec3 & p3d) const;
            inline HPoint2 screenProjectionInHPoint(const Vec3 & p3d) const { return HPoint2(screenProjection(p3d), 1.0); }
            Vec3 spatialDirection(const Vec2 & p2d) const;
            inline Vec3 spatialDirection(const PixelLoc & p) const { return spatialDirection(Vec2(p.x, p.y)); }

            PanoramicCamera toPanoramic() const { return PanoramicCamera(_focal, _eye, _center, _up); }

        private:
            double _screenW, _screenH;
            double _focal;
            Vec3 _eye, _center, _up;
            Vec3 _xaxis, _yaxis, _zaxis;

            template <class Archive> inline void serialize(Archive & ar) {
                ar(_screenW, _screenH);
                ar(_focal);
                ar(_eye, _center, _up);
                ar(_xaxis, _yaxis, _zaxis);
            }
            friend class cereal::access;
        };


        // sample image from image using camera conversion
        template <class OutCameraT, class InCameraT>
        class CameraSampler {
        public:
            CameraSampler(const OutCameraT & outCam, const InCameraT & inCam)
                : _outCam(outCam), _inCam(inCam) {
                assert(outCam.eye() == inCam.eye());
                auto outCamSize = _outCam.screenSize();
                _mapx = cv::Mat::zeros(outCamSize, CV_32FC1);
                _mapy = cv::Mat::zeros(outCamSize, CV_32FC1);
                for (int j = 0; j < outCamSize.height; j++) {
                    for (int i = 0; i < outCamSize.width; i++) {
                        Vec2 screenp(i, j);
                        Vec3 p3 = _outCam.spatialDirection(screenp);
                        if (!_inCam.isVisibleOnScreen(p3)){
                            _mapx.at<float>(j, i) = -1;
                            _mapy.at<float>(j, i) = -1;
                            continue;
                        }
                        Vec2 screenpOnInCam = _inCam.screenProjection(p3);
                        _mapx.at<float>(j, i) = static_cast<float>(screenpOnInCam(0));
                        _mapy.at<float>(j, i) = static_cast<float>(screenpOnInCam(1));
                    }
                }
            }

            Image operator() (const Image & inputIm,
                int borderMode = cv::BORDER_REPLICATE,
                const cv::Scalar & borderValue = cv::Scalar(0, 0, 0, 0)) const {
                Image outputIm;
                cv::remap(inputIm, outputIm, _mapx, _mapy,
                    inputIm.channels() <= 4 ? cv::INTER_LINEAR : cv::INTER_NEAREST, borderMode, borderValue);
                return outputIm;
            }

            template <class T>
            ImageOfType<T> operator() (const ImageOfType<T> & inputIm,
                int borderMode = cv::BORDER_REPLICATE,
                const T & borderValue = 0) const {
                ImageOfType<T> outputIm;
                cv::remap(inputIm, outputIm, _mapx, _mapy,
                    cv::INTER_NEAREST, borderMode, borderValue);
                return outputIm;
            }

        private:
            OutCameraT _outCam;
            InCameraT _inCam;
            cv::Mat _mapx, _mapy;

            template <class Archive>
            inline void serialize(Archive & ar) {
                ar(_outCam, _inCam, _mapx, _mapy);
            }
            friend class cereal::access;
        };


        template <class OutCameraT, class InCameraT>
        inline CameraSampler<OutCameraT, InCameraT> MakeCameraSampler(const OutCameraT & outCam, const InCameraT & inCam){
            return CameraSampler<OutCameraT, InCameraT>(outCam, inCam);
        }




      
        namespace {
            template <class T>
            struct IsCameraImpl {
                template <class TT>
                static auto test(int) -> decltype(
                    std::declval<TT>().eye(),
                    std::declval<TT>().screenProjection(std::declval<core::Vec3>()),
                    std::declval<TT>().screenProjectionInHPoint(std::declval<core::Vec3>()),
                    std::declval<TT>().isVisibleOnScreen(std::declval<core::Vec3>()),
                    std::declval<TT>().spatialDirection(std::declval<core::Vec2>()),
                    std::declval<TT>().spatialDirection(std::declval<core::PixelLoc>()),
                    std::declval<TT>().screenSize(),
                    std::true_type()
                    );
                template <class>
                static std::false_type test(...);
                static const bool value = std::is_same<decltype(test<T>(0)), std::true_type>::value;
            };
        }

        // judge whether T is a camera type
        template <class T>
        struct IsCamera : std::integral_constant<bool, IsCameraImpl<T>::value> {};



        // create horizontal cameras
        std::vector<PerspectiveCamera> CreateHorizontalPerspectiveCameras(const PanoramicCamera & panoCam,
            int num = 16, int width = 500, int height = 500, double focal = 250.0);

        std::vector<PerspectiveCamera> CreateHorizontalPerspectiveCameras(const PanoramicCamera & panoCam, 
            const std::vector<Vec3> & dirs,
            int width = 500, int height = 500, double focal = 250.0, double angleThreshold = 0.1);

        // create cameras toward cubic faces
        std::vector<PerspectiveCamera> CreateCubicFacedCameras(const PanoramicCamera & panoCam,
            int width = 500, int height = 500, double focal = 250.0);



        // view class
        template <class CameraT, class ImageT = Image, class = std::enable_if_t<IsCamera<CameraT>::value>>
        struct View {
            ImageT image;
            CameraT camera;

            template <class AnotherCameraT, class = std::enable_if_t<IsCamera<std::decay_t<AnotherCameraT>>::value>>
            inline View<std::decay_t<AnotherCameraT>, ImageT> sampled(AnotherCameraT && cam) const {
                View<std::decay_t<AnotherCameraT>, ImageT> v;
                v.image = MakeCameraSampler(cam, camera)(image);
                v.camera = std::forward<AnotherCameraT>(cam);
                return v;
            }

            template <class Archiver>
            void serialize(Archiver & ar) {
                ar(image, camera);
            }
        };

        using PerspectiveView = View<PerspectiveCamera>;
        using PanoramicView = View<PanoramicCamera>;
        using PartialPanoramicView = View<PartialPanoramicCamera>;




        // create panoramic view
        PanoramicView CreatePanoramicView(const Image & panorama,
            const Point3 & eye = Point3(0, 0, 0),
            const Point3 & center = Point3(1, 0, 0),
            const Vec3 & up = Vec3(0, 0, 1));

        // create perspective view
        PerspectiveView CreatePerspectiveView(const Image & perspectiveImage,
            const Point3 & eye = Point3(0, 0, 0),
            const Point3 & center = Point3(1, 0, 0),
            const Vec3 & up = Vec3(0, 0, -1),
            const LineSegmentExtractor & lse = LineSegmentExtractor(),
            const VanishingPointsDetector & vpd = VanishingPointsDetector(),
            std::vector<Classified<Line3>> * line3s = nullptr,
            std::vector<Classified<Line2>> * line2s = nullptr,
            std::vector<Vec3> * vps = nullptr,
            double * focal = nullptr);


        template <class CameraT, class = std::enable_if_t<IsCamera<CameraT>::value>>
        View<CameraT> CreateView(const Image & image,
            const Point3 & eye = Point3(0, 0, 0),
            const Point3 & center = Point3(1, 0, 0),
            const Vec3 & up = Vec3(0, 0, -1));


        template <>
        inline View<PanoramicCamera> CreateView(const Image & image,
            const Point3 & eye, const Point3 & center, const Vec3 & up){
            return CreatePanoramicView(image, eye, center, up);
        }

        template <>
        inline View<PerspectiveCamera> CreateView(const Image & image,
            const Point3 & eye, const Point3 & center, const Vec3 & up){
            return CreatePerspectiveView(image, eye, center, up);
        }

    }
}
 
#endif