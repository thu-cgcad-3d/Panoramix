#ifndef PANORAMIX_CORE_FEATURE_HPP
#define PANORAMIX_CORE_FEATURE_HPP

#include <opencv2/features2d/features2d.hpp>
#include <opencv2/stitching/detail/matchers.hpp>
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

            template <class Archive> inline void serialize(Archive & ar) { ar(_outCam, _inCam); }

        private:
            OutCameraT _outCam;
            InCameraT _inCam;
            cv::Mat _mapx, _mapy;
        };






        

        /// feature extractors

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
                template <class Archive> inline void serialize(Archive & ar) { 
                    ar(minLength, xBorderWidth, yBorderWidth, numDirs, useExperimentalAlgorithm); 
                }
            };
        public:
            inline explicit LineSegmentExtractor(const Params & params = Params()) : _params(params){}
            Feature operator() (const Image & im) const;
            template <class Archive> inline void serialize(Archive & ar) { ar(_params); }
        private:
            Params _params;
        };



        // compute line intersections
        std::vector<HPoint2> ComputeLineIntersections(const std::vector<Line2> & lines,
            std::vector<std::pair<int, int>> * lineids = nullptr,
            bool suppresscross = true);


        // compute folding difficulty along vanishing point
        double ComputeFoldingDifficultyAlongVanishingPoint(const std::vector<Point2> & points, const HPoint2 & vp);



        // point feature extractor
        template <class CVFeatureExtractorT, class FeatureT = std::vector<cv::KeyPoint>>
        class CVFeatureExtractor {
        public:
            using Feature = FeatureT;
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
            CVFeatureExtractorT _feature2D;
        };




        // segmentation
        class SegmentationExtractor {
        public:
            using Feature = Image; // CV_32SC1, from 0 to numRegion, use at<int32_t> to extract
            struct Params {
                inline Params() : sigma(0.8f), c(100.0f), minSize(200) {}
                float sigma; // for smoothing
                float c; // threshold function
                int minSize; // min component size
                template <class Archive> inline void serialize(Archive & ar) { ar(sigma, c, minSize); }
            };
        public:
            inline explicit SegmentationExtractor(const Params & params = Params()) : _params(params){}
            Feature operator() (const Image & im, bool forVisualization = false) const;
            template <class Archive> inline void serialize(Archive & ar) { ar(_params); }
        private:
            Params _params;
        };




        // vanishing point extractor
        class VanishingPointsExtractor {
        public:
            using Feature = std::array < HPoint2, 3 > ;
            struct Params {
                inline Params() {}

                template <class Archive> inline void serialize(Archive & ar) {}
            };
        public:
            inline VanishingPointsExtractor(const Params & params = Params()) : _params(params) {}
            // find vanishing points with known focal length
            Feature operator() (const LineSegmentExtractor::Feature & lines, const Point2 & projCenter, double focalLength) const;
            // find vanishing points without known focal length
            std::pair<Feature, double> operator() (const LineSegmentExtractor::Feature & lines, const Point2 & projCenter) const;

            template <class Archive> inline void serialize(Archive & ar) { ar(_params); }
        private:
            Params _params;
        };





        // geometric context
        class GeometricContextExtractor {
        public:
            using Feature = ImageWithType < int > ;
            struct Params {
                inline Params() {}

                template <class Archive> inline void serialize(Archive & ar) {}
            };
        public:
            inline GeometricContextExtractor(const Params & params = Params()) : _params(params) {}
            Feature operator() (const Image & image) const;
            template <class Archive> inline void serialize(Archive & ar) { ar(_params); }
        private:
            Params _params;
        };



        /// homography estimation
        std::pair<double, double> ComputeFocalsFromHomography(const Mat3 & H, std::pair<bool, bool> * ok = nullptr);









    }
}
 
#endif