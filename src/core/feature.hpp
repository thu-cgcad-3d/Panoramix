#ifndef PANORAMIX_CORE_FEATURE_HPP
#define PANORAMIX_CORE_FEATURE_HPP

#include <opencv2/features2d/features2d.hpp>
#include <opencv2/stitching/detail/matchers.hpp>
#include <opencv2/nonfree/features2d.hpp>

#include "basic_types.hpp"
 
namespace panoramix {
    namespace core { 


        // non maxima suppression
        void NonMaximaSuppression(const Image & src, Image & dst, int sz = 50,
            std::vector<PixelLoc> * pixels = nullptr,
            const Imageb & mask = Imageb());




        // line extractor
        class LineSegmentExtractor {
        public:
            using Feature = std::vector<Line2>;
            enum Algorithm {
                GradientGrouping,
                LSD
            };
            struct Params {
                inline Params() : minLength(15), xBorderWidth(1), yBorderWidth(1), numDirs(8), algorithm(LSD) {}
                int minLength;
                int xBorderWidth, yBorderWidth;
                int numDirs;
                Algorithm algorithm;
                template <class Archive> 
                inline void serialize(Archive & ar) { 
                    ar(minLength, xBorderWidth, yBorderWidth, numDirs, algorithm);
                }
            };
        public:
            inline explicit LineSegmentExtractor(const Params & params = Params()) : _params(params){}
            const Params & params() const { return _params; }
            Params & params() { return _params; }
            Feature operator() (const Image & im) const;
            template <class Archive> inline void serialize(Archive & ar) { ar(_params); }
        private:
            Params _params;
        };






        // compute line intersections
        std::vector<HPoint2> ComputeLineIntersections(const std::vector<Line2> & lines,
            std::vector<std::pair<int, int>> * lineids = nullptr,
            bool suppresscross = true, 
            double minDistanceBetweenLinePairs = std::numeric_limits<double>::max());


        // classify lines in 2d
        void ClassifyLines(std::vector<Classified<Line2>> &lines, const std::vector<HPoint2> & vps,
            double angleThreshold, double sigma, double scoreThreshold = 0.8);

        // classify lines in 3d
        void ClassifyLines(std::vector<Classified<Line3>> &lines, const std::vector<Vec3> & vps,
            double angleThreshold, double sigma, double scoreThreshold = 0.8);

        // compute straightness of points
        std::pair<double, InfiniteLine2> ComputeStraightness(const std::vector<std::vector<PixelLoc>> & edges,
            double * interArea = nullptr, double * interLen = nullptr);



        // find 3 orthogonal directions
        std::vector<Vec3> FindThreeOrthogonalPrinicipleDirections(const std::vector<Vec3> & directions,
            int longitudeDivideNum = 1000, int latitudeDivideNum = 500);


        // 2d vanishing point detection
        class VanishingPointsDetector {
        public:
            enum Algorithm {
                Hedau,
                MSAC
            };
            struct Params {
                inline Params(double maxPPOffset = 80, double minFocal = 40, double maxFocal = 1e5)
                    : maxPrinciplePointOffset(maxPPOffset), minFocalLength(minFocal), maxFocalLength(maxFocal), algorithm(Hedau) {}
                double maxPrinciplePointOffset;
                double minFocalLength, maxFocalLength;
                Algorithm algorithm;
                template <class Archive> inline void serialize(Archive & ar) { 
                    ar(maxPrinciplePointOffset, minFocalLength, maxFocalLength, algorithm);
                }
            };

        public:
            inline explicit VanishingPointsDetector(const Params & params = Params()) : _params(params) {}
            const Params & params() const { return _params; }
            Params & params() { return _params; }
            
            // accepts (lines, projection center)
            // returns (3 vanishing points, the focal length, line classes)
            std::tuple<std::array<HPoint2, 3>, double, std::vector<int>>
                operator() (const std::vector<Line2> & lines, const Point2 & projCenter) const;            
            template <class Archive> inline void serialize(Archive & ar) { ar(_params); }
        
        private:
            std::tuple<std::array<HPoint2, 3>, double, std::vector<int>>
                estimateWithProjectionCenterAtOrigin (const std::vector<Line2> & lines) const;
            Params _params;
        };









        // vanishing point detection for local manhattan scenes
        class LocalManhattanVanishingPointsDetector {
        public:
            struct Params {
				inline Params(double maxPPOffset = 200, double minFocal = 40, double maxFocal = 1e5) 
					: verticalVPAngleRange(M_PI_4 / 6.0), 
					verticalVPMinDistanceRatioToCenter(2.0), 
					maxPrinciplePointOffset(maxPPOffset), minFocalLength(minFocal), maxFocalLength(maxFocal) {
				}
				// the angle between {the line connecting verticalVP and image center} and {the vertical line} 
				// should be within [-verticalVPAngleRange, +verticalVPAngleRange] 
				double verticalVPAngleRange; 
				double verticalVPMinDistanceRatioToCenter;
                double maxPrinciplePointOffset;
                double minFocalLength, maxFocalLength;
                template <class Archive> inline void serialize(Archive & ar) {
                    ar(verticalVPAngleRange, verticalVPMinDistanceRatioToCenter, 
						maxPrinciplePointOffset, minFocalLength, maxFocalLength);
                }

                Image image;
            };
            struct Result {
                std::vector<HPoint2> vanishingPoints;
                int verticalVanishingPointId;
                std::vector<std::pair<int, int>> horizontalVanishingPointIds;
                double focalLength;
                InfiniteLine2 horizon;
                Point2 principlePoint;
                std::vector<int> lineClasses;
                template <class Archive> inline void serialize(Archive & ar) {
                    ar(vanishingPoints, verticalVanishingPointId,
                        horizontalVanishingPointIds, focalLength, horizon, principlePoint,
                        lineClasses);
                }

                std::vector<InfiniteLine2> hlineCands;
            };
        public:
            inline explicit LocalManhattanVanishingPointsDetector(const Params & params = Params()) : _params(params) {}
            const Params & params() const { return _params; }
            Params & params() { return _params; }
            Result operator() (const std::vector<Line2> & lines, const Point2 & projCenter) const;
            template <class Archive> inline void serialize(Archive & ar) { ar(_params); }
        private:
            //Result estimateWithProjectionCenterAtOrigin(const std::vector<Line2> & lines) const;
            //Result estimateWithProjectionCenterAtOriginII(const std::vector<Line2> & lines) const;
            Result estimateWithProjectionCenterAtOriginIII(const std::vector<Line2> & lines) const;

            Params _params;
        };




        /// homography estimation
        std::pair<double, double> ComputeFocalsFromHomography(const Mat3 & H, std::pair<bool, bool> * ok = nullptr);





        // segmentation
        class SegmentationExtractor {
        public:
            using Feature = Imagei; // CV_32SC1
            enum Algorithm {
                GraphCut,
                SLIC,
                QuickShiftCPU,
                QuickShiftGPU
            };
            struct Params {
                inline Params() : sigma(0.8f), c(100.0f), minSize(200), algorithm(GraphCut),
                    superpixelSizeSuggestion(1000), superpixelNumberSuggestion(100) {
                }
                float sigma; // for smoothing
                float c; // threshold function
                int minSize; // min component size
                Algorithm algorithm;
                int superpixelSizeSuggestion; // use superpixel size suggestion if [superpixelSizeSuggestion > 0]
                int superpixelNumberSuggestion; // use superpixel number suggestion if [superpixelSizeSuggestion < 0]
                template <class Archive> inline void serialize(Archive & ar) {
                    ar(sigma, c, minSize, algorithm, superpixelSizeSuggestion, superpixelNumberSuggestion);
                }
            };
        public:
            inline explicit SegmentationExtractor(const Params & params = Params()) : _params(params){}
            const Params & params() const { return _params; }
            Params & params() { return _params; }
            std::pair<Feature, int> operator() (const Image & im) const;
            std::pair<Feature, int> operator() (const Image & im, const std::vector<Line2> & lines) const;
            template <class Archive> inline void serialize(Archive & ar) { ar(_params); }
        private:
            Params _params;
        };




        /// geometric context estimator
        class GeometricContextEstimator {
        public:
            struct Params {
                Params();
                bool useMatlab;
                std::string matlabCodeFolder;
            };

        public:
            inline explicit GeometricContextEstimator(const Params & params = Params()) : _params(params) {}
            const Params & params() const { return _params; }
            Params & params() { return _params; }
            Image operator() (const Image & im) const; // CV_64FCX

        private:
            Params _params;
        };


    }
}
 
#endif