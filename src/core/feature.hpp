#ifndef PANORAMIX_CORE_FEATURE_HPP
#define PANORAMIX_CORE_FEATURE_HPP

#include <opencv2/features2d/features2d.hpp>
#include <opencv2/stitching/detail/matchers.hpp>
#include <opencv2/nonfree/features2d.hpp>

#include "basic_types.hpp"
 
namespace panoramix {
    namespace core { 


        class PerspectiveCamera;
        class PanoramicCamera;


        // non maxima suppression
        void NonMaximaSuppression(const Image & src, Image & dst, int sz = 50,
            std::vector<Pixel> * pixels = nullptr,
            const Imageb & mask = Imageb());

        std::vector<KeyPoint> SIFT(const Image & im, cv::OutputArray descriptors = cv::noArray(),
            int nfeatures = 0, int nOctaveLayers = 3,
            double contrastThreshold = 0.04, double edgeThreshold = 10,
            double sigma = 1.6);
        std::vector<KeyPoint> SURF(const Image & im, cv::OutputArray descriptors,
            double hessianThreshold,
            int nOctaves = 4, int nOctaveLayers = 2,
            bool extended = true, bool upright = false);
        std::vector<KeyPoint> SURF(const Image & im, cv::OutputArray descriptors = cv::noArray());

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
            Feature operator() (const Image & im, int pyramidHeight, int minSize = 100) const;
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
            double angleThreshold, double sigma, double scoreThreshold = 0.8,
            double avoidVPDistanceThreshold = -1.0);

        // classify lines in 3d
        void ClassifyLines(std::vector<Classified<Line3>> &lines, const std::vector<Vec3> & vps,
            double angleThreshold, double sigma, double scoreThreshold = 0.8, 
            double avoidVPAngleThreshold = M_PI / 18.0);



        // compute straightness of points
        std::pair<double, Ray2> ComputeStraightness(const std::vector<std::vector<Pixel>> & edges,
            double * interArea = nullptr, double * interLen = nullptr);



        // find 3 orthogonal directions
        Failable<std::vector<Vec3>> FindOrthogonalPrinicipleDirections(const std::vector<Vec3> & directions,
            int longitudeDivideNum = 1000, int latitudeDivideNum = 500, 
            bool allowMoreThan2HorizontalVPs = false, const Vec3 & verticalSeed = Vec3(0, 0, 1));

        int NearestDirectionId(const std::vector<Vec3> & directions,
            const Vec3 & verticalSeed = Vec3(0, 0, 1));


        // estimate vanishing points
        std::vector<Vec3> EstimateVanishingPointsAndClassifyLines(const PerspectiveCamera & cams,
            std::vector<Classified<Line2>> & lineSegments);
        std::vector<Vec3> EstimateVanishingPointsAndClassifyLines(const std::vector<PerspectiveCamera> & cams,
            std::vector<std::vector<Classified<Line2>>> & lineSegments);



        // compute pp and focal from 3 orthogonal vps
        std::pair<Point2, double> ComputePrinciplePointAndFocalLength(const Point2 & vp1, const Point2 & vp2, const Point2 & vp3);


        // 2d vanishing point detection
        class VanishingPointsDetector {
        public:
            enum Algorithm {
                Naive,
                TardifSimplified,
                MATLAB_PanoContext,
                MATLAB_Tardif
            };
            struct Params {
                inline Params(Algorithm algo = Naive, double maxPPOffsetRatio = 2.0, 
                    double minFocalRatio = 0.05, double maxFocalRatio = 20.0)
                    : maxPrinciplePointOffsetRatio(maxPPOffsetRatio),
                    minFocalLengthRatio(minFocalRatio), 
                    maxFocalLengthRatio(maxFocalRatio), 
                    algorithm(algo) {
                }

                double maxPrinciplePointOffsetRatio;
                double minFocalLengthRatio, maxFocalLengthRatio;
                Algorithm algorithm;
                template <class Archive> inline void serialize(Archive & ar) { 
                    ar(maxPrinciplePointOffsetRatio, minFocalLengthRatio, maxFocalLengthRatio, algorithm);
                }
            };

        public:
            inline explicit VanishingPointsDetector(const Params & params = Params()) : _params(params) {}
            const Params & params() const { return _params; }
            Params & params() { return _params; }
            
            // accepts (lines, projection center)
            // returns (>= 3 vanishing points (the first 3 vanishing points should be the Manhattan VPs), the focal length, line classes)
            Failable<std::tuple<std::vector<HPoint2>, double, std::vector<int>>> operator() (
                const std::vector<Line2> & lines, const Sizei & imSize) const;  
            Failable<std::tuple<std::vector<HPoint2>, double>> operator() (
                std::vector<Classified<Line2>> & lines, const Sizei & imSize) const;

            template <class Archive> inline void serialize(Archive & ar) { ar(_params); }
        
        private:
            Params _params;
        };





        /// homography estimation
        std::pair<Failable<double>, Failable<double>> ComputeFocalsFromHomography(const Mat3 & H);





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
                    superpixelSizeSuggestion(1000), superpixelNumberSuggestion(100), useYUVColorSpace(false) {
                }
                float sigma; // for smoothing
                float c; // threshold function
                int minSize; // min component size
                Algorithm algorithm;
                int superpixelSizeSuggestion; // use superpixel size suggestion if [superpixelSizeSuggestion > 0]
                int superpixelNumberSuggestion; // use superpixel number suggestion if [superpixelSizeSuggestion < 0]
                bool useYUVColorSpace;
                template <class Archive> inline void serialize(Archive & ar) {
                    ar(sigma, c, minSize, algorithm, superpixelSizeSuggestion, superpixelNumberSuggestion, useYUVColorSpace);
                }
            };
        public:
            inline explicit SegmentationExtractor(const Params & params = Params()) : _params(params){}
            const Params & params() const { return _params; }
            Params & params() { return _params; }
            std::pair<Feature, int> operator() (const Image & im, bool isPanorama = false) const;
            std::pair<Feature, int> operator() (const Image & im, const std::vector<Line2> & lines,
                double extensionLength = 0.0) const;
            std::pair<Feature, int> operator() (const Image & im, const std::vector<Line3> & lines,
                const PanoramicCamera & cam, double extensionAngle = 0.0) const;
            template <class Archive> inline void serialize(Archive & ar) { ar(_params); }
        private:
            Params _params;
        };


        bool IsDenseSegmentation(const Imagei & segRegions);
        std::map<std::pair<int, int>, std::vector<std::vector<Pixel>>> FindContoursOfRegionsAndBoundaries(
            const Imagei & segRegions, int connectionExtendSize, bool simplifyStraightEdgePixels = true);
        std::map<std::set<int>, std::vector<Pixel>> ExtractBoundaryJunctions(const Imagei & regions);



        Imagei ComputeOrientationMaps(const std::vector<Classified<Line2>> & lines,
            const std::vector<HPoint2> & vps, const Sizei & imSize);



        /// geometric context estimator
        Image7d ComputeRawGeometricContext(const Image & im, bool outdoor, bool useHedauForIndoor);

        enum class GeometricContextIndex : size_t {
            FloorOrGround = 0,
            CeilingOrSky = 1,
            Vertical = 2,
            ClutterOrPorous = 3,
            Other = 4
        };

        Image5d MergeGeometricContextLabelsHoiem(const Image7d & rawgc);
        Image5d MergeGeometricContextLabelsHedau(const Image7d & rawgc);

        Image5d ComputeGeometricContext(const Image & im, bool outdoor, bool useHedauForIndoor = false);


    }
}
 
#endif