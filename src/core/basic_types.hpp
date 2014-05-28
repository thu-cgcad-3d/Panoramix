#ifndef PANORAMIX_CORE_BASIC_TYPES_HPP
#define PANORAMIX_CORE_BASIC_TYPES_HPP

#include <vector>
#include <list>
#include <deque>
#include <set>
#include <unordered_set>
#include <forward_list>
#include <array>
#include <memory>
#include <cstdint>
#include <complex>
#include <functional>

#include <opencv2/opencv.hpp>
#include <opencv2/stitching/detail/matchers.hpp>

#include "version.hpp"

namespace panoramix {
    namespace core {
        
        // vectors/points
        template <class T, int N> using Vec = cv::Vec<T, N>;
        using Vec2 = Vec<double, 2>;
        using Vec3 = Vec<double, 3>;
        using Vec4 = Vec<double, 4>;
        template <class T, int N> using Point = cv::Vec<T, N>;
        using Point2 = Point<double, 2>;
        using Point3 = Point<double, 3>;
        using Point4 = Point<double, 4>;
        template <class T, int M, int N> using Mat = cv::Matx<T, M, N>;
        using Mat3 = Mat<double, 3, 3>;
        using Mat4 = Mat<double, 4, 4>;

        using cv::norm;


        // homogeneous point
        template <class T, int N>
        struct HPoint {
            inline HPoint() : coord(), scalar(1){}
            inline HPoint(const Point<T, N> & c, T s = 1) : coord(c), scalar(s) {}
            inline Point<T, N> toPoint() const { return coord / scalar; }
            Vec<T, N + 1> toVector() const { 
                Vec<T, N + 1> v; 
                std::copy(coord.val, coord.val + N, v.val);
                v[N] = scalar;
                return v;
            }
            Point<T, N> coord;
            T scalar;
        };
        template <class T, int N>
        HPoint<T, N-1> HPointFromVector(const Vec<T, N> & v){
            HPoint<T, N-1> hp;
            std::copy(v.val, v.val + N - 1, hp.coord.val);
            hp.scalar = v[N-1];
            return hp;
        }
        using HPoint2 = HPoint<double, 2>;
        using HPoint3 = HPoint<double, 3>;
        using HPoint4 = HPoint<double, 4>;


        // geographic coordinate
        struct GeoCoord {
            inline explicit GeoCoord(double longi = 0.0, double lati = 0.0)
                : longitude(longi), latitude(lati) {}
            template <class T>
            inline GeoCoord(const Vec<T, 3> & d)
                : longitude(atan2(d(1), d(0))),
                latitude(atan(d(2) / sqrt((d(1)*d(1)) + (d(0)*d(0)))))
            {}
            template <class T = double>
            inline Vec<T, 3> toVector() const {
                return Vec<T, 3>(static_cast<T>(cos(longitude)*cos(latitude)),
                    static_cast<T>(sin(longitude)*cos(latitude)),
                    static_cast<T>(sin(latitude)));
            }
            double longitude; // - M_PI ~ + M_PI
            double latitude; // - M_PI_2 ~ + M_PI_2
        };



        // key point
        using KeyPoint = cv::KeyPoint;
        

        // size
        using Size = cv::Size2f;
        using SizeI = cv::Size2i;


        // line
        template <class T, int N> 
        struct Line {
            inline Line(){}
            inline Line(const Point<T, N> & f, const Point<T, N> & s) 
                : first(f), second(s) {}
            inline Point<T, N> center() const { return (first + second) / 2.0; }
            inline T length() const { return norm(first - second); }
            inline Vec<T, N> direction() const { return second - first; }
            inline Line reversed() const { return Line(second, first); }
            Point<T, N> first, second;
        };
        using Line2 = Line<double, 2>;
        using Line3 = Line<double, 3>;

        
        // position on line
        template <class T, int N>
        struct PositionOnLine {
            inline PositionOnLine(){}
            inline PositionOnLine(const Line<T, N> & line, const T & r)
                : ratio(r), position(line.first + (line.second - line.first) * ratio) {}
            T ratio; // [0 ~ 1]: on line
            Point<T, N> position; // position = line.first + (line.second - line.fist) * ratio
        };


        // homogeneous line
        template <class T, int N>
        struct HLine {
            HPoint<T, N> first, second;
            inline Line<T, N> toLine() const {
                return Line<T, N>{first.toPoint(), second.toPoint()};
            }
        };
        using HLine2 = HLine<double, 2>;
        using HLine3 = HLine<double, 3>;


        // image
        using Image = cv::Mat;
        using PixelLoc = cv::Point;        
        std::pair<PixelLoc, PixelLoc> MinMaxLocOfImage(const Image & im);
        std::pair<double, double> MinMaxValOfImage(const Image & im);


        // color
        using Color = cv::Scalar;
        enum class ColorTag {
            Transparent,

            White,
            Black,

            DimGray,
            Gray,
            DarkGray,
            Silver,
            LightGray,

            Red,
            Green,
            Blue,

            Yellow,
            Magenta,
            Cyan,
            Orange
        };
        Color ColorFromTag(ColorTag t);
        enum class ColorTableDescriptor : int8_t {
            RGB,
            WRGB,
            AllColors
        };
        const std::vector<Color> & 
            PredefinedColorTable(ColorTableDescriptor descriptor = ColorTableDescriptor::AllColors);


        // somthing classified
        template <class T>
        struct Classified {
            int claz;
            T component;
        }; 


        // somthing transformed in N-d space
        template <class T, int N>
        struct Transformed {
            Mat<double, N+1, N+1> transform;
            T component;
        };

 
    }
}
 
#endif