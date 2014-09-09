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
#include "macros.hpp"
#include "serialization.hpp"

namespace panoramix {
    namespace core {
        
        // vectors/points
        template <class T, int N> using Vec = cv::Vec<T, N>;
        using Vec2 = Vec<double, 2>;
        using Vec3 = Vec<double, 3>;
        using Vec4 = Vec<double, 4>;
        using Vec2i = Vec<int, 2>;
        using Vec3i = Vec<int, 3>;
        using Vec4i = Vec<int, 4>;
        using Vec2b = Vec<uint8_t, 2>;
        using Vec3b = Vec<uint8_t, 3>;
        using Vec4b = Vec<uint8_t, 4>;
        template <class T, int N> using Point = cv::Vec<T, N>;
        using Point2 = Point<double, 2>;
        using Point3 = Point<double, 3>;
        using Point4 = Point<double, 4>;
        using Point2i = Point<int, 2>;
        using Point3i = Point<int, 3>;
        using Point4i = Point<int, 4>;

        template <class T, int M, int N> using Mat = cv::Matx<T, M, N>;
        using Mat3 = Mat<double, 3, 3>;
        using Mat4 = Mat<double, 4, 4>;

        using cv::norm;
        template <class T>
        inline T normalize(const T & d) { return d / norm(d); }

        template <class To, class From, int N>
        inline Vec<To, N> ConvertTo(const Vec<To, N> & v) {
            Vec<To, N> out;
            for (int i = 0; i < N; i++) {
                out[i] = static_cast<To>(v[i]);
            }
            return out;
        }


        // ratio
        template <class T, class S>
        struct Ratio {
            inline Ratio() : denominator(1) {}
            inline Ratio(const T & c) : numerator(c), denominator(1) {}
            inline Ratio(const T & c, const S & s) : numerator(c), denominator(s) {}
            inline T value() const { return numerator / denominator; }
            T numerator;
            S denominator;
        };
        template <class Archive, class T, class S>
        inline void serialize(Archive & ar, Ratio<T, S> & r) {
            ar(r.numerator, r.denominator);
        }
        template <class T, class S>
        inline Ratio<T, S> MakeRatio(const T & c, const S & s) { return Ratio<T, S>(c, s); }

        namespace {
            enum RatioOperationType { Plus, Minus, Mult, Div };
            template <class T1, class S1, class T2, class S2, RatioOperationType>
            struct RatioOperationResult {};
            template <class T1, class S1, class T2, class S2>
            struct RatioOperationResult<T1, S1, T2, S2, Plus> {
                using type = Ratio<
                decltype(std::declval<T1>() * std::declval<S2>() + std::declval<T2>() * std::declval<S1>()),
                decltype(std::declval<S1>() * std::declval<S2>())
                >;
            };
            template <class T1, class S1, class T2, class S2>
            struct RatioOperationResult<T1, S1, T2, S2, Minus> {
                using type = Ratio<
                decltype(std::declval<T1>() * std::declval<S2>() - std::declval<T2>() * std::declval<S1>()),
                decltype(std::declval<S1>()  *std::declval<S2>())
                >;
            };
            template <class T1, class S1, class T2, class S2>
            struct RatioOperationResult<T1, S1, T2, S2, Mult> {
                using type = Ratio<
                decltype(std::declval<T1>() * std::declval<T2>()),
                decltype(std::declval<S1>() * std::declval<S2>())
                >;
            };
            template <class T1, class S1, class T2, class S2>
            struct RatioOperationResult<T1, S1, T2, S2, Div> {
                using type = Ratio<
                decltype(std::declval<T1>() * std::declval<S2>()),
                decltype(std::declval<S1>() * std::declval<T2>())
                >;
            };
        }

        template <class T1, class S1, class T2, class S2>
        inline typename RatioOperationResult<T1, S1, T2, S2, Plus>::type
            operator + (const Ratio<T1, S1> & a, const Ratio<T2, S2> & b) {
                return MakeRatio(a.numerator * b.denominator + b.numerator * a.denominator, a.denominator * b.denominator);
            }
        template <class T1, class S1, class T2, class S2>
        inline typename RatioOperationResult<T1, S1, T2, S2, Minus>::type
            operator - (const Ratio<T1, S1> & a, const Ratio<T2, S2> & b) {
                return MakeRatio(a.numerator * b.denominator - b.numerator * a.denominator, a.denominator * b.denominator);
            }
        template <class T1, class S1, class T2, class S2>
        inline typename RatioOperationResult<T1, S1, T2, S2, Mult>::type
            operator * (const Ratio<T1, S1> & a, const Ratio<T2, S2> & b) {
                return MakeRatio(a.numerator * b.numerator, a.denominator * b.denominator);
            }
        template <class T1, class S1, class T2, class S2>
        inline typename RatioOperationResult<T1, S1, T2, S2, Div>::type
            operator / (const Ratio<T1, S1> & a, const Ratio<T2, S2> & b) {
                return MakeRatio(a.numerator * b.denominator, a.denominator * b.numerator);
            }
        template <class T, class S>
        inline Ratio<T, S> operator - (const Ratio<T, S> & r) {
            return Ratio<T, S>(-r.numerator, r.denominator);
        }
        template <class T, class S>
        inline bool operator == (const Ratio<T, S> & a, const Ratio<T, S> & b) {
            return a.numerator == b.numerator && a.denominator == b.denominator;
        }

        using Rational = Ratio<double, double>;
        template <class T, class = std::enable_if_t<std::is_floating_point<T>::value>>
        inline Ratio<T, T> MakeRational(const T & v) {
            return std::isinf(v) ? (v > 0 ? Ratio<T, T>(1.0, 0.0) : Ratio<T, T>(-1.0, 0.0))
                : (std::isnan(v) ? Ratio<T, T>(0.0, 0.0) : Ratio<T, T>(v, 1.0));
        }


        // homogeneous point
        template <class T, int N>
        using HPoint = Ratio<Point<T, N>, T>;
        template <class T, int N>
        Vec<T, N + 1> VectorFromHPoint(const HPoint<T, N> & p, const T & scale = 1.0) {
            Vec<T, N + 1> v;
            std::copy(p.numerator.val, p.numerator.val + N, v.val);
            v[N] = p.denominator * scale;
            return v;
        }
        template <class T, int N>
        HPoint<T, N - 1> HPointFromVector(const Vec<T, N> & v) {
            HPoint<T, N - 1> hp;
            std::copy(v.val, v.val + N - 1, hp.numerator.val);
            hp.denominator = v[N - 1];
            return hp;
        }
        template <class T, int N>
        inline Ratio<T, T> norm(const HPoint<T, N> & p) {
            return Ratio<T, T>(norm(p.numerator), p.denominator);
        }
        template <class T, int N>
        inline Ratio<T, T> dot(const HPoint<T, N> & a, const HPoint<T, N> &b) {
            return Ratio<T, T>(a.numerator.dot(b.numerator), a.denominator * b.denominator);
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
        inline bool operator == (const GeoCoord & a, const GeoCoord & b) {
            return a.longitude == b.longitude && a.latitude == b.latitude;
        }
        template <class Archive>
        inline void serialize(Archive & ar, GeoCoord & gc) {
            ar(gc.longitude, gc.latitude);
        }


        // key point
        using KeyPoint = cv::KeyPoint;
        

        // size
        using Size = cv::Size2f;
        using SizeI = cv::Size2i;



        // infinite line
        template <class T, int N>
        struct InfiniteLine {
            inline InfiniteLine() {}
            inline InfiniteLine(const Point<T, N> & a, const Vec<T, N> & d) : anchor(a), direction(d) {}
            Point<T, N> anchor;
            Vec<T, N> direction;
        };
        template <class T, int N>
        inline bool operator == (const InfiniteLine<T, N> & a, const InfiniteLine<T, N> & b) {
            return a.anchor == b.anchor && a.direction == b.direction;
        }
        template <class Archive, class T, int N>
        inline void serialize(Archive & ar, InfiniteLine<T, N> & p) {
            ar(p.anchor, p.direction);
        }
        using InfiniteLine2 = InfiniteLine<double, 2>;
        using InfiniteLine3 = InfiniteLine<double, 3>;
        template <class T>
        inline Vec<T, 3> GetLine2Coeffs(const InfiniteLine<T, 2> & line) {
            return Vec<T, 3>{line.direction[1], 
                -line.direction[0], 
                -(line.direction[1] * line.anchor[0] - line.direction[0] * line.anchor[1])
            };
        }


        // plane
        template <class T, int N>
        struct Plane {
            Point<T, N> anchor;
            Vec<T, N> normal;
        };
        template <class T, int N>
        inline bool operator == (const Plane<T, N> & a, const Plane<T, N> & b) {
            return a.anchor == b.anchor && a.normal == b.normal;
        }
        template <class Archive, class T, int N>
        inline void serialize(Archive & ar, Plane<T, N> & p) {
            ar(p.anchor, p.normal);
        }
        using Plane3 = Plane<double, 3>;



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
            inline InfiniteLine<T, N> infinieLine() const { return InfiniteLine<T, N>{ first, second - first }; }
            Point<T, N> first, second;
        };
        template <class T, int N>
        inline bool operator == (const Line<T, N> & a, const Line<T, N> & b) {
            return a.first == b.first && a.second == b.second;
        }
        template <class Archive, class T, int N>
        inline void serialize(Archive & ar, Line<T, N> & l) {
            ar(l.first, l.second);
        }
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
        template <class T, int N>
        inline bool operator == (const PositionOnLine<T, N> & a, const PositionOnLine<T, N> & b) {
            return a.ratio == b.ratio && a.position == b.position;
        }
        template <class Archive, class T, int N>
        inline void serialize(Archive & ar, PositionOnLine<T, N> & p) {
            ar(p.ratio, p.position);
        }
        using PositionOnLine2 = PositionOnLine<double, 2>;
        using PositionOnLine3 = PositionOnLine<double, 3>;


        // homogeneous line
        template <class T, int N>
        struct HLine {
            HPoint<T, N> first, second;
            inline Line<T, N> toLine() const {
                return Line<T, N>{first.value(), second.value()};
            }
        };
        template <class T, int N>
        inline bool operator ==  (const HLine<T, N> & a, const HLine<T, N> & b) {
            return a.first == b.first && a.second == b.second;
        }
        template <class Archive, class T, int N>
        inline void serialize(Archive & ar, HLine<T, N> & l) {
            ar(l.first, l.second);
        }
        using HLine2 = HLine<double, 2>;
        using HLine3 = HLine<double, 3>;






        // image
        using Image = cv::Mat;
        template <class T>
        using ImageWithType = cv::Mat_<T> ;
        using PixelLoc = cv::Point;        
        std::pair<PixelLoc, PixelLoc> MinMaxLocOfImage(const Image & im);
        std::pair<double, double> MinMaxValOfImage(const Image & im);



        // somthing classified
        template <class T>
        struct Classified {
            int claz;
            T component;
        };
        template <class T>
        inline bool operator == (const Classified<T> & a, const Classified<T> & b) {
            return a.claz == b.claz && a.component == b.component;
        }
        template <class Archive, class T>
        inline void serialize(Archive & ar, Classified<T> & c) {
            ar(c.claz, c.component);
        }



        namespace {
            template <class T, int N>
            Vec<T, N> MakeMin(const Vec<T, N>& v1, const Vec<T, N>& v2) {
                Vec<T, N> v;
                for (int i = 0; i < N; i++)
                    v[i] = v1[i] < v2[i] ? v1[i] : v2[i];
                return v;
            }
            template <class T, int N>
            Vec<T, N> MakeMax(const Vec<T, N>& v1, const Vec<T, N>& v2) {
                Vec<T, N> v;
                for (int i = 0; i < N; i++)
                    v[i] = v1[i] < v2[i] ? v2[i] : v1[i];
                return v;
            }
        }

        // box
        template <class T, int N>
        struct Box {
            using Type = T;
            static const int Dimension = N;

            Point<T, N> minCorner, maxCorner;
            bool isNull;

            inline Box() : isNull(true) {}
            inline Box(const Point<T, N> & c1, const Point<T, N> & c2)
                : minCorner(MakeMin(c1, c2)), maxCorner(MakeMax(c1, c2)), isNull(false) {}

            inline Vec<T, N> size() const { return maxCorner - minCorner; }
            inline Point<T, N> center() const { return (maxCorner + minCorner) * (0.5); }

            inline bool contains(const Point<T, N> & p) const {
                if (isNull)
                    return false;
                for (int i = 0; i < N; i++) {
                    if (minCorner[i] > p[i] || maxCorner[i] < p[i])
                        return false;
                }
                return true;
            }
            inline bool contains(const Box & b) const {
                return b.isNull ? true : contains(b.minCorner) && contains(b.maxCorner);
            }
            inline bool operator != (const Box & b) const { return !(*this == b); }

            inline Box & operator |= (const Box & b) {
                if (isNull) {
                    *this = b;
                    return *this;
                }
                if (b.isNull)
                    return *this;
                minCorner = MakeMin(minCorner, b.minCorner);
                maxCorner = MakeMax(maxCorner, b.maxCorner);
                return *this;
            }
        };
        template <class T, int N>
        inline bool operator == (const Box<T, N> & a, const Box<T, N> & b) {
            return a.isNull ? b.isNull : (!b.isNull && a.minCorner == b.minCorner && a.maxCorner == b.maxCorner);
        }
        template <class Archive, class T, int N>
        inline void serialize(Archive & ar, Box<T, N> & b) {
            ar(b.isNull, b.minCorner, b.maxCorner);
        }
        template <class T, int N>
        inline Box<T, N> operator | (const Box<T, N> & b1, const Box<T, N> & b2) {
            Box<T, N> b12 = b1;
            return b12 |= b2;
        }
        using Box2 = Box<double, 2>;
        using Box3 = Box<double, 3>;

 
    }
}


 
#endif