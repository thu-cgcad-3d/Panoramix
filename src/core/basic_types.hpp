#ifndef PANORAMIX_CORE_BASIC_TYPES_HPP
#define PANORAMIX_CORE_BASIC_TYPES_HPP

#include <iterator>
#include <vector>
#include <list>
#include <deque>
#include <forward_list>
#include <array>
#include <set>
#include <map>
#include <unordered_set>
#include <unordered_map>
#include <memory>
#include <numeric>
#include <cstdint>
#include <complex>
#include <functional>

#include <opencv2/opencv.hpp>

#include "version.hpp"
#include "macros.hpp"
#include "serialization.hpp"
#include "ratio.hpp"
#include "any.hpp"
#include "failable.hpp"

namespace panoramix {
    namespace core {
        
        // vectors/points
        template <class T, int N> using Vec = cv::Vec<T, N>;
        using Vec2 = Vec<double, 2>;
        using Vec3 = Vec<double, 3>;
        using Vec4 = Vec<double, 4>;
        using Vec2f = Vec<float, 2>;
        using Vec3f = Vec<float, 3>;
        using Vec4f = Vec<float, 4>;
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
        using Point2f = Point<float, 2>;
        using Point3f = Point<float, 3>;
        using Point4f = Point<float, 4>;
        using Point2i = Point<int, 2>;
        using Point3i = Point<int, 3>;
        using Point4i = Point<int, 4>;

        template <class T> struct IsVecOrPoint : std::false_type {};
        template <class T, int N> struct IsVecOrPoint<Point<T, N>> : std::true_type{};

        // matrix
        template <class T, int M, int N> using Mat = cv::Matx<T, M, N>;
        using Mat3 = Mat<double, 3, 3>;
        using Mat4 = Mat<double, 4, 4>;
        using Mat3f = Mat<float, 3, 3>;
        using Mat4f = Mat<float, 4, 4>;

        using cv::norm;
        template <class T>
        inline T normalize(const T & d) { return d / norm(d); }

        namespace {
            template <class To, class From, int N>
            inline Vec<To, N> VecCastPrivate(const Vec<From, N> & v, std::false_type) {
                Vec<To, N> out;
                for (int i = 0; i < N; i++) {
                    out[i] = static_cast<To>(v[i]);
                }
                return out;
            }
            template <class To, int N>
            inline const Vec<To, N> & VecCastPrivate(const Vec<To, N> & v, std::true_type) {
                return v;
            }
        }

        template <class To, class From, int N>
        inline Vec<To, N> vec_cast(const Vec<From, N> & v) {
            return VecCastPrivate<To>(v, std::integral_constant<bool, std::is_same<To, From>::value>());
        }

        template <class T, int M, int N>
        inline Vec<T, M + N> Concat(const Vec<T, M> & a, const Vec<T, N> & b) {
            Vec<T, M + N> ab;
            std::copy(a.val, a.val + M, ab.val);
            std::copy(b.val, b.val + N, ab.val + M);
            return ab;
        }

        template <class T, int M>
        inline Vec<T, M + 1> Concat(const Vec<T, M> & a, const T & b) {
            Vec<T, M + 1> ab;
            std::copy(a.val, a.val + M, ab.val);
            ab[M] = b;
            return ab;
        }

        template <class T, int M>
        inline Vec<T, M + 1> Concat(const T & a, const Vec<T, M> & b) {
            Vec<T, M + 1> ab;
            ab.val[0] = a;
            std::copy(b.val, b.val + M, ab.val + 1);
            return ab;
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
        HPoint<T, N - 1> HPointFromVector(const Vec<T, N> & v, const T & scale = 1.0) {
            HPoint<T, N - 1> hp;
            std::copy(v.val, v.val + N - 1, hp.numerator.val);
            hp.denominator = v[N - 1] / scale;
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
                : longitude(std::atan2(d(1), d(0))),
                latitude(std::atan(d(2) / std::sqrt((d(1)*d(1)) + (d(0)*d(0)))))
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
        struct Ray {
            inline Ray() {}
            inline Ray(const Point<T, N> & a, const Vec<T, N> & d) : anchor(a), direction(d) {}
            Point<T, N> anchor;
            Vec<T, N> direction;
        };
        template <class T, int N>
        inline bool operator == (const Ray<T, N> & a, const Ray<T, N> & b) {
            return a.anchor == b.anchor && a.direction == b.direction;
        }
        template <class Archive, class T, int N>
        inline void serialize(Archive & ar, Ray<T, N> & p) {
            ar(p.anchor, p.direction);
        }
        using Ray2 = Ray<double, 2>;
        using Ray3 = Ray<double, 3>;
        template <class T>
        inline Vec<T, 3> GetCoeffs(const Ray<T, 2> & line) {
            return Vec<T, 3>{line.direction[1], 
                -line.direction[0], 
                -(line.direction[1] * line.anchor[0] - line.direction[0] * line.anchor[1])
            };
        }
        template <class T>
        inline Ray<T, 2> Ray2FromCoeffs(const Vec<T, 3> & c) {
            T d = c[0] * c[0] + c[1] * c[1];
            Point<T, 2> anchor(-c[2] * c[0] / d, -c[2] * c[1] / d);
            Vec<T, 2> dir(c[1], -c[0]);
            return Ray<T, 2>(anchor, dir);
        }



        // plane
        template <class T, int N>
        struct Plane {
            inline Plane() {}
            inline Plane(const Point<T, N> & a, const Vec<T, N> & n) : anchor(a), normal(n) {}
            inline Point<T, N> root() const { return normal * (anchor.dot(normal)) / norm(normal) / norm(normal); }
            inline T signedDistanceTo(const Point<T, N> & p) const { return (p - anchor).dot(normalize(normal)); }
            inline T distanceTo(const Point<T, N> & p) const { return abs((p - anchor).dot(normalize(normal))); }
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
        template <class T>
        inline Plane<T, 3> Plane3From3Points(const Point<T, 3> & a, const Point<T, 3> & b, const Point<T, 3> & c){
            return Plane<T, 3>(a, normalize((b - a).cross(c - a)));
        }
        // ax + by + cz = 1
        template <class T>
        inline Plane<T, 3> Plane3FromEquation(T a, T b, T c){
            T k = (a * a + b * b + c * c);
            return Plane<T, 3>(Point<T, 3>(a, b, c) / k, normalize(Vec<T, 3>(a, b, c)));
        }
        template <class T>
        inline Vec<T, 3> Plane3ToEquation(const Plane<T, 3> & p){
            return p.normal / p.anchor.dot(p.normal);
        }



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
            inline Ray<T, N> infiniteLine() const { return Ray<T, N>{ first, second - first }; }
            inline Line & operator *= (const T & factor) { first *= factor; second *= factor; return *this; }
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
        template <class T, int N>
        inline Line<T, N> operator * (const Line<T, N> & line, const T & factor){
            return Line<T, N>(line.first * factor, line.second * factor);
        }
        template <class T, int N>
        inline Line<T, N> operator * (const T & factor, const Line<T, N> & line){
            return Line<T, N>(line.first * factor, line.second * factor);
        }
        template <class T, int N>
        inline Line<T, N> operator + (const Line<T, N> & line, const Vec<T, N> & trans){
            return Line<T, N>(line.first + trans, line.second + trans);
        }
        template <class T, int N>
        inline Line<T, N> operator - (const Line<T, N> & line, const Vec<T, N> & trans){
            return Line<T, N>(line.first - trans, line.second - trans);
        }
        using Line2 = Line<double, 2>;
        using Line3 = Line<double, 3>;

        template <class T, int N>
        inline Line<T, N> normalize(const Line<T, N> & line) { 
            return Line<T, N>(normalize(line.first), normalize(line.second)); 
        }


        
        // position on line/infline
        template <class T, int N>
        struct PositionOnLine {
            inline PositionOnLine(){}
            inline PositionOnLine(const Line<T, N> & line, const T & r)
                : ratio(r), position(line.first + (line.second - line.first) * ratio) {}
            inline PositionOnLine(const Ray<T, N> & line, const T & r)
                : ratio(r), position(line.anchor + line.direction * ratio) {}
            T ratio; // [0 ~ 1] on line, or [-inf, +inf] on infinite line
            // position = line.first + (line.second - line.fist) * ratio
            // or       = line.anchor + line.direction * ratio
            Point<T, N> position; 
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
        template <class T> using ImageOfType = cv::Mat_<T> ;
        using Imageb = ImageOfType<bool>;
        using Imageb3 = ImageOfType<Vec<bool, 3>>;
        using Imageub = ImageOfType<uint8_t>;
        using Imageub3 = ImageOfType<Vec<uint8_t, 3>>;
        using Imagei = ImageOfType<int>;
        using Imagei3 = ImageOfType<Vec<int, 3>>;
        using Imagef = ImageOfType<float>;
        using Imagef3 = ImageOfType<Vec<float, 3>>;
        using Imaged = ImageOfType<double>;
        using Imaged3 = ImageOfType<Vec<double, 3>>;


        namespace {
            template <class To, class From>
            inline ImageOfType<To> VecCastPrivate(const ImageOfType<From> & im, std::false_type) {
                ImageOfType<To> cim(im.size());
                for (auto it = im.begin(); it != im.end(); ++it){
                    cim(it.pos()) = static_cast<To>(*it);
                }
                return cim;
            }
            template <class To>
            inline ImageOfType<To> VecCastPrivate(const ImageOfType<To> & v, std::true_type) {
                return v.clone();
            }

            template <class To, class From, int N>
            inline ImageOfType<Vec<To, N>> VecCastPrivate(const ImageOfType<Vec<From, N>> & im, std::false_type) {
                ImageOfType<Vec<To, N>> cim(im.size());
                for (auto it = im.begin(); it != im.end(); ++it){
                    cim(it.pos()) = vec_cast<To>(*it);
                }
                return cim;
            }
            template <class To, int N>
            inline ImageOfType<Vec<To, N>> VecCastPrivate(const ImageOfType<Vec<To, N>> & v, std::true_type) {
                return v.clone();
            }
        }

        template <class To, class From>
        inline ImageOfType<To> vec_cast(const ImageOfType<From> & v) {
            return VecCastPrivate<To>(v, std::integral_constant<bool, std::is_same<To, From>::value>());
        }

        template <class To, class From, int N>
        inline ImageOfType<Vec<To, N>> vec_cast(const ImageOfType<Vec<From, N>> & v) {
            return VecCastPrivate<To>(v, std::integral_constant<bool, std::is_same<To, From>::value>());
        }


        using PixelLoc = cv::Point;
        template <class T>
        inline Vec<T, 2> vec_cast(const PixelLoc & p) { 
            return Vec<T, 2>(static_cast<T>(p.x), static_cast<T>(p.y)); 
        }

        inline int AreaOfImage(const Image & im) { return im.cols * im.rows; }
        void ResizeToWidth(Image & im, int width);
        void ResizeToHeight(Image & im, int height);
		void ResizeToMakeWidthUnder(Image & im, int widthUpperBound);
        void ResizeToMakeHeightUnder(Image & im, int heightUpperBound);
        bool MayBeAPanorama(const Image & im);
        bool MakePanorama(Image & im);
        
        std::pair<PixelLoc, PixelLoc> MinMaxLocOfImage(const Image & im);
        std::pair<double, double> MinMaxValOfImage(const Image & im);

        template <class T> 
        inline PixelLoc ToPixelLoc(const Point<T, 2> & p){
            return PixelLoc(static_cast<int>(p[0]), static_cast<int>(p[1]));
        }

        PixelLoc PixelLocFromGeoCoord(const GeoCoord & p, int longidiv, int latidiv);
        GeoCoord GeoCoordFromPixelLoc(const PixelLoc & pixel, int longidiv, int latidiv);


        // dense mat
        template <class T>
        using DenseMat = cv::Mat_<T>;
        using DenseMatd = DenseMat<double>;

        // sparse mat
        template <class T>
        using SparseMat = cv::SparseMat_<T>;
        using SparseMatd = SparseMat<double>;
        template <class T>
        struct SparseMatElement {
            using ValueType = T;
            int row, col;
            T value;
            inline SparseMatElement() : row(-1), col(-1) {}
            inline SparseMatElement(int r, int c, T v) : row(r), col(c), value(v) {}
        };
        template <class SparseMatElementIteratorT, class T = typename std::iterator_traits<SparseMatElementIteratorT>::value_type::ValueType>
        inline SparseMat<T> MakeSparseMatFromElements(int row, int col, 
            SparseMatElementIteratorT && begin, SparseMatElementIteratorT && end){
            int dims[] = { row, col };
            SparseMat<T> mat(2, dims);
            while (begin != end){
                mat.ref(begin->row, begin->col) = begin->value;
                ++begin;
            }
            return mat;
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


        // sphere
        template <class T, int N>
        struct Sphere {
            Point<T, N> center;
            T radius;
        };
        template <class T, int N>
        inline bool operator == (const Sphere<T, N> & a, const Sphere<T, N> & b){
            return a.center == b.center && a.radius == b.radius;
        }
        template <class Archive, class T, int N>
        inline void serialize(Archive & ar, Sphere<T, N> & s){
            ar(s.center, s.radius);
        }
        using Sphere2 = Sphere<double, 2>;
        using Sphere3 = Sphere<double, 3>;


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
            inline Sphere<T, N> outerSphere() const { 
                return Sphere<T, N>{center(), static_cast<T>(norm(maxCorner - minCorner) / 2.0)}; 
            }
            inline Box & expand(const Vec<T, N> & s) {
                minCorner -= s;
                maxCorner += s;
                return *this;
            }
            inline Box & expand(const T & s) {
                for (int i = 0; i < N; i++){
                    minCorner[i] -= s;
                    maxCorner[i] += s;
                }
                return *this;
            }

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

 
        // polygon
        template <class T, int N>
        struct Polygon {
            std::vector<core::Point<T, N>> corners;
            core::Vec<T, N> normal;
            inline Plane<T, N> plane() const { return Plane<T, N>(corners.front(), normal); }
        };

        template <class Archive, class T, int N>
        inline void serialize(Archive & ar, Polygon<T, N> & p) {
            ar(p.corners, p.normal);
        }

        using Polygon3 = Polygon<double, 3>;
        using Polygon3f = Polygon<float, 3>;

        template <class T>
        inline Polygon<T, 3> MakeTriangle(const Point<T, 3> & p1, const Point<T, 3> & p2, const Point<T, 3> & p3){
            return Polygon<T, 3>{{ p1, p2, p3 }, (p2 - p1).cross(p3 - p1)};
        }





        // something classified
        template <class T, class ClassT = int>
        struct Classified {
            ClassT claz;
            T component;
        };
        template <class T, class ClassT>
        inline Classified<std::decay_t<T>, std::decay_t<ClassT>> ClassifyAs(T && comp, ClassT && claz){
            return Classified<std::decay_t<T>, std::decay_t<ClassT>>{std::forward<ClassT>(claz), std::forward<T>(comp)};
        }
        template <class T, class ClassT>
        inline std::vector<Classified<T, ClassT>> ClassifyEachAs(const std::vector<T> & comps, const ClassT & claz){
            std::vector<Classified<T, ClassT>> cs(comps.size());
            for (int i = 0; i < comps.size(); i++)
                cs[i] = ClassifyAs(comps[i], claz);
            return cs;
        }
        template <class T, class ClassT>
        inline std::vector<Classified<T, ClassT>> ClassifyEachAs(std::vector<T> && comps, const ClassT & claz){
            std::vector<Classified<T, ClassT>> cs(comps.size());
            for (int i = 0; i < comps.size(); i++)
                cs[i] = ClassifyAs(std::move(comps[i]), claz);
            return cs;
        }
        template <class T, class ClassT>
        inline bool operator == (const Classified<T, ClassT> & a, const Classified<T, ClassT> & b) {
            return a.claz == b.claz && a.component == b.component;
        }
        template <class Archive, class T, class ClassT>
        inline void serialize(Archive & ar, Classified<T, ClassT> & c) {
            ar(c.claz, c.component);
        }


        // something noted
        template <class T>
        struct Noted {
            std::string note;
            T component;
        };
        template <class T, class StringT>
        inline Noted<std::decay_t<T>> NoteAs(T && comp, StringT && note){
            return Noted<std::decay_t<T>>{std::forward<StringT>(note), std::forward<T>(comp)};
        }
        template <class T>
        inline bool operator == (const Noted<T> & a, const Noted<T> & b) {
            return a.note == b.note && a.component == b.component;
        }
        template <class Archive, class T>
        inline void serialize(Archive & ar, Noted<T> & c) {
            ar(c.note, c.component);
        }


        // something scored for sorting
        template <class T, class S = double>
        struct Scored {
            S score;
            T component;
        };
        template <class T, class S>
        inline Scored<std::decay_t<T>, std::decay_t<S>> ScoreAs(T && comp, S && score){
            return Scored<std::decay_t<T>, std::decay_t<S>>{std::forward<S>(score), std::forward<T>(comp)};
        }
        // note this!!!
        template <class T, class S>
        inline bool operator == (const Scored<T, S> & a, const Scored<T, S> & b) {
            return a.score == b.score && a.component == b.component;
        }

        // score comparison
        template <class T, class S>
        inline bool operator < (const Scored<T, S> & a, const Scored<T, S> & b){
            return a.score < b.score;
        }
        template <class T, class S>
        inline bool operator > (const Scored<T, S> & a, const Scored<T, S> & b){
            return a.score > b.score;
        }
        template <class T, class S>
        inline bool operator <= (const Scored<T, S> & a, const Scored<T, S> & b){
            return a.score <= b.score;
        }
        template <class T, class S>
        inline bool operator >= (const Scored<T, S> & a, const Scored<T, S> & b){
            return a.score >= b.score;
        }
        template <class Archive, class T, class S>
        inline void serialize(Archive & ar, Scored<T, S> & c) {
            ar(c.score, c.component);
        }


        // something enabled
        template <class T>
        struct Enabled {
            bool enabled;
            T component;
        };
        template <class T>
        inline Enabled<std::decay_t<T>> EnableAs(T && comp, bool e = true){
            return Enabled<std::decay_t<T>>{e, std::forward<T>(comp)};
        }
        template <class T>
        inline bool operator == (const Enabled<T> & a, const Enabled<T> & b) {
            return a.enabled == b.enabled && a.component == b.component;
        }
        template <class Archive, class T>
        inline void serialize(Archive & ar, Enabled<T> & c) {
            ar(c.enabled, c.component);
        }


    }
}


 
#endif