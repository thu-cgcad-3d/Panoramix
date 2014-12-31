#ifndef PANORAMIX_CORE_UTILITIES_HPP
#define PANORAMIX_CORE_UTILITIES_HPP

#include <iterator>
#include <Eigen/Dense>
#include "rtree.h"

#include "basic_types.hpp"
 
namespace panoramix {
    namespace core {

        // test value
        // can be used to check whether NaN exists by invoking: HasValue(a, std::isnan)
        template <class T, class TesterT, class = std::enable_if_t<std::is_floating_point<T>::value>>
        inline bool HasValue(const T & v, const TesterT & tester) {
            return tester(v);
        }

        template <class T, int N, class TesterT, class = std::enable_if_t<std::is_floating_point<T>::value>>
        inline bool HasValue(const Vec<T, N> & v, const TesterT & tester){
            for (int i = 0; i < N; i++){
                if (tester(v[i]))
                    return true;
            }
            return false;
        }

        template <class T, class S, class TesterT>
        inline bool HasValue(const Ratio<T, S> & r, const TesterT & tester){
            return HasValue(r.numerator, tester) || HasValue(r.denominator, tester);
        }

        template <class T, int N, class TesterT, class = std::enable_if_t<std::is_floating_point<T>::value>>
        inline bool HasValue(const Line<T, N> & v, const TesterT & tester){
            return HasValue(v.first, tester) || HasValue(v.second, tester);
        }

        template <class T, int N, class TesterT, class = std::enable_if_t<std::is_floating_point<T>::value>>
        inline bool HasValue(const InfiniteLine<T, N> & v, const TesterT & tester){
            return HasValue(v.anchor, tester) || HasValue(v.direction, tester);
        }

        template <class T, int N, class TesterT, class = std::enable_if_t<std::is_floating_point<T>::value>>
        inline bool HasValue(const Plane<T, N> & p, const TesterT & tester){
            return HasValue(p.anchor, tester) || HasValue(p.normal, tester);
        }

        template <class T, int N, class TesterT, class = std::enable_if_t<std::is_floating_point<T>::value>>
        inline bool HasValue(const Sphere<T, N> & s, const TesterT & tester){
            return HasValue(s.center) || tester(s.radius);
        }

        template <class T, int N, class TesterT, class = std::enable_if_t<std::is_floating_point<T>::value>>
        inline bool HasValue(const Box<T, N> & b, const TesterT & tester){
            return HasValue(b.minCorner, tester) || HasValue(b.maxCorner, tester);
        }

        template <class T1, class T2, class TesterT>
        inline bool HasValue(const std::pair<T1, T2> & p, const TesterT & tester){
            return HasValue(p.first) || HasValue(p.second);
        }

        template <class T, class AllocT, class TesterT>
        inline bool HasValue(const std::vector<T, AllocT> & v, const TesterT & tester) {
            for (auto & e : v){
                if (HasValue(e, tester))
                    return true;
            }
            return false;
        }




        template <class T, class = std::enable_if_t<std::is_floating_point<T>::value>>
        inline bool IsInfOrNaN(const T & v){
            return std::isinf(v) || std::isnan(v);
        }



        // contains
        template <class ContainerT, class ValueT>
        inline bool Contains(const ContainerT & c, const ValueT & v) {
            return std::find(c.cbegin(), c.cend(), v) != c.cend();
        }

        template <class KeyT, class ValueT, class PredT, class AllocT>
        inline bool Contains(const std::map<KeyT, ValueT, PredT, AllocT> & m, const KeyT & k) {
            return m.find(k) != m.end();
        }

        template <class KeyT, class PredT, class AllocT>
        inline bool Contains(const std::set<KeyT, PredT, AllocT> & m, const KeyT & k) {
            return m.find(k) != m.end();
        }

        template <class KeyT, class ValueT, class HasherT, class KeyeqT, class AllocT>
        inline bool Contains(const std::unordered_map<KeyT, ValueT, HasherT, KeyeqT, AllocT> & m, const KeyT & k) {
            return m.find(k) != m.end();
        }

        template <class KeyT, class HasherT, class KeyeqT, class AllocT>
        inline bool Contains(const std::unordered_set<KeyT, HasherT, KeyeqT, AllocT> & m, const KeyT & k) {
            return m.find(k) != m.end();
        }


        // all same
        template <class IteratorT, class EqualT = std::equal_to<void>>
        inline bool AllSameInRange(IteratorT begin, IteratorT end, EqualT eq = EqualT()){
            for (auto i = begin; i != end; ++i){
                if (!eq(*i, *begin))
                    return false;
            }
            return true;
        }

        template <class ContainerT, class EqualT = std::equal_to<void>>
        inline bool AllSameInContainer(const ContainerT & c, EqualT eq = EqualT()){
            return AllSameInRange(std::begin(c), std::end(c), eq);
        }



        // zero
        template <class T>
        INLINE bool IsFuzzyZero(const T & t, const T & epsilon = 1e-8){
            return t < epsilon && t > -epsilon;
        }


        // squared
        template <class T>
        INLINE T Square(const T & v) {
            return v * v;
        }

        // gaussian
        template <class T, class K>
        INLINE T Gaussian(const T & x, const K & sigma) {
            return std::exp(- Square(x / sigma) / 2.0);
        }

        // pitfall
        template <class T, class K>
        INLINE T Pitfall(const T & x, const K & sigma) {
            return abs(x) <= sigma ? Square(x / sigma) : 1;
        }

        // entropy [-factor, 0]
        template <class IteratorT, class T = double>
        inline T EntropyOfRange(IteratorT begin, IteratorT end, const T & factor = 1.0){
            T e = 0;
            while (begin != end){
                auto & v = *begin;
                auto ve = (v * log2(v));
                e -= (std::isnan(ve) ? 0.0 : ve);
                ++begin;
            }
            return e * factor;
        }

        template <class ContainerT, class T = double>
        inline T EntropyOfContainer(const ContainerT & cont, const T & factor = 1.0){
            return EntropyOfRange(std::begin(cont), std::end(cont), factor);
        }




        /// distance functions
        // for real numbers
        template <class T, class = std::enable_if_t<std::is_arithmetic<T>::value>>
        inline T Distance(const T & a, const T & b) {
            return std::abs(a - b);
        }

        // for complex numbers
        template <class T>
        inline T Distance(const std::complex<T> & a, const std::complex<T> & b) {
            return std::abs(a - b);
        }

        template <class T, int N>
        inline T Distance(const Point<T, N> & a, const Point<T, N> & b) {
            return norm(a - b);
        }

        inline double Distance(const PixelLoc & a, const PixelLoc & b) {
            return sqrt(Square(a.x - b.x) + Square(a.y - b.y));
        }

        inline double Distance(const KeyPoint & a, const KeyPoint & b) {
            return sqrt(Square(a.pt.x - b.pt.x) + Square(a.pt.y - b.pt.y));
        }

        template <class T, int N>
        inline T Distance(const HPoint<T, N> & a, const HPoint<T, N> & b) {
            return norm(a.value() - b.value());
        }

        template <class T, int N>
        inline T Distance(const PositionOnLine<T, N> & a, const PositionOnLine<T, N> & b) {
            return norm(a.position - b.position);
        }

        inline double Distance(const cv::Scalar & a, const cv::Scalar & b) {
            return norm(a - b);
        }


        // standard distance functor
        template <class T>
        struct DefaultDistanceFunctor {
            using DistanceType = decltype(Distance(std::declval<T>(), std::declval<T>()));
            inline DistanceType operator()(const T & a, const T & b) const {
                return Distance(a, b);
            }
        };






        /// bounding box functions for basic types

        // for scalars
        template <class T, class = std::enable_if_t<std::is_arithmetic<T>::value>>
        inline Box<T, 1> BoundingBox(const T & t) {
            return Box<T, 1>(Point<T, 1>(t), Point<T, 1>(t));
        }

        template <class T>
        inline Box<T, 2> BoundingBox(const std::complex<T> & c) {
            return Box<T, 2>(Point<T, 2>(c.real(), c.imag()), Point<T, 2>(c.real(), c.imag()));
        }


        template <class T, int N>
        inline Box<T, N> BoundingBox(const Box<T, N> & b){
            return b;
        }

        template <class T, int N>
        inline Box<T, N> BoundingBox(const Point<T, N> & p){
            return Box<T, N>(p, p);
        }

        template <class T, int N>
        inline Box<T, N> BoundingBox(const HPoint<T, N> & hp) {
            return BoundingBox(hp.value());
        }

        inline Box2 BoundingBox(const PixelLoc & p) {
            return Box2(
                Point2(static_cast<double>(p.x), static_cast<double>(p.y)),
                Point2(static_cast<double>(p.x), static_cast<double>(p.y)));
        }

        inline Box2 BoundingBox(const KeyPoint & p) {
            return Box2(Point2(p.pt.x, p.pt.y), Point2(p.pt.x, p.pt.y));
        }

        template <class T, int N>
        inline Box<T, N> BoundingBox(const Line<T, N> & l) {
            return Box<T, N>(l.first, l.second);
        }

        template <class T, int N>
        inline Box<T, N> BoundingBox(const HLine<T, N> & l) {
            return BoundingBox(l.toLine());
        }

        template <class T, int N>
        inline Box<T, N> BoundingBox(const PositionOnLine<T, N> & p) {
            return BoundingBox(p.position);
        }

        inline Box2 BoundingBox(const Image & im){
            return Box2(Point2(0, 0), Point2(im.cols, im.rows));
        }

        template <class T, int N>
        inline Box<T, N> BoundingBox(const Sphere<T, N> & s){
            return Box<T, N>(s.center).expand(s.radius);
        }

        template <class T, int N>
        inline Box<T, N> BoundingBox(const Polygon<T, N> & p){
            Box<T, N> b;
            for (auto & c : p.corners){
                b |= BoundingBox(c);
            }
            return b;
        }

        // pointers
        template <class T>
        inline auto BoundingBox(T const * const ptr) -> decltype(BoundingBox(*ptr)) {
            return BoundingBox(*ptr);
        }

        template <class T>
        inline auto BoundingBox(std::shared_ptr<T> ptr) -> decltype(BoundingBox(*ptr)) {
            return BoundingBox(*ptr);
        }

        template <class T>
        inline auto BoundingBox(std::unique_ptr<T> ptr) -> decltype(BoundingBox(*ptr)) {
            return BoundingBox(*ptr);
        }

        template <class T>
        inline auto BoundingBox(std::weak_ptr<T> ptr) -> decltype(BoundingBox(*ptr)) {
            return BoundingBox(*ptr);
        }

        // decorators
        template <class T>
        inline auto BoundingBox(const Classified<T> & c) -> decltype(BoundingBox(c.component)) {
            return BoundingBox(c.component);
        }

        template <class T>
        inline auto BoundingBox(const Noted<T> & n) -> decltype(BoundingBox(n.component)){
            return BoundingBox(n.component);
        }

        template <class T>
        inline auto BoundingBox(const Scored<T> & s) -> decltype(BoundingBox(s.component)){
            return BoundingBox(s.component);
        }

        // return null box if s.enabled == false
        template <class T>
        inline auto BoundingBox(const Enabled<T> & s) -> decltype(BoundingBox(s.component)){
            using BoxType = decltype(BoundingBox(s.component));
            return s.enabled ? BoundingBox(s.component) : BoxType();
        }
       

        // bounding box of range
        template <class IteratorT>
        auto BoundingBoxOfRange(IteratorT begin, IteratorT end) -> decltype(BoundingBox(*begin)) {
            using BoxType = decltype(BoundingBox(*begin));
            BoxType box; // a null box
            while (begin != end){
                auto b = BoundingBox(*begin);
                box |= b;
                ++begin;
            }
            return box;
        }

        // bounding box of container
        template <class ContainerT>
        inline auto BoundingBoxOfContainer(const ContainerT & cont)
            -> decltype(BoundingBoxOfRange(std::begin(cont), std::end(cont))) {
            return BoundingBoxOfRange(std::begin(cont), std::end(cont));
        }

        // bounding box of pair-range
        template <class PairIteratorT>
        auto BoundingBoxOfPairRange(PairIteratorT begin, PairIteratorT end) -> decltype(BoundingBox((*begin).second)) {
            using BoxType = decltype(BoundingBox((*begin).second));
            BoxType box;
            while (begin != end){
                auto b = BoundingBox((*begin).second);
                box |= b;
                ++begin;
            }
            return box;
        }


        // the default bounding box functor
        template <class T>
        struct DefaultBoundingBoxFunctor {
            using BoxType = decltype(BoundingBox(std::declval<T>()));
            inline BoxType operator()(const T & t) const {
                return BoundingBox(t);
            }
        };


        // the default influence box functor
        template <class T>
        struct DefaultInfluenceBoxFunctor {
            using BoxType = decltype(BoundingBox(std::declval<T>()));
            using ValueType = typename BoxType::Type;

            inline explicit DefaultInfluenceBoxFunctor(const ValueType & extSz = 0) : extendedSize(extSz){}
            inline BoxType operator()(const T & t) const {
                auto box = BoundingBox(t);
                for (int i = 0; i < BoxType::Dimension; i++){
                    box.minCorner[i] -= extendedSize;
                    box.maxCorner[i] += extendedSize;
                }
                return box;
            }
            const ValueType extendedSize;
        };







        // tools
        template <class T, class K>
        inline bool FuzzyEquals(const T & a, const T & b, const K & epsilon){
            return Distance(a, b) <= epsilon; // not only numerics
        }

        template <class T>
        inline int DiracDelta(const T & v) {
            return v == 0 ? 1 : 0;
        }

        template <class T, class K1, class K2>
        inline bool IsBetween(const T & v, const K1 & low, const K2 & high) {
            return !(v < low) && !(high < v);
        }

        template <class T, class K1, class K2>
        inline T BoundBetween(const T & v, const K1 & low, const K2 & high) {
            if (v < low)
                return low;
            return v < high ? v : high;
        }

        namespace {

            template <class T>
            T WrapBetweenPrivate(const T & input, const T & low, const T & high, const std::false_type &) {
                if (low >= high)
                    return input;
                if (low <= input && input < high)
                    return input;
                const auto sz = high - low;
                auto result = input - int((input - low) / sz) * sz + (input < low ? sz : 0);
                return result == high ? low : result;
            }

            template <class T>
            T WrapBetweenPrivate(const T & input, const T & low, const T & high, const std::true_type &) {
                if (low >= high)
                    return input;
                if (low <= input && input < high)
                    return input;
                const auto sz = high - low;
                auto result = (input - low) % sz + low + (input < low ? sz : 0);
                return result == high ? low : result;
            }

        }

        // returns [low, high)
        template <class T>
        inline T WrapBetween(const T & input, const T & low, const T & high) {
            return WrapBetweenPrivate(input, low, high, std::is_integral<T>());
        }


        template <class T, int N>
        T EncodeSubscriptToIndex(const Point<T, N> & subscript, const Vec<T, N> & dimension) {
            T index = subscript[0];
            for (int i = 1; i < N; i++){
                index = index * dimension[i] + subscript[i];
            }
            return index;
        }

        namespace {
            template <class T, int N>
            Point<T, N> DecodeIndexToSubscriptPrivate(T index, const Vec<T, N> & dimension, const std::false_type &) {
                Point<T, N> subscript;
                for (int i = N - 1; i >= 0; i--){
                    subscript[i] = WrapBetween(index, T(0.0), dimension[i]);
                    index = (index - subscript[i]) / dimension[i];
                }
                return subscript;
            }

            template <class T, int N>
            Point<T, N> DecodeIndexToSubscriptPrivate(T index, const Vec<T, N> & dimension, const std::true_type &) {
                Point<T, N> subscript;
                for (int i = N - 1; i >= 0; i--){
                    subscript[i] = index % dimension[i];
                    index = index / dimension[i];
                }
                return subscript;
            }
        }

        template <class T, int N>
        inline Point<T, N> DecodeIndexToSubscript(T index, const Vec<T, N> & dimension) {
            return DecodeIndexToSubscriptPrivate(index, dimension, std::is_integral<T>());
        }







        // for vectors

        template <class T>
        inline std::pair<Vec<T, 3>, Vec<T, 3>> ProposeXYDirectionsFromZDirection(const Vec<T, 3> & z) {
            auto y = z.cross(Vec<T, 3>(1, 0, 0));
            if (y(0) == 0 && y(1) == 0 && y(2) == 0) {
                y = z.cross(Vec<T, 3>(0, 1, 0));
            }
            auto x = y.cross(z);
            return std::make_pair(normalize(x), normalize(y));
        }

        // returns [0, pi]
        template <class T, int N>
        inline T AngleBetweenDirections(const Vec<T, N> & v1, const Vec<T, N> & v2) {
            auto s = v1.dot(v2) / norm(v1) / norm(v2);
            return s >= 1.0 - 1e-9 ? 0.0 : (s <= -1.0 + 1e-9 ? M_PI : acos(s));
        }

        // returns [0, pi/2]
        template <class T, int N>
        inline T AngleBetweenUndirectedVectors(const Vec<T, N> & v1, const Vec<T, N> & v2) {
            auto s = abs(v1.dot(v2) / norm(v1) / norm(v2));
            return s >= 1.0 ? 0.0 : acos(s);
        }

        template <class T, int N>
        inline bool IsFuzzyParallel(const Vec<T, N> & v1, const Vec<T, N> & v2, const T & epsilon = 0.1) {
            auto s = v1.dot(v2) / norm(v1) / norm(v2);
            return s >= 1.0 - epsilon || s <= -1.0 + epsilon;
        }

        template <class T, int N>
        inline bool IsFuzzyPerpenducular(const Vec<T, N> & v1, const Vec<T, N> & v2, const T & epsilon = 0.1) {
            auto s = v1.dot(v2) / norm(v1) / norm(v2);
            return abs(s) <= epsilon;
        }

        template <class T>
        inline T SignedAngleBetweenDirections(const Vec<T, 2> & from, const Vec<T, 2> & to,
            bool clockwiseAsPositive = true) {
            double angle = atan2(-from(0)*to(1) + to(0)*from(1), from(1)*to(1) + from(0)*to(0));
            return clockwiseAsPositive ? angle : -angle;
        }

        // for lines and points
        // returns projection position
        template <class T, int N>
        PositionOnLine<T, N> ProjectionOfPointOnLine(const Point<T, N> & p, const Line<T, N> & line) {
            Vec<T, N> lineDir = line.direction();
            lineDir /= norm(lineDir);
            T projRatio = (p - line.first).dot(lineDir) / line.length();
            return PositionOnLine<T, N>(line, projRatio);
        }

        // returns (distance, nearest point)
        template <class T, int N>
        std::pair<T, Point<T, N>> DistanceFromPointToLine(const Point<T, N> & p, const InfiniteLine<T, N> & line) {
            Vec<T, N> lineDir = line.direction;
            lineDir /= norm(lineDir);
            auto root = (p - line.anchor).dot(lineDir) * lineDir + line.anchor;
            return std::make_pair(norm(p - root), root);
        }

        // returns signed distance
        template <class T>
        T SignedDistanceFromPointToLine(const Point<T, 2> & p, const InfiniteLine<T, 2> & line) {
            auto coeffs = GetCoeffs(line);
            return (coeffs[0] * p[0] + coeffs[1] * p[1] + coeffs[2]) / sqrt(Square(coeffs[0]) + Square(coeffs[1]));
        }

        // returns (distance, nearest position)
        template <class T, int N>
        std::pair<T, PositionOnLine<T, N>> DistanceFromPointToLine(const Point<T, N> & p, const Line<T, N> & line) {
            Vec<T, N> lineDir = line.direction();
            lineDir /= norm(lineDir);
            T projRatio = (p - line.first).dot(lineDir) / line.length();
            projRatio = BoundBetween(projRatio, 0, 1);
            PositionOnLine<T, N> pos(line, projRatio);
            return std::make_pair(norm(p - pos.position), pos);
        }

        // returns (distance, (nearest point on line1, nearest point on line2))
        // see http://geomalgorithms.com/a07-_distance.html for explainations
        template <class T, int N>
        std::pair<T, std::pair<Point<T, N>, Point<T, N>>> DistanceBetweenTwoLines(
            const InfiniteLine<T, N> & line1, const InfiniteLine<T, N> & line2) {
            
            auto u = normalize(line1.direction);
            auto v = normalize(line2.direction);
            auto w = line1.anchor - line2.anchor;
            auto a = u.dot(u);         // always >= 0
            auto b = u.dot(v);
            auto c = v.dot(v);         // always >= 0
            auto d = u.dot(w);
            auto e = v.dot(w);

            T P = a * c - b * b;
            if (P == 0){
                auto ppair = std::make_pair(line1.anchor, 
                    DistanceFromPointToLine(line1.anchor, line2).second);
                return std::make_pair(Distance(ppair.first, ppair.second), ppair);
            }

            T sc = (b*e - c*d) / P;
            T tc = (a*e - b*d) / P;
            
            auto p1 = line1.anchor + sc * u;
            auto p2 = line2.anchor + tc * v;
            return std::make_pair(Distance(p1, p2), std::make_pair(p1, p2));
        }

        // returns (distance, (nearest position on line1, nearest position on line2))
        // see http://geomalgorithms.com/a07-_distance.html for explainations
        template <class T, int N>
        std::pair<T, std::pair<PositionOnLine<T, N>, PositionOnLine<T, N>>> DistanceBetweenTwoLines(
            const Line<T, N> & line1, const Line<T, N> & line2) {

            auto u = line1.direction();
            auto v = line2.direction();
            auto w = line1.first - line2.first;
            auto a = u.dot(u);         // always >= 0
            auto b = u.dot(v);
            auto c = v.dot(v);         // always >= 0
            auto d = u.dot(w);
            auto e = v.dot(w);
            auto D = a*c - b*b;        // always >= 0
            T    sc, sN, sD = D;       // sc = sN / sD, default sD = D >= 0
            T    tc, tN, tD = D;       // tc = tN / tD, default tD = D >= 0

            static const double SMALL_NUM = 0;
            // compute the line parameters of the two closest points
            if (D <= SMALL_NUM) { // the lines are almost parallel
                sN = 0.0;         // force using point P0 on segment S1
                sD = 1.0;         // to prevent possible division by 0.0 later
                tN = e;
                tD = c;
            } else {                 // get the closest points on the infinite lines
                sN = (b*e - c*d);
                tN = (a*e - b*d);
                if (sN < 0.0) {        // sc < 0 => the s=0 edge is visible
                    sN = 0.0;
                    tN = e;
                    tD = c;
                } else if (sN > sD) {  // sc > 1  => the s=1 edge is visible
                    sN = sD;
                    tN = e + b;
                    tD = c;
                }
            }

            if (tN < 0.0) {            // tc < 0 => the t=0 edge is visible
                tN = 0.0;
                // recompute sc for this edge
                if (-d < 0.0)
                    sN = 0.0;
                else if (-d > a)
                    sN = sD;
                else {
                    sN = -d;
                    sD = a;
                }
            } else if (tN > tD) {      // tc > 1  => the t=1 edge is visible
                tN = tD;
                // recompute sc for this edge
                if ((-d + b) < 0.0)
                    sN = 0;
                else if ((-d + b) > a)
                    sN = sD;
                else {
                    sN = (-d + b);
                    sD = a;
                }
            }

            // finally do the division to get sc and tc
            sc = (IsBetween(sN, -SMALL_NUM, +SMALL_NUM) ? 0.0 : sN / sD);
            tc = (IsBetween(tN, -SMALL_NUM, +SMALL_NUM) ? 0.0 : tN / tD);

            PositionOnLine<T, N> pos1(line1, sc);
            PositionOnLine<T, N> pos2(line2, tc);
            auto dist = norm(pos1.position - pos2.position);
            return std::make_pair(dist, std::make_pair(pos1, pos2));
        }

        // intersecton between line and plane
        template <class T, int N>
        PositionOnLine<T, N> IntersectionOfLineAndPlane(const InfiniteLine<T, N> & line, const Plane<T, N> & plane){
            T lambda = (plane.anchor - line.anchor).dot(plane.normal) / line.direction.dot(plane.normal);
            return PositionOnLine<T, N>(line, lambda);
        }

        // distance from point to plane
        template <class T, int N>
        inline std::pair<T, Point<T, N>> DistanceFromPointToPlane(const Point<T, N> & p, const Plane<T, N> & plane) {
            auto signedDist = plane.signedDistanceTo(p);
            return std::make_pair(abs(signedDist), p - signedDist * normalize(plane.normal));
        }


        template <class T, int N>
        Vec<T, N> BarycentricCoordinatesOfLineAndPlaneUnitIntersection(
            const InfiniteLine<T, N> & line, const Point<T, N> * cornersData){
            using namespace Eigen;
            Map<const Matrix<T, N, N, Eigen::ColMajor>> corners((const T*)cornersData, N, N);
            Map<const Matrix<T, N, 1>> D(line.direction.val, N, 1);
            Map<const Matrix<T, N, 1>> A(line.anchor.val, N, 1);
            
            Matrix<T, N, N> M = - (corners.colwise() - corners.col(0)); // 0 : p1-p2 : p1-p3
            M.col(0) = D; // d : p1-p2 : p1-p3
            
            Vec<T, N> coord;
            Map<Matrix<T, N, 1>> X(coord.val, N, 1);
            X = M.fullPivLu().solve(corners.col(0) - A);
            T depthOfLineIntersection = X(0);
            X(0) = 1.0f;
            for (int i = 1; i < N; i++){
                X(0) -= X(i);
            }
            return coord;
        }

        template <class T, int N>
        inline Vec<T, N> BarycentricCoordinatesOfLineAndPlaneUnitIntersection(
            const InfiniteLine<T, N> & line, const std::array<Point<T, N>, N> & corners){
            return BarycentricCoordinatesOfLineAndPlaneUnitIntersection(line, corners.data());
        }




        // eigen vectors and eigen values from points
        template <class T, int N>
        std::array<Scored<Vec<T, N>, T>, N> EigenVectorAndValuesFromPoints(const Point<T, N> * ptsData, size_t n) {
            using namespace Eigen;
            Map<const Matrix<T, Dynamic, N, RowMajor>> mat((const T*)ptsData, n, N);
            //std::cout << "points data mat: " << std::endl << mat << std::endl;
            auto centered = (mat.rowwise() - mat.colwise().mean()).eval();
            auto cov = ((centered.adjoint() * centered) / double(mat.rows() - 1)).eval();
            EigenSolver<decltype(cov)> es(cov);
            std::array<Scored<Vec<T, N>, T>, N> result;
            for (int i = 0; i < N; i++){
                for (int k = 0; k < N; k++){
                    result[i].component(k) = es.eigenvectors()(i, k).real();
                }
                result[i].score = es.eigenvalues()(i).real();
            }
            return result;
        }

        template <class T, int N, int M>
        inline std::array<Scored<Vec<T, N>, T>, N> EigenVectorAndValuesFromPoints(const Point<T, N>(&ptsData)[M]) {
            return EigenVectorAndValuesFromPoints((const Point<T, N>*)ptsData, M);
        }

        template <class T, int N>
        inline std::array<Scored<Vec<T, N>, T>, N> EigenVectorAndValuesFromPoints(const std::vector<Point<T, N>> & pts) {
            return EigenVectorAndValuesFromPoints((const Point<T, N>*)pts.data(), pts.size());
        }


        // on left
        template <class T>
        inline bool IsOnLeftSide(const Point<T, 2> & p, const Point<T, 2> & a, const Point<T, 2> & b){
            T data[9] = { a[0], a[1], 1, b[0], b[1], 1, p[0], p[1], 1 };
            T tmp1 = data[0 * 3 + 0] * (data[1 * 3 + 1] * data[2 * 3 + 2] - data[1 * 3 + 2] * data[2 * 3 + 1]);
            T tmp2 = data[0 * 3 + 1] * (data[1 * 3 + 0] * data[2 * 3 + 2] - data[1 * 3 + 2] * data[2 * 3 + 0]);
            T tmp3 = data[0 * 3 + 2] * (data[1 * 3 + 0] * data[2 * 3 + 1] - data[1 * 3 + 1] * data[2 * 3 + 0]);
            return tmp1 - tmp2 + tmp3 > 0;
        }

        // in triangle
        template <class T>
        inline bool IsInTriangle(const Point<T, 2> & p, const Point<T, 2> & a, const Point<T, 2> & b, const Point<T, 2> & c){
            bool lab = IsOnLeftSide(p, a, b);
            bool lbc = IsOnLeftSide(p, b, c);
            bool lca = IsOnLeftSide(p, c, a);
            return lab == lbc && lbc == lca;
        }



        // camera functions with matrix
        // make a lookat view matrix
        template <class T>
        Mat<T, 4, 4> MakeMat4LookAt(const Vec<T, 3> & eye, const Vec<T, 3> & center,
            const Vec<T, 3> & up) {
            Vec<T, 3> zaxis = (center - eye); zaxis /= core::norm(zaxis);
            Vec<T, 3> xaxis = up.cross(zaxis); xaxis /= core::norm(xaxis);
            Vec<T, 3> yaxis = zaxis.cross(xaxis);
            Mat<T, 4, 4> m(
                -xaxis(0), yaxis(0), -zaxis(0), 0,
                -xaxis(1), yaxis(1), -zaxis(1), 0,
                -xaxis(2), yaxis(2), -zaxis(2), 0,
                xaxis.dot(eye), -yaxis.dot(eye), zaxis.dot(eye), 1);
            return m.t();
        }

        // make a perspective projection matrix
        template <class T>
        Mat<T, 4, 4> MakeMat4Perspective(const T & fovyRadians, const T & aspect,
            const T & nearZ, const T & farZ) {
            T cotan = T(1.0) / std::tan(fovyRadians / 2.0);
            Mat<T, 4, 4> m(
                cotan / aspect, 0, 0, 0,
                0, cotan, 0, 0,
                0, 0, (farZ + nearZ) / (nearZ - farZ), -1,
                0, 0, (2 * farZ * nearZ) / (nearZ - farZ), 0);
            return m.t();
        }



        // generate Fibonacci Directions
        template <class AddVec3FunT>
        inline void GenerateFibonacciDirections(int n, AddVec3FunT addVec3){
            static const float ga = M_PI * (- 1.0f + std::sqrt(5.0f)) / 2.0f;
            for (int i = 0; i < n; i++){
                auto d = GeoCoord(ga * (float)i, asin(-1.0f + 2.0f * float(i) / n)).toVector();
                addVec3(d[0], d[1], d[2]);
            }
        }



        // mesh makers
        // make tetrahedron
        template <class AddVertex3FunT, class AddTriFaceFunT>
        void MakeTetrahedron(AddVertex3FunT addVertex, AddTriFaceFunT addFace) {
            auto v1 = addVertex(0.0f, 0.0f, 0.0f);
            auto v2 = addVertex(0.0f, 0.0f, 1.0f);
            auto v3 = addVertex(0.0f, 1.0f, 0.0f);
            auto v4 = addVertex(1.0f, 0.0f, 0.0f);

            addFace(v1, v2, v3);
            addFace(v1, v4, v2);
            addFace(v1, v3, v4);
            addFace(v2, v4, v3);
        }

        // make quad faced cube
        template <class AddVertex3FunT, class AddQuadFaceFunT>
        void MakeQuadFacedCube(AddVertex3FunT addVertex, AddQuadFaceFunT addFace) {
            /*
                4 ----- 5
               /       /|
              0 ----- 1 |
              |	      | |
              | 7	  | 6  -- x
              |	      |/
              3 ----- 2
             /
            y

            */
            auto v1 = addVertex(0.0f, 1.0f, 1.0f);
            auto v2 = addVertex(1.0f, 1.0f, 1.0f);
            auto v3 = addVertex(1.0f, 1.0f, 0.0f);
            auto v4 = addVertex(0.0f, 1.0f, 0.0f);

            auto v5 = addVertex(0.0f, 0.0f, 1.0f);
            auto v6 = addVertex(1.0f, 0.0f, 1.0f);
            auto v7 = addVertex(1.0f, 0.0f, 0.0f);
            auto v8 = addVertex(0.0f, 0.0f, 0.0f);

            addFace(v1, v2, v3, v4);
            addFace(v2, v6, v7, v3);
            addFace(v6, v5, v8, v7);
            addFace(v5, v1, v4, v8);
            addFace(v5, v6, v2, v1);
            addFace(v4, v3, v7, v8);
        }

        // make tri faced cube
        template <class AddVertex3FunT, class AddTriFaceFunT>
        void MakeTriFacedCube(AddVertex3FunT addVertex, AddTriFaceFunT addFace) {
            auto v1 = addVertex(0.0f, 1.0f, 1.0f);
            auto v2 = addVertex(1.0f, 1.0f, 1.0f);
            auto v3 = addVertex(1.0f, 1.0f, 0.0f);
            auto v4 = addVertex(0.0f, 1.0f, 0.0f);

            auto v5 = addVertex(0.0f, 0.0f, 1.0f);
            auto v6 = addVertex(1.0f, 0.0f, 1.0f);
            auto v7 = addVertex(1.0f, 0.0f, 0.0f);
            auto v8 = addVertex(0.0f, 0.0f, 0.0f);

            addFace(v1, v2, v3);
            addFace(v1, v3, v4);

            addFace(v2, v6, v7);
            addFace(v2, v7, v3);
            
            addFace(v6, v5, v8);
            addFace(v6, v8, v7);
            
            addFace(v5, v1, v4);
            addFace(v5, v4, v8);
            
            addFace(v5, v6, v2);
            addFace(v5, v2, v1);

            addFace(v4, v3, v7);
            addFace(v4, v7, v8);
        }

        // make quad faced sphere
        template <class AddVertex3FunT, class AddQuadFaceFunT>
        void MakeQuadFacedSphere(AddVertex3FunT addVertex, AddQuadFaceFunT addFace, int m, int n) {
            using ThisVertHandle = decltype(addVertex(0.0f, 0.0f, 0.0f));
            std::vector<std::vector<ThisVertHandle>> vhs(m, std::vector<ThisVertHandle>(n - 1));
            for (int i = 0; i < m; i++){
                for (int j = 0; j < n - 1; j++){
                    double xratio = 1.0f - 1.0f / (n - 1) * j;
                    double yratio = 1.0f / (m - 1) * i;
                    double xangle = M_PI * 2 * xratio;
                    double yangle = M_PI * yratio - M_PI_2;
                    vhs[i][j] = addVertex(sin(xangle - M_PI_2)*cos(yangle), cos(xangle - M_PI_2)*cos(yangle), sin(yangle));
                }
            }
            for (int i = 1; i < m; i++){
                for (int j = 1; j < n - 1; j++){
                    addFace(vhs[i][j], vhs[i][j - 1], vhs[i - 1][j - 1], vhs[i - 1][j]);
                }
                addFace(vhs[i][0], vhs[i][n - 2], vhs[i - 1][n - 2], vhs[i - 1][0]);
            }
        }

        // make tri faced sphere
        template <class AddVertex3FunT, class AddTriFaceFunT>
        void MakeTriFacedSphere(AddVertex3FunT addVertex, AddTriFaceFunT addFace, int m, int n) {
            using VertHandle = decltype(addVertex(0.0f, 0.0f, 0.0f));
            std::vector<std::vector<VertHandle>> vhs(m, std::vector<VertHandle>(n - 1));
            for (int i = 0; i < m; i++){
                for (int j = 0; j < n - 1; j++){
                    double xratio = 1.0f - 1.0f / (n - 1) * j;
                    double yratio = 1.0f / (m - 1) * i;
                    double xangle = M_PI * 2 * xratio;
                    double yangle = M_PI * yratio - M_PI_2;
                    vhs[i][j] = addVertex(sin(xangle - M_PI_2)*cos(yangle), cos(xangle - M_PI_2)*cos(yangle), sin(yangle));
                }
            }
            for (int i = 1; i < m; i++){
                for (int j = 1; j < n - 1; j++){
                    addFace(vhs[i][j], vhs[i][j - 1], vhs[i - 1][j - 1]);
                    addFace(vhs[i][j], vhs[i - 1][j - 1], vhs[i - 1][j]);
                }
                addFace(vhs[i][0], vhs[i][n - 2], vhs[i - 1][n - 2]);
                addFace(vhs[i][j], vhs[i - 1][j - 1], vhs[i - 1][j]);
            }
        }

        // make an icosahedron
        template <class AddVertex3FunT, class AddTriFaceFunT>
        void MakeIcosahedron(AddVertex3FunT addVertex, AddTriFaceFunT addFace) {
            using VertHandle = decltype(addVertex(0.0f, 0.0f, 0.0f));
            // create a basic icosahedron
            static const float t = (1.0f + std::sqrt(5.0f)) / 2.0f;
            VertHandle vhs[] = {
                addVertex(-1, t, 0.0f),
                addVertex(1, t, 0.0f),
                addVertex(-1, -t, 0.0f),
                addVertex(1, -t, 0.0f),
                addVertex(0.0f, -1, t),
                addVertex(0.0f, 1, t),
                addVertex(0.0f, -1, -t),
                addVertex(0.0f, 1, -t),
                addVertex(t, 0.0f, -1),
                addVertex(t, 0.0f, 1),
                addVertex(-t, 0.0f, -1),
                addVertex(-t, 0.0f, 1)
            };
            addFace(vhs[0], vhs[11], vhs[5]);
            addFace(vhs[0], vhs[5], vhs[1]);
            addFace(vhs[0], vhs[1], vhs[7]);
            addFace(vhs[0], vhs[7], vhs[10]);
            addFace(vhs[0], vhs[10], vhs[11]);
            addFace(vhs[1], vhs[5], vhs[9]);
            addFace(vhs[5], vhs[11], vhs[4]);
            addFace(vhs[11], vhs[10], vhs[2]);
            addFace(vhs[10], vhs[7], vhs[6]);
            addFace(vhs[7], vhs[1], vhs[8]);
            addFace(vhs[3], vhs[9], vhs[4]);
            addFace(vhs[3], vhs[4], vhs[2]);
            addFace(vhs[3], vhs[2], vhs[6]);
            addFace(vhs[3], vhs[6], vhs[8]);
            addFace(vhs[3], vhs[8], vhs[9]);
            addFace(vhs[4], vhs[9], vhs[5]);
            addFace(vhs[2], vhs[4], vhs[11]);
            addFace(vhs[6], vhs[2], vhs[10]);
            addFace(vhs[8], vhs[6], vhs[7]);
            addFace(vhs[9], vhs[8], vhs[1]);
        }


        // make mesh from simple mesh file
        template <class AddVertex3FunT, class AddTriFaceFunT>
        void MakeTriMeshFromSMFFile(AddVertex3FunT addVertex, AddTriFaceFunT addFace, const std::string & fileName) {
            using VertHandle = decltype(addVertex(0.0f, 0.0f, 0.0f));
            std::vector<VertHandle> vhs;
            std::vector<std::array<int, 3>> fvids;
            int minFvid = std::numeric_limits<int>::max();
            std::ifstream ifs(fileName);
            std::string line;
            while (std::getline(ifs, line)){
                std::istringstream iss(line);
                char tag;
                iss >> tag;
                if (tag == 'v'){
                    float x, y, z;
                    iss >> x >> y >> z;
                    vhs.push_back(addVertex(x, y, z));
                }
                else if (tag == 'f'){
                    int a, b, c;
                    iss >> a >> b >> c;
                    fvids.push_back({ { a, b, c } });
                    minFvid = std::min({ minFvid, a, b, c });
                }
            }
            // install faces
            for (auto & fvs : fvids){
                addFace(vhs[fvs[0] - minFvid], vhs[fvs[1] - minFvid], vhs[fvs[2] - minFvid]);
            }
        }


       

    }
}
 
#endif