#ifndef PANORAMIX_CORE_UTILITIES_HPP
#define PANORAMIX_CORE_UTILITIES_HPP

#include <iterator>
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
        inline bool HasValue(const Box<T, N> & b, const TesterT & tester){
            return HasValue(b.minCorner, tester) || HasValue(b.maxCorner, tester);
        }

        template <class T, class AllocT, class TesterT>
        inline bool HasValue(const std::vector<T, AllocT> & v, const TesterT & tester) {
            for (auto & e : v){
                if (HasValue(e, tester))
                    return true;
            }
            return false;
        }



        // squared
        template <class T>
        inline T Square(const T & v) {
            return v * v;
        }

        // gaussian
        template <class T, class K>
        inline T Gaussian(const T & x, const K & sigma) {
            return std::exp(- Square(x / sigma) / 2.0);
        }

        // pitfall
        template <class T, class K>
        inline T Pitfall(const T & x, const K & sigma) {
            return abs(x) <= sigma ? Square(x / sigma) : 1;
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



        /// bounding box functions

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
            return Box2(Point2(p.x, p.y), Point2(p.x, p.y));
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

        template <class T>
        inline auto BoundingBox(const Classified<T> & c) -> decltype(BoundingBox(c.component)) {
            return BoundingBox(c.component);
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




        // RTree Wrapper
        template <class T, class BoundingBoxFunctorT = DefaultBoundingBoxFunctor<T>>
        class RTreeWrapper {
        public:
            using BoxType = decltype(std::declval<BoundingBoxFunctorT>()(std::declval<T>()));
            using ValueType = typename BoxType::Type;
            static const int Dimension = BoxType::Dimension;

            inline explicit RTreeWrapper(BoundingBoxFunctorT getBB = BoundingBoxFunctorT())
                : _rtree(std::make_shared<third_party::RTree<T, ValueType, Dimension>>()),
                _getBoundingBox(getBB) {}

            template <class IteratorT>
            inline RTreeWrapper(IteratorT begin, IteratorT end, BoundingBoxFunctorT getBB = BoundingBoxFunctorT())
                : _rtree(std::make_shared<third_party::RTree<T, ValueType, Dimension>>()),
                _getBoundingBox(getBB){
                insert(begin, end);
            }

        public:
            inline size_t size() const { return _rtree->Count(); }
            inline bool empty() const { return size() == 0; }

            inline void clear() { return _rtree->RemoveAll(); }
            inline const BoundingBoxFunctorT & getBoundingBox() const { return _getBoundingBox; }

            inline void insert(const T & t) {
                BoxType box = _getBoundingBox(t);
                for (int i = 0; i < Dimension; i++){
                    if (isnan(box.minCorner[i]) || isnan(box.maxCorner[i])){
#ifdef _DEBUG
                        std::cout << "invalid box type (NaN value), ignore this element" << std::endl;
#endif
                        return;
                    }
                    if (!(box.minCorner[i] <= box.maxCorner[i])) {
#ifdef _DEBUG
                        std::cout << "invalid box type (minCorner[i] > maxCorner[i]), ignore this element" << std::endl;
#endif
                        return;
                    }
                }
                _rtree->Insert(box.minCorner.val, box.maxCorner.val, t);
            }

            template <class IteratorT>
            void insert(IteratorT begin, IteratorT end) {
                while (begin != end){
                    insert(*begin);
                    ++begin;
                }
            }

            template <class CallbackFunctorT>
            inline int search(const BoxType & b, CallbackFunctorT callback) const {
                return _rtree->Search(b.minCorner.val, b.maxCorner.val, callback);
            }

            template <class CallbackFunctorT>
            inline int searchNear(const T & t, CallbackFunctorT callback) const {
                auto b = _getBoundingBox(t);
                return _rtree->Search(b.minCorner.val, b.maxCorner.val, callback);
            }

            template <class IsEqualT = std::equal_to<T>>
            bool contains(const T & t, IsEqualT && cmp = IsEqualT()) const {
                bool exists = false;
                search(_getBoundingBox(t), [&exists, &t, &cmp](const T & ele) {
                    if (cmp(ele, t)) {
                        exists = true;
                        return false;
                    }
                    return true;
                });
                return exists;
            }

        private:
            std::shared_ptr<third_party::RTree<T, ValueType, Dimension>> _rtree;
            BoundingBoxFunctorT _getBoundingBox;
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

        // returns [0, pi]
        template <class T, int N>
        inline T AngleBetweenDirections(const Vec<T, N> & v1, const Vec<T, N> & v2) {
            auto s = v1.dot(v2) / norm(v1) / norm(v2);
            return s >= 1.0 ? 0.0 : (s <= -1.0 ? M_PI : acos(s));
        }

        // returns [0, pi/2]
        template <class T, int N>
        inline T AngleBetweenUndirectedVectors(const Vec<T, N> & v1, const Vec<T, N> & v2) {
            auto s = abs(v1.dot(v2) / norm(v1) / norm(v2));
            return s >= 1.0 ? 0.0 : acos(s);
        }

        template <class T, int N>
        inline bool IsApproxParallel(const Vec<T, N> & v1, const Vec<T, N> & v2, const T & epsilon = 0.1) {
            auto s = v1.dot(v2) / norm(v1) / norm(v2);
            return s >= 1.0 - epsilon || s <= -1.0 + epsilon;
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





        // generic algorithms

        // fill the container with linear sequence
        template <class IteratorT>
        void CreateLinearSequence(IteratorT begin, IteratorT end,
            const typename std::iterator_traits<IteratorT>::value_type& low,
            const typename std::iterator_traits<IteratorT>::value_type& high){
            auto dist = std::distance(begin, end);
            auto w = high - low;
            int id = 0;
            while (begin != end){
                *begin = (id++) * w / dist + low;
                ++begin;
            }
        }


        // merge, rearrange the input array
        // DistanceFunctorT(a, b) -> DistanceT : compute the distance from a to b
        // returns the begin iterators of merged groups
        template < class IteratorT, class IterOutIteratorT, class DistanceT,
        class DistanceFunctorT = DefaultDistanceFunctor < typename std::iterator_traits<IteratorT>::value_type >
        >
        IterOutIteratorT MergeNearNaive(IteratorT begin, IteratorT end, IterOutIteratorT itersOut, std::true_type,
            DistanceT thres, DistanceFunctorT distFun = DistanceFunctorT()) {
            if (begin == end)
                return itersOut;

            std::vector<IteratorT> gBegins(1, begin);
            for (auto i = std::next(begin); i != end; ++i){
                DistanceT minDist = std::numeric_limits<DistanceT>::max();
                auto nearestGBeginIter = gBegins.end();
                for (auto giter = gBegins.begin(); giter != gBegins.end(); ++giter) {
                    auto gBegin = *giter;
                    DistanceT dist = distFun(*gBegin, *i);
                    if (dist <= thres && dist < minDist){
                        minDist = dist;
                        nearestGBeginIter = giter;
                    }
                }
                if (nearestGBeginIter != gBegins.end()){ // found group
                    if (std::next(nearestGBeginIter) != gBegins.end()){
                        auto nextGBegin = *std::next(nearestGBeginIter);
                        std::rotate(nextGBegin, i, std::next(i));
                        for (auto j = std::next(nearestGBeginIter); j != gBegins.end(); ++j)
                            ++ (*j);
                    }
                }
                else { // add new group
                    gBegins.push_back(i);
                }
            }
            return std::copy(gBegins.begin(), gBegins.end(), itersOut);
        }

        // merge, without rearrangement
        // DistanceFunctorT(a, b) -> DistanceT : compute the distance from a to b
        // returns the iterators pointing to group leaders
        template < class IteratorT, class IterOutIteratorT, class DistanceT,
        class DistanceFunctorT = DefaultDistanceFunctor < typename std::iterator_traits<IteratorT>::value_type >
        >
        IterOutIteratorT MergeNearNaive(IteratorT begin, IteratorT end, IterOutIteratorT itersOut, std::false_type,
            DistanceT thres, DistanceFunctorT distFun = DistanceFunctorT()) {
            if (begin == end)
                return itersOut;

            *(itersOut++) = begin;
            std::vector<IteratorT> gBegins(1, begin);
            for (auto i = std::next(begin); i != end; ++i){
                auto giter = gBegins.begin();
                for (; giter != gBegins.end(); ++giter) {
                    auto gBegin = *giter;
                    auto dist = distFun(*gBegin, *i);
                    if (dist <= thres){
                        break;
                    }
                }
                if (giter == gBegins.end()){ // add new group
                    gBegins.push_back(i);
                    *(itersOut++) = i;
                }
            }

            return itersOut;
        }



        // merge using RTree, without rearrangement
        // DistanceFunctorT(a, b) -> ? : compute the distance from a to b
        // BoundingBoxFunctorT(a) -> Box<?,?> : compute the bounding box of a
        // returns the iterators pointing to group leaders
        template <class IteratorT, class IterOutIteratorT, class DistanceT,
        class DistanceFunctorT = DefaultDistanceFunctor<typename std::iterator_traits<IteratorT>::value_type>,
        class BoundingBoxFunctorT = DefaultBoundingBoxFunctor<typename std::iterator_traits<IteratorT>::value_type>
        >
        IterOutIteratorT MergeNearRTree(IteratorT begin, IteratorT end, IterOutIteratorT itersOut, std::false_type,
            DistanceT thres, DistanceFunctorT distFun = DistanceFunctorT(),
            BoundingBoxFunctorT getBoundingBox = BoundingBoxFunctorT()) {

            if (begin == end)
                return itersOut;

            using BoxType = decltype(getBoundingBox(*begin));
            using T = typename BoxType::Type;
            static const int N = BoxType::Dimension;

            third_party::RTree<IteratorT, T, N> rtree;
            for (auto i = begin; i != end; ++i){
                Box<T, N> box = getBoundingBox(*i);
                for (int k = 0; k < N; k++){ // extend the box
                    box.minCorner[k] -= thres * 2;
                    box.maxCorner[k] += thres * 2;
                }
                // search in RTree
                int foundCount = 0;
                rtree.Search(box.minCorner.val, box.maxCorner.val,
                    [distFun, i, thres, &foundCount](IteratorT it){
                    if (distFun(*i, *it) <= thres){
                        foundCount++;
                        return false;
                    }
                    return true;
                });
                if (foundCount == 0){
                    rtree.Insert(box.minCorner.val, box.maxCorner.val, i);
                    *(itersOut)++ = i;
                }
            }

            return itersOut;
        }



        // Minimum Spanning Tree
        // EdgeVertsGetterT(Edge e)->std::pair<Vert,Vert>
        // EdgeCompareOnWeightT(Edge e1, Edge e2)->bool 
        //     determins whether weight of e1 is lower than weight of e2
        // VertCompareT(Vert v1, Vert v2)->bool
        //     used in std::map to register set id of vertices
        template < class VertIteratorT, class EdgeIteratorT,
        class EdgeVertsGetterT,
        class EdgeOutputIteratorT,
        class EdgeCompareOnWeightT,
        class VertCompareT = std::less < typename std::iterator_traits<VertIteratorT>::value_type >
        >
        void MinimumSpanningTree(
            VertIteratorT vertsBegin, VertIteratorT vertsEnd,
            EdgeIteratorT edgesBegin, EdgeIteratorT edgesEnd,
            EdgeOutputIteratorT MSTedges,
            EdgeVertsGetterT vertsGetter,
            EdgeCompareOnWeightT edgeCompareOnWeight,
            VertCompareT vertCompare = VertCompareT()
        ) {

            using Edge = typename std::iterator_traits<typename EdgeIteratorT>::value_type;
            using Vert = typename std::iterator_traits<typename VertIteratorT>::value_type;
            static_assert(typename std::is_same<std::pair<Vert, Vert>, decltype(vertsGetter(*edgesBegin))>::value,
                "result of EdgeVertsGetterT must be std::pair<Vert, Vert>!");

            std::vector<Edge> edges(edgesBegin, edgesEnd);
            std::sort(edges.begin(), edges.end(), edgeCompareOnWeight);

            std::map<Vert, int, VertCompareT> vertSetIds(vertCompare);
            int idx = 0;
            for (auto i = vertsBegin; i != vertsEnd; ++i)
                vertSetIds.insert(std::make_pair((*i), idx++));

            auto remainedEdgesBegin = edges.begin();
            while (remainedEdgesBegin != edges.end()){
                Edge e = *remainedEdgesBegin;
                auto verts = vertsGetter(e);
                int fromid = vertSetIds[verts.first];
                int toid = vertSetIds[verts.second];
                if (fromid != toid){
                    *MSTedges++ = e;
                    for (auto & vtoid : vertSetIds){
                        if (vtoid.second == toid){
                            vtoid.second = fromid;
                        }
                    }
                }
                ++remainedEdgesBegin;
            }

        }


        // DepthFirstSearch
        template < class VertIteratorT, class NeighborVertsContainerGetterT, class VertCallbackT,
        class VertCompareT = std::less <typename std::iterator_traits<VertIteratorT>::value_type>
        >
        void DepthFirstSearch(VertIteratorT vertsBegin, VertIteratorT vertsEnd,
            NeighborVertsContainerGetterT neighborVertsContainerGetter,
            VertCallbackT vertCallback,
            VertCompareT vertCompare = VertCompareT()
        ) {

            using Vert = typename std::iterator_traits<typename VertIteratorT>::value_type;
            static_assert(std::is_same<Vert, 
                std::decay_t<decltype(*std::begin(neighborVertsContainerGetter(std::declval<Vert>())))>>::value,
                "NeighborVertsContainerGetterT should returns a container of Vert");
            static_assert(std::is_same<Vert, 
                std::decay_t<decltype(*std::end(neighborVertsContainerGetter(std::declval<Vert>())))>>::value,
                "NeighborVertsContainerGetterT should returns a container of Vert");

            struct {
                bool operator()(Vert root, std::map<Vert, bool, VertCompareT>& vVisited, 
                    NeighborVertsContainerGetterT vNeighborsGetter, VertCallbackT vCallback) {
                    if (vVisited[root])
                        return true;
                    if (!vCallback(root))
                        return false;

                    vVisited[root] = true;
                    auto vNeighborsContainer = vNeighborsGetter(root);
                    for (const auto & v : vNeighborsContainer) {
                        if (!(*this)(v, vVisited, vNeighborsGetter, vCallback))
                            return false;
                    }
                    return true;
                }
            } depthFirstSearchOneTree;

            std::map<Vert, bool, VertCompareT> visited(vertCompare);
            for (auto i = vertsBegin; i != vertsEnd; ++i)
                visited[*i] = false;
            while (true) {
                auto rootIter = vertsBegin;
                while (rootIter != vertsEnd && visited[*rootIter]){
                    ++rootIter;
                }
                if (rootIter == vertsEnd)
                    break;
                if (!depthFirstSearchOneTree(*rootIter, visited, neighborVertsContainerGetter, vertCallback))
                    break;
            }
        }



        // Topological Sort (using Depth First Search)
        template <class VertIteratorT, class VertOutIteratorT, class PredecessorVertsContainerGetterT,
        class VertCompareT = std::less<typename std::iterator_traits<VertIteratorT>::value_type>
        >
        void TopologicalSort(VertIteratorT vertsBegin, VertIteratorT vertsEnd,
            VertOutIteratorT sortedVertsBegin,
            PredecessorVertsContainerGetterT predecessorVertsContainerGetter,
            VertCompareT vertCompare = VertCompareT()){

            using Vert = typename std::iterator_traits<typename VertIteratorT>::value_type;
            static_assert(std::is_same<Vert,
                std::decay_t<decltype(*std::begin(predecessorVertsContainerGetter(std::declval<Vert>())))>> ::value,
                "PrecidentVertsContainerGetterT should returns a container of Vert");
            static_assert(std::is_same<Vert,
                std::decay_t<decltype(*std::end(predecessorVertsContainerGetter(std::declval<Vert>())))>> ::value,
                "PrecidentVertsContainerGetterT should returns a container of Vert");

            struct {
                void operator()(Vert root, std::map<Vert, bool, VertCompareT>& vVisited,
                    VertOutIteratorT sortedVertsOut,
                    PredecessorVertsContainerGetterT predecessorVertsContainerGetter) {
                    if (vVisited[root])
                        return;

                    vVisited[root] = true;
                    auto vPredecessors = predecessorVertsContainerGetter(root);
                    for (const auto & v : vPredecessors) {
                        (*this)(v, vVisited, sortedVertsOut, predecessorVertsContainerGetter);
                    }

                    *sortedVertsOut = root;
                    ++sortedVertsOut;
                }
            } depthFirstSearchOneTree;

            std::map<Vert, bool, VertCompareT> visited(vertCompare);
            for (auto i = vertsBegin; i != vertsEnd; ++i)
                visited[*i] = false;
            while (true) {
                auto rootIter = vertsBegin;
                while (rootIter != vertsEnd && visited[*rootIter]){
                    ++rootIter;
                }
                if (rootIter == vertsEnd)
                    break;
                depthFirstSearchOneTree(*rootIter, visited, sortedVertsBegin, predecessorVertsContainerGetter);
            }

        }


        // Connected Components
        template <class VertIteratorT, class NeighborVertsContainerGetterT, class VertexTypeRecorderT,
        class VertCompareT = std::less<typename std::iterator_traits<VertIteratorT>::value_type>
        >
        int ConnectedComponents(VertIteratorT vertsBegin, VertIteratorT vertsEnd,
            NeighborVertsContainerGetterT neighborVertsContainerGetter,
            VertexTypeRecorderT vertTypeRecorder,
            VertCompareT vertCompare = VertCompareT()) {

            using Vert = typename std::iterator_traits<typename VertIteratorT>::value_type;
            static_assert(std::is_same<Vert, std::decay_t<decltype(*std::begin(neighborVertsContainerGetter(std::declval<Vert>())))>>::value,
                "NeighborVertsContainerGetterT should returns a container of Vert");
            static_assert(std::is_same<Vert, std::decay_t<decltype(*std::end(neighborVertsContainerGetter(std::declval<Vert>())))>>::value,
                "NeighborVertsContainerGetterT should returns a container of Vert");

            struct {
                void operator()(const Vert & root, std::map<Vert, bool, VertCompareT>& vVisited, 
                    const NeighborVertsContainerGetterT & vNeighborsGetter, 
                    const VertexTypeRecorderT & vTypeRecorder, int cid) {

                    if (vVisited[root])
                        return;
                    vTypeRecorder(root, cid);
                    vVisited[root] = true;
                    auto vNeighborsContainer = vNeighborsGetter(root);
                    for (const auto & v : vNeighborsContainer) {
                        (*this)(v, vVisited, vNeighborsGetter, vTypeRecorder, cid);
                    }
                }
            } depthFirstSearchOneTree;


            std::map<Vert, bool, VertCompareT> visited(vertCompare);
            for (auto i = vertsBegin; i != vertsEnd; ++i)
                visited[*i] = false;

            int cid = 0;
            while (true) {
                auto rootIter = vertsBegin;
                while (rootIter != vertsEnd && visited[*rootIter]){
                    ++rootIter;
                }
                if (rootIter == vertsEnd)
                    break;
                depthFirstSearchOneTree(*rootIter, visited, neighborVertsContainerGetter, vertTypeRecorder, cid);
                cid++;
            }

            return cid;
        }


        // graph cut



    }
}
 
#endif