#ifndef PANORAMIX_CORE_UTILITIES_HPP
#define PANORAMIX_CORE_UTILITIES_HPP

#include <iterator>
#include "rtree.h"

#include "basic_types.hpp"
 
namespace panoramix {
    namespace core {

        // squared
        template <class T>
        inline T Square(const T & v) {
            return v * v;
        }


        /// distance functions
        // for real numbers
        template <class T>
        inline std::enable_if_t<std::is_arithmetic<T>::value, T> Distance(const T & a, const T & b) {
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
            return norm(a.toPoint() - b.toPoint());
        }

        template <class T, int N>
        inline T Distance(const PositionOnLine<T, N> & a, const PositionOnLine<T, N> & b) {
            return norm(a.position - b.position);
        }


        // standard distance functor
        template <class T>
        struct DefaultDistanceFunctor {
            using DistanceType = decltype(Distance(std::declval<T>(), std::declval<T>()));
            inline DistanceType operator()(const T & a, const T & b) const {
                return Distance(a, b);
            }
        };



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


        /// bounding box functions

        // box
        template <class T, int N>
        struct Box {
            using Type = T;
            static const int Dimension = N;

            Point<T, N> minCorner, maxCorner;
            bool isNull;

            inline Box() : isNull(true){}
            inline Box(const Point<T, N> & c1, const Point<T, N> & c2)
                : minCorner(MakeMin(c1, c2)), maxCorner(MakeMax(c1, c2)), isNull(false) {}

            inline Vec<T, N> size() const { return maxCorner - minCorner; }
            inline Point<T, N> center() const { return (maxCorner + minCorner) * (0.5); }
            inline bool contains(const Point<T, N> & p) const {
                if (isNull)
                    return false;
                for (int i = 0; i < N; i++){
                    if (minCorner[i] > p[i] || maxCorner[i] < p[i])
                        return false;
                }
                return true;
            }
            inline bool contains(const Box & b) const {
                return b.isNull ? true : contains(b.minCorner) && contains(b.maxCorner);
            }
            inline bool operator == (const Box & b) const {
                return isNull ? b.isNull : (!b.isNull && minCorner == b.minCorner && maxCorner == b.maxCorner);
            }
            inline bool operator != (const Box & b) const { return !(*this == b); }

            inline Box & operator |= (const Box & b) {
                if (isNull){
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
        inline Box<T, N> operator | (const Box<T, N> & b1, const Box<T, N> & b2){
            Box<T, N> b12 = b1;
            return b12 |= b2;
        }
        using Box2 = Box<double, 2>;
        using Box3 = Box<double, 3>;

        // for scalars
        template <class T>
        inline std::enable_if_t<std::is_arithmetic<T>::value, Box<T, 1>> BoundingBox(const T & t) {
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
            return BoundingBox(hp.toPoint());
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



        // the standard bounding box functor
        template <class T>
        struct DefaultBoundingBoxFunctor {
            using BoxType = decltype(BoundingBox(std::declval<T>()));
            inline BoxType operator()(const T & t) const {
                return BoundingBox(t);
            }
        };


        // the standard influence box functor
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

            inline void insert(const T & t) {
                BoxType box = _getBoundingBox(t);
                for (int i = 0; i < Dimension; i++){
                    if (!(box.minCorner[i] <= box.maxCorner[i])){
                        std::cerr << "invalid box type, replaced with a null box" << std::endl;
                        box = BoxType();
                        //return;
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

        private:
            std::shared_ptr<third_party::RTree<T, ValueType, Dimension>> _rtree;
            BoundingBoxFunctorT _getBoundingBox;
        };




        // spheres
        template <class T, int N>
        struct Sphere {
            Point<T, N> center;
            T radius;
        };
        using Sphere2 = Sphere<double, 2>;
        using Sphere3 = Sphere<double, 3>;





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

        template <class T>
        T WrapBetween(const T & input, const T & low, const T & high) {
            if (low >= high)
                return input;
            if (low <= input && input < high)
                return input;
            const auto sz = high - low;
            auto result = input - int((input - low) / sz) * sz + (input < low ? sz : 0);
            return result == high ? low : result;
        }

        template <class T, int N>
        T EncodeSubscriptToIndex(const Point<T, N> & subscript, const Vec<T, N> & dimension) {
            T index = subscript[0];
            for (int i = 1; i < N; i++){
                index = index * dimension[i] + subscript[i];
            }
            return index;
        }

        template <class T, int N>
        Point<T, N> DecodeIndexToSubscript(T index, const Vec<T, N> & dimension) {
            Point<T, N> subscript;
            for (int i = N - 1; i >= 0; i--){
                subscript[i] = WrapBetween(index, 0, dimension[i]);
                index = (index - subscript[i]) / dimension[i];
            }
            return subscript;
        }







        // for vectors
        template <class T, int N>
        inline T AngleBetweenDirections(const Vec<T, N> & v1, const Vec<T, N> & v2) {
            return acos(v1.dot(v2) / norm(v1) / norm(v2));
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
            Vec<T, N> lineDir = line.second - line.first;
            lineDir /= norm(lineDir);
            T projRatio = (p - line.first).dot(lineDir) / line.length();
            return PositionOnLine<T, N>(line, projRatio);
        }

        // returns (distance, nearest position)
        template <class T, int N>
        std::pair<T, PositionOnLine<T, N>> DistanceFromPointToLine(const Point<T, N> & p, const Line<T, N> & line) {
            Vec<T, N> lineDir = line.second - line.first;
            lineDir /= norm(lineDir);
            T projRatio = (p - line.first).dot(lineDir) / line.length();
            projRatio = BoundBetween(projRatio, 0, 1);
            PositionOnLine<T, N> pos(line, projRatio);
            return std::make_pair(norm(p - pos.position), pos);
        }

        // returns (distance, (nearest position on line1, nearest position on line2))
        template <class T, int N>
        std::pair<T, std::pair<PositionOnLine<T, N>, PositionOnLine<T, N>>> DistanceBetweenTwoLines(
            const Line<T, N> & line1, const Line<T, N> & line2) {
            auto u = line1.direction();
            auto v = line2.direction();
            auto w0 = line1.first - line2.first;
            auto a = u.dot(u), b = u.dot(v), c = v.dot(v), d = u.dot(w0), e = v.dot(w0);
            auto p = (a*c - b*b);
            auto t1 = (b*e - c*d) / p;
            auto t2 = (a*e - b*d) / p;
            if (p == 0){ // is parallel
                t1 = -e / b;
                t2 = e / c;
            }
            t1 = BoundBetween(t1, 0, 1);
            t2 = BoundBetween(t2, 0, 1);
            PositionOnLine<T, N> pos1(line1, t1);
            PositionOnLine<T, N> pos2(line2, t2);
            return std::make_pair(norm(pos1.position - pos2.position), std::make_pair(pos1, pos2));
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
        class DistanceFunctorT = DefaultDistanceFunctor < typename std::iterator_traits<IteratorT>::value_type >>
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
        class DistanceFunctorT = DefaultDistanceFunctor < typename std::iterator_traits<IteratorT>::value_type >>
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
        class VertCompareT = std::less < typename std::iterator_traits<VertIteratorT>::value_type >>
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
        template < class VertIteratorT, class NeighborVertsRangeGetterT, class VertCallbackT,
        class VertCompareT = std::less < typename std::iterator_traits<VertIteratorT>::value_type >>
            void DepthFirstSearch(VertIteratorT vertsBegin, VertIteratorT vertsEnd,
            NeighborVertsRangeGetterT neighborVertsRangeGetter,
            VertCallbackT vertCallback,
            VertCompareT vertCompare = VertCompareT()) {

                using Vert = typename std::iterator_traits<typename VertIteratorT>::value_type;
                static_assert(std::is_same<Vert, std::decay_t<decltype(*std::get<0>(neighborVertsRangeGetter(std::declval<Vert>())))>>::value,
                    "NeighborVertsRangeGetterT should returns (VertIterator vbegin, VertIterator vend)");
                static_assert(std::is_same<Vert, std::decay_t<decltype(*std::get<1>(neighborVertsRangeGetter(std::declval<Vert>())))>>::value,
                    "NeighborVertsRangeGetterT should returns (VertIterator vbegin, VertIterator vend)");

                struct {
                    bool operator()(Vert root, std::map<Vert, bool, VertCompareT>& vVisited, NeighborVertsRangeGetterT vNeighborsGetter, VertCallbackT vCallback) {
                        if (vVisited[root])
                            return true;
                        if (!vCallback(root))
                            return false;

                        vVisited[root] = true;
                        auto vNeighborsRange = vNeighborsGetter(root);
                        for (auto nextv = std::get<0>(vNeighborsRange); nextv != std::get<1>(vNeighborsRange); ++nextv) {
                            //static_assert(typename std::is_same<Vert, decltype(*nextv)>>::value, "*nextv should be of type Vert");
                            if (!(*this)(*nextv, vVisited, vNeighborsGetter, vCallback))
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
                    if (!depthFirstSearchOneTree(*rootIter, visited, neighborVertsRangeGetter, vertCallback))
                        break;
                }
            }


        // Connected Components
        template < class VertIteratorT, class NeighborVertsRangeGetterT, class VertexTypeRecorderT,
        class VertCompareT = std::less < typename std::iterator_traits<VertIteratorT>::value_type >
        >
        int ConnectedComponents(VertIteratorT vertsBegin, VertIteratorT vertsEnd,
            NeighborVertsRangeGetterT neighborVertsRangeGetter,
            VertexTypeRecorderT vertTypeRecorder,
            VertCompareT vertCompare = VertCompareT()) {

            using Vert = typename std::iterator_traits<typename VertIteratorT>::value_type;
            static_assert(std::is_same<Vert, std::decay_t<decltype(*std::get<0>(neighborVertsRangeGetter(std::declval<Vert>())))>>::value,
                "NeighborVertsRangeGetterT should returns (VertIterator vbegin, VertIterator vend)");
            static_assert(std::is_same<Vert, std::decay_t<decltype(*std::get<1>(neighborVertsRangeGetter(std::declval<Vert>())))>>::value,
                "NeighborVertsRangeGetterT should returns (VertIterator vbegin, VertIterator vend)");

            struct {
                void operator()(Vert root, std::map<Vert, bool, VertCompareT>& vVisited, NeighborVertsRangeGetterT vNeighborsGetter, VertexTypeRecorderT vTypeRecorder, int cid) {
                    if (vVisited[root])
                        return;
                    vTypeRecorder(root, cid);
                    vVisited[root] = true;
                    auto vNeighborsRange = vNeighborsGetter(root);
                    for (auto nextv = std::get<0>(vNeighborsRange); nextv != std::get<1>(vNeighborsRange); ++nextv) {
                        //static_assert(typename std::is_same<Vert, decltype(*nextv)>>::value, "*nextv should be of type Vert");
                        (*this)(*nextv, vVisited, vNeighborsGetter, vTypeRecorder, cid);
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
                depthFirstSearchOneTree(*rootIter, visited, neighborVertsRangeGetter, vertTypeRecorder, cid);
                cid++;
            }

            return cid;
        }

    }
}
 
#endif