#ifndef PANORAMIX_CORE_UTILITIES_HPP
#define PANORAMIX_CORE_UTILITIES_HPP

#include <iterator>
#include "basic_types.hpp"
 
namespace panoramix {
    namespace core {
        

        // for numerics
        template <class T>
        inline bool FuzzyEquals(const T & a, const T & b, const T & epsilon){
            return std::abs(a - b) <= epsilon;
        }

        template <class T>
        inline int DiracDelta(const T & v) {
            return v == 0 ? 1 : 0;
        }

        template <class T>
        inline T Square(const T & v) {
            return v * v;
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
            for (int i = N-1; i >=0; i--){
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
            double angle = atan2(- from(0)*to(1) + to(0)*from(1), from(1)*to(1) + from(0)*to(0));
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












        // generic algorithms

        // fill the container with linear sequence
        template <class IteratorT>
        void CreateLinearSequence(IteratorT begin, IteratorT end, 
            const typename std::iterator_traits<IteratorT>::value_type& low, 
            const typename std::iterator_traits<IteratorT>::value_type& high){
            auto dist = std::distance(begin, end);
            auto w = high - low;
            int id = 0;
            while(begin != end){
                *begin = (id++) * w / dist + low;
                ++begin;
            }
        }        


        // merge, rearrange the input array
        // DistanceFunctorT(a, b) -> DistanceT : compute the distance from a to b
        // returns the begin iterators of merged groups
        template <class IteratorT, class DistanceT, class DistanceFunctorT = std::minus<DistanceT>>
        std::vector<IteratorT> NaiveMergeNear(IteratorT begin, IteratorT end, std::true_type,
            DistanceT thres, DistanceFunctorT distFun = DistanceFunctorT()) {
            if (begin == end)
                return std::vector<IteratorT>();
            
            std::vector<IteratorT> gBegins(1, begin);
            for (auto i = std::next(begin); i != end; ++i){
                DistanceT minDist = std::numeric_limits<DistanceT>::max();
                auto nearestGBeginIter = gBegins.end();
                for (auto giter = gBegins.begin(); giter != gBegins.end(); ++ giter) {
                    auto gBegin = *giter;
                    DistanceT dist = std::abs(distFun(*gBegin, *i));
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
                } else { // add new group
                    gBegins.push_back(i);
                }
            }

            return gBegins;
        }

        // merge, without rearrangement
        // DistanceFunctorT(a, b) -> DistanceT : compute the distance from a to b
        // returns the iterators pointing to group leaders
        template <class IteratorT, class DistanceT, class DistanceFunctorT = std::minus<DistanceT>>
        std::vector<IteratorT> NaiveMergeNear(IteratorT begin, IteratorT end, std::false_type,
            DistanceT thres, DistanceFunctorT distFun = DistanceFunctorT()) {
            if (begin == end)
                return std::vector<IteratorT>();

            std::vector<IteratorT> gBegins(1, begin);
            for (auto i = std::next(begin); i != end; ++i){
                auto giter = gBegins.begin();
                for (; giter != gBegins.end(); ++giter) {
                    auto gBegin = *giter;
                    DistanceT dist = std::abs(distFun(*gBegin, *i));
                    if (dist <= thres){
                        break;
                    }
                }
                if (giter == gBegins.end()){ // add new group
                    gBegins.push_back(i);
                }
            }

            return gBegins;
        }


        
        // Minimum Spanning Tree
        // EdgeVertsGetterT(Edge e)->std::pair<Vert,Vert>
        //      bool operator==(Vert,Vert) must be available
        // EdgeWeightGetterT(Edge e)->Scalar
        template <class VertIteratorT, class EdgeIteratorT, 
        class EdgeVertsGetterT,
        class EdgeOutputIteratorT, 
        class EdgeCompareOnWeightT,
        class VertCompareT = std::less<typename std::iterator_traits<VertIteratorT>::value_type>>
        void MinimumSpanningTree(
            VertIteratorT vertsBegin, VertIteratorT vertsEnd,
            EdgeIteratorT edgesBegin, EdgeIteratorT edgesEnd,
            EdgeOutputIteratorT MSTedges,
            EdgeVertsGetterT vertsGetter, 
            EdgeCompareOnWeightT edgeCompareOnWeight = EdgeCompareOnWeightT(),
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

        
        


    }
}
 
#endif