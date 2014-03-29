#ifndef PANORAMIX_CORE_UTILITIES_HPP
#define PANORAMIX_CORE_UTILITIES_HPP

#include "basic_types.hpp"
 
namespace panoramix {
    namespace core {

        // for numerics
        template <class T>
        inline T Square(const T & v) {
            return v * v;
        }

        template <class T, class K>
        inline T WrapBetween(const T& input, const K& low, const K& high) {
            if (low >= high)
                return input;
            if (low <= input && input < high)
                return input;
            const K sz = high - low;
            return input - int((input - low) / sz) * sz + (input < low ? sz : 0);
        }

        // for vectors
        template <class T, int N>
        inline T AngleBetweenDirections(const Vec<T, N> & v1, const Vec<T, N> & v2) {
            return acos(v1.dot(v2) / norm(v1) / norm(v2));
        }

        template <class T>
        inline T SignedAngleBetweenDirections(const Vec<T, 2> & from, const Vec<T, 2> & to, 
            bool closewiseIsPositive = true) {
            double angle = atan2(- from(0)*to(1) + to(0)*from(1), from(1)*to(1) + from(0)*to(0));
            return closewiseIsPositive ? angle : -angle;
        }

        // generic algorithms

        // merge, rearrange the input array
        // DistanceFunctorT(a, b) -> DistanceT : compute the distance from a to b
        // returns the begin iterators of merged groups
        template <class IteratorT, class DistanceT, class DistanceFunctorT = std::minus<DistanceT>>
        std::vector<IteratorT> MergeNear(IteratorT begin, IteratorT end, std::true_type,
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
        std::vector<IteratorT> MergeNear(IteratorT begin, IteratorT end, std::false_type,
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

        /// optimizations

        // graph topology
        template <class VertHandleIteratorT, class EdgeHandleIteratorT,
        class EdgeFromGetterT, class EdgeToGetterT, 
        class VertEdgeBeginGetterT, class VertEdgeEndGetterT
        >
        struct DirectedGraphTopology {
            using VertHandle = typename std::iterator_traits<VertHandleIteratorT>::value_type;
            using EdgeHandle = typename std::iterator_traits<EdgeHandleIteratorT>::value_type;
            
            VertHandleIteratorT vertsBegin, vertsEnd;
            EdgeHandleIteratorT edgesBegin, edgesEnd;
            EdgeFromGetterT fromOfEdge; // EdgeFromGetterT(EdgeHandle h) -> VertHandle
            EdgeToGetterT toOfEdge; 
            VertEdgeBeginGetterT edgesBeginOfVert; // VertEdgeBeginGetterT(VertHandle h) -> EdgeHandleIteratorT
            VertEdgeEndGetterT edgesEndOfVert; // VertEdgeEndGetterT(VertHandle h) -> EdgeHandleIteratorT
        };
        template <class VertHandleIteratorT, class EdgeHandleIteratorT,
        class EdgeFromGetterT, class EdgeToGetterT,
        class VertEdgeBeginGetterT, class VertEdgeEndGetterT
        >
        inline DirectedGraphTopology <VertHandleIteratorT, EdgeHandleIteratorT, 
                                        EdgeFromGetterT, EdgeToGetterT, 
                                        VertEdgeBeginGetterT, VertEdgeEndGetterT>
        MakeDirectedGraphTopology(VertHandleIteratorT vertsBegin, VertHandleIteratorT vertsEnd, 
            EdgeHandleIteratorT edgesBegin, EdgeHandleIteratorT edgesEnd,
            EdgeFromGetterT fromOfEdge, EdgeToGetterT toOfEdge, 
            VertEdgeBeginGetterT edgesBeginOfVert, VertEdgeEndGetterT edgesEndOfVert) {
            return { vertsBegin, vertsEnd, edgesBegin, edgesEnd, fromOfEdge, toOfEdge, edgesBeginOfVert, edgesEndOfVert };
        }


        // graph cut
        template <class VertHandleIteratorT, class EdgeHandleIteratorT, 
        class EdgeFromGetterT, class EdgeToGetterT, 
        class VertEdgeBeginGetterT, class VertEdgeEndGetterT,
        class VertEnergyGetterT, class EdgeEnergyGetterT,
        class VertLabelGetterSetterT
        >
        void GraphCut(
            const DirectedGraphTopology<VertHandleIteratorT, EdgeHandleIteratorT,
                EdgeFromGetterT, EdgeToGetterT, VertEdgeBeginGetterT, VertEdgeEndGetterT> & graphTopo,
            VertEnergyGetterT energyOfVert, EdgeEnergyGetterT energyOfEdge, VertLabelGetterSetterT labelOfVert){



        }

        


    }
}
 
#endif