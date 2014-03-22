#ifndef PANORAMIX_CORE_UTILITIES_HPP
#define PANORAMIX_CORE_UTILITIES_HPP

#include "basic_types.hpp"
 
namespace panoramix {
    namespace core {

        template <class T>
        inline T Square(const T & v) {
            return v * v;
        }

        template <class Vec3T1, class Vec3T2, class ValueT>
        inline void AngleBetweenDirections(const Vec3T1 & v1, const Vec3T2 & v2, ValueT & angle) {
            angle = acos(v1.dot(v2) / v1.norm() / v2.norm());
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
 
        template <class Mat4T, class Vec3T>
        Mat4T Matrix4MakeLookAt(const Vec3T & eye, const Vec3T & center, 
            const Vec3T & up, const Mat4T & base) {
            Vec3T zaxis = (center - eye).normalized();
            Vec3T xaxis = up.cross(zaxis).normalized();
            Vec3T yaxis = zaxis.cross(xaxis);
            Mat4T m;
            m <<
                xaxis(0), yaxis(0), zaxis(0), 0,
                xaxis(1), yaxis(1), zaxis(1), 0,
                xaxis(2), yaxis(2), zaxis(2), 0,
                -xaxis.dot(eye), -yaxis.dot(eye), -zaxis.dot(eye), 1;
            return m.transpose() * base;
        }

        template <class Mat4T, class ValueT>
        Mat4T Matrix4MakePerspective(const ValueT & fovyRadians, const ValueT & aspect, 
            const ValueT & nearZ, const ValueT & farZ, const Mat4T & base) {
            ValueT cotan = ValueT(1.0) / std::tan(fovyRadians / 2.0);
            Mat4T m;
            m <<
                cotan / aspect, 0, 0, 0,
                0, cotan, 0, 0,
                0, 0, (farZ + nearZ) / (nearZ - farZ), -1,
                0, 0, (2 * farZ * nearZ) / (nearZ - farZ), 0;
            return m.transpose() * base;
        }

        // merge, rearrange the input array
        // DistanceFunctorT(a, b) -> DistanceT
        // returns the begin iterators of each group
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
        // DistanceFunctorT(a, b) -> DistanceT
        // returns the begin iterators of each group
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


    }
}
 
#endif