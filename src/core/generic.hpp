#ifndef PANORAMIX_CORE_GENERIC_HPP
#define PANORAMIX_CORE_GENERIC_HPP

#include "basic_types.hpp"

namespace panoramix {
    namespace core {

        /// distance functions
        template <class T>
        inline std::enable_if_t<std::is_arithmetic<T>::value, T> 
            Distance(const T & a, const T & b) {
            return std::abs(a - b);
        }

        template <class T, int N>
        inline T Distance(const Point<T, N> & a, const Point<T, N> & b) {
            return norm(a - b);
        }

        template <class T, int N>
        inline T Distance(const HPoint<T, N> & a, const HPoint<T, N> & b) {
            return norm(a.toPoint() - b.toPoint());
        }

        // distance functor
        template <class T>
        struct DistanceFunctor {
            using DistanceType = decltype(Distance(T(), T()));
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
        inline std::enable_if_t<std::is_arithmetic<T>::value, Box<T, 1>> 
            BoundingBox(const T & t) {
            return Box<T, 1>(Point<T, 1>(t), Point<T, 1>(t));
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
       


        // bounding box functor
        template <class T>
        struct BoundingBoxFunctor {
            using BoxType = decltype(BoundingBox(T()));
            inline BoxType operator()(const T & t) const {
                return BoundingBox(t);
            }
        };


        // influence box functor
        template <class T>
        struct InfluenceBoxFunctor {
            using BoxType = decltype(BoundingBox(T()));
            using ValueType = typename BoxType::Type;

            inline explicit InfluenceBoxFunctor(const ValueType & extSz = 0) : extendedSize(extSz){}
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


        // spheres
        template <class T, int N>
        struct Sphere {
            Point<T, N> center;
            T radius;
        };
        using Sphere2 = Sphere<double, 2>;
        using Sphere3 = Sphere<double, 3>;





    }
}
 
#endif