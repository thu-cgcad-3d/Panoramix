#ifndef PANORAMIX_CORE_P3GRAPH_HPP
#define PANORAMIX_CORE_P3GRAPH_HPP
 
#include "view.hpp"

// perspective projection planer recovery graph

namespace panoramix {
    namespace core {


        struct P3Unary;
        struct P3Binary;

        // environment configurations
        struct P3Environment {
            std::vector<Vec3> vanishingPoints;
            template <class Archive> void serialize(Archive & ar) {
                ar(vanishingPoints);
            }
        };


        // unary variable
        class P3UnaryVariable {
        public:
            double norm(int k = 2) const;
            
            template <class T> T interpretAs(const P3Unary & unary, const P3Environment & env) const { NOT_IMPLEMENTED_YET(); }
            template <> Plane3 interpretAs<Plane3>(const P3Unary & unary, const P3Environment & env) const;
            template <> Line3 interpretAs<Line3>(const P3Unary & unary, const P3Environment & env) const;
            template <> Sphere3 interpretAs<Sphere3>(const P3Unary & unary, const P3Environment & env) const;

            std::vector<double> variableCoeffsForInverseDepthAtDirection(const Vec3 & direction,
                const P3Unary & unary, const P3Environment & env) const;
            double inverseDepthAtDirection(const Vec3 & direction, 
                const P3Unary & unary, const P3Environment & env) const;
            double depthAtCenter(const P3Unary & unary, const P3Environment & env) const;
                        
            template <class Archive> void serialize(Archive & ar) {ar(_variables);}

        private:
            // (a, b, c) for region plane ax+by+c=1,
            // 1/centerDepth for classified line, 
            // (1/firstCornerDepth, 1/secondCornerDepth) for unclassified line
            std::vector<double> _variables;
        };


        // binary variable
        class P3BinaryVariable {
        public:
            
            
            

            template <class Archive> void serialize(Archive & ar) {
                ar(_linearCoefficients, _quadraticCoefficients);
            }

        private:
            std::vector<double> _linearCoefficients;
            Imaged _quadraticCoefficients;
        };


        // unary entity
        struct MGUnary {
            enum Type { Region, Line } type;
            std::vector<Vec3> normalizedCorners;
            Vec3 normalizedCenter;
            int orientationClaz; // -1 means undetermined

            template <class Archive> void serialize(Archive & ar) {
                ar(type, normalizedCorners, normalizedCenter, orientationClaz);
            }
        };




        // binary entity
        struct MGBinary {
            enum Type {
                RegionRegionConnection,
                RegionRegionOverlapping,
                RegionLineConnection,
                LineLineIntersection,
                LineLineIncidence
            } type;
            double weight;
            std::vector<Vec3> normalizedAnchors;
            std::array<double, 2> importanceRatioInRelatedUnaries;
            template <class Archive> void serialize(Archive & ar) {
                ar(type, weight, normalizedAnchors, importanceRatioInRelatedUnaries);
            }
        };

    }
}
 
#endif