#ifndef PANORAMIX_REC_LINES_NET_HPP
#define PANORAMIX_REC_LINES_NET_HPP

#include "../core/graphical_model.hpp"
#include "../core/utilities.hpp"
#include "../core/feature.hpp"
 
namespace panoramix {
    namespace rec {

        using namespace core;

        // net of line segments
        class LinesNet {
        public:
            struct Params {
                Params();
                LineSegmentExtractor lineSegmentExtractor;
                double intersectionDistanceThreshold;
                double incidenceDistanceAlongDirectionThreshold; // distance along the direction with lines
                double incidenceDistanceVerticalDirectionThreshold; // distance along the vertical direction
                template <class Archiver>
                void serialize(Archiver & ar) {
                    ar(lineSegmentExtractor, intersectionDistanceThreshold,
                        incidenceDistanceAlongDirectionThreshold, incidenceDistanceVerticalDirectionThreshold);
                }
            };

            struct LineData {
                Classified<Line2> line;
                template <class Archiver>
                void serialize(Archiver & ar) {
                    ar(line);
                }
            };
            struct LineRelationData {
                Point2 relationCenter;
                float junctionWeight;
                enum Type {
                    Incidence,
                    Intersection
                } type;
                template <class Archiver>
                void serialize(Archiver & ar) {
                    ar(relationCenter, junctionWeight, type);
                }
            };
            using LinesGraph = GraphicalModel02 < LineData, LineRelationData > ;
            using LineHandle = HandleAtLevel < 0 > ;
            using LineRelationHandle = HandleAtLevel < 1 > ;
        
        public:
            inline LinesNet() {}
            explicit LinesNet(const Image & image, const Params & params = Params());
            void buildNetAndComputeFeaturesUsingVanishingPoints(const std::array<HPoint2, 3> & vps, 
                const std::vector<int> & lineSegmentClasses = std::vector<int>());

            inline const LineSegmentExtractor::Feature & lineSegments() const { return _lineSegments; }
            inline const std::vector<HPoint2> & lineSegmentIntersections() const { return _lineSegmentIntersections; }
            inline const LinesGraph & lines() const { return _lines; }
            inline const Image & image() const { return _image; }

            // 32FC(3x2)
            inline const ImageWithType<Mat<float, 3, 2>> & lineVotingDistribution() const { 
                return _lineVotingDistribution; 
            }

            inline const Params & params() const { return _params; }

            enum LineVotingDirection : int {
                TowardsVanishingPoint = 0,
                TowardsOppositeOfVanishingPoint = 1
            };
            enum class JunctionType : int {
                L, T, Y, W, X
            };

        private:
            Image _image;
            ImageWithType<Mat<float, 3, 2>> _lineVotingDistribution; // 32FC(3x2)
            std::map<JunctionType, ImageWithType<float>> _junctionDistributions;
            LinesGraph _lines;
            LineSegmentExtractor::Feature _lineSegments;
            std::vector<HPoint2> _lineSegmentIntersections;
            std::vector<std::pair<int, int>> _lineSegmentIntersectionIds;
            Params _params;

            template <class Archiver>
            void serialize(Archiver & ar) {
                ar(_image, _lineVotingDistribution, _junctionDistributions, 
                    _lines, _lineSegments, 
                    _lineSegmentIntersections, _lineSegmentIntersectionIds, _params);
            }
            friend class cereal::access;
        };
    
    }
}
 
#endif