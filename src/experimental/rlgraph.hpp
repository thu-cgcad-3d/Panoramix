#ifndef PANORAMIX_EXPERIMENTAL_RL_GRAPH_HPP
#define PANORAMIX_EXPERIMENTAL_RL_GRAPH_HPP

#include "../core/basic_types.hpp"
#include "../core/cons_graph.hpp"

namespace panoramix {
    namespace experimental {

        using namespace panoramix::core;

        // rl graph
        struct RData;
        struct RRData;
        struct LData;
        struct LLData;
        struct RLData;
        struct RRLData;
        struct RRRData;
        struct RRRRData;

        using RLGraph = ConstraintGraph<std::tuple<RData, LData>,
            std::tuple<
            ConstraintConfig<RRData, RData, RData>,
            ConstraintConfig<LLData, LData, LData>,
            ConstraintConfig<RLData, RData, LData>,
            ConstraintConfig<RRLData, RData, RData, LData>,
            ConstraintConfig<RRRData, RData, RData, RData>,
            ConstraintConfig<RRRRData, RData, RData, RData, RData>
            >
        >;

        using RHandle = ComponentHandle<RData>;
        using RRHandle = ConstraintHandle<RRData>;
        using LHandle = ComponentHandle<LData>;
        using LLHandle = ConstraintHandle<LLData>;
        using RLHandle = ConstraintHandle<RLData>;
        using RRLHandle = ConstraintHandle<RRLData>;
        using RRRHandle = ConstraintHandle<RRRData>;
        using RRRRHandle = ConstraintHandle<RRRRData>;

        struct RData {
            Vec3 normalizedCenter;
            double area;
            std::vector<std::vector<Vec3>> normalizedContours;
            Vec5 gcMean;
            template <class Archiver>
            void serialize(Archiver & ar) {
                ar(normalizedCenter, area, normalizedContours, gcMean);
            }
        };

        struct RRData {
            std::vector<std::vector<Vec3>> normalizedEdges;
            double length;
            std::vector<std::vector<Vec3>> normalizedSampledPoints;
            enum BoundaryType {Connected = 0, NotConnected = 1} occDetectionResult;
            template <class Archiver>
            void serialize(Archiver & ar) {
                ar(normalizedEdges, length, normalizedSampledPoints, occDetectionResult);
            }
        };

        struct LData {
            Line3 normalizedLine;
            std::vector<double> vpScores;
            int initialClaz;
            template <class Archiver>
            void serialize(Archiver & ar) {
                ar(normalizedLine, vpScores, initialClaz);
            }
        };

        struct LLData {
            enum ManhattanJunctionType {  
                /*intersection:*/ Y = 0, W, L, K, X, T, Hex, Outsided,  
                /*incidence:*/ I, IAcrossView,
                ManhattanJunctionTypeNum
            } mjType;
            Vec3 normalizedRelationCenter;
            template <class Archiver>
            void serialize(Archiver & ar) {
                ar(mjType, normalizedRelationCenter);
            }
        };

        namespace {
            // connection data
            struct RL_RRL_DataBase {
                std::vector<Vec3> normalizedAnchors;
                double length;
                template <class Archiver>
                void serialize(Archiver & ar) {
                    ar(normalizedAnchors, length);
                }
            };
        }
        struct RLData : RL_RRL_DataBase {};
        struct RRLData : RL_RRL_DataBase {};

        namespace {
            // boundary junction data
            struct RRR_RRRR_DataBase { 
                Vec3 normalizedCenter;
                template <class Archiver>
                void serialize(Archiver & ar) {
                    ar(normalizedCenter);
                }
            };
        }
        struct RRRData : RRR_RRRR_DataBase {};
        struct RRRRData : RRR_RRRR_DataBase {};


        
        


        // graph builder
        class RLGraphBuilder {
        public:
            struct Params {
                double intersectionAngleThresholdForLL;
                double incidenceParaAngleThresholdForLL;
                double incidenceVertAngleThresholdForLL;
                double incidenceParaAngleThresholdForLLAcrossViews;
                double incidenceVertAngleThresholdForLLAcrossViews;
                double samplingStepAngleOnRR;
                double samplingStepAngleOnRL;
                double samplingStepAngleOnOcc;
                int samplingStepSizeOnRR;
                int samplingStepSizeOnRL;
                int samplingStepSizeOnOcc;
                
                Params();

                template <class Archiver>
                void serialize(Archiver & ar) {
                    ar(intersectionAngleThresholdForLL,
                        incidenceParaAngleThresholdForLL, incidenceVertAngleThresholdForLL,
                        incidenceParaAngleThresholdForLLAcrossViews, incidenceVertAngleThresholdForLLAcrossViews,
                        samplingStepAngleOnRR, samplingStepAngleOnRL, samplingStepAngleOnOcc,
                        samplingStepSizeOnRR, samplingStepSizeOnRL, samplingStepSizeOnOcc);
                }
            };

        public:
            explicit RLGraphBuilder(const Params & params = Params()) : _params(params) {}
            const Params & params() const { return _params; }
            Params & params() { return _params; }

            RLGraph operator()(const std::vector<Vec3> & vps, 
                const std::vector<std::vector<Classified<Line2>>> & lines, 
                const std::vector<DenseMatd> & lineVPScores,
                const std::vector<PerspectiveCamera> & cams,
                const Imagei & segmentedRegions, 
                const Image5d & gc,
                const std::vector<Chain3> & occ,
                const PerspectiveCamera & cam, std::vector<RHandle> * segId2Rhs = nullptr) const;

            RLGraph operator()(const std::vector<Vec3> & vps, 
                const std::vector<std::vector<Classified<Line2>>> & lines, 
                const std::vector<DenseMatd> & lineVPScores, 
                const std::vector<PerspectiveCamera> & cams,
                const Imagei & segmentedRegions, 
                const Image5d & gc, 
                const std::vector<Chain3> & occ,
                const PanoramicCamera & cam, std::vector<RHandle> * segId2Rhs = nullptr) const;

            template <class Archiver>
            void serialize(Archiver & ar) {
                ar(_params);
            }

        private:
            Params _params;
        };






    }

}

#endif