#pragma once



#include "../core/basic_types.hpp"
#include "../core/utility.hpp"
#include "../core/generic_topo.hpp"
#include "../core/cons_graph.hpp"
#include "../core/cameras.hpp"



namespace pano {
    namespace experimental {

        using namespace core;

        // the region-line graph definition
        struct RegionData;
        struct RegionBoundaryData;
        struct LineData;
        struct LineRelationData;
        struct RegionLineConnectionData;


        using RLGraph = ConstraintGraph<std::tuple<RegionData, LineData>,
            std::tuple<
            ConstraintConfig<RegionBoundaryData, RegionData, RegionData>,
            ConstraintConfig<LineRelationData, LineData, LineData>,
            ConstraintConfig<RegionLineConnectionData, RegionData, LineData>
            >
        >;   

        template <class T>
        using RLGraphComponentTable = typename ComponentHandledTableFromConstraintGraph<T, RLGraph>::type;
        template <class T>
        using RLGraphConstraintTable = typename ConstraintHandledTableFromConstraintGraph<T, RLGraph>::type;



        // regions
        struct RegionData {
            Vec3 normalizedCenter;
            double area; // projected onto the plane rooted at normalizedCenter
            std::vector<std::vector<Vec3>> normalizedContours;
            template <class Archiver>
            void serialize(Archiver & ar) {
                ar(normalizedCenter, area, normalizedContours);
            }
        };
        using RegionHandle = ComponentHandle<RegionData>;
        struct RegionBoundaryData {
            std::vector<std::vector<Vec3>> normalizedEdges;
            double length;
            std::vector<std::vector<Vec3>> normalizedSampledPoints;
            template <class Archiver>
            void serialize(Archiver & ar) {
                ar(normalizedEdges, length, normalizedSampledPoints);
            }
        };
        using RegionBoundaryHandle = ConstraintHandle<RegionBoundaryData>;



        // lines
        struct LineData {
            Line3 line;
            int initialClaz;
            template <class Archiver>
            void serialize(Archiver & ar) {
                ar(line, initialClaz);
            }
        };
        using LineHandle = ComponentHandle<LineData>;
        struct LineRelationData {
            Vec3 normalizedRelationCenter;
            float junctionWeight;
            enum Type {
                Incidence,
                Intersection
            } type;
            template <class Archiver>
            void serialize(Archiver & ar) {
                ar(normalizedRelationCenter, junctionWeight, type);
            }
        };
        using LineRelationHandle = ConstraintHandle<LineRelationData>;


        // region and line
        struct RegionLineConnectionData {
            std::vector<Vec3> normalizedAnchors;
            double length;
            bool detachable; // if the line lies on the edge of a region, it may be detachable from the region
            template <class Archiver>
            void serialize(Archiver & ar) {
                ar(normalizedAnchors, length, detachable);
            }
        };
        using RegionLineConnectionHandle = ConstraintHandle<RegionLineConnectionData>;




        // junction weight 
        // 7.0, 10.0
        float IncidenceJunctionWeight(bool acrossViews);
        float OutsiderIntersectionJunctionWeight();
        // [0.0 ~ 5.0]
        float ComputeIntersectionJunctionWeightWithLinesVotes(const Mat<float, 3, 2> & votes);




        // add lines to rl graph from classified line segments
        void AppendLines(RLGraph & mg, const std::vector<Classified<Line2>> & lineSegments,
            const PerspectiveCamera & cam,
            const std::vector<Vec3> & vps,
            double intersectionAngleThreshold = 0.04,
            double incidenceAngleAlongDirectionThreshold = 0.1,
            double incidenceAngleVerticalDirectionThreshold = 0.02,
            double interViewIncidenceAngleAlongDirectionThreshold = 0.15, // for new line-line incidence recognition
            double interViewIncidenceAngleVerticalDirectionThreshold = 0.03);

        std::vector<LineHandle> AppendLines2(RLGraph & mg, const std::vector<Classified<Line3>> & lineSegments,
            const std::vector<Vec3> & vps,
            double intersectionAngleThreshold = 0.04,
            double incidenceAngleAlongDirectionThreshold = M_PI_2,
            double incidenceAngleVerticalDirectionThreshold = 0.02);





        // add more regions and related constraints to rl graph
        // AppendRegions
        std::vector<RegionHandle> AppendRegions(RLGraph & mg, const Imagei & segmentedRegions, const PerspectiveCamera & cam,
            double samplingStepAngleOnBoundary, double samplingStepAngleOnLine,
            int samplerSizeOnBoundary = 3, int samplerSizeOnLine = 3, bool noBoundaryUnderLines = false);
        std::vector<RegionHandle> AppendRegions(RLGraph & mg, const Imagei & segmentedRegions, const PanoramicCamera & cam,
            double samplingStepAngleOnBoundary, double samplingStepAngleOnLine, 
            int samplerSizeOnBoundary = 3, int samplerSizeOnLine = 3, bool noBoundaryUnderLines = false);

        // AppendRegions
        std::pair<std::vector<RegionHandle>, std::vector<RegionBoundaryHandle>> 
            AppendRegions(RLGraph & mg, const Imagei & segmentedRegions, 
            const std::vector<std::vector<Pixel>> & bndpixels, 
            const std::vector<std::pair<int, int>> & bnd2segs,
            const PanoramicCamera & cam,
            double samplingStepAngleOnBoundary, double samplingStepAngleOnLine,
            int samplerSizeOnBoundary = 3, int samplerSizeOnLine = 3, bool noBoundaryUnderLines = false);

        std::pair<std::vector<RegionHandle>, std::vector<RegionBoundaryHandle>>
            AppendRegions2(RLGraph & mg, const Imagei & segmentedRegions,
            const std::vector<std::vector<Pixel>> & bndpixels,
            const std::vector<std::pair<int, int>> & bnd2segs,
            const PanoramicCamera & cam,
            double samplingStepAngleOnBoundary, double samplingStepAngleOnLine);


        // get a perfect mask view for a region
        View<PartialPanoramicCamera, Imageub> PerfectRegionMaskView(const RLGraph & mg, RegionHandle rh, double focal = 100.0);






        template <class T, int N, class CameraT>
        HandledTable<RegionHandle, Vec<T, N>> CollectFeatureMeanOnRegions(const RLGraph & mg,
            const CameraT & pcam, const ImageOf<Vec<T, N>> & feature) {
            HandledTable<RegionHandle, Vec<T, N>> featureMeanTable = mg.createComponentTable<RegionData, Vec<T, N>>();
            for (auto & r : mg.components<RegionData>()) {
                auto rh = r.topo.hd;
                auto regionMaskView = PerfectRegionMaskView(mg, rh);
                auto sampler = MakeCameraSampler(regionMaskView.camera, pcam);
                auto featureOnRegion = sampler(feature);
                int votes = 0;
                Vec<T, N> featureSum;
                for (auto it = regionMaskView.image.begin(); it != regionMaskView.image.end(); ++it) {
                    if (!*it) {
                        continue;
                    }
                    featureSum += featureOnRegion(it.pos());
                    votes += 1;
                }
                auto featureMean = featureSum / std::max(votes, 1);
                featureMeanTable[rh] = featureMean;
            }
            return featureMeanTable;
        }



        HandledTable<RegionBoundaryHandle, int> CollectOcclusionResponseOnBoundaries(const RLGraph & mg,
            const std::vector<Chain3> & occ, const PanoramicCamera & cam,
            double samplingStepAngleOnOcc = 1 / 500.0, int samplerSizeOnOcc = 1);

        HandledTable<RegionBoundaryHandle, int> CollectOcclusionResponseOnBoundaries(const RLGraph & mg,
            const std::vector<Scored<Chain3>> & occ, const PanoramicCamera & cam,
            double samplingStepAngleOnOcc = 1 / 500.0, int samplerSizeOnOcc = 1);

    }
}


