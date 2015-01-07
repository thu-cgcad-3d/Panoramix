#ifndef PANORAMIX_CORE_VIEW_HPP
#define PANORAMIX_CORE_VIEW_HPP

#include "basic_types.hpp"
#include "graph.hpp"
#include "cameras.hpp"
 
namespace panoramix {
    namespace core {

        // view class
        template <class CameraT, class = std::enable_if_t<IsCamera<CameraT>::value>>
        struct View {
            Image image;
            CameraT camera;
            template <class Archiver>
            void serialize(Archiver & ar) {
                ar(image, camera);
            }
        };

        // create panoramic view
        View<PanoramicCamera> CreatePanoramicView(const Image & panorama);

        // view sampling
        template <class OriginalCameraT, class TargetCameraIteratorT, class ViewOutIteratorT>
        inline void SampleViews(const View<OriginalCameraT> & originalView,
            TargetCameraIteratorT camBegin, TargetCameraIteratorT camEnd, ViewOutIteratorT outputViewIter){
            using TargetCameraType = typename std::iterator_traits<TargetCameraIteratorT>::value_type;
            for (; camBegin != camEnd; ++camBegin, ++outputViewIter){
                CameraSampler<TargetCameraType, OriginalCameraT> camSampler(*camBegin, originalView.camera);
                *outputViewIter = View<TargetCameraType>{camSampler(originalView.image), *camBegin};
            }
        }


        // regions graph
        struct RegionData {
            Vec2 center;
            double area;
            std::vector<std::vector<PixelLoc>> contours; 
            std::vector<std::vector<PixelLoc>> dilatedContours;
            Box2 boundingBox;
            template <class Archiver>
            void serialize(Archiver & ar) {
                ar(center, area, contours, dilatedContours, boundingBox);
            }
        };
        struct RegionBoundaryData {
            std::vector<std::vector<PixelLoc>> edges;
            double length;
            InfiniteLine2 fittedLine;
            double tjunctionLikelihood;
            double straightness;
            std::vector<std::vector<Point2>> sampledPoints;
            template <class Archiver>
            void serialize(Archiver & ar) {
                ar(edges, length, fittedLine, tjunctionLikelihood, straightness, sampledPoints);
            }
        };
        using RegionsGraph = HomogeneousGraph02<RegionData, RegionBoundaryData>;
        using RegionHandle = HandleAtLevel<0>;
        using RegionBoundaryHandle = HandleAtLevel<1>;

        // create a compressed regions graph from segmented regions
        // ensurance: RegionHandle(i) always represents the region data for region mask: segmentedRegions == i
        RegionsGraph CreateRegionsGraph(const Imagei & segmentedRegions,
            double samplingStepLengthOnBoundary = 15.0, int dilationSize = 3);



        // junction weight 
        float IncidenceJunctionWeight(bool acrossViews);
        float OutsiderIntersectionJunctionWeight();
        float ComputeIntersectionJunctionWeightWithLinesVotes(const Mat<float, 3, 2> & votes);


        // lines graph
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
        using LinesGraph = HomogeneousGraph02<LineData, LineRelationData>;
        using LineHandle = HandleAtLevel<0>;
        using LineRelationHandle = HandleAtLevel<1>;


        // create lines graph from classified line segments
        // ignore claz == -1 line segments
        LinesGraph CreateLinesGraph(const std::vector<Classified<Line2>> & lineSegments,
            const std::vector<HPoint2> & vps,
            double intersectionDistanceThreshold = 8,
            double incidenceDistanceAlongDirectionThreshold = 15,
            double incidenceDistanceVerticalDirectionThreshold = 3);


        // recognize region-line connections
        // ensurance: RegionHandle(i) always represents the region data for region mask: segmentedRegions == i
        std::map<std::pair<RegionHandle, LineHandle>, std::vector<Point2>>
            RecognizeRegionLineConnections(const Imagei & segmentedRegions,
            const LinesGraph & lines, double samplingStepLengthOnLines);






        // component indices classes
        template <class HandleT>
        struct ComponentIndexInView {
            static_assert(IsHandle<HandleT>::value, "HandleT must be core::Handle<T>!");
            int viewId;
            HandleT handle;
            inline explicit ComponentIndexInView(int vid = -1, HandleT hd = HandleT()) 
                : viewId(vid), handle(hd) {}
            template <class Archiver>
            void serialize(Archiver & ar) {
                ar(viewId, handle);
            }
        };

        // helper functions
        template <class HandleT>
        inline bool operator == (const ComponentIndexInView<HandleT> & a, const ComponentIndexInView<HandleT> & b) {
            return a.viewId == b.viewId && a.handle == b.handle;
        }

        template <class HandleT>
        inline bool operator != (const ComponentIndexInView<HandleT> & a, const ComponentIndexInView<HandleT> & b) {
            return !(a == b);
        }

        template <class HandleT>
        inline bool operator < (const ComponentIndexInView<HandleT> & a, const ComponentIndexInView<HandleT> & b) {
            return std::tie(a.viewId, a.handle) < std::tie(b.viewId, b.handle);
        }

        template <class IndexT>
        struct ComponentIndexInViewHasher {
            inline size_t operator()(const IndexT & idx) const {
                return ((idx.viewId) << 4) + (idx.handle.id);
            }
        };

        template <class IndexT1, class IndexT2>
        struct ComponentIndexInViewHasher<std::pair<IndexT1, IndexT2>> {
            inline size_t operator()(const std::pair<IndexT1, IndexT2> & idx) const {
                return (((idx.first.viewId << 3) + idx.first.handle.id) << 10) +
                    (idx.second.viewId << 3) + idx.second.handle.id;
            }
        };

        // containers
        template <class T> using ComponentIndexHashSet = std::unordered_set<T, ComponentIndexInViewHasher<T>>;
        template <class KeyT, class ValueT> using ComponentIndexHashMap = std::unordered_map<KeyT, ValueT, ComponentIndexInViewHasher<KeyT>>;

        // instantial types
        using RegionIndex = ComponentIndexInView<RegionHandle>;
        using RegionBoundaryIndex = ComponentIndexInView<RegionBoundaryHandle>;
        using LineIndex = ComponentIndexInView<LineHandle>;
        using LineRelationIndex = ComponentIndexInView<LineRelationHandle>;


        // region/line data getter
        template <class VT, class ET, class HandleT>
        inline auto GetData(const ComponentIndexInView<HandleT> & i, const std::vector<HomogeneousGraph02<VT, ET>> & nets)
            -> decltype(std::declval<const HomogeneousGraph02<VT, ET>>().data(std::declval<HandleT>())) {
            return nets[i.viewId].data(i.handle);
        }

        template <class VT, class ET, class HandleT>
        inline auto GetData(const ComponentIndexInView<HandleT> & i, std::vector<HomogeneousGraph02<VT, ET>> & nets)
            -> decltype(std::declval<HomogeneousGraph02<VT, ET>>().data(std::declval<HandleT>())) {
            return nets[i.viewId].data(i.handle);
        }




        // estimate vanishing points and build lines graphs
        void EstimateVanishingPointsAndBuildLinesGraphs(const std::vector<View<PerspectiveCamera>> & views,
            std::vector<Vec3> & vanishingPoints,
            std::vector<LinesGraph> & linesGraphs,
            double intersectionDistanceThreshold,
            double incidenceDistanceAlongDirectionThreshold,
            double incidenceDistanceVerticalDirectionThreshold,
            bool considerNonManhattanVPs = false);


        // recognize region overlappings across views
        ComponentIndexHashMap<std::pair<RegionIndex, RegionIndex>, double>
            RecognizeRegionOverlappingsAcrossViews(const std::vector<View<PerspectiveCamera>> & views,
            const std::vector<RegionsGraph> & regionsGraphs);


        // recognize region overlappings across views
        ComponentIndexHashMap<std::pair<RegionIndex, RegionIndex>, double>
            RecognizeRegionOverlappingsAcrossViews(const std::vector<View<PanoramicCamera>> & views,
            const std::vector<RegionsGraph> & regionsGraphs);


        // recognize line incidences across views
        ComponentIndexHashMap<std::pair<LineIndex, LineIndex>, Vec3>
            RecognizeLineIncidencesAcrossViews(const std::vector<View<PerspectiveCamera>> & views,
            const std::vector<LinesGraph> & linesGraphs,
            double interViewIncidenceAngleAlongDirectionThreshold,
            double interViewIncidenceAngleVerticalDirectionThreshold);




      

    }
}
 
#endif