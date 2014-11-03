#ifndef PANORAMIX_REC_RECONSTRUCTION_HPP
#define PANORAMIX_REC_RECONSTRUCTION_HPP

#include "../core/basic_types.hpp"
#include "../core/cameras.hpp"

#include "regions_net.hpp"
#include "lines_net.hpp"

namespace panoramix {
    namespace rec {

        using namespace core;

        // view class
        template <class CameraT>
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

        // perspective sampling
        std::vector<View<PerspectiveCamera>> PerspectiveSampling(const View<PanoramicCamera> & panoView,
            const std::vector<PerspectiveCamera> & cameras);

        // initialize regions net and lines net
        std::pair<RegionsNet, LinesNet> InitializeFeatureNets(const View<PerspectiveCamera> & view,
            double samplingStepLengthOnRegionBoundaries,
            double intersectionDistanceThreshold, 
            double incidenceDistanceVerticalDirectionThreshold, 
            double incidenceDistanceAlongDirectionThreshold);

        // estimate vanishing points and classify lines
        std::array<Vec3, 3> EstimateVanishingPointsAndClassifyLines(const std::vector<View<PerspectiveCamera>> & views,
            std::vector<LinesNet> & linesNets);



        // component indices classes
        template <class HandleT>
        struct ComponentIndexInView {
            int viewId;
            HandleT handle;
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
            if (a.viewId != b.viewId)
                return a.viewId < b.viewId;
            return a.handle.id < b.handle.id;
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
        using RegionIndex = ComponentIndexInView<RegionsNet::RegionHandle>;
        using RegionBoundaryIndex = ComponentIndexInView<RegionsNet::BoundaryHandle>;
        using LineIndex = ComponentIndexInView<LinesNet::LineHandle>;
        using LineRelationIndex = ComponentIndexInView<LinesNet::LineRelationHandle>;


        void RecognizeRegionLineConstraints(const std::vector<View<PerspectiveCamera>> & views,
            const std::vector<RegionsNet> & regionsNets, const std::vector<LinesNet> & linesNets,
            ComponentIndexHashMap<std::pair<RegionIndex, RegionIndex>, double> & regionOverlappings,
            ComponentIndexHashMap<std::pair<RegionIndex, LineIndex>, std::vector<Vec3>> & regionLineConnections,
            ComponentIndexHashMap<std::pair<LineIndex, LineIndex>, Vec3> & interViewLineIncidences,
            double interViewIncidenceAngleAlongDirectionThreshold,
            double samplingStepLengthOnLines);


        void ComputeConnectedComponentsUsingRegionLineConstraints(const std::vector<View<PerspectiveCamera>> & views,
            const std::vector<RegionsNet> & regionsNets, const std::vector<LinesNet> & linesNets,
            const ComponentIndexHashMap<std::pair<RegionIndex, RegionIndex>, double> & regionOverlappings,
            const ComponentIndexHashMap<std::pair<RegionIndex, LineIndex>, std::vector<Vec3>> & regionLineConnections,
            const ComponentIndexHashMap<std::pair<LineIndex, LineIndex>, Vec3> & interViewLineIncidences,
            int & regionConnectedComponentsNum, ComponentIndexHashMap<RegionIndex, int> & regionConnectedComponentIds,
            int & lineConnectedComponentsNum, ComponentIndexHashMap<LineIndex, int> & lineConnectedComponentIds);


        void EstimateSpatialLineDepths(const std::vector<View<PerspectiveCamera>> & views,
            const std::vector<LinesNet> & linesNets, 
            const std::array<Vec3, 3> & vanishingPoints,
            const ComponentIndexHashMap<std::pair<LineIndex, LineIndex>, Vec3> & interViewLineIncidences, 
            int lineConnectedComponentsNum, const ComponentIndexHashMap<LineIndex, int> & lineConnectedComponentIds,
            ComponentIndexHashMap<LineIndex, Line3> & reconstructedLines,
            double constantEtaForFirstLineInEachConnectedComponent);


        


    }
}
 
#endif