#ifndef PANORAMIX_CORE_RECONSTRUCTION_ENGINE_HPP
#define PANORAMIX_CORE_RECONSTRUCTION_ENGINE_HPP

#include "../core/basic_types.hpp"
#include "../core/feature.hpp"
#include "../core/utilities.hpp"

#include "optimization.hpp"

#include "regions_net.hpp"
#include "lines_net.hpp"

namespace panoramix {
    namespace rec { 

        using namespace core;
       
        // engine
        // ReconstructionEngine
        class ReconstructionEngine {
        public:
            struct Params {
                Params();
                PanoramicCamera camera; // camera for generating the panorama

                // angle scalar to judge whether two views may share certain common features
                double cameraAngleScaler; 
                // angle scalar to judge whether two views are too close
                double smallCameraAngleScalar; 

                double samplingStepLengthOnRegionBoundaries;
                double samplingStepLengthOnLines;

                // for line connections
                double intersectionDistanceThreshold;
                double incidenceDistanceAlongDirectionThreshold;
                double incidenceDistanceVerticalDirectionThreshold;
            };

            struct ViewData;
            struct ViewConnectionData;
            using ViewsGraph = GraphicalModel02<ViewData, ViewConnectionData>;
            using ViewHandle = HandleAtLevel<0>;
            using ViewConnectionHandle = HandleAtLevel<1>;

        private:            
            template <class HandleT>
            struct IndexOfSubStructureInView {
                ViewHandle viewHandle;
                HandleT handle;
            };

            template <class IndexT>
            struct IndexOfSubStructureInViewHasher {
                inline size_t operator()(const IndexT & idx) const {
                    return ((idx.viewHandle.id) << 4) + (idx.handle.id);
                }
            };

            template <class IndexT1, class IndexT2>
            struct IndexOfSubStructureInViewHasher < std::pair<IndexT1, IndexT2> > {
                inline size_t operator()(const std::pair<IndexT1, IndexT2> & idx) const {
                    return (((idx.first.viewHandle.id << 3) + idx.first.handle.id) << 10) +
                        (idx.second.viewHandle.id << 3) + idx.second.handle.id;
                }
            };

            using RegionIndex = IndexOfSubStructureInView <RegionsNet::RegionHandle>;
            using RegionBoundaryIndex = IndexOfSubStructureInView < RegionsNet::BoundaryHandle >;
            using LineIndex = IndexOfSubStructureInView <LinesNet::LineHandle>;
            using LineRelationIndex = IndexOfSubStructureInView < LinesNet::LineRelationHandle >;

            inline const RegionsNet::RegionData & regionData(const RegionIndex & ri) const { 
                return _views.data(ri.viewHandle).regionNet->regions().data(ri.handle); 
            }
            inline const RegionsNet::BoundaryData & regionBoundaryData(const RegionBoundaryIndex & rbi) const {
                return _views.data(rbi.viewHandle).regionNet->regions().data(rbi.handle);
            }

            inline const LinesNet::LineData & lineData(const LineIndex & li) const {
                return _views.data(li.viewHandle).lineNet->lines().data(li.handle);
            }
            inline const LinesNet::LineRelationData & lineRelationData(const LineRelationIndex & lri) const {
                return _views.data(lri.viewHandle).lineNet->lines().data(lri.handle);
            }

            template <class T>
            using IndexHashSet = std::unordered_set<T, IndexOfSubStructureInViewHasher<T>>;

            template <class KeyT, class ValueT>
            using IndexHashMap = std::unordered_map<KeyT, ValueT, IndexOfSubStructureInViewHasher<KeyT>>;

            template <class HandleT>
            friend bool operator == (const IndexOfSubStructureInView<HandleT> & a, 
                const IndexOfSubStructureInView<HandleT> & b);

            template <class HandleT>
            friend bool operator < (const IndexOfSubStructureInView<HandleT> & a,
                const IndexOfSubStructureInView<HandleT> & b);

        public:
            inline explicit ReconstructionEngine(const Params params = Params()) : _params(params) {}
            inline const Params & params() const { return _params; }
            
            inline ViewHandle insertView(const ViewData & vd) { return _views.add(vd); }

            // insert a new photo, with known parameters
            ViewHandle insertPhoto(const Image & im, const PerspectiveCamera & cam, 
                double cameraDirectionErrorScale = 0.0);

            // insert a new panoramic image and connect close views and set panorama
            void insertPanorama(const Image & panorama, const std::vector<PerspectiveCamera> & viewCams, 
                const PanoramicCamera & panCam);

            // compute features for a single view
            void computeFeatures(ViewHandle h);

            // connect this view to neighbor views who may overlap with h
            // we don't use delaunay triangulation here because overlapping view connections may be 'crossed' 
            size_t updateConnections(ViewHandle h);

            // whether this view overlaps some existing views a lot, measured by the smallCameraAngleScalar parameter
            ViewHandle isTooCloseToAnyExistingView(ViewHandle h) const;
            
            // estimate vanishing points using lines extract from all views, classify this lines and lift them all to space
            void estimateVanishingPointsAndClassifyLines();

            // construct region-line net
            void recognizeRegionLineRelations();



            // estimate spatial line depths
            void estimateSpatialLineDepths();

            // classify regions using classfied line labels
            void initializeRegionOrientations();

                
        public:

            //// views
            // view data
            struct ViewData {
                // cameras
                PerspectiveCamera originalCamera, camera;
                double cameraDirectionErrorScale;
                
                // image and 2d image features
                Image image;
                std::shared_ptr<RegionsNet> regionNet;
                std::shared_ptr<LinesNet> lineNet;
            };

            // view connection data
            struct ViewConnectionData { };

            // global data
            struct GlobalData {
                Image panorama;
                std::array<Vec3, 3> vanishingPoints;
                IndexHashMap<std::pair<RegionIndex, RegionIndex>, double> overlappedRegionIndexPairs;
                IndexHashMap<std::pair<LineIndex, LineIndex>, Vec3> lineIncidenceRelationsAcrossViews;
                IndexHashMap<std::pair<RegionIndex, LineIndex>, std::vector<Vec3>> regionLineIntersectionSampledPoints;
            };

            inline const ViewsGraph & views() const { return _views; }
            inline const GlobalData & globalData() const { return _globalData; }

        private:
            ViewsGraph _views;
            Params _params;
            GlobalData _globalData;
        };




        template <class HandleT>
        inline bool operator == (const ReconstructionEngine::IndexOfSubStructureInView<HandleT> & a,
            const ReconstructionEngine::IndexOfSubStructureInView<HandleT> & b) {
            return a.viewHandle == b.viewHandle && a.handle == b.handle;
        }

        template <class HandleT>
        inline bool operator < (const ReconstructionEngine::IndexOfSubStructureInView<HandleT> & a,
            const ReconstructionEngine::IndexOfSubStructureInView<HandleT> & b) {
            if (a.viewHandle.id != b.viewHandle.id)
                return a.viewHandle.id < b.viewHandle.id;
            return a.handle.id < b.handle.id;
        }

 
    }
}
 
#endif