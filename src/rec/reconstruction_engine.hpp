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

            struct ComponentData;
            struct ConstraintData;
            using ConstraintGraph = GraphicalModel02<ComponentData, ConstraintData>;
            using ComponentHandle = HandleAtLevel<0>;
            using ConstraintHandle = HandleAtLevel<1>;

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
            size_t updateConnections(ViewHandle h);

            // connect all current views using delauny triangulation 
            void updateConnections();

            // whether this view overlaps some existing views a lot, measured by the smallCameraAngleScalar parameter
            ViewHandle isTooCloseToAnyExistingView(ViewHandle h) const;
            
            // estimate vanishing points using lines extract from all views, classify this lines and lift them all to space
            void estimateVanishingPointsAndClassifyLines();

            // reconstruct regions
            void reconstructLinesAndFaces();
                
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
            struct ViewConnectionData {
                cv::detail::MatchesInfo matchInfo;
            };



            //// components and constraints
            // component data
            struct LineComponentData {
                OptimizibleExpression<double> etaExpr;
                ReconstructionEngine::ViewHandle viewHandle;
                LinesNet::LineHandle lineHandle;
                int connectedComponentId;
            };

            struct RegionComponentData {
                OptimizibleExpression<Eigen::Vector3d> thetaExpr;
                DisableableExpression<double> manhattanEnergyExpr;
                ReconstructionEngine::ViewHandle viewHandle;
                RegionsNet::RegionHandle regionHandle;
            };

            struct ComponentData {
                enum class Type {
                    UnInitialized,
                    Line,
                    Region
                };
                explicit ComponentData(Type t = Type::UnInitialized);
                Type type;
                LineComponentData asLine;
                RegionComponentData asRegion;
                DisableableExpression<double> reserveScaleEnergyExpr;
            };

            // constraint data
            struct RegionOverlapConstraintData {
                double overlapRatio;
            };

            struct RegionConnectivityConstraintData {
                ReconstructionEngine::ViewHandle viewHandle;
                RegionsNet::BoundaryHandle boundaryHandle;
            };

            struct LineInterViewIncidenceConstraintData {
                Vec3 relationCenter;
            };

            struct LineConnectivityConstraintData {
                ReconstructionEngine::ViewHandle viewHandle;
                LinesNet::LineRelationHandle lineRelationHandle;
            };

            struct RegionLineConnectivityConstraintData {
                std::vector<Vec3> sampledPoints;
            };

            struct ConstraintData {
                enum class Type {
                    UnInitialized,
                    RegionOverlap,
                    RegionConnectivity,
                    LineConnectivity,
                    LineInterViewIncidence,
                    RegionLineConnectivity
                };
                explicit ConstraintData(Type t = Type::UnInitialized);

                Type type;
                RegionOverlapConstraintData asRegionOverlap;
                RegionConnectivityConstraintData asRegionConnectivity;
                LineConnectivityConstraintData asLineConnectivity;
                LineInterViewIncidenceConstraintData asLineInterViewIncidence;
                RegionLineConnectivityConstraintData asRegionLineConnectivity;

                DisableableExpression<double> constraintEnergyExpr;
                deriv::Expression<double> invalidityExpr; // used to select necessary constraints
            };


            // global data
            struct GlobalData {
                Image panorama;
                std::array<Vec3, 3> vanishingPoints;
            };

            inline const ViewsGraph & views() const { return _views; }
            inline const ConstraintGraph & constraints() const { return _constraints; }
            inline const GlobalData & globalData() const { return _globalData; }

        private:
            ViewsGraph _views;
            ConstraintGraph _constraints;
            Params _params;
            GlobalData _globalData;
        };

 
    }
}
 
#endif