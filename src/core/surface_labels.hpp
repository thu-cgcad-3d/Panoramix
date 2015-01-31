#ifndef PANORAMIX_CORE_SURFACE_LABELS_HPP
#define PANORAMIX_CORE_SURFACE_LABELS_HPP

#include "feature.hpp"
#include "cameras.hpp"
 
namespace panoramix {
    namespace core { 

        enum class SurfaceLabelNames : int {
            None = -1,
            Vertical,
            Horizontal,
            Planar,
            NonPlanar,
            Void,
            Count
        };

        using SurfaceLabelDistribution = std::array<Imaged, (unsigned)SurfaceLabelNames::Count>;
        using SurfaceLabelPixelVec = Vec<double, (unsigned)SurfaceLabelNames::Count>;

        SurfaceLabelPixelVec MeanPixelVecInRegion(const SurfaceLabelDistribution & distribution,
            const std::vector<PixelLoc> & contour);
        SurfaceLabelPixelVec PixelVecAtPosition(const SurfaceLabelDistribution & distribution,
            const PixelLoc & p);

        template <class CameraT, class = std::enable_if_t<IsCamera<CameraT>::value>>
        class SurfaceLabels {
        public:
            inline SurfaceLabels() {}
            inline explicit SurfaceLabels(const CameraT cam) : _camera(cam) {
                for (auto & prob : _distributions){
                    prob = Imaged::zeros(_camera.screenSize());
                }
            }
            inline SurfaceLabels(const CameraT cam, const std::array<Imaged, (unsigned)SurfaceLabelNames::Count> & dists)
                : _camera(cam), _distributions(dists) {}

            const CameraT & camera() const { return _camera; }
            const SurfaceLabelDistribution & distributions() const { return _distributions; }
            SurfaceLabelDistribution & distributions() { return _distributions; }

            const Imaged & distribution(SurfaceLabelNames label) const { return _distributions[(int)(label)]; }
            Imaged & distribution(SurfaceLabelNames label) { return _distributions[(int)(label)]; }


            // pixel vector at ...
            SurfaceLabelPixelVec pixelVecAt(const Vec3 & dir) const {
                return PixelVecAtPosition(_distributions, ToPixelLoc(_camera.screenProjection(dir)));
            }

            SurfaceLabelPixelVec pixelVecAt(const std::vector<Vec3> & regionContour) const {
                std::vector<PixelLoc> contour(regionContour.size());
                for (int i = 0; i < regionContour.size(); i++){
                    contour[i] = ToPixelLoc(_camera.screenProjection(regionContour[i]));
                }
                return PixelVecAtPosition(_distributions, contour);
            }

            SurfaceLabelPixelVec pixelVecAt(const Line3 & lineProj, double sampleStepAngle = 0.005) const {
                double angle = AngleBetweenDirections(lineProj.first, lineProj.second);
                SurfaceLabelPixelVec pv;
                int count = 0;
                for (double a = 0.0; a <= angle; a += sampleStepAngle){
                    Vec3 dir = RotateDirection(lineProj.first, lineProj.second, a);
                    Point2 p = _camera.screenProjection(dir);
                    pv += PixelVecAtPosition(_distributions, ToPixelLoc(p));
                    count++;
                }
                return pv / (double)count;
            }

        public:
            template <class Archiver> void serialize(Archiver & ar) {
                ar(_camera, _distributions);
            }

        private:
            CameraT _camera;
            SurfaceLabelDistribution _distributions;
        };


        SurfaceLabels<PanoramicCamera> ComputeSurfaceLabels(const View<PanoramicCamera> & view,
            const std::vector<View<PerspectiveCamera, GeometricContextEstimator::Feature>> & perspectiveGCs);

        SurfaceLabels<PanoramicCamera> ComputeSurfaceLabels(const View<PanoramicCamera> & view, SceneClass sceneClaz,
            const GeometricContextEstimator & gcEstimator = GeometricContextEstimator(),
            int num = 16, int width = 400, int height = 700, double focal = 250.0);
        SurfaceLabels<PerspectiveCamera> ComputeSurfaceLabels(const View<PerspectiveCamera> & view, SceneClass sceneClaz,
            const GeometricContextEstimator & gcEstimator = GeometricContextEstimator());
        

    }
}
 
#endif