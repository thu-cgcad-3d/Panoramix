#ifndef PANORAMIX_CORE_SURFACE_LABELS_HPP
#define PANORAMIX_CORE_SURFACE_LABELS_HPP

#include "feature.hpp"
#include "view.hpp"
 
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

        SurfaceLabelNames MostLikelySurfaceLabelAtPosition(const SurfaceLabelDistribution & distribution, 
            const PixelLoc & p);
        SurfaceLabelNames MostLikelySurfaceLabelInRegion(const SurfaceLabelDistribution & distribution, 
            const std::vector<PixelLoc> & contour);


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

            SurfaceLabelNames mostLikelyLabelAt(const Vec3 & dir) const {
                return MostLikelyLabelAtPosition(_distributions, ToPixelLoc(_camera.screenProjection(dir)));
            }

            SurfaceLabelNames mostLikelyLabelAt(const std::vector<Vec3> & regionContour) const {
                std::vector<PixelLoc> contour(regionContour.size());
                for (int i = 0; i < regionContour.size(); i++){
                    contour[i] = ToPixelLoc(_camera.screenProjection(regionContour[i]));
                }
                return MostLikelySurfaceLabelInRegion(_distributions, contour);
            }

            SurfaceLabelNames mostLikelyLabelAt(const Line3 & lineProj, double thickness = 1.0) const {
                //// FIXME: not correct using this method!!!!
                Point2 p1 = _camera.screenProjection(lineProj.first);
                Point2 p2 = _camera.screenProjection(lineProj.second);
                Vec2 y = normalize(p2 - p1);
                Vec2 x(y[1], - y[0]);                
                std::vector<PixelLoc> contour = {
                    ToPixelLoc(p1 - x * thickness),
                    ToPixelLoc(p1 - y * thickness),
                    ToPixelLoc(p1 + x * thickness),
                    ToPixelLoc(p2 + x * thickness),
                    ToPixelLoc(p2 + y * thickness),
                    ToPixelLoc(p2 - x * thickness)
                };
                return MostLikelySurfaceLabelInRegion(_distributions, contour);
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