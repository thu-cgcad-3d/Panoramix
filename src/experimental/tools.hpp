#ifndef PANORAMIX_EXPERIMENTAL_TOOLS_HPP
#define PANORAMIX_EXPERIMENTAL_TOOLS_HPP

#include "../gui/basic_types.hpp"
#include "../gui/scene.hpp"
#include "../gui/canvas.hpp"

#include "rl_graph.hpp"

namespace panoramix {
    namespace experimental {



        // output a 2d image of this rl graph
        template <
            class CameraT,
            class RhColorerT = core::ConstantFunctor<gui::Color>,
            class LhColorerT = core::ConstantFunctor<gui::Color>,
            class RBhColorerT = core::ConstantFunctor<gui::Color>
        >
        Image Print(const RLGraph & mg, const Imagei & segmentedRegions, const CameraT & camera, const std::vector<RegionHandle> & rh2rids,
        RhColorerT && rhColor = gui::Transparent,
        LhColorerT && lhColor = gui::Transparent,
        RBhColorerT && rbColor = gui::Transparent,
        int boundaryWidth = 1, int lineWidth = 2){
            Image3f rendered = Image3f::zeros(segmentedRegions.size());
            // regions
            for (auto it = rendered.begin(); it != rendered.end(); ++it){
                RegionHandle rh = rh2rids[segmentedRegions(it.pos())];
                if (rh.invalid())
                    continue;
                gui::Color color = rhColor(rh);
                *it = Vec3f(color.bluef(), color.greenf(), color.redf());
            }
            // lines
            if (lineWidth > 0){
                for (auto & ld : mg.components<LineData>()){
                    gui::Color color = lhColor(ld.topo.hd);
                    if (color.isTransparent())
                        continue;
                    static const double sampleAngle = M_PI / 100.0;
                    auto & line = ld.data.line;
                    double spanAngle = AngleBetweenDirections(line.first, line.second);
                    std::vector<PixelLoc> ps; ps.reserve(spanAngle / sampleAngle);
                    for (double angle = 0.0; angle <= spanAngle; angle += sampleAngle){
                        Vec3 dir = RotateDirection(line.first, line.second, angle);
                        ps.push_back(ToPixelLoc(camera.toScreen(dir)));
                    }
                    for (int i = 1; i < ps.size(); i++){
                        auto & p1 = ps[i - 1];
                        auto & p2 = ps[i];
                        cv::clipLine(cv::Rect(0, 0, rendered.cols, rendered.rows), p1, p2);
                        cv::line(rendered, p1, p2, (cv::Scalar)color, lineWidth);
                    }
                }
            }
            // region boundary
            if (boundaryWidth > 0){
                for (auto & b : mg.constraints<RegionBoundaryData>()){
                    gui::Color color = rbColor(b.topo.hd);
                    if (color.isTransparent())
                        continue;
                    for (auto & e : b.data.normalizedEdges){
                        if (e.size() <= 1) continue;
                        for (int i = 1; i < e.size(); i++){
                            auto p1 = core::ToPixelLoc(camera.toScreen(e[i - 1]));
                            auto p2 = core::ToPixelLoc(camera.toScreen(e[i]));
                            cv::clipLine(cv::Rect(0, 0, rendered.cols, rendered.rows), p1, p2);
                            cv::line(rendered, p1, p2, (cv::Scalar)color, boundaryWidth);
                        }
                    }
                }
            }
            return rendered;
        }



        // visualize current mixed graph
        void Visualize(gui::SceneBuilder & vis, const View<PanoramicCamera> & texture,
            const RLGraph & mg, const RLGraphControls & controls, const RLGraphVars & vars);
        void Visualize(gui::SceneBuilder & vis, const View<PerspectiveCamera> & texture,
            const RLGraph & mg, const RLGraphControls & controls, const RLGraphVars & vars);






        // cut loop
        struct SectionalPiece {
            RegionHandle rh;
            std::pair<Point3, Point3> range;
        };
        std::vector<SectionalPiece> MakeSectionalPieces(const HandledTable<RegionHandle, std::vector<Polygon3>> & polygons,
            const Plane3 & cutplane);

        Chain3 MakeChain(const std::vector<SectionalPiece> & pieces, bool closed = true);
        inline double Area(const std::vector<SectionalPiece> & loop) { return Area(Polygon3(MakeChain(loop))); }



        // estimate 
        std::pair<double, double> EstimateEffectiveRangeAlongDirection(
            const HandledTable<RegionHandle, std::vector<Polygon3>> & polygons,
            const Vec3 & direction, double stepLen, double minEffectiveAreaRatio = 0.6,
            double gamma1 = 0.05, double gamma2 = 0.05);


        // smooth
        void SmoothInstances(const RLGraph & mg, const InstanceTable<RegionData> & planes, const InstanceTable<LineData> & lines);


    }
}

#endif