#ifndef PANORAMIX_EXPERIMENTAL_TOOLS_HPP
#define PANORAMIX_EXPERIMENTAL_TOOLS_HPP

#include "../gui/basic_types.hpp"
#include "../gui/visualize2d.hpp"
#include "../gui/visualizers.hpp"

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
        RBhColorerT && rbColor = gui::Transparent){
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
                    ps.push_back(ToPixelLoc(camera.screenProjection(dir)));
                }
                for (int i = 1; i < ps.size(); i++){
                    auto & p1 = ps[i - 1];
                    auto & p2 = ps[i];
                    cv::clipLine(cv::Rect(0, 0, rendered.cols, rendered.rows), p1, p2);
                    cv::line(rendered, p1, p2, (cv::Scalar)color, 2);
                }
            }
            // region boundary
            for (auto & b : mg.constraints<RegionBoundaryData>()){
                gui::Color color = rbColor(b.topo.hd);
                if (color.isTransparent())
                    continue;
                for (auto & e : b.data.normalizedEdges){
                    if (e.size() <= 1) continue;
                    for (int i = 1; i < e.size(); i++){
                        auto p1 = core::ToPixelLoc(camera.screenProjection(e[i - 1]));
                        auto p2 = core::ToPixelLoc(camera.screenProjection(e[i]));
                        cv::clipLine(cv::Rect(0, 0, rendered.cols, rendered.rows), p1, p2);
                        cv::line(rendered, p1, p2, (cv::Scalar)color, 1);
                    }
                }
            }
            return rendered;
        }



        // visualize current mixed graph
        void Visualize(gui::Visualizer & vis, const View<PanoramicCamera> & texture,
            const RLGraph & mg, const RLGraphControls & controls, const RLGraphVars & vars);
        void Visualize(gui::Visualizer & vis, const View<PerspectiveCamera> & texture,
            const RLGraph & mg, const RLGraphControls & controls, const RLGraphVars & vars);



    }
}

#endif