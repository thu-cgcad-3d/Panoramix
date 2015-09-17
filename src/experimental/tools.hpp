#pragma once


#include "../gui/basic_types.hpp"
#include "../gui/scene.hpp"
#include "../gui/canvas.hpp"

#include "rl_graph_solver.hpp"

namespace pano {
    namespace experimental {



        // output a 2d image of this rl graph
        template <
            class CameraT,
            class RhColorerT = core::ConstantFunctor<gui::ColorTag>,
            class LhColorerT = core::ConstantFunctor<gui::ColorTag>,
            class RBhColorerT = core::ConstantFunctor<gui::ColorTag>
        >
        Image3f Print(const RLGraph & mg, const Imagei & segmentedRegions, const CameraT & camera, const std::vector<RegionHandle> & rh2rids,
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
                    std::vector<Pixel> ps; ps.reserve(spanAngle / sampleAngle);
                    for (double angle = 0.0; angle <= spanAngle; angle += sampleAngle){
                        Vec3 dir = RotateDirection(line.first, line.second, angle);
                        ps.push_back(ToPixel(camera.toScreen(dir)));
                    }
                    for (int i = 1; i < ps.size(); i++){
                        auto & p1 = ps[i - 1];
                        auto & p2 = ps[i];
                        if (Distance(p1, p2) >= rendered.cols / 2) {
                            continue;
                        }
                        cv::clipLine(cv::Rect(0, 0, rendered.cols, rendered.rows), p1, p2);
                        cv::line(rendered, p1, p2, (cv::Scalar)color / 255.0, lineWidth);
                        //rendered(ps[i]) = color;
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
                            auto p1 = core::ToPixel(camera.toScreen(e[i - 1]));
                            auto p2 = core::ToPixel(camera.toScreen(e[i]));
                            cv::clipLine(cv::Rect(0, 0, rendered.cols, rendered.rows), p1, p2);
                            cv::line(rendered, p1, p2, (cv::Scalar)color / 255.0, boundaryWidth);
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






    


    }
}
