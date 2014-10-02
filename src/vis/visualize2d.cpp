#include "../core/utilities.hpp"

#include "visualize2d.hpp"

namespace panoramix {
    namespace vis {

        Visualizer2D::Visualizer2D(const Image & im, const Params & p) : params(p){
            setImage(im);
        }

        void Visualizer2D::setImage(const Image & im) {
            _image = im.clone();
        }


        Visualizer2D operator << (Visualizer2D viz, const Image & im) {
            assert(!im.empty());
            cv::addWeighted(viz.image(), (1.0f - viz.params.alphaForNewImage), im, viz.params.alphaForNewImage, 0.0, viz.image());
            return viz;
        }

        namespace manip2d {

            Manipulator<int> Show(int delay) {
                return Manipulator<int>([](Visualizer2D & viz, int d) {
                    static int id = 0;
                    cv::imshow(viz.params.winName, viz.image());
                    cv::imwrite("./vis_outputs/" + std::to_string(id++) + viz.params.winName + ".png", viz.image());
                    cv::waitKey(d);
                }, delay);
            }

        }


        inline PixelLoc ToPixelLoc(const core::Point2 & p){
            return PixelLoc(static_cast<int>(p[0]), static_cast<int>(p[1]));
        }

        Visualizer2D operator << (Visualizer2D viz, const InfiniteLine<double, 2> & line){
            Point2 imCenter(viz.image().cols / 2.0, viz.image().rows / 2.0);
            auto d = DistanceFromPointToLine(imCenter, line);
            const Point2 & root = d.second;
            auto dir = normalize(line.direction);
            double scale = norm(imCenter) * 2;
            PixelLoc p1 = ToPixelLoc(root - dir * scale);
            PixelLoc p2 = ToPixelLoc(root + dir * scale);
            cv::clipLine(cv::Rect(0, 0, viz.image().cols, viz.image().rows), p1, p2);
            cv::line(viz.image(), p1, p2, viz.params.color, viz.params.thickness, viz.params.lineType, viz.params.shift);
            return viz;
        }

    }
}