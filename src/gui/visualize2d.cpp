#include "../core/utilities.hpp"

#include "visualize2d.hpp"

namespace panoramix {
    namespace gui {

        inline core::Vec3b ToVec3b(const Color & color) {
            return core::Vec3b(color.red(), color.green(), color.blue());
        }

        Visualizer2D::Visualizer2D(const Image & im, const Params & p) : params(p){
            setImage(im);
        }

        void Visualizer2D::setImage(const Image & im) {
            if (im.type() == CV_32SC1){
                ImageOfType<Vec3b> imv;
                imv.create(im.size());
                for (int x = 0; x < im.cols; x++){
                    for (int y = 0; y < im.rows; y++){
                        imv(y, x) = ToVec3b(params.colorTable.roundedAt(im.at<int32_t>(y, x)));
                    }
                }
                _image = imv;
            }
            else{
                _image = im.clone();
            }
        }


        Visualizer2D operator << (Visualizer2D viz, const Image & im) {
            assert(!im.empty());
            cv::addWeighted(viz.image(), (1.0f - viz.params.alphaForNewImage), im, viz.params.alphaForNewImage, 0.0, viz.image());
            return viz;
        }

        Visualizer2D operator << (Visualizer2D viz, const ImageOfType<int32_t> & im) {
            ImageOfType<Vec3b> imv;
            imv.create(im.size());
            for (int x = 0; x < im.cols; x++){
                for (int y = 0; y < im.rows; y++){
                    imv(y, x) = ToVec3b(viz.params.colorTable[im(y, x)]);
                }
            }
            return viz << (const Image&)imv;
        }

        namespace manip2d {

            Manipulator<std::pair<int, bool>> Show(int delay, bool asLabels) {
                return Manipulator<std::pair<int, bool>>([](Visualizer2D & viz, std::pair<int, bool> d) {
                    static int id = 0;
                    core::Image im = viz.image().clone();
                    if (viz.image().channels() > 3){
                        std::vector<cv::Mat> cs(3);
                        for (int i = 0; i < 3; i++)
                            cv::extractChannel(im, cs[i], i);
                        cv::merge(cs, im);
                    }
                    cv::imshow(viz.params.winName, im);
                    cv::imwrite("./vis_outputs/" + std::to_string(id++) + viz.params.winName + ".png", im);
                    cv::waitKey(d.first);
                }, std::make_pair(delay, asLabels));
            }

        }


        inline PixelLoc ToPixelLoc(const core::Point2 & p){
            return PixelLoc(static_cast<int>(p[0]), static_cast<int>(p[1]));
        }

        Visualizer2D operator << (Visualizer2D viz, const Ray<double, 2> & line){
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