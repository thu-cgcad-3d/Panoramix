#include "visualize2d.hpp"

namespace panoramix {
    namespace vis {

        Visualizer2D::Visualizer2D(const Image & im, const Params & p) : params(p) {
            setImage(im);
        }

        namespace {
            inline Image ConvertImage(const Image & im) {
                Image image;
                im.copyTo(image);
                return image;
            }
        }

        void Visualizer2D::setImage(const Image & im) {
            _image = ConvertImage(im);
        }


        Visualizer2D operator << (Visualizer2D viz, const Image & im) {
            assert(!im.empty());
            cv::addWeighted(viz.image(), (1.0f - viz.params.alphaForNewImage), im, viz.params.alphaForNewImage, 0.0, viz.image());
            return viz;
        }

        namespace manip {

            Manipulator<int> Show(int delay) {
                return Manipulator<int>([](Visualizer2D & viz, int d) {
                    cv::imshow(viz.params.winName, viz.image());
                    cv::waitKey(d);
                }, delay);
            }

        }


    }
}