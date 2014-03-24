#include "feature_visualize.hpp"

namespace panoramix {
    namespace core {

        ImageFeatureVisualizer::ImageFeatureVisualizer(const Image & im, const Params & p) : params(p) {
            setImage(im);
        }

        namespace {
            inline Image ConvertImage(const Image & im) {
                Image image;
                im.copyTo(image);
                return image;
            }
        }

        void ImageFeatureVisualizer::setImage(const Image & im) {
            _image = ConvertImage(im);
        }


        ImageFeatureVisualizer operator << (ImageFeatureVisualizer viz, const Image & im) {
            assert(!im.empty());
            cv::addWeighted(viz.image(), (1.0f - viz.params.alphaForNewImage), im, viz.params.alphaForNewImage, 0.0, viz.image());
            return viz;
        }

        namespace manip {

            Manipulator<int> Show(int delay) {
                return Manipulator<int>([](ImageFeatureVisualizer & viz, int d) {
                    cv::imshow(viz.params.winName, viz.image());
                    cv::waitKey(d);
                }, delay);
            }

        }


    }
}