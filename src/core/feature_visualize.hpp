#ifndef PANORAMIX_CORE_FEATURE_VISUALIZE_HPP
#define PANORAMIX_CORE_FEATURE_VISUALIZE_HPP

#include "basic_types.hpp"

namespace panoramix {
    namespace core {
 
        class ImageFeatureVisualizer {
        public:
            using Image = cv::Mat;
            struct Params {
                inline Params() 
                    : winName("Image Feature Visualizer"), 
                      color(255, 255, 255), thickness(1), 
                      lineType(8), shift(0), alphaForNewImage(0.5f), 
                      colorTableDescriptor(ColorTableDescriptor::AllColors) {}

                std::string winName;
                Color color;
                int thickness;
                int lineType;
                int shift;
                float alphaForNewImage;
                ColorTableDescriptor colorTableDescriptor;
            };

        public:
            explicit ImageFeatureVisualizer(const Image & im = Image(), const Params & p = Params());

            void setImage(const Image & im);
            inline Image & image() { return _image; }

        public:
            Params params;
        private:
            Image _image;
        };


        // manipulators
        namespace manip{
            template <class ArgT>
            struct Manipulator {
                inline Manipulator(void(f)(ImageFeatureVisualizer &, ArgT), ArgT a) : func(f), arg(a){}
                void(*func)(ImageFeatureVisualizer &, ArgT);	// the function pointer
                ArgT arg;	// the argument value
            };

            inline Manipulator<const std::string &> SetWindowName(const std::string & name) {
                return Manipulator<const std::string &>(
                    [](ImageFeatureVisualizer & viz, const std::string & name){
                    viz.params.winName = name; },
                        name);
            }

            inline Manipulator<const Color &> SetColor(const Color & color) {
                return Manipulator<const Color &>(
                    [](ImageFeatureVisualizer & viz, const Color & c){
                    viz.params.color = c; },
                        color);
            }

            inline Manipulator<int> SetThickness(int thickness) {
                return Manipulator<int>(
                    [](ImageFeatureVisualizer & viz, int t){
                    viz.params.thickness = t; },
                        thickness);
            }

            inline Manipulator<int> SetLineType(int type) {
                return Manipulator<int>(
                    [](ImageFeatureVisualizer & viz, int t){
                    viz.params.lineType = t; },
                        type);
            }

            inline Manipulator<float> SetAlphaForNewImage(float a) {
                return Manipulator<float>(
                    [](ImageFeatureVisualizer & viz, float alpha){
                    viz.params.alphaForNewImage = alpha; },
                        a);
            }

            inline Manipulator<ColorTableDescriptor> SetColorTableDescriptor(ColorTableDescriptor descriptor) {
                return Manipulator<ColorTableDescriptor>(
                    [](ImageFeatureVisualizer & viz, ColorTableDescriptor d){
                    viz.params.colorTableDescriptor = d; },
                        descriptor);
            }
            
            Manipulator<int> Show(int delay = 0);
        }

        inline ImageFeatureVisualizer operator << (ImageFeatureVisualizer viz, void (*func)(ImageFeatureVisualizer&)) {
            func(viz);
            return viz;
        }

        template <class ArgT>
        inline ImageFeatureVisualizer operator << (ImageFeatureVisualizer viz, manip::Manipulator<ArgT> smanip) {
            smanip.func(viz, smanip.arg);
            return viz;
        }

        


        // points
        template <class T, int N>
        inline ImageFeatureVisualizer operator << (ImageFeatureVisualizer viz, const Point<T, N> & p) {
            int x = static_cast<int>(std::round(p[0]));
            int y = N < 2 ? 0 : static_cast<int>(std::round(p[1]));
            cv::circle(viz.image(), cv::Point(x, y), 1, viz.params.color, viz.params.thickness, viz.params.lineType, viz.params.shift);
            return viz;
        }

        template <class T, int N>
        inline ImageFeatureVisualizer operator << (ImageFeatureVisualizer viz, const HPoint<T, N> & p) {
            return viz << p.toPoint();
        }

        // lines
        template <class T>
        inline ImageFeatureVisualizer operator << (ImageFeatureVisualizer viz, const Line<T, 2> & line) {
            cv::line(viz.image(), cv::Point(line.first(0), line.first(1)),
                cv::Point(line.second(0), line.second(1)),
                viz.params.color, viz.params.thickness, viz.params.lineType, viz.params.shift);
            return viz;
        }

        template <class T>
        inline ImageFeatureVisualizer operator << (ImageFeatureVisualizer viz, const HLine<T, 2> & line) {
            return viz << line.toLine();
        }

        // keypoints
        inline ImageFeatureVisualizer operator << (ImageFeatureVisualizer viz, const KeyPoint & p) {
            cv::drawKeypoints(viz.image(), std::vector<KeyPoint>(1, p), viz.image(), viz.params.color);
            return viz;
        }

        inline ImageFeatureVisualizer operator << (ImageFeatureVisualizer viz, const std::vector<KeyPoint> & ps) {
            cv::drawKeypoints(viz.image(), ps, viz.image(), viz.params.color);
            return viz;
        }

        // classified thing
        template <class T>
        inline ImageFeatureVisualizer operator << (ImageFeatureVisualizer viz, const Classified<T> & thing) {
            static const auto WhiteColor = ColorFromTag(ColorTag::White);
            auto & predefinedColorTable = PredefinedColorTable(viz.params.colorTableDescriptor);
            viz.params.color = thing.claz < 0 ? WhiteColor : predefinedColorTable[thing.claz % predefinedColorTable.size()];
            return viz << thing.component;
        }


        // image
        ImageFeatureVisualizer operator << (ImageFeatureVisualizer viz, const Image & im);


        // containers
        template <class T, int N>
        inline ImageFeatureVisualizer operator << (ImageFeatureVisualizer viz, const std::array<T, N> & a) {
            for (auto & e : a)
                viz << e;
            return viz;
        }

        template <class ContainerT>
        inline ImageFeatureVisualizer VisualizeAllInContainer(ImageFeatureVisualizer viz, const ContainerT & c) {
            for (auto & e : c)
                viz << e;
            return viz;
        }

        #define VISUALIZE_AS_CONTAINER(claz) \
            template <class T> \
            inline ImageFeatureVisualizer operator << (ImageFeatureVisualizer viz, const claz<T> & c) { \
                return VisualizeAllInContainer(viz, c); \
            }

        VISUALIZE_AS_CONTAINER(std::list)
        VISUALIZE_AS_CONTAINER(std::vector)
        VISUALIZE_AS_CONTAINER(std::deque)
        VISUALIZE_AS_CONTAINER(std::set)
        VISUALIZE_AS_CONTAINER(std::unordered_set)
        VISUALIZE_AS_CONTAINER(std::forward_list)

 
    }
}
 
#endif