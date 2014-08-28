#ifndef PANORAMIX_VIS_VISUALIZE2D_HPP
#define PANORAMIX_VIS_VISUALIZE2D_HPP

#include "../core/basic_types.hpp"
#include "basic_types.hpp"

namespace panoramix {
    namespace vis {

        using namespace core;

        class Visualizer2D {
        public:
            struct Params {
                inline Params() 
                    : winName("2D Visualizer"), 
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
            explicit Visualizer2D(const Image & im = Image(), const Params & p = Params());

            void setImage(const Image & im);
            inline Image & image() { return _image; }

        public:
            Params params;

        private:
            Image _image;
        };


        // manipulators
        namespace manip2d{
            template <class ArgT>
            struct Manipulator {
                inline Manipulator(void(f)(Visualizer2D &, ArgT), ArgT a) : func(f), arg(a){}
                void(*func)(Visualizer2D &, ArgT);    // the function pointer
                ArgT arg;   // the argument value
            };

            inline Manipulator<const std::string &> SetWindowName(const std::string & name) {
                return Manipulator<const std::string &>(
                    [](Visualizer2D & viz, const std::string & name){
                    viz.params.winName = name; },
                        name);
            }

            inline Manipulator<Color> SetColor(const Color & color) {
                return Manipulator<Color>(
                    [](Visualizer2D & viz, Color c){
                    viz.params.color = c; },
                        color);
            }

            inline Manipulator<Color> SetColor(const ColorTag & color) {
                return SetColor(ColorFromTag(color));
            }

            inline Manipulator<int> SetThickness(int thickness) {
                return Manipulator<int>(
                    [](Visualizer2D & viz, int t){
                    viz.params.thickness = t; },
                        thickness);
            }

            inline Manipulator<int> SetLineType(int type) {
                return Manipulator<int>(
                    [](Visualizer2D & viz, int t){
                    viz.params.lineType = t; },
                        type);
            }

            inline Manipulator<float> SetAlphaForNewImage(float a) {
                return Manipulator<float>(
                    [](Visualizer2D & viz, float alpha){
                    viz.params.alphaForNewImage = alpha; },
                        a);
            }

            inline Manipulator<ColorTableDescriptor> SetColorTableDescriptor(ColorTableDescriptor descriptor) {
                return Manipulator<ColorTableDescriptor>(
                    [](Visualizer2D & viz, ColorTableDescriptor d){
                    viz.params.colorTableDescriptor = d; },
                        descriptor);
            }
            
            Manipulator<int> Show(int delay = 1);
        }

        inline Visualizer2D operator << (Visualizer2D viz, void (*func)(Visualizer2D&)) {
            func(viz);
            return viz;
        }

        template <class ArgT>
        inline Visualizer2D operator << (Visualizer2D viz, manip2d::Manipulator<ArgT> smanip) {
            smanip.func(viz, smanip.arg);
            return viz;
        }

        


        // points
        template <class T>
        inline Visualizer2D operator << (Visualizer2D viz, const Point<T, 2> & p) {
            int x = static_cast<int>(std::round(p[0]));
            int y = static_cast<int>(std::round(p[1]));
            cv::circle(viz.image(), cv::Point(x, y), 1, viz.params.color, viz.params.thickness, viz.params.lineType, viz.params.shift);
            return viz;
        }

        template <class T>
        inline Visualizer2D operator << (Visualizer2D viz, const HPoint<T, 2> & p) {
            return viz << p.toPoint();
        }

        // lines
        template <class T>
        inline Visualizer2D operator << (Visualizer2D viz, const Line<T, 2> & line) {
            cv::line(viz.image(), cv::Point(static_cast<int>(line.first(0)), static_cast<int>(line.first(1))),
                cv::Point(static_cast<int>(line.second(0)), static_cast<int>(line.second(1))),
                viz.params.color, viz.params.thickness, viz.params.lineType, viz.params.shift);
            return viz;
        }

        template <class T>
        inline Visualizer2D operator << (Visualizer2D viz, const HLine<T, 2> & line) {
            return viz << line.toLine();
        }

        // keypoints
        inline Visualizer2D operator << (Visualizer2D viz, const KeyPoint & p) {
            cv::drawKeypoints(viz.image(), std::vector<KeyPoint>(1, p), viz.image(), viz.params.color);
            return viz;
        }

        inline Visualizer2D operator << (Visualizer2D viz, const std::vector<KeyPoint> & ps) {
            cv::drawKeypoints(viz.image(), ps, viz.image(), viz.params.color);
            return viz;
        }

        // classified thing
        template <class T>
        inline Visualizer2D operator << (Visualizer2D viz, const Classified<T> & thing) {
            static const auto WhiteColor = ColorFromTag(ColorTag::White);
            auto & predefinedColorTable = PredefinedColorTable(viz.params.colorTableDescriptor);
            viz.params.color = thing.claz < 0 ? WhiteColor : predefinedColorTable[thing.claz % predefinedColorTable.size()];
            return viz << thing.component;
        }


        // image
        Visualizer2D operator << (Visualizer2D viz, const Image & im);


        // containers
        template <class T, int N>
        inline Visualizer2D operator << (Visualizer2D viz, const std::array<T, N> & a) {
            for (auto & e : a)
                viz << e;
            return viz;
        }

        template <class ContainerT>
        inline Visualizer2D VisualizeAllInContainer(Visualizer2D viz, const ContainerT & c) {
            for (auto & e : c)
                viz << e;
            return viz;
        }

        #define VISUALIZE2D_AS_CONTAINER(claz) \
            template <class T> \
            inline Visualizer2D operator << (Visualizer2D viz, const claz<T> & c) { \
                return VisualizeAllInContainer(viz, c); \
            }

        VISUALIZE2D_AS_CONTAINER(std::list)
        VISUALIZE2D_AS_CONTAINER(std::vector)
        VISUALIZE2D_AS_CONTAINER(std::deque)
        VISUALIZE2D_AS_CONTAINER(std::set)
        VISUALIZE2D_AS_CONTAINER(std::unordered_set)
        VISUALIZE2D_AS_CONTAINER(std::forward_list)


    }
}
 
#endif