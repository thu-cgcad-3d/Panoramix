#ifndef PANORAMIX_GUI_VISUALIZE2D_HPP
#define PANORAMIX_GUI_VISUALIZE2D_HPP

#include "../core/basic_types.hpp"
#include "basic_types.hpp"

namespace panoramix {
    namespace gui {

        using namespace core;

        class Visualizer2D {
        public:
            struct Params {
                inline Params() 
                    : winName("2D Visualizer"), 
                      color(255, 255, 255), thickness(1), 
                      lineType(8), shift(0), alphaForNewImage(0.5f), 
                      colorTable(ColorTableDescriptor::AllColors) {}

                std::string winName;
                Color color;
                int thickness;
                int lineType;
                int shift;
                float alphaForNewImage;
                ColorTable colorTable;
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

            inline Manipulator<Color> SetColor(Color color) {
                return Manipulator<Color>(
                    [](Visualizer2D & viz, Color c){
                    viz.params.color = c; },
                        color);
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

            inline Manipulator<ColorTable> SetColorTable(const ColorTable & colorTable) {
                return Manipulator<ColorTable>(
                    [](Visualizer2D & viz, ColorTable ct) {
                    viz.params.colorTable = ct; },
                        colorTable);
            }

            inline Manipulator<ColorTable> SetColorTable(const ColorTableDescriptor & d) {
                return SetColorTable(PredefinedColorTable(d));
            }
            
            Manipulator<std::pair<int, bool>> Show(int delay = 0, bool asLabels = false);
        }

        inline Visualizer2D & operator << (Visualizer2D & viz, void (*func)(Visualizer2D&)) {
            func(viz);
            return viz;
        }

        template <class ArgT>
        inline Visualizer2D & operator << (Visualizer2D & viz, manip2d::Manipulator<ArgT> smanip) {
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
            return viz << p.value();
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

        Visualizer2D operator << (Visualizer2D viz, const Ray<double, 2> & line);
        
        template <class T>
        inline Visualizer2D operator << (Visualizer2D viz, const Ray<T, 2> & line){
            return viz << Ray2(core::vec_cast<double>(line.anchor), core::vec_cast<double>(line.direction));
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
            viz.params.color = viz.params.colorTable[thing.claz];
            return viz << thing.component;
        }

        // noted thing
        template <class T>
        inline Visualizer2D operator << (Visualizer2D viz, const Noted<T> & thing){
            auto center = core::BoundingBox(thing.component).center();
            cv::putText(viz.image(), thing.note, cv::Point(static_cast<int>(center[0]), static_cast<int>(center[1])), 1, 1, viz.params.color);
            return viz << thing.component;
        }

        // enabled thing
        template <class T>
        inline Visualizer2D operator << (Visualizer2D viz, const Enabled<T> & thing){
            if (thing.enabled){
                return viz << thing.component;
            }
            return viz;
        }


        // image
        Visualizer2D operator << (Visualizer2D viz, const Image & im);
        Visualizer2D operator << (Visualizer2D viz, const ImageOfType<int32_t> & im);


        // integer image
        template <class T, class = std::enable_if_t<std::is_integral<T>::value>>
        Visualizer2D operator << (Visualizer2D viz, const ImageOfType<T> & im) {
            NOT_IMPLEMENTED_YET();
        }


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