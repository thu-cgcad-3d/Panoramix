#ifndef PANORAMIX_VIS_VISUALIZE3D_HPP
#define PANORAMIX_VIS_VISUALIZE3D_HPP

#include "../core/basic_types.hpp"
#include "../core/feature.hpp"

#include "basic_types.hpp"

namespace panoramix {
    namespace vis {

        class Visualizer3D {
        public:
            struct Params { // global parameters
                Params();
                std::string winName;
                core::Color backgroundColor;
                core::PerspectiveCamera camera;
                core::Color defaultColor;
                float pointSize;
                float lineWidth;
                core::ColorTableDescriptor colorTableDescriptor;
                RenderModeFlags renderMode;
                core::Mat4 modelMatrix;
            };
            struct VisualData;
            struct Widgets;
            using VisualDataPtr = std::shared_ptr<VisualData>;
            using WidgetsPtr = std::shared_ptr<Widgets>;
            
            explicit Visualizer3D(const Params & p = Params());
            ~Visualizer3D();

            inline VisualDataPtr data() { return _data; }
            inline const VisualDataPtr & data() const { return _data; }
            inline WidgetsPtr widgets() { return _widgets; }
            inline const WidgetsPtr & widgets() const{ return _widgets; }

        public:
            Params & params() const;

        private:
            VisualDataPtr _data;
            WidgetsPtr _widgets;
        };


        namespace manip3d {

            template <class ArgT>
            struct Manipulator {
                inline Manipulator(void(f)(Visualizer3D &, ArgT), ArgT a) : func(f), arg(a){}
                void(*func)(Visualizer3D &, ArgT);    // the function pointer
                ArgT arg;   // the argument value
            };

            Manipulator<std::string> SetWindowName(std::string name);
            Manipulator<core::Color> SetDefaultColor(core::Color color);
            inline Manipulator<core::Color> SetDefaultColor(core::ColorTag tag){ return SetDefaultColor(core::ColorFromTag(tag)); }
            Manipulator<core::Color> SetBackgroundColor(core::Color color);
            inline Manipulator<core::Color> SetBackgroundColor(core::ColorTag tag){ return SetBackgroundColor(core::ColorFromTag(tag)); }
            Manipulator<core::PerspectiveCamera> SetCamera(core::PerspectiveCamera camera);
            Manipulator<float> SetPointSize(float pointSize);
            Manipulator<float> SetLineWidth(float lineWidth);
            Manipulator<core::ColorTableDescriptor> SetColorTableDescriptor(core::ColorTableDescriptor descriptor);
            Manipulator<RenderModeFlags> SetRenderMode(RenderModeFlags mode);
            Manipulator<core::Mat4> SetModelMatrix(core::Mat4 mat);
            Manipulator<bool> Show(bool block = true);

            void AutoSetCamera(Visualizer3D & viz);
        }


        inline Visualizer3D operator << (Visualizer3D viz, void(*func)(Visualizer3D&)) {
            func(viz);
            return viz;
        }

        template <class ArgT>
        inline Visualizer3D operator << (Visualizer3D viz, manip3d::Manipulator<ArgT> smanip) {
            smanip.func(viz, smanip.arg);
            return viz;
        }





        // points
        Visualizer3D operator << (Visualizer3D viz, const core::Point3 & p);
        inline Visualizer3D operator << (Visualizer3D viz, const core::HPoint3 & p) {
            return viz << p.toPoint();
        }

        // lines
        Visualizer3D operator << (Visualizer3D viz, const core::Line3 & l);
        inline Visualizer3D operator << (Visualizer3D viz, const core::HLine3 & l) {
            return viz << l.toLine();
        }

        // classified thing
        template <class T>
        inline Visualizer3D operator << (Visualizer3D viz, const core::Classified<T> & thing) {
            static const auto WhiteColor = core::ColorFromTag(core::ColorTag::White);
            auto & predefinedColorTable = core::PredefinedColorTable(viz.params().colorTableDescriptor);
            viz.params().defaultColor = thing.claz < 0 ? WhiteColor : predefinedColorTable[thing.claz % predefinedColorTable.size()];
            return viz << thing.component;
        }


        // containers
        template <class T, int N>
        inline Visualizer3D operator << (Visualizer3D viz, const std::array<T, N> & a) {
            for (auto & e : a)
                viz = viz << e;
            return viz;
        }

        template <class ContainerT>
        inline Visualizer3D VisualizeAllInContainer(Visualizer3D viz, const ContainerT & c) {
            for (auto & e : c)
                viz = viz << e;
            return viz;
        }


#define VISUALIZE3D_AS_CONTAINER(claz) \
    template <class T> \
    inline Visualizer3D operator << (Visualizer3D viz, const claz<T> & c) { \
        return VisualizeAllInContainer(viz, c); \
    }

        VISUALIZE3D_AS_CONTAINER(std::list)
        VISUALIZE3D_AS_CONTAINER(std::vector)
        VISUALIZE3D_AS_CONTAINER(std::deque)
        VISUALIZE3D_AS_CONTAINER(std::set)
        VISUALIZE3D_AS_CONTAINER(std::unordered_set)
        VISUALIZE3D_AS_CONTAINER(std::forward_list)
            
    }
}
 
#endif