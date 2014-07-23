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
                vis::Color backgroundColor;
                core::PerspectiveCamera camera;
                vis::Color defaultColor;
                float pointSize;
                float lineWidth;
                vis::ColorTableDescriptor colorTableDescriptor;
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
            Manipulator<vis::Color> SetDefaultColor(vis::Color color);
            inline Manipulator<vis::Color> SetDefaultColor(vis::ColorTag tag){ return SetDefaultColor(vis::ColorFromTag(tag)); }
            Manipulator<vis::Color> SetBackgroundColor(vis::Color color);
            inline Manipulator<vis::Color> SetBackgroundColor(vis::ColorTag tag){ return SetBackgroundColor(vis::ColorFromTag(tag)); }
            Manipulator<core::PerspectiveCamera> SetCamera(core::PerspectiveCamera camera);
            Manipulator<float> SetPointSize(float pointSize);
            Manipulator<float> SetLineWidth(float lineWidth);
            Manipulator<vis::ColorTableDescriptor> SetColorTableDescriptor(vis::ColorTableDescriptor descriptor);
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

        // image as texture
        Visualizer3D operator << (Visualizer3D viz, const core::Image & tex);

        // polygons
        Visualizer3D operator << (Visualizer3D viz, const std::vector<std::pair<core::Point3, core::Point2>> & polygonWithTexCoords);


        // classified thing
        template <class T>
        inline Visualizer3D operator << (Visualizer3D viz, const core::Classified<T> & thing) {
            auto oldDefaultColor = viz.params().defaultColor;
            auto & predefinedColorTable = PredefinedColorTable(viz.params().colorTableDescriptor);
            if (thing.claz >= 0){
                viz.params().defaultColor = predefinedColorTable[thing.claz % predefinedColorTable.size()];
            }
            viz << thing.component;
            viz.params().defaultColor = oldDefaultColor;
            return viz;
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
           



        class AdvancedVisualizer3D {
        public:
            struct Params { // global parameters
                Params();
                std::string winName;
                vis::Color backgroundColor;
                core::PerspectiveCamera camera;
                vis::ColorTableDescriptor colorTableDescriptor;
                core::Mat4 modelMatrix;
            };

            struct VisualData;
            struct Widgets;
            using VisualDataPtr = std::shared_ptr<VisualData>;
            using WidgetsPtr = std::shared_ptr<Widgets>;

            explicit AdvancedVisualizer3D(const Params & params = Params());
            ~AdvancedVisualizer3D();

        public:
            Params & params() const;

        private:
            VisualDataPtr _data;
            WidgetsPtr _widgets;
        };

    }
}
 
#endif