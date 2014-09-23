#ifndef PANORAMIX_VIS_VISUALIZE3D_HPP
#define PANORAMIX_VIS_VISUALIZE3D_HPP

#include "../core/feature.hpp"
#include "basic_types.hpp"
#include "renderable_object_tree.hpp"

namespace panoramix {
    namespace vis {
            
        class Visualizer3D {
            struct PrivateData;
        public:
            struct Params { // predefined parameters for widgets
                Params();
                std::string winName;
                vis::Color backgroundColor;
                core::PerspectiveCamera camera;
                RenderModeFlags renderMode;
            };
            
            Visualizer3D();
            explicit Visualizer3D(const Params & p, const DefaultRenderState & s);

        public:
            template <class T>
            inline void addAndActivate(const T & data) {
                RenderableObject * newRo = MakeRenderable(data, _defaultRenderState, _activeObject);
                _activeObject = newRo;
            }
            void deactivateLast();

            std::shared_ptr<RenderableObject> root() const { return _root; }
            RenderableObject * activeObject() const { return _activeObject; }
            std::shared_ptr<PrivateData> data() const { return _data; }

        public:
            Params params;
            DefaultRenderState defaultRenderState;

        private:
            DefaultRenderState _defaultRenderState;
            std::shared_ptr<RenderableObject> _root; // never empty
            RenderableObject * _activeObject; // never empty
            std::shared_ptr<PrivateData> _data; // GUI data
        };


        namespace manip3d {

            template <class ArgT>
            struct Manipulator {
                inline Manipulator(void(f)(Visualizer3D &, ArgT), ArgT a) : func(f), arg(a){}
                void(*func)(Visualizer3D &, ArgT);    // the function pointer
                ArgT arg;   // the argument value
            };

            Manipulator<const std::string &> SetWindowName(const std::string & name);
            Manipulator<vis::Color> SetDefaultForegroundColor(vis::Color color);
            inline Manipulator<vis::Color> SetDefaultForegroundColor(vis::ColorTag tag) { return SetDefaultForegroundColor(vis::ColorFromTag(tag)); }
            Manipulator<vis::Color> SetBackgroundColor(vis::Color color);
            inline Manipulator<vis::Color> SetBackgroundColor(vis::ColorTag tag){ return SetBackgroundColor(vis::ColorFromTag(tag)); }
            Manipulator<const core::PerspectiveCamera &> SetCamera(const core::PerspectiveCamera & camera);
            Manipulator<float> SetDefaultPointSize(float pointSize);
            Manipulator<float> SetDefaultLineWidth(float lineWidth);
            Manipulator<const vis::ColorTable &> SetDefaultColorTable(const vis::ColorTable & colorTable);
            inline Manipulator<const vis::ColorTable &> SetDefaultColorTable(vis::ColorTableDescriptor d) { return SetDefaultColorTable(vis::ColorTable(d)); }
            Manipulator<RenderModeFlags> SetRenderMode(RenderModeFlags mode);
            Manipulator<std::pair<bool, bool>> Show(bool doModel = true, bool autoSetCamera = true);

            template <class T>
            Manipulator<const T &> Begin(const T & t) {
                return Manipulator<const T &>([](Visualizer3D & viz, const T & t) {
                    viz.addAndActivate(t);
                }, t);
            }

            void End(Visualizer3D & viz) {
                viz.deactivateLast();
            }

            Manipulator<const core::Mat4 &> SetModelMatrix(const core::Mat4 & mat);
            Manipulator<const core::Image &> SetTexture(const core::Image & tex);

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
            return viz << p.value();
        }

        // lines
        Visualizer3D operator << (Visualizer3D viz, const core::Line3 & l);
        inline Visualizer3D operator << (Visualizer3D viz, const core::HLine3 & l) {
            return viz << l.toLine();
        }


        // classified thing
        template <class T>
        inline Visualizer3D operator << (Visualizer3D viz, const core::Classified<T> & thing) {
            viz.params().defaultColor = viz.params().colorTable[thing.claz];
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