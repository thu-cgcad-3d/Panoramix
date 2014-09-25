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

            inline void End(Visualizer3D & viz) {
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

        template <class T, class = std::enable_if_t<CanMakeRenderable<T>::value>>
        inline Visualizer3D operator << (Visualizer3D viz, const T & v) {
            return viz << manip3d::Begin(v) << manip3d::End;
        }

    }
}
 
#endif