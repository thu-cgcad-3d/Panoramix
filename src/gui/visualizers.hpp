#ifndef PANORAMIX_GUI_VISUALIZERS_HPP
#define PANORAMIX_GUI_VISUALIZERS_HPP

#include "../core/feature.hpp"
#include "../core/cameras.hpp"
#include "../core/iterators.hpp"
#include "../core/containers.hpp"
#include "../core/generic_topo.hpp"
#include "basic_types.hpp"
#include "scene.hpp"

class QWidget;
class QAction;

namespace panoramix {
    namespace gui {     


        class Visualizer {
        public:
            explicit Visualizer(const std::string & winName = "panoramix::vis::Visualizer") {
                _activeOH = _tree.addRoot(std::make_shared<VisualObject>());

                renderOptions.winName = winName;
                renderOptions.backgroundColor = ColorTag::White;
                renderOptions.renderMode = RenderModeFlag::All;
                renderOptions.camera = core::PerspectiveCamera(500, 500, 250, { 1.0, 1.0, 1.0 }, { 0.0, 0.0, 0.0 }, { 0.0, 0.0, 1.0 });
                renderOptions.bwColor = 0.3;
                renderOptions.bwTexColor = 0.7;
                renderOptions.showInside = true;

                installingOptions.discretizeOptions.color = ColorTag::Black;
                installingOptions.discretizeOptions.colorTable = PredefinedColorTable(ColorTableDescriptor::AllColors);
                installingOptions.discretizeOptions.isolatedTriangles = false;
                installingOptions.discretizeOptions.subdivisionNums[0] = 32;
                installingOptions.discretizeOptions.subdivisionNums[1] = 64;
                installingOptions.defaultShaderSource = PredefinedShaderSource(OpenGLShaderSourceDescriptor::XTriangles);
                installingOptions.pointSize = 10.0;
                installingOptions.lineWidth = 5.0;
            }
                        
            VisualObjectInstallingOptions installingOptions;
            RenderOptions renderOptions;

        public:
            const VisualObjectTree & tree() const { return _tree; }

            VisualObjectHandle activeObjectHandle() const { return _activeOH; }
            VisualObject & activeObject() const { return *_tree.data(_activeOH); }

            const TriMesh & activeMesh() const { return _tree.data(_activeOH)->mesh(); }
            TriMesh & activeMesh() { return _tree.data(_activeOH)->mesh(); }

        public:

            template <class T>
            inline Visualizer & add(const T & data) {
                _tree.add(_activeOH, Visualize(data, installingOptions));
                return *this;
            }
            template <class T, class FunT>
            inline Visualizer & add(T & data, const FunT & fun) {
                _tree.add(_activeOH, Visualize(data, fun, doptions,
                    std::integral_constant<CallbackFunctionType, CallbackFunctionTraits<FunT, T>::value>()));
                return *this;
            }


            template <class T>
            inline Visualizer & begin(const T & data) {
                _activeOH = _tree.add(_activeOH, Visualize<T>(data, installingOptions));
                return *this;
            }

            template <class T>
            inline Visualizer & begin(const std::vector<T> & data) {
                _activeOH = _tree.add(_activeOH, VisualizeCollection<T>(data, installingOptions));
                return *this;
            }

            template <class T, class FunT>
            inline Visualizer & begin(T & data, const FunT & fun) { 
                _activeOH = _tree.add(_activeOH, Visualize<T, FunT>(data, fun, installingOptions,
                    std::integral_constant<CallbackFunctionType, CallbackFunctionTraits<FunT, T>::value>()));
                return *this;
            }

            template <class T, class FunT>
            inline Visualizer & begin(std::vector<T> & data, const FunT & fun) {
                _activeOH = _tree.add(_activeOH, VisualizeCollection<T, FunT>(data, fun, installingOptions,
                    std::integral_constant<CallbackFunctionType, CallbackFunctionTraits<FunT, T>::value>()));
                return *this;
            }


            inline Visualizer & shaderSource(const OpenGLShaderSource & ss) { 
                activeObject().setShaderSource(ss);
                return *this; 
            }
            inline Visualizer & resource(const std::string resourceName) {
                activeObject().resources().push_back(ResourceStore::get(resourceName));
                return *this;
            }


            inline Visualizer & end() { 
                if (!_tree.isRoot(_activeOH)){
                    _activeOH = _tree.parent(_activeOH);
                }
                return *this; 
            }

            inline Visualizer & renderMode(RenderModeFlags flags) { renderOptions.renderMode = flags; return *this; }
            inline Visualizer & camera(const core::PerspectiveCamera & cam) { renderOptions.camera = cam; return *this; }

            enum CameraScalePolicy {
                WatchAtMedianScale,
                WatchAtMeanScale,
                WatchAtMaxScale
            };

            QWidget * createWidget(bool autoSetCamera, QWidget * parent);
            void show(bool doModal = true, bool autoSetCamera = true, CameraScalePolicy csp = WatchAtMedianScale);

        private:
            VisualObjectTree _tree;
            VisualObjectHandle _activeOH;
        };

       





    }

}
 
#endif