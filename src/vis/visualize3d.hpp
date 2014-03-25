#ifndef PANORAMIX_VIS_VISUALIZE3D_HPP
#define PANORAMIX_VIS_VISUALIZE3D_HPP

#include "../core/basic_types.hpp"
#include "../core/feature.hpp"

class QWidget;

namespace panoramix {
    namespace vis {

        using namespace core;

        class VisualizerGUIData;
        class Visualizer3D {
        public:
            struct Params {
                Params();
                std::string winName;
                Color color, backgroundColor;
                PerspectiveCamera camera;
                int thickness;
                int lineType;
                ColorTableDescriptor colorTableDescriptor;
            };

        public:
            explicit Visualizer3D(const Params & p = Params(), QWidget * parent = nullptr);

            const Params & params() const;
            Params & params();

        private:
            std::unique_ptr<VisualizerGUIData> _gui;
        };

    }
}
 
#endif