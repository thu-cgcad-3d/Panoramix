#ifndef CONFIGURATION_HPP
#define CONFIGURATION_HPP

#include "../../src/core/basic_types.hpp"
#include "../../src/core/cameras.hpp"
#include "../../src/core/mixed_graph.hpp"

#include "steps.hpp"

struct PanoView {
    panoramix::core::View<panoramix::core::PanoramicCamera> view;
};

StepWidgetInterface * CreateBindingWidgetAndActions(DataOfType<PanoView> & pv,
    QList<QAction*> & actions, QWidget * parent);




struct Segmentation {
    panoramix::core::Imagei segmentation;
    int segmentsNum;
};

StepWidgetInterface * CreateBindingWidgetAndActions(DataOfType<Segmentation> & segs,
    QList<QAction*> & actions, QWidget * parent);



struct LinesAndVPs {
    std::vector<panoramix::core::PerspectiveCamera> cams;
    std::vector<std::vector<panoramix::core::Classified<panoramix::core::Line2>>> lines;
    std::vector<panoramix::core::Line3> line3ds;
    std::vector<panoramix::core::Vec3> vps;
};

StepWidgetInterface * CreateBindingWidgetAndActions(DataOfType<LinesAndVPs> & segs,
    QList<QAction*> & actions, QWidget * parent);




struct ReconstructionSetup {
    panoramix::core::View<panoramix::core::PanoramicCamera> view;
    panoramix::core::MixedGraph mg;
    panoramix::core::MixedGraphPropertyTable props;
    panoramix::core::Imagei segmentation;
};

StepWidgetInterface * CreateBindingWidgetAndActions(DataOfType<ReconstructionSetup> & segs,
    QList<QAction*> & actions, QWidget * parent);





struct Reconstruction {
    panoramix::core::View<panoramix::core::PanoramicCamera> view;
    panoramix::core::MixedGraph mg;
    panoramix::core::MixedGraphPropertyTable props;
};

StepWidgetInterface * CreateBindingWidgetAndActions(DataOfType<Reconstruction> & segs,
    QList<QAction*> & actions, QWidget * parent);



struct ReconstructionRefinement {
    panoramix::core::View<panoramix::core::PanoramicCamera> view;
    panoramix::core::MixedGraph mg;
    panoramix::core::MixedGraphPropertyTable props;
    
};




#endif