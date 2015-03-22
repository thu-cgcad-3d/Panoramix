#ifndef CONFIGURATION_HPP
#define CONFIGURATION_HPP

#include "../../src/core/basic_types.hpp"
#include "../../src/core/cameras.hpp"
#include "../../src/core/mixed_graph.hpp"

#include "stepsdag.hpp"


StepWidgetInterface * CreateBindingWidgetAndActions(
    DataOfType<panoramix::core::PanoramicView> & pv,
    QList<QAction*> & actions, QWidget * parent);


struct Segmentation {
    panoramix::core::Imagei segmentation;
    int segmentsNum;
};

StepWidgetInterface * CreateBindingWidgetAndActions(
    DataOfType<Segmentation> & segs,
    QList<QAction*> & actions, QWidget * parent);


struct LinesAndVPs {
    std::vector<panoramix::core::PerspectiveView> perspectiveViews;
    std::vector<std::vector<panoramix::core::Classified<panoramix::core::Line2>>> lines;
    std::vector<panoramix::core::Line3> line3ds;
    std::vector<panoramix::core::Vec3> vps;
};

StepWidgetInterface * CreateBindingWidgetAndActions(
    DataOfType<LinesAndVPs> & segs,
    QList<QAction*> & actions, QWidget * parent);



template <class CameraT>
struct ReconstructionSetup {
    panoramix::core::View<CameraT> view;
    panoramix::core::MixedGraph mg;
    panoramix::core::MixedGraphPropertyTable props;
    panoramix::core::Imagei segmentation;
};

StepWidgetInterface * CreateBindingWidgetAndActions(
    DataOfType<ReconstructionSetup<panoramix::core::PanoramicCamera>> & rec,
    QList<QAction*> & actions, QWidget * parent);

StepWidgetInterface * CreateBindingWidgetAndActions(
    DataOfType<ReconstructionSetup<panoramix::core::PerspectiveCamera>> & rec,
    QList<QAction*> & actions, QWidget * parent);


template <class CameraT>
struct Reconstruction {
    panoramix::core::View<CameraT> view;
    panoramix::core::MixedGraph mg;
    panoramix::core::MixedGraphPropertyTable props;
};

StepWidgetInterface * CreateBindingWidgetAndActions(
    DataOfType<Reconstruction<panoramix::core::PanoramicCamera>> & rec,
    QList<QAction*> & actions, QWidget * parent);

StepWidgetInterface * CreateBindingWidgetAndActions(
    DataOfType<Reconstruction<panoramix::core::PerspectiveCamera>> & rec,
    QList<QAction*> & actions, QWidget * parent);



struct ReconstructionRefinement {
    panoramix::core::View<panoramix::core::PanoramicCamera> view;
    panoramix::core::MixedGraph mg;
    panoramix::core::MixedGraphPropertyTable props;
    
};




#endif