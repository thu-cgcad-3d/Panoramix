#pragma once

#include "../../src/core/basic_types.hpp"
#include "../../src/core/cameras.hpp"
#include "../../src/experimental/rl_graph.hpp"

#include "stepsdag.hpp"

using namespace panoramix;

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
    panoramix::experimental::RLGraph mg;
    panoramix::experimental::RLGraphControls controls;
    panoramix::core::Imagei segmentation;
    std::vector<panoramix::experimental::RegionHandle> segId2RegionHandles;
};

StepWidgetInterface * CreateBindingWidgetAndActions(
    DataOfType<ReconstructionSetup<panoramix::core::PanoramicCamera>> & rec,
    QList<QAction*> & actions, QWidget * parent);

StepWidgetInterface * CreateBindingWidgetAndActions(
    DataOfType<ReconstructionSetup<panoramix::core::PerspectiveCamera>> & rec,
    QList<QAction*> & actions, QWidget * parent);


struct RLGraphAndControl {
    panoramix::experimental::RLGraph mg;
    panoramix::experimental::RLGraphControls controls;
    panoramix::experimental::RLGraphVars vars;
};

template <class CameraT>
struct Reconstruction {
    panoramix::core::View<CameraT> view;
    std::vector<RLGraphAndControl> gs;
};

StepWidgetInterface * CreateBindingWidgetAndActions(
    DataOfType<Reconstruction<panoramix::core::PanoramicCamera>> & rec,
    QList<QAction*> & actions, QWidget * parent);

StepWidgetInterface * CreateBindingWidgetAndActions(
    DataOfType<Reconstruction<panoramix::core::PerspectiveCamera>> & rec,
    QList<QAction*> & actions, QWidget * parent);


template <class CameraT>
struct ReconstructionRefinement {
    panoramix::core::View<CameraT> view;
    std::vector<panoramix::core::LayeredShape3> shapes;
};

StepWidgetInterface * CreateBindingWidgetAndActions(
    DataOfType<ReconstructionRefinement<panoramix::core::PanoramicCamera>> & rec,
    QList<QAction*> & actions, QWidget * parent);

