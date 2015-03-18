#ifndef CONFIGURATION_HPP
#define CONFIGURATION_HPP

#include "../../src/core/basic_types.hpp"
#include "../../src/core/cameras.hpp"
#include "../../src/core/mixed_graph.hpp"

#include "steps.hpp"

struct PanoView {
    panoramix::core::View<panoramix::core::PanoramicCamera> view;
};

WidgetLike * CreateBindingWidgetAndActions(DataOfType<PanoView> & pv,
    QList<QAction*> & actions, QWidget * parent);




struct Segmentation {
    panoramix::core::Imagei segmentation;
    int segmentsNum;
};

WidgetLike * CreateBindingWidgetAndActions(DataOfType<Segmentation> & segs,
    QList<QAction*> & actions, QWidget * parent);





//struct RegionHints {
//    enum RegionLabel {
//        TowardVP1, TowardVP2, TowardVP3,
//        AlongVP1, AlongVP2, AlongVP3,
//        Planar, Clutter, Void
//    };
//    enum BoundaryLabel {
//        Connected,
//        Disconnected
//    };
//    std::vector<RegionLabel> regionIdToLabels;
//    std::vector<panoramix::core::Classified<std::vector<panoramix::core::Point2>, BoundaryLabel>> labeledBoundaries;
//    panoramix::core::View<panoramix::core::PanoramicCamera> view;
//    panoramix::core::Imagei segmentation;
//    std::vector<panoramix::core::Vec3> vps;
//};
//
//WidgetLike * CreateBindingWidgetAndActions(DataOfType<RegionHints> & pv,
//    QList<QAction*> & actions, QWidget * parent);





struct LinesAndVPs {
    std::vector<panoramix::core::PerspectiveCamera> cams;
    std::vector<std::vector<panoramix::core::Classified<panoramix::core::Line2>>> lines;
    std::vector<panoramix::core::Line3> line3ds;
    std::vector<panoramix::core::Vec3> vps;
};

WidgetLike * CreateBindingWidgetAndActions(DataOfType<LinesAndVPs> & segs,
    QList<QAction*> & actions, QWidget * parent);




struct ReconstructionSetup {
    panoramix::core::View<panoramix::core::PanoramicCamera> view;
    panoramix::core::MixedGraph mg;
    panoramix::core::MixedGraphPropertyTable props;
    panoramix::core::Imagei segmentation;
};

WidgetLike * CreateBindingWidgetAndActions(DataOfType<ReconstructionSetup> & segs,
    QList<QAction*> & actions, QWidget * parent);





struct Reconstruction {
    panoramix::core::View<panoramix::core::PanoramicCamera> view;
    panoramix::core::MixedGraph mg;
    panoramix::core::MixedGraphPropertyTable props;
};

WidgetLike * CreateBindingWidgetAndActions(DataOfType<Reconstruction> & segs,
    QList<QAction*> & actions, QWidget * parent);







#endif