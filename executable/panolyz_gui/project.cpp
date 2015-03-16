#include <QtGui>
#include <QtWidgets>

#include "../../src/core/basic_types.hpp"
#include "../../src/core/mixed_graph.hpp"
#include "../../src/vis/qt_glue.hpp"
#include "../../src/vis/project.hpp"

#include "widgets.hpp"
#include "project.hpp"

namespace px = panoramix;

namespace data {

    struct Input {
        QImage qimage;
        px::core::Image image;
    };

    struct Segmentation {
        QImage qimage;
        px::core::Imagei segmentation;
    };

    struct Lines {

    };

}



Project::Project(QObject * parent) : QObject(parent) {
}

Project::Project(const QString & image, bool isPano, QObject * parent) : QObject(parent) {
    initialize(image, isPano);
}

Project::Project(const QString & projFile, QObject * parent) : QObject(parent) {
    loadFromDisk(projFile);
}

Project::~Project() {}

void Project::initialize(const QString & image, bool isPano) {
    _projectFileInfo = QFileInfo(tr("Unnamed"));
    _imageFileInfo = QFileInfo(image);
    _isPanorama = isPano;

    _core = std::make_unique<px::vis::ProjectCore>();

    // append steps
    int stepLoadImage = _core->addStep("Load Input", [this](){
        QImage im;
        im.load(_imageFileInfo.absoluteFilePath());
        auto pim = px::vis::MakeCVMat(im);
        return pim;
    });

    
    
    if (isPano){
        //int stepPerspSampling = 
    }

}


void Project::saveToDisk(const QString & filename) const {

}

void Project::loadFromDisk(const QString & filename){
    _projectFileInfo = QFileInfo(filename);
    
}




void Project::update() {

    _core->updateAll([this](int stepId){
        

        
    });

}



//struct InputImageItem : ProjectItem {
//    virtual QString name() { return fileInfo.fileName(); }
//    virtual QWidget * popWidget() { NOT_IMPLEMENTED_YET(); }
//    virtual void update(const QList<ProjectItem*> & dependencies) {
//        qimage.load(fileInfo.absoluteFilePath());
//        image = px::vis::MakeCVMat(qimage, true);
//    }
//
//    bool isPanorama;
//    QFileInfo fileInfo;
//    QImage qimage;
//
//    px::core::Image image;
//
//
//};
//
//struct LineSegmentsItem : ProjectItem {
//
//};
//
//struct ImageSegmentationItem : ProjectItem {
//    virtual QString name() { return QObject::tr("Segmented Image"); }
//    virtual QWidget * popWidget() { NOT_IMPLEMENTED_YET(); }
//    virtual void update(const QList<ProjectItem*> & dependencies) {
//        Q_ASSERT(dependencies.size() == 1 && IsA(dependencies.first(), InputImageItem*));
//        auto inputImageItem = dynamic_cast<InputImageItem*>(dependencies.first());
//        //segmentedImage = segmenter()
//        NOT_IMPLEMENTED_YET();
//    }
//
//    px::core::SegmentationExtractor::Feature segmentedImage;
//    px::core::SegmentationExtractor segmenter;
//};