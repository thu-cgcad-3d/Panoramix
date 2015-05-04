#include "../../src/gui/singleton.hpp"
#include "../../src/gui/qt_glue.hpp"
#include "../../src/core/feature.hpp"
#include "../../src/misc/matlab_api.hpp"

#include "mainwin.hpp"


Widget::Widget(QWidget * parent) : QWidget(parent), scale(1.0) {
    setMouseTracking(true);
    setMinimumSize(500, 500);
    displayOptions = DisplayAll;
}

void Widget::loadImage(const QImage & im){
    lock.lockForWrite();
    image = gui::MakeCVMat(im);
    SegmentationExtractor segExtractor;
    segExtractor.params().algorithm = SegmentationExtractor::GraphCut;
    std::tie(segs, nsegs) = segExtractor(image, false);
    LineSegmentExtractor lineExtractor;
    lineExtractor.params().algorithm = LineSegmentExtractor::LSD;
    lines = lineExtractor(image);
    lineClaz = std::vector<int>(lines.size(), -1);
    // todo
    lock.unlock();
}

void Widget::loadMAT(const QString & filename){
    misc::MAT matfile(filename.toStdString(), misc::MAT::Update);
    misc::FromMXArray(matfile.var("im"), image);
    misc::FromMXArray(matfile.var("segs"), segs);
    nsegs = matfile.var("nsegs").scalar();

}

void Widget::paintEvent(QPaintEvent * e){
    QPainter painter(this);
    painter.setBackground(Qt::white);
    painter.eraseRect(rect());

    if (image.empty())
        return;

    lock.lockForRead();

    QImage im = gui::MakeQImage(image);
    painter.resetTransform();
    painter.translate(rect().center() + translate);
    painter.scale(scale, scale);
    painter.translate(-im.rect().center());
    painter.fillRect(im.rect().translated(QPoint(5, 5) / scale), QColor(35, 30, 30, 100));

    if (displayOptions & DisplayImage){
        painter.drawImage(im.rect(), im);
    }


    lock.unlock();
}



MainWin::MainWin(QWidget *parent) : QMainWindow(parent) {
   
    //// menus
    auto menuFile = menuBar()->addMenu(tr("&File"));
    auto menuView = menuBar()->addMenu(tr("&View"));
    auto menuSettings = menuBar()->addMenu(tr("&Settings"));
    auto menuTools = menuBar()->addMenu(tr("&Tools"));
    auto menuHelp = menuBar()->addMenu(tr("Help"));

    w = new Widget(this);
    setCentralWidget(w);

    // open image
    auto actLoadImage = menuFile->addAction(tr("Load Image"));
    actLoadImage->setShortcut(QKeySequence::Open);
    connect(actLoadImage, &QAction::triggered, [this](){
        QString filename = QFileDialog::getOpenFileName(this, tr("Select an Image"),
            tr(PROJECT_TEST_DATA_DIR_STR),
            tr("Image Files (*.png;*.jpg;*.jpeg);;All Files (*.*)"));
        if (filename.isEmpty())
            return;
        QImage im(filename);
        w->loadImage(im);
        update();
    });

    auto actLoadMat = menuFile->addAction(tr("Load MAT"));
    connect(actLoadMat, &QAction::triggered, [this](){
        QString filename = QFileDialog::getOpenFileName(this, tr("Select a MAT file"),
            tr(PROJECT_TEST_DATA_DIR_STR),
            tr("MAT File (*.mat);;All Files (*.*)"));
        if (filename.isEmpty())
            return;
        w->loadMAT(filename);
        update();
    });

}

