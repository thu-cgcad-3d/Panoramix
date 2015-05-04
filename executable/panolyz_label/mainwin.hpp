#pragma once

#include <functional>
#include <QtWidgets>

#include "../../src/core/basic_types.hpp"

using namespace panoramix;
using namespace panoramix::core;


class Widget : public QWidget {
    Q_OBJECT

public:
    Widget(QWidget *parent = 0);
    ~Widget(){}

    void loadImage(const QImage & im);
    void loadMAT(const QString & filename);

protected:
    void paintEvent(QPaintEvent * e) override;

public:
    Image image;
    Imagei segs;
    int nsegs;
    std::vector<Line2> lines;
    std::vector<int> lineClaz;
    std::vector<Point2> vps;
    Imagei segClaz;

private:
    double scale;
    QPoint translate;
    QPoint lastPos;
    QReadWriteLock lock;

    enum DisplayOption : uint32_t  {
        DisplayImage = 0x01,
        DisplaySegments = 0x02,
        DisplaySegmentLabels = 0x04,
        DisplayLines = 0x08,
        DisplayVPRays = 0x10,
        DisplayAll = DisplayImage | DisplaySegments | DisplaySegmentLabels | DisplayLines | DisplayVPRays
    };
    uint32_t displayOptions;
};


class MainWin : public QMainWindow {
	Q_OBJECT
	
public:
    MainWin(QWidget *parent = 0);
    ~MainWin(){}

private:
    Widget * w;
};
 