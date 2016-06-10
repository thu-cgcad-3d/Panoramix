#pragma once

#include <QtOpenGL>
#include <QtWidgets>
#include <functional>

#include "line_drawing.hpp"

namespace pano {
namespace experimental {



// LineDrawingWidget
class LineDrawingWidget : public QWidget {
public:
  explicit LineDrawingWidget(const QImage &im, QWidget *parent = nullptr);

protected:
  virtual void paintEvent(QPaintEvent *e) override;
  virtual void mousePressEvent(QMouseEvent *e) override;
  virtual void mouseMoveEvent(QMouseEvent *e) override;
  virtual void mouseReleaseEvent(QMouseEvent *e) override;
  virtual void wheelEvent(QWheelEvent *e) override;
  virtual void keyPressEvent(QKeyEvent *e) override;

private:
  int selectPoint(const QPointF & p, double dist_thres = 0.1) const;
  std::pair<int, int> selectEdge(const QPointF &p,
                                 double dist_thres = 0.1) const;

public:
  QImage image;
  std::vector<Point2> points;
  std::vector<std::vector<int>> faces;
  std::vector<std::vector<int>> coplanar_points;
  std::vector<std::vector<int>> colinear_points;
};
}
}
