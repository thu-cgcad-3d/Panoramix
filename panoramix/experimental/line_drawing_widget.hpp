#pragma once

#include <QtOpenGL>
#include <QtWidgets>
#include <functional>

#include "line_drawing.hpp"

namespace pano {
namespace experimental {

struct LineDrawingAnnotation {
  std::vector<Point2> points;
  std::vector<std::vector<int>> coplanar_points;
  std::vector<std::vector<int>> colinear_points;
  template <class ArchiveT>
  void serialize(ArchiveT & ar) {
    ar(points, coplanar_points, colinear_points);
  }
};

// LineDrawingWidget
class LineDrawingWidget : public QWidget {
public:
  explicit LineDrawingWidget(const Image &im, QWidget *parent = nullptr);

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
  LineDrawingAnnotation anno;

private:
  QPointF last_point_;
  int last_touched_point_;
  std::pair<int, int> last_touched_edge_;
  std::vector<int> selected_points_;
  std::set<std::pair<int, int>> selected_edges_;
  enum Mode { AddMode, SelectOrMoveMode } mode_;

  bool show_points_;
  bool show_edges_;
  bool show_coplanarities_;
  bool show_colinearities_;
};
}
}
