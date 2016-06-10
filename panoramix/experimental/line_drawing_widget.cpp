#include "line_drawing_widget.hpp"
#include "qttools.hpp"

namespace pano {
namespace experimental {
LineDrawingWidget::LineDrawingWidget(const QImage &im, QWidget *parent)
    : QWidget(parent), image(im) {
  setMinimumSize(im.size());
  setMaximumSize(im.size());
}

void LineDrawingWidget::paintEvent(QPaintEvent *e) {
  QPainter painter(this);
  painter.drawImage(rect(), image);
  // draw faces loops
  painter.setPen(QPen(Qt::black, 2.0, Qt::SolidLine));
  for (auto &f : faces) {
    if (f.size() < 2) {
      continue;
    }
    for (int i = 0; i < f.size(); i++) {
      int p1 = f[i];
      int p2 = f[(i + 1) % f.size()];
      painter.drawLine(gui::MakeQPointF(points[p1]),
                       gui::MakeQPointF(points[p2]));
    }
  }
  // draw points
  for (auto &p : points) {
    painter.drawEllipse(gui::MakeQPointF(p), 4.0, 4.0);
  }
  // draw colinear points
  painter.setPen(QPen(Qt::red, 1.0, Qt::DashLine));
  for (auto &ps : colinear_points) {
    for (int i = 0; i + 1 < ps.size(); i++) {
      int p1 = ps[i];
      int p2 = ps[i + 1];
      painter.drawLine(gui::MakeQPointF(points[p1]),
                       gui::MakeQPointF(points[p2]));
    }
  }
  // draw coplanar points
  painter.setPen(QPen(Qt::blue, 0.6, Qt::DotLine));
  for (auto &ps : coplanar_points) {
    if (ps.size() < 2) {
      continue;
    }
    for (int i = 0; i < ps.size(); i++) {
      int p1 = ps[i];
      int p2 = ps[(i + 1) % f.size()];
      painter.drawLine(gui::MakeQPointF(points[p1]),
                       gui::MakeQPointF(points[p2]));
    }
  }
}

void LineDrawingWidget::mousePressEvent(QMouseEvent *e) {
  
}

void LineDrawingWidget::mouseMoveEvent(QMouseEvent *e) {}

void LineDrawingWidget::mouseReleaseEvent(QMouseEvent *e) {}

void LineDrawingWidget::wheelEvent(QWheelEvent *e) {}

void LineDrawingWidget::keyPressEvent(QKeyEvent *e) {}

int LineDrawingWidget::selectPoint(const QPointF &p, double dist_thres) const {
  double min_dist = dist_thres;
  int selected = -1;
  for (int i = 0; i < points.size(); i++) {
    double dist = Distance(gui::MakeCorePoint(p), points[i]);
    if (dist < min_dist) {
      selected = i;
      min_dist = dist;
    }
  }
  return selected;
}

std::pair<int, int> LineDrawingWidget::selectEdge(const QPointF &p,
                                                  double dist_thres) const {
  std::pair<int, int> selected{-1, -1};
  double min_dist = dist_thres;
  for (auto &f : faces) {
    if (f.size() < 2) {
      continue;
    }
    for (int i = 0; i < f.size(); i++) {
      int p1 = f[i];
      int p2 = f[(i + 1) % f.size()];
      double dist =
          Distance(gui::MakeCorePoint(p), Line2(points[p1], points[p2]));
      if (dist < min_dist) {
        selected = std::make_pair(p1, p2);
        min_dist = dist;
      }
    }
  }
  return selected;
}
}
}
