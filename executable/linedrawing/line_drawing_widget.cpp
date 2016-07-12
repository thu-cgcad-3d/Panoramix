#include "pch.hpp"

#include "cache.hpp"
#include "gui_util.hpp"
#include "iterators.hpp"
#include "line_drawing_widget.hpp"
#include "qttools.hpp"
#include "utility.hpp"

namespace pano {
namespace experimental {
LineDrawingWidget::LineDrawingWidget(const Image &im, QWidget *parent)
    : QWidget(parent), image(gui::MakeQImage(im)), mode_(SelectOrMoveMode),
      last_touched_point_(-1), last_touched_edge_{-1, -1} {
  setMinimumSize(image.size());
  setMaximumSize(image.size());
  setCursor(Qt::PointingHandCursor);

  show_edges_ = true;
  show_coplanarities_ = true;
  show_colinearities_ = true;
}

static const double kPointRadius = 6.0;
static const double kLineWidth = 3.0;
static const gui::ColorTable kColorTable =
    gui::CreateRandomColorTableWithSize(1024);

void LineDrawingWidget::paintEvent(QPaintEvent *e) {
  QPainter painter(this);
  painter.setRenderHint(QPainter::Antialiasing, true);
  painter.setBackground(Qt::black);
  painter.eraseRect(rect());
  painter.drawImage(rect(), image);

  // draw points
  painter.setBrush(Qt::black);
  if (show_points_) {
    for (auto &p : projection.line_drawing.points) {
      painter.drawEllipse(gui::MakeQPointF(p), kPointRadius, kPointRadius);
    }
  }

  // draw edges
  painter.setPen(QPen(Qt::black, kLineWidth, Qt::SolidLine));
  if (show_edges_) {
    for (auto &e : projection.line_drawing.topo.edges) {
      painter.drawLine(
          gui::MakeQPointF(projection.line_drawing.points[e.first]),
          gui::MakeQPointF(projection.line_drawing.points[e.second]));
    }
  }

  // draw coplanar points
  painter.setPen(Qt::NoPen);
  if (show_coplanarities_) {
    for (int k = 0; k < projection.line_drawing.topo.coplanar_points.size();
         k++) {
      auto &ps = projection.line_drawing.topo.coplanar_points[k];
      assert(ps.size() >= 3);
      auto color = gui::MakeQColor(kColorTable[k % kColorTable.size()]);
      color.setAlpha(100);
      painter.setBrush(color);
      QPolygonF polygon(ps.size());
      for (int i = 0; i < ps.size(); i++) {
        polygon[i] = gui::MakeQPointF(projection.line_drawing.points[ps[i]]);
      }
      painter.drawPolygon(polygon);
    }
  }

  // draw colinear points
  painter.setBrush(Qt::NoBrush);
  if (show_colinearities_) {
    painter.setPen(QPen(Qt::blue, 2.0, Qt::DashLine));
    for (auto &ps : projection.line_drawing.topo.colinear_points) {
      assert(ps.size() >= 3);
      for (int i = 0; i < ps.size(); i++) {
        int p1 = ps[i];
        int p2 = ps[(i + 1) % ps.size()];
        painter.drawLine(gui::MakeQPointF(projection.line_drawing.points[p1]),
                         gui::MakeQPointF(projection.line_drawing.points[p2]));
      }
    }
  }

  // draw selected points
  painter.setBrush(Qt::NoBrush);
  for (int p : selected_points_) {
    painter.setPen(QPen(Qt::red, 3.0, Qt::DashLine));
    painter.drawEllipse(gui::MakeQPointF(projection.line_drawing.points[p]),
                        kPointRadius * 2, kPointRadius * 2);
  }
  // draw selected edges
  for (auto &e : selected_edges_) {
    painter.setPen(QPen(Qt::red, kLineWidth + 2, Qt::SolidLine));
    painter.drawLine(
        gui::MakeQPointF(projection.line_drawing.points[e.first]),
        gui::MakeQPointF(projection.line_drawing.points[e.second]));
  }
}

void LineDrawingWidget::mousePressEvent(QMouseEvent *e) {
  last_point_ = e->pos();
  last_touched_point_ = selectPoint(e->pos(), kPointRadius * 2);
  last_touched_edge_ = selectEdge(e->pos(), kLineWidth);

  if (e->buttons() & Qt::LeftButton) {
    if (mode_ == SelectOrMoveMode) {
      if (last_touched_point_ != -1) {
        selected_edges_.clear();
        if (e->modifiers() & Qt::ShiftModifier) {
          if (!Contains(selected_points_, last_touched_point_)) {
            selected_points_.push_back(last_touched_point_);
          }
        } else if (e->modifiers() & Qt::ControlModifier) {
          if (Contains(selected_points_, last_touched_point_)) {
            selected_points_.erase(std::remove(selected_points_.begin(),
                                               selected_points_.end(),
                                               last_touched_point_),
                                   selected_points_.end());
          } else {
            selected_points_.push_back(last_touched_point_);
          }
        } else {
          selected_points_ = {last_touched_point_};
        }
      } else if (last_touched_edge_.first != -1) {
        selected_points_.clear();
        if (e->modifiers() & Qt::ShiftModifier) {
          selected_edges_.insert(last_touched_edge_);
        } else if (e->modifiers() & Qt::ControlModifier) {
          if (Contains(selected_edges_, last_touched_edge_)) {
            selected_edges_.erase(last_touched_edge_);
          } else {
            selected_edges_.insert(last_touched_edge_);
          }
        } else {
          selected_edges_ = {last_touched_edge_};
        }
      }
    } else if (mode_ == AddMode) {
      if (last_touched_point_ == -1 && last_touched_edge_.first == -1) {
        // add a point
        projection.line_drawing.points.push_back(gui::MakeCorePoint(e->pos()));
        selected_points_ = {
            static_cast<int>(projection.line_drawing.points.size()) - 1};
      }
    }
  }

  update();
}

void LineDrawingWidget::mouseMoveEvent(QMouseEvent *e) {
  if (e->buttons() & Qt::LeftButton) {
    if (mode_ == SelectOrMoveMode) {
      for (int p : selected_points_) {
        projection.line_drawing.points[p] +=
            gui::MakeCoreVec(e->pos() - last_point_);
      }
    }
  }
  last_point_ = e->pos();
  update();
}

void LineDrawingWidget::mouseReleaseEvent(QMouseEvent *e) {}

void LineDrawingWidget::wheelEvent(QWheelEvent *e) {}

void LineDrawingWidget::keyPressEvent(QKeyEvent *e) {
  if (e->key() == Qt::Key_C) { // switch edit mode
    Println("C pressed");
    if (mode_ == AddMode) {
      mode_ = SelectOrMoveMode;
      setCursor(Qt::PointingHandCursor);
      core::Println("switched to select/move mode");
    } else if (mode_ == SelectOrMoveMode) {
      mode_ = AddMode;
      setCursor(Qt::ArrowCursor);
      core::Println("switched to add mode");
    }
  } else if (e->key() == Qt::Key_Escape) { // clear all selection
    core::Println("Escape pressed");
    selected_points_.clear();
    selected_edges_.clear();
  } else if (e->key() == Qt::Key_A) { // select all/none points
    core::Println("A pressed");
    if (selected_points_.size() < projection.line_drawing.points.size()) {
      selected_points_ = std::vector<int>(
          MakeIotaIterator<int>(0),
          MakeIotaIterator<int>(projection.line_drawing.points.size()));
    } else {
      selected_points_.clear();
    }
  } else if (e->key() == Qt::Key_E) { // make edges
    core::Println("E pressed");
    if (selected_points_.size() >= 2) {
      for (int i = 0; i + 1 < selected_points_.size(); i++) {
        auto e = MakeOrderedPair(selected_points_[i], selected_points_[i + 1]);
        if (e.first != e.second &&
            !Contains(projection.line_drawing.topo.edges, e)) {
          projection.line_drawing.topo.edges.push_back(e);
        } else {
          core::Println("two corners of this edge are same!");
        }
      }
      selected_points_.clear();
    } else {
      core::Println("two few points selected");
    }
  } else if (e->key() == Qt::Key_F) { // make coplanar constraints (faces)
    core::Println("F pressed");
    if (selected_points_.size() >= 3) {
      projection.line_drawing.topo.coplanar_points.push_back(
          std::move(selected_points_));
      selected_points_.clear();
    } else {
      core::Println("too few points selected");
    }
  } else if (e->key() == Qt::Key_L) { // make colinear constraints
    core::Println("L pressed");
    if (selected_points_.size() >= 3) {
      projection.line_drawing.topo.colinear_points.push_back(
          std::move(selected_points_));
      selected_points_.clear();
    } else {
      core::Println("too few points selected");
    }
  } else if (e->key() ==
             Qt::Key_X) { // create a intersection point from two selected edges
    core::Println("X pressed");
    if (selected_edges_.size() == 2) {
      auto it = selected_edges_.begin();
      Line2 line1(projection.line_drawing.points[it->first],
                  projection.line_drawing.points[it->second]);
      ++it;
      Line2 line2(projection.line_drawing.points[it->first],
                  projection.line_drawing.points[it->second]);
      Point2 inter =
          DistanceBetweenTwoLines(line1.ray(), line2.ray()).second.first;
      if (line1.first == line2.first || line1.first == line2.second ||
          line1.second == line2.first || line1.second == line2.second) {
        core::Println("the two edges share same end points and cannot create "
                      "intersection");
      } else if (Distance(inter, gui::MakeCorePoint(this->rect().center())) >
                 1000) {
        core::Println("the intersection lies too far away, are you sure the "
                      "edges you selected should intersect, or is it necessary "
                      "to create such an intersection?");
      } else {
        projection.line_drawing.points.push_back(inter);
        int inter_id = projection.line_drawing.points.size() - 1;
        auto it = selected_edges_.begin();
        projection.line_drawing.topo.colinear_points.push_back(
            {it->first, it->second, inter_id});
        ++it;
        projection.line_drawing.topo.colinear_points.push_back(
            {it->first, it->second, inter_id});
        selected_edges_.clear();
      }
    } else {
      core::Println("select exactly two edges to make intersection!");
    }
  } else if (e->key() == Qt::Key_Delete) { // delete selected points
    core::Println("Delete pressed");
    std::vector<Point2> new_points;
    new_points.reserve(projection.line_drawing.points.size());
    for (int i = 0; i < projection.line_drawing.points.size(); i++) {
      if (!Contains(selected_points_, i)) {
        new_points.push_back(projection.line_drawing.points[i]);
      }
    }
    projection.line_drawing.points = std::move(new_points);
    selected_points_.clear();
    selected_edges_.clear();
  } else if (e->key() == Qt::Key_S && (e->modifiers() & Qt::ControlModifier) &&
             (e->modifiers() & Qt::ShiftModifier)) { // ctrl + shift + s
    QString filename = QFileDialog::getSaveFileName(
        this, tr("Save current line drawing"), QString(), tr("*.cereal"));
    if (!filename.isEmpty()) {
      SaveToDisk(filename.toStdString(), projection);
    }
  } else if (e->key() == Qt::Key_0) {
    show_edges_ = !show_edges_;
  } else if (e->key() == Qt::Key_1) {
    show_points_ = !show_points_;
  } else if (e->key() == Qt::Key_2) {
    show_coplanarities_ = !show_coplanarities_;
  } else if (e->key() == Qt::Key_3) {
    show_colinearities_ = !show_colinearities_;
  }
  update();
}

int LineDrawingWidget::selectPoint(const QPointF &p, double dist_thres) const {
  double min_dist = dist_thres;
  int selected = -1;
  for (int i = 0; i < projection.line_drawing.points.size(); i++) {
    double dist =
        Distance(gui::MakeCorePoint(p), projection.line_drawing.points[i]);
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
  for (auto &f : projection.line_drawing.topo.coplanar_points) {
    assert(f.size() >= 3);
    for (int i = 0; i < f.size(); i++) {
      int p1 = f[i];
      int p2 = f[(i + 1) % f.size()];
      double dist = Distance(gui::MakeCorePoint(p),
                             Line2(projection.line_drawing.points[p1],
                                   projection.line_drawing.points[p2]));
      if (dist < min_dist) {
        selected = std::make_pair(p1, p2);
        min_dist = dist;
      }
    }
  }
  return MakeOrderedPair(selected);
}
}
}
