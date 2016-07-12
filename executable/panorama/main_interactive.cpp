#include <QtWidgets>

#include "parallel.hpp"
#include "cache.hpp"
#include "clock.hpp"

#include "canvas.hpp"
#include "qttools.hpp"
#include "singleton.hpp"
#include "gui_util.hpp"

#include "pi_graph_annotation.hpp"
#include "pi_graph_cg.hpp"
#include "pi_graph_control.hpp"
#include "pi_graph_occlusion.hpp"
#include "pi_graph_optimize.hpp"
#include "pi_graph_vis.hpp"

#define DISABLED_main MACRO_CONCAT(main_, __COUNTER__)

using namespace pano;
using namespace pano::core;
using namespace pano::experimental;

void DrawOn(const Image &im, std::vector<Classified<Line2>> &lines) {
  class Widget : public QWidget {
  public:
    Widget(const Image &im, std::vector<Classified<Line2>> &lines,
           QWidget *p = nullptr)
        : QWidget(p), _im(gui::MakeQImage(im)), _lines(lines),
          _isDrawing(false) {
      setMouseTracking(true);
      setMinimumSize(_im.size());
      setMaximumSize(_im.size());
      estimateVPs();
    }

  protected:
    virtual void paintEvent(QPaintEvent *e) override {
      QPainter painter(this);
      painter.setRenderHint(QPainter::Antialiasing);
      painter.drawImage(QPoint(), _im);

      if (!_vpFailed) {
        static const QColor _colors[] = {Qt::red, Qt::green, Qt::blue};
        for (int i = 0; i < _lines.size(); i++) {
          int claz = _lines[i].claz;
          if (claz == -1 || claz >= 3) {
            continue;
          }
          auto p1 = ToPixel(_vps[claz].value());
          auto p2 = ToPixel(_lines[i].component.center());
          cv::clipLine(cv::Rect(Point2i(), Point2i(width(), height())), p1, p2);
          painter.setPen(QPen(_colors[claz], 1, Qt::DotLine));
          painter.drawLine(gui::MakeQPoint(p1), gui::MakeQPoint(p2));
        }
      }

      painter.setPen(QPen(Qt::black, 2));
      std::vector<QLineF> lines(_lines.size());
      for (int i = 0; i < _lines.size(); i++) {
        lines[i] = gui::MakeQLineF(_lines[i].component);
      }
      painter.drawLines(lines.data(), lines.size());
      if (_isDrawing) {
        painter.setPen(QPen(Qt::black, 2, Qt::DashLine));
        painter.drawLine(_lastP, _curP);
      }
    }
    virtual void mousePressEvent(QMouseEvent *e) override {
      _curP = e->pos();
      if (_isDrawing) {
        _lines.emplace_back(ClassifyAs(
            Line2(gui::MakeCorePoint(_lastP), gui::MakeCorePoint(_curP)), -1));
        estimateVPs();
        _isDrawing = false;
      } else {
        _lastP = e->pos();
        _isDrawing = true;
      }
      update();
    }
    virtual void mouseMoveEvent(QMouseEvent *e) override {
      _curP = e->pos();
      update();
    }
    virtual void mouseReleaseEvent(QMouseEvent *e) override {}

  private:
    void estimateVPs() {
      if (_lines.size() <= 3) {
        return;
      }
      auto result = _vpd(_lines, gui::MakeCoreSize(_im.size()));
      _vpFailed = result.failed();
      if (!_vpFailed) {
        std::tie(_vps, _focal) = result.unwrap();
      }
    }

  private:
    QImage _im;
    std::vector<Classified<Line2>> &_lines;
    QPointF _lastP;
    QPointF _curP;
    bool _isDrawing;

    bool _vpFailed;
    std::vector<HPoint2> _vps;
    double _focal;
    VanishingPointsDetector _vpd;
  };

  auto temp = im.clone();
  Widget w(temp, lines);
  w.show();
  gui::Singleton::ContinueGui();
}

int main_interactive(int argc, char **argv) {
  gui::Singleton::InitGui(argc, argv);
  misc::SetCachePath("D:\\PanoramaReconstruction\\Interactive\\");
  std::string path;
  auto im = gui::FileDialog::PickAnImage(PANORAMIX_TEST_DATA_DIR_STR "\\house\\", &path);
  if (im.empty()) {
    return 0;
  }

  ResizeToArea(im, 500 * 500);

  std::vector<Classified<Line2>> lines;
  DrawOn(im, lines);
  misc::SaveCache(path, "im_lines", im, lines);

  return 0;
}

int DISABLED_main(int argc, char **argv) {
  gui::Singleton::InitGui(argc, argv);
  misc::SetCachePath("D:\\PanoramaReconstruction\\Interactive\\");
  std::string path;
  gui::FileDialog::PickAnImage(PANORAMIX_TEST_DATA_DIR_STR "\\house\\", &path);
  if (path.empty()) {
    return 0;
  }

  core::Image im;
  std::vector<Classified<Line2>> lines;
  misc::LoadCache(path, "im_lines", im, lines);

  DrawOn(im, lines);

  // split regions by lines

  return 0;
}