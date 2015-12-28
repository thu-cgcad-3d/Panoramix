#include <QtWidgets>

#include "../../src/core/parallel.hpp"
#include "../../src/misc/cache.hpp"
#include "../../src/misc/clock.hpp"

#include "../../src/gui/canvas.hpp"
#include "../../src/gui/qttools.hpp"
#include "../../src/gui/scene.hpp"
#include "../../src/gui/singleton.hpp"
#include "../../src/gui/utility.hpp"

#include "../../src/experimental/pi_graph_annotation.hpp"
#include "../../src/experimental/pi_graph_cg.hpp"
#include "../../src/experimental/pi_graph_control.hpp"
#include "../../src/experimental/pi_graph_occlusion.hpp"
#include "../../src/experimental/pi_graph_optimize.hpp"
#include "../../src/experimental/pi_graph_vis.hpp"

using namespace pano;
using namespace pano::core;
using namespace pano::experimental;

#define CONCAT_IMPL(x, y) x##y
#define MACRO_CONCAT(x, y) CONCAT_IMPL(x, y)
#define DISABLED_main MACRO_CONCAT(main_, __COUNTER__)

template <class PointT> struct LineDrawing {
  std::vector<PointT> vertPositions;
  std::vector<std::pair<int, int>> line2verts;
  std::vector<std::vector<int>> face2verts;
  DenseMatd vertAdjMat;
};

LineDrawing<Point3> LoadLineDrawing(const std::string &filename,
                                    const std::string &gtfilename) {
  LineDrawing<Point3> drawing;

  std::ifstream ifs(filename);
  if (!ifs.is_open()) {
    return drawing;
  }

  int lineNum = 0;
  ifs >> lineNum;
  std::vector<Line2> lineMatrix(lineNum);
  for (int i = 0; i < lineNum; i++) {
    ifs >> lineMatrix[i].first[0] >> lineMatrix[i].first[1] >>
        lineMatrix[i].second[0] >> lineMatrix[i].second[1];
  }

  int vertexNum = 0;
  ifs >> vertexNum;
  drawing.vertAdjMat = DenseMat<bool>::zeros(vertexNum, vertexNum);
  for (int i = 0; i < vertexNum; i++) {
    for (int j = 0; j < vertexNum; j++) {
      ifs >> drawing.vertAdjMat(i, j);
    }
  }

  int faceNum = 0;
  ifs >> faceNum;
  drawing.face2verts.resize(faceNum);
  std::string line;
  std::getline(ifs, line);
  for (int i = 0; i < faceNum; i++) {
    std::getline(ifs, line);
    std::istringstream ss(line);
    int vnum = 0;
    ss >> vnum;
    if (vnum == 0) {
      continue;
    }
    int id = 0;
    for (int k = 0; k < vnum; k++) {
      ss >> id;
      drawing.face2verts[i].push_back(id);
    }
  }

  int dummy = 0;
  ifs >> dummy >> dummy;
  drawing.line2verts.resize(lineNum);
  for (int i = 0; i < lineNum; i++) {
    ifs >> drawing.line2verts[i].first >> drawing.line2verts[i].second;
  }

  ifs.close();
  ifs.open(gtfilename);
  if (!ifs.is_open()) {
    return drawing;
  }
  drawing.vertPositions.resize(vertexNum);
  for (int k = 0; k < 3; k++) {
    for (int i = 0; i < vertexNum; i++) {
      ifs >> drawing.vertPositions[i][k];
    }
  }

  return drawing;
}

template <class FunT>
LineDrawing<Point2> Transform(const LineDrawing<Point3> &in, const FunT &fun) {
  LineDrawing<Point2> out;
  out.vertPositions.resize(in.vertPositions.size());
  std::transform(in.vertPositions.begin(), in.vertPositions.end(),
                 out.vertPositions.begin(), fun);
  out.line2verts = in.line2verts;
  out.face2verts = in.face2verts;
  out.vertAdjMat = in.vertAdjMat.clone();
  return out;
}

int DISABLED_main(int argc, char **argv) {
  gui::Singleton::InitGui(argc, argv);
  misc::SetCachePath("D:\\Cache\\Panoramix\\LineDrawing\\");
  misc::Matlab matlab;

  std::string folder = "F:\\DataSets\\linedrawing_"
                       "dataset\\DatabaseFile\\bigHouse\\";
  auto drawing =
      LoadLineDrawing(folder + "bigHouse_cy.3dr", folder + "bigHouse.gt");

  gui::SceneBuilder sb;
  sb.installingOptions().defaultShaderSource =
      gui::OpenGLShaderSourceDescriptor::XLines;
  sb.installingOptions().discretizeOptions.color = gui::Black;
  sb.installingOptions().lineWidth = 3.0;
  double scale =
      BoundingBoxOfContainer(drawing.vertPositions).outerSphere().radius;
  for (auto &l : drawing.line2verts) {
    Line3 line3(drawing.vertPositions[l.first],
                drawing.vertPositions[l.second]);
    sb.add(line3 / scale);
  }
  sb.show(
      true, true,
      gui::RenderOptions().backgroundColor(gui::White).renderMode(gui::Lines));

  VanishingPointsDetector vpd;
  return 0;
}

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

int DISABLED_main(int argc, char **argv) {
  gui::Singleton::InitGui(argc, argv);
  misc::SetCachePath("D:\\Panoramix\\LineDrawing\\");
  std::string path;
  auto im =
      gui::PickAnImage(PANORAMIX_TEST_DATA_DIR_STR "\\localmanh\\", &path);
  ResizeToArea(im, 500 * 500);
  // Imageb edge;
  // cv::Canny(im, edge, 100, 300);
  // cv::imshow("im", im);
  // cv::imshow("edge", edge);
  // cv::waitKey();

  std::vector<Classified<Line2>> lines;
  DrawOn(im, lines);

  misc::SaveCache(path, "im_lines", im, lines);

  return 0;
}

int main(int argc, char **argv) {
  gui::Singleton::InitGui(argc, argv);
  misc::SetCachePath("D:\\Panoramix\\LineDrawing\\");
  std::string path;
  gui::PickAnImage(PANORAMIX_TEST_DATA_DIR_STR "\\localmanh\\", &path);
  core::Image im;
  std::vector<Classified<Line2>> lines;
  misc::LoadCache(path, "im_lines", im, lines);

  DrawOn(im, lines);

  return 0;
}