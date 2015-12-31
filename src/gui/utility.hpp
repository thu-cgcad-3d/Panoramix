#pragma once

#include "../core/cameras.hpp"

#include "color.hpp"

class QWidget;

namespace pano {
namespace gui {

class Scene;
class RenderOptions;

// spatial projected polygon for panorama reconstruction
struct SpatialProjectedPolygon {
  std::vector<Vec3> corners;
  Point3 projectionCenter;
  Plane3 plane;
};

int SelectFrom(const std::vector<std::string> &strs,
               const std::string &title = std::string(),
               const std::string &text = std::string(), int acceptId = -1,
               int rejectId = -1);

Image PickAnImage(const std::string &dir = std::string(),
                        std::string *picked = nullptr);

std::vector<Image> PickImages(const std::string &dir = std::string(),
                                    std::vector<std::string> *picked = nullptr);

std::vector<Image>
PickAllImagesFromAFolder(const std::string &dir = std::string(),
                         std::vector<std::string> *picked = nullptr);
void ForEachImageFromAFolder(
    const std::string &dir,
    const std::function<bool(const std::string &impath)> &fun);

bool MakePanoramaByHand(Image &im, bool *extendedOnTop = nullptr,
                        bool *extendedOnBottom = nullptr,
                        bool *topIsPlanar = nullptr,
                        bool *bottomIsPlanar = nullptr);

void PaintWith(
    const std::function<Image()> &updater,
    const std::vector<PenConfig> &penConfigs,
    const std::function<bool(const std::vector<Point2> &polyline,
                             int penId)> &callback);

void VisualizeWithPanoramicOperation(const Scene &scene,
                                     const RenderOptions &options);

void VisualizeAll(const View<PanoramicCamera, Image3ub> &view,
                  const std::vector<Classified<Line3>> &lines,
                  const Imagei &segs, int nsegs, const Image5d &gc);

void DrawChainsInPanorama(const PanoramicView &view,
                          const std::vector<PenConfig> &penConfigs,
                          std::vector<Chain3> &chains);
}

namespace core {
Box3 BoundingBox(const gui::SpatialProjectedPolygon &spp);
}
}
