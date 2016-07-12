#pragma once

#include "algorithms.hpp"
#include "cache.hpp"
#include "cameras.hpp"
#include "canvas.hpp"
#include "clock.hpp"
#include "line_drawing.hpp"
#include "line_drawing_widget.hpp"
#include "optimization.hpp"
#include "parallel.hpp"
#include "scene.hpp"
#include "single_view.hpp"
#include "singleton.hpp"
#include "step.hpp"

namespace pano {
namespace experimental {
// LineDrawingGroundTruth
struct LineDrawingGroundTruth {
  LineDrawing<Point3> line_drawing;
  PerspectiveCamera camera;
  template <class Archive> void serialize(Archive &ar) {
    ar(line_drawing, camera);
  }
};

// LineDrawingInput
class LineDrawingInput {
public:
  std::string model_name, cam_name;
  LineDrawingImageProjection projection;
  std::shared_ptr<LineDrawingGroundTruth> groundtruth;

public:
  template <class Archive> void serialize(Archive &ar) {
    ar(model_name, cam_name, projection, groundtruth);
  }

  // FromObjFile
  static LineDrawingInput FromObjFile(const std::string &model_name,
                                      const std::string &cam_name);
  // FromImageAnnotation
  static LineDrawingInput FromImageAnnotation(const std::string &image_name);
};

// LineDrawingReconstructionSteps
struct LineDrawingReconstructionSteps {
  misc::StepConfig preprocess = "preprocess", find_vps = "find_vps",
                   estimate_orientations = "estimate_orientations",
                   reconstruction = "reconstruction";
  template <class Archive> void serialize(Archive &ar) {
    ar(preprocess, find_vps, estimate_orientations, reconstruction);
  }
};

// GUILevel
enum GUILevel { GUILevelNone = 0, GUILevelNecessary = 1, GUILevelAll = 2 };

// RunLineDrawingReconstruction
void RunLineDrawingReconstruction(const LineDrawingInput &input,
                                  LineDrawingReconstructionSteps &steps,
                                  GUILevel gui_level = GUILevelAll);
}
}