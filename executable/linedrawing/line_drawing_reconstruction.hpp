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

// LineDrawingInput
struct LineDrawingInput {
  std::string model_name, cam_name;
  std::string folder;
  LineDrawingImageProjection projection;
  std::shared_ptr<LineDrawingReconstruction> groundtruth;
  template <class Archive> void serialize(Archive &ar) {
    ar(model_name, cam_name, folder, projection, groundtruth);
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
LineDrawing3 RunLineDrawingReconstruction(const LineDrawingInput &input,
                                          LineDrawingReconstructionSteps &steps,
                                          GUILevel gui_level = GUILevelAll);
}
}