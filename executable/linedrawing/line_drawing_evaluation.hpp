#pragma once

#include <functional>
#include <map>

#include "cameras.hpp"
#include "line_drawing.hpp"

namespace pano {
namespace experimental {
// MakeLineDrawingFeatureExtractor
using LineDrawingFeatureExtractor =
    std::function<std::vector<double>(const std::vector<double> &point_depths)>;
LineDrawingFeatureExtractor MakeLineDrawingFeatureExtractor(
    const LineDrawing2 &line_drawing, const AuxiliaryData<LineDrawing2> &aux,
    const std::map<std::pair<int, int>, bool> &faces_overlap,
    const std::vector<std::set<int>> &face_sets, const PerspectiveCamera &cam);

//std::vector<double> ExtractFeaturesFromLineDrawing()
}
}
