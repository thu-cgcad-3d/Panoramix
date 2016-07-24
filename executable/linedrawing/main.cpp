#include "cache.hpp"
#include "matlab_api.hpp"
#include "singleton.hpp"

#include "line_drawing_evaluation.hpp"
#include "line_drawing_reconstruction.hpp"
#include "line_drawing_tools.hpp"

using namespace pano;
using namespace pano::core;
using namespace pano::experimental;

int main(int argc, char **argv, char **env) {
  gui::Singleton::SetCmdArgs(argc, argv, env);
  gui::Singleton::InitGui(argc, argv);
  misc::SetCachePath(PANORAMIX_CACHE_DATA_DIR_STR "\\LineDrawing\\");

  enum TaskName { EditImageAnnotation, TrainEnergyWeights, RunReconstruction };

  TaskName task = RunReconstruction;

  if (task == EditImageAnnotation) { // edit annotation on image
    std::string im_file_name;
    Image3ub im = gui::FileDialog::PickAnImage(PANORAMIX_TEST_DATA_DIR_STR,
                                               &im_file_name);
    ResizeToMakeHeightUnder(im, 600);
    std::string folder = misc::FolderOfFile(im_file_name);
    std::string name = misc::NameOfFile(im_file_name);
    std::string anno_file = folder + "\\" + name + ".ldanno.cereal";
    while (true) {
      LineDrawingWidget w(im);
      LoadFromDisk(anno_file, w.projection);
      w.projection.image = im;
      w.show();

      gui::Singleton::ContinueGui();
      int choice = gui::SelectFrom({"save", "abandon"}, "make your choice",
                                   "your choice?", 0, 1);
      if (choice == 0) {
        SaveToDisk(anno_file, w.projection);
      }
    }
  } else if (task == TrainEnergyWeights) {

    std::default_random_engine rng;
    std::vector<std::string> model_names = {
        "hex",       "triangle", "twotriangles", "stool", "plane",
        "towers",    "tower",    "towerx",       "car",   "gate",
        "smallgate", "castle",   "desk"};

    std::vector<std::vector<double>> feature_vecs;
    std::vector<double> losses;

    for (auto &model_name : model_names) {
      Println("model_name = ", model_name);
      std::string folder = PANORAMIX_TEST_DATA_DIR_STR "\\linedrawing\\objs\\" +
                           model_name + "\\";
      auto line_drawing =
          LineDrawing3FromObjFile(folder + model_name + "_w_intf.obj");

      // randomly generate cameras
      ForEachSampledCamera(
          6,
          {DegreesToRadians(60), DegreesToRadians(90), DegreesToRadians(120)},
          3, {0.1}, BoundingBox(line_drawing), Sizei(500, 500),
          [&](const PerspectiveCamera &cam) {
            bool has_parallel_edge = false;
            for (auto &e : line_drawing.topo.edges) {
              if (IsFuzzyPerpendicular(cam.forward(),
                                       normalize(line_drawing.points[e.first] -
                                                 line_drawing.points[e.second]),
                                       1e-8)) {
                has_parallel_edge = true;
                break;
              }
            }
            if (has_parallel_edge) {
              return;
            }

            if (false) {
              gui::SceneBuilder sb;
              sb.installingOptions().defaultShaderSource =
                  gui::OpenGLShaderSourceDescriptor::XLines;
              sb.add(line_drawing);
              sb.show(
                  true, false,
                  gui::RenderOptions().backgroundColor(gui::White).camera(cam));
            }

            LineDrawing2 projection;
            std::vector<double> groundtruth_depths;
            std::tie(projection, groundtruth_depths) =
                DecomposeProjectionAndDepths(line_drawing, cam);
            Normalize(groundtruth_depths.begin(), groundtruth_depths.end());
            auto aux = MakeAuxiliary(projection);

            auto face_sets = DecomposeFaces(line_drawing.topo.coplanar_points,
                                            projection.points);
            auto faces_overlap = ComputeFacesOverlap(projection, aux);
            auto feature_extractor = MakeLineDrawingFeatureExtractor(
                projection, aux, faces_overlap, face_sets, cam);

            auto groundtruth_feature_vec = feature_extractor(groundtruth_depths);
            feature_vecs.push_back(std::move(groundtruth_feature_vec));
            losses.push_back(0.0);

            // add noise to depths
            for (int k = 0; k < 10; k++) {
              std::vector<double> noised_depths = groundtruth_depths;
              std::normal_distribution<> dist(0.0, 0.05);
              for (double &d : noised_depths) {
                d += dist(rng);
                d = std::max(d, 0.0001);
              }
              Normalize(noised_depths.begin(), noised_depths.end());

              auto feature_vec = feature_extractor(noised_depths);
              double loss =
                  L2Distance(groundtruth_depths.begin(),
                             groundtruth_depths.end(), noised_depths.begin(),
                             noised_depths.end()) /
                  groundtruth_depths.size();
              feature_vecs.push_back(std::move(feature_vec));
              losses.push_back(loss);
            }
          });
    }

    Println("nsamples = ", feature_vecs.size());

    // save all feature vecs and losses as matlab file
    misc::MAT mat_file(PANORAMIX_TEST_DATA_DIR_STR
                       "\\linedrawing\\line_drawing_features_losses.mat",
                       misc::MAT::Write);
    DenseMatd features_mat(feature_vecs.size(), feature_vecs.front().size(), 0.0);
    for (int i = 0; i < feature_vecs.size(); i++) {
      auto & feature_vec = feature_vecs[i];
      for (int j = 0; j < feature_vec.size(); j++) {
        features_mat(i, j) = feature_vec[j];
      }
    }
    DenseMatd losses_vec(losses);
    mat_file.setVar("features", features_mat);
    mat_file.setVar("losses", losses_vec);

  } else if (task == RunReconstruction) {
    LineDrawingReconstructionSteps steps;
    steps.preprocess.rerun = false;
    steps.find_vps.rerun = steps.preprocess.rerun || false;
    steps.estimate_orientations.rerun =
        steps.find_vps.rerun || false;
    steps.reconstruction.rerun =
        steps.estimate_orientations.rerun || false;

    std::vector<LineDrawingInput> inputs = {
        LineDrawingInput::FromObjFile("hex", "cam1"),          //
        LineDrawingInput::FromObjFile("triangle", "cam1"),     //
        LineDrawingInput::FromObjFile("twotriangles", "cam1"), //
        LineDrawingInput::FromObjFile("stool", "cam1"),        //
        LineDrawingInput::FromObjFile("plane", "cam1"),        //
        LineDrawingInput::FromObjFile("towers", "cam1"),       // skewed a bit
        LineDrawingInput::FromObjFile("tower", "cam1"),        //
        LineDrawingInput::FromObjFile("towerx", "cam1"),       //
        LineDrawingInput::FromObjFile("car", "cam1"),          //
        LineDrawingInput::FromObjFile("gate", "cam1"),         // skewed a bit
        LineDrawingInput::FromObjFile("smallgate", "cam1"),    //
        LineDrawingInput::FromObjFile("castle", "cam1"),       //
        LineDrawingInput::FromObjFile("desk", "cam1")          //
    };

    for (auto &input : inputs) {
      auto result =
          RunLineDrawingReconstruction(input, steps, GUILevelNecessary);
      SaveLineDrawing3ToObjFile(input.folder + "\\" + input.model_name + "_" +
                                    input.cam_name + "_result.obj",
                                result);
    }
  }

  return 0;
}