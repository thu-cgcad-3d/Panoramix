#include "../../src/core/basic_types.hpp"
#include "../../src/core/containers.hpp"
#include "../../src/experimental/pi_graph.hpp"
#include "../../src/experimental/pi_graph_control.hpp"
#include "../../src/experimental/pi_graph_occlusion.hpp"
#include "../../src/experimental/pi_graph_solve.hpp"
#include "../../src/experimental/pi_graph_vis.hpp"
#include "../../src/experimental/pi_graph_annotation.hpp"
#include "../../src/experimental/tools.hpp"
#include "../../src/misc/matlab_api.hpp"
#include "../../src/gui/scene.hpp"
#include "../../src/gui/canvas.hpp"
#include "../../src/gui/utility.hpp"

#include "routines.hpp"


namespace panolyz {

    namespace PanoramaIndoorTest {

        static const bool refresh = true;

        void Run() {

            using namespace pano;
            using namespace core;
            using namespace experimental;

            misc::Matlab matlab;

            std::vector<std::string> paths;
            std::vector<Image> images = gui::PickImages("H:\\DataSet\\pi\\", &paths);
            assert(paths.size() == images.size());

            for (int k = 0; k < images.size(); k++) {
                // get annotation
                const std::string & path = paths[k];
                PILayoutAnnotation anno = LoadOrInitializeNewLayoutAnnotation(path);

            }

        }
    }
}