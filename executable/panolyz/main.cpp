#include "../../src/misc/cmd_tools.hpp"
#include "routines.hpp"

using namespace panoramix;
using namespace panolyz;

int main(int argc, char ** argv) {
    
    std::string defaultFileName;
    //defaultFileName = PROJECT_TEST_DATA_DIR_STR"/panorama/indoor/13.jpg";
    //defaultFileName = PROJECT_TEST_DATA_DIR_STR"/panorama/indoor/14.jpg";
    //defaultFileName = PROJECT_TEST_DATA_DIR_STR"/panorama/indoor/x3.jpg";
    //defaultFileName = PROJECT_TEST_DATA_DIR_STR"/panorama/indoor/45.jpg";
    //defaultFileName = PROJECT_TEST_DATA_DIR_STR"/panorama/indoor/x2.jpg";
    //defaultFileName = PROJECT_TEST_DATA_DIR_STR"/panorama/outdoor/univ1.jpg";
    //defaultFileName = PROJECT_TEST_DATA_DIR_STR"/panorama/indoor/k (9).jpg";
    //defaultFileName = PROJECT_TEST_DATA_DIR_STR"/panorama/outdoor/yard.jpg";
    //defaultFileName = PROJECT_TEST_DATA_DIR_STR"/panorama/indoor/k (11).jpg";
    //defaultFileName = PROJECT_TEST_DATA_DIR_STR"/panorama/indoor/k (10).jpg";
    //defaultFileName = PROJECT_TEST_DATA_DIR_STR"/panorama/indoor/k (7).jpg";
    //defaultFileName = PROJECT_TEST_DATA_DIR_STR"/panorama/indoor/univlab.jpg";
    //defaultFileName = PROJECT_TEST_DATA_DIR_STR"/panorama/indoor/univlab2.jpg";
    //defaultFileName = PROJECT_TEST_DATA_DIR_STR"/panorama/indoor/univlab3.jpg";
    //defaultFileName = PROJECT_TEST_DATA_DIR_STR"/panorama/indoor/k (2).jpg";
    //defaultFileName = PROJECT_TEST_DATA_DIR_STR"/panorama/indoor/x5.jpg";
    //defaultFileName = PROJECT_TEST_DATA_DIR_STR"/panorama/indoor/x6.jpg";
    //defaultFileName = PROJECT_TEST_DATA_DIR_STR"/panorama/indoor/x7.jpg";
    //defaultFileName = PROJECT_TEST_DATA_DIR_STR"/panorama/indoor/x8.jpg";
    //defaultFileName = PROJECT_TEST_DATA_DIR_STR"/panorama/indoor/x9.jpg";
    //defaultFileName = PROJECT_TEST_DATA_DIR_STR"/panorama/indoor/google_chinese.png";

    //defaultFileName = PROJECT_TEST_DATA_DIR_STR"/normal/ng1.png";
    //defaultFileName = PROJECT_TEST_DATA_DIR_STR"/normal/room.png";
    //defaultFileName = PROJECT_TEST_DATA_DIR_STR"/normal/room1.jpg";
    //defaultFileName = PROJECT_TEST_DATA_DIR_STR"/normal/room2e.jpg";
    //defaultFileName = PROJECT_TEST_DATA_DIR_STR"/normal/room3.jpg";
    //defaultFileName = PROJECT_TEST_DATA_DIR_STR"/normal/room4.jpg";
    //defaultFileName = PROJECT_TEST_DATA_DIR_STR"/normal/room5.jpg"; //!
    //defaultFileName = PROJECT_TEST_DATA_DIR_STR"/normal/room6.jpg";
    //defaultFileName = PROJECT_TEST_DATA_DIR_STR"/normal/room7.jpg";
    //defaultFileName = PROJECT_TEST_DATA_DIR_STR"/normal/room8.jpg";
    //defaultFileName = PROJECT_TEST_DATA_DIR_STR"/normal/room9.jpg";
    //defaultFileName = PROJECT_TEST_DATA_DIR_STR"/normal/room10.jpg";
    //defaultFileName = PROJECT_TEST_DATA_DIR_STR"/normal/room11.jpg";
    //defaultFileName = PROJECT_TEST_DATA_DIR_STR"/normal/room12.png";
    //defaultFileName = PROJECT_TEST_DATA_DIR_STR"/normal/room13.jpg";
    //defaultFileName = PROJECT_TEST_DATA_DIR_STR"/normal/room14.jpg";
    //defaultFileName = PROJECT_TEST_DATA_DIR_STR"/normal/room15.jpg";
    //defaultFileName = PROJECT_TEST_DATA_DIR_STR"/normal/room16.jpg";
    //defaultFileName = PROJECT_TEST_DATA_DIR_STR"/normal/room17.jpg";
    //defaultFileName = PROJECT_TEST_DATA_DIR_STR"/normal/room18.jpg";
    defaultFileName = PROJECT_TEST_DATA_DIR_STR"/normal/room19.jpg";
    //defaultFileName = PROJECT_TEST_DATA_DIR_STR"/normal/room20.jpg";
    //defaultFileName = PROJECT_TEST_DATA_DIR_STR"/normal/room21.jpg";
    //defaultFileName = PROJECT_TEST_DATA_DIR_STR"/normal/room22.jpg";

    misc::CmdOptions cmdOptions = {
        {"f", defaultFileName, "input image file path"},
        {"p", false, "whether the input image is a panorama"},
        {"i", true, "whether the scene is indoor"}
    };
    
    if (!cmdOptions.parseArguments(argc, argv)){
        return 0;
    }

    core::Image image = cv::imread(cmdOptions.value<std::string>("f"));
    core::ResizeToMakeHeightUnder(image, cmdOptions.value<bool>("p") ? 900 : 400);

    if (cmdOptions.value<bool>("i")){
        if (cmdOptions.value<bool>("p")){
            Routine<PanoramaIndoor_v1>(image);
        }
        else {
            Routine<NormalIndoor_v1>(image);
        }
    }
    else{
        if (cmdOptions.value<bool>("p")){
            Routine<PanoramaOutdoor_v1>(image);
        }
        else {
            Routine<NormalOutdoor_v1>(image);
        }
    }

    return 0;
}