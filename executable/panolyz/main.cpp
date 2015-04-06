#include "../../src/misc/cmd_tools.hpp"
#include "routines.hpp"

using namespace panoramix;
using namespace panolyz;

int main(int argc, char ** argv) {
    
    std::string defaultFileName;
    //defaultFileName = PROJECT_TEST_DATA_DIR_STR"/panorama/indoor/13.jpg";
    //defaultFileName = PROJECT_TEST_DATA_DIR_STR"/panorama/indoor/14.jpg";
    //defaultFileName = PROJECT_TEST_DATA_DIR_STR"/panorama/indoor/x3.jpg";
    defaultFileName = PROJECT_TEST_DATA_DIR_STR"/panorama/indoor/45.jpg";
    //defaultFileName = PROJECT_TEST_DATA_DIR_STR"/panorama/indoor/x2.jpg";
    //defaultFileName = PROJECT_TEST_DATA_DIR_STR"/panorama/outdoor/univ1.jpg";
    //defaultFileName = PROJECT_TEST_DATA_DIR_STR"/panorama/indoor/k (9).jpg";// too small
    //defaultFileName = PROJECT_TEST_DATA_DIR_STR"/panorama/outdoor/yard.jpg";
    //defaultFileName = PROJECT_TEST_DATA_DIR_STR"/panorama/indoor/k (11).jpg"; // too small
    //efaultFileName = PROJECT_TEST_DATA_DIR_STR"/panorama/indoor/k (10).jpg"; // too small
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


    //defaultFileName = PROJECT_TEST_DATA_DIR_STR"/normal/room.png"; // bingo
    //defaultFileName = PROJECT_TEST_DATA_DIR_STR"/normal/room4.jpg"; // bingo
    //defaultFileName = PROJECT_TEST_DATA_DIR_STR"/normal/room7.jpg"; // bingo!!!
    //defaultFileName = PROJECT_TEST_DATA_DIR_STR"/normal/room10.jpg"; // bingo!
    //defaultFileName = PROJECT_TEST_DATA_DIR_STR"/normal/room12.png"; // bingo!
    //defaultFileName = PROJECT_TEST_DATA_DIR_STR"/normal/room13.jpg"; // bingo
    //defaultFileName = PROJECT_TEST_DATA_DIR_STR"/normal/room14.jpg"; // bingo
    //defaultFileName = PROJECT_TEST_DATA_DIR_STR"/normal/room21.jpg"; // bingo! ?
    //defaultFileName = PROJECT_TEST_DATA_DIR_STR"/normal/room22.jpg"; // bingo
    //defaultFileName = PROJECT_TEST_DATA_DIR_STR"/normal/room23.jpg"; // bingo



    //defaultFileName = PROJECT_TEST_DATA_DIR_STR"/normal/ng1.png";
    //defaultFileName = PROJECT_TEST_DATA_DIR_STR"/normal/room2e.jpg";
    //defaultFileName = PROJECT_TEST_DATA_DIR_STR"/normal/room3.jpg";
    //defaultFileName = PROJECT_TEST_DATA_DIR_STR"/normal/room5.jpg"; //!
    //defaultFileName = PROJECT_TEST_DATA_DIR_STR"/normal/room6.jpg";
    //defaultFileName = PROJECT_TEST_DATA_DIR_STR"/normal/room8.jpg"; //!
    //defaultFileName = PROJECT_TEST_DATA_DIR_STR"/normal/room11.jpg"; // vp failed
    //defaultFileName = PROJECT_TEST_DATA_DIR_STR"/normal/room15.jpg"; // vp failed
    //defaultFileName = PROJECT_TEST_DATA_DIR_STR"/normal/room16.jpg";
    //defaultFileName = PROJECT_TEST_DATA_DIR_STR"/normal/room17.jpg";
    //defaultFileName = PROJECT_TEST_DATA_DIR_STR"/normal/room18.jpg";
    //defaultFileName = PROJECT_TEST_DATA_DIR_STR"/normal/room19.jpg"; // flat...
    //defaultFileName = PROJECT_TEST_DATA_DIR_STR"/normal/room20.jpg";

    misc::CmdOptions cmdOptions = {
        {"f", defaultFileName, "input image file path"},
        {"p", true, "whether the input image is a panorama"},
        {"i", true, "whether the scene is indoor"}
    };
    
    if (!cmdOptions.parseArguments(argc, argv)){
        return 0;
    }

    core::Image image = cv::imread(cmdOptions.value<std::string>("f"));
    core::ResizeToHeight(image, cmdOptions.value<bool>("p") ? 700 : 400);

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