#include <QApplication>
#include <QFileDialog>
#include <QtCore>

#include "eval.hpp"

int main(int argc, char ** argv) {
    
    pano::gui::Singleton::InitGui();

    pano::misc::Matlab matlab;
    
    std::vector<std::string> impaths;
    pano::gui::PickImages("H:\\DataSet\\PanoContext\\bedroom\\", &impaths);
    auto directions = panolyz::FibonacciDirections(1000);

    for (auto & impath : impaths) {
        auto anno = pano::experimental::LoadOrInitializeNewLayoutAnnotation(impath);

        auto modelGT = panolyz::PredictionOfGT(impath);
        //modelGT->visualize(directions);
        
      /*  auto modelPanoContext = panolyz::PredictionOfPanoContext(impath, matlab);
        modelPanoContext->visualize(directions);*/

        panolyz::PredictOptions options;

      /*  options.useGroundTruthOcclusions = false;
        auto modelPxWithoutGTOcc = panolyz::PredictionOfPanoramix(impath, options, matlab);
        modelPxWithoutGTOcc->visualize();*/

        options.useGroundTruthOcclusions = true;
        auto modelPxWithGTOcc = panolyz::PredictionOfPanoramix(impath, options, matlab);
        modelPxWithGTOcc->visualize();

     /*   std::cout << "average variation of PanoContext is "
            << panolyz::AverageVariation(*modelGT, *modelPanoContext, directions)
            << std::endl;
        std::cout << "average variation of Panoramix (without groundtruth occlusion) is "
            << panolyz::AverageVariation(*modelGT, *modelPxWithoutGTOcc, directions)
            << std::endl;*/
        std::cout << "average variation of Panoramix (with groundtruth occlusion) is "
            << panolyz::AverageVariation(*modelGT, *modelPxWithGTOcc, directions)
            << std::endl;
    }

    return 0;

}