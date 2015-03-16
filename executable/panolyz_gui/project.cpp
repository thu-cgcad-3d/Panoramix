#include <QtGui>
#include <QtWidgets>

#include "../../src/core/basic_types.hpp"
#include "../../src/core/mixed_graph.hpp"
#include "../../src/vis/qt_glue.hpp"
#include "../../src/vis/project.hpp"

#include "widgets.hpp"
#include "project.hpp"

using namespace panoramix;

namespace data {

    struct PanoImage {
        core::View<core::PanoramicCamera> view;
       
    };

    struct Segmentation {
        core::Imagei segmentation;
        int segmentsNum;
    };

    struct LinesAndVPs {
        std::vector<core::PerspectiveCamera> cams;
        std::vector<std::vector<core::Classified<core::Line2>>> lines;
        std::vector<core::Line3> line3ds;
        std::vector<core::Vec3> vps;
    };

    struct Reconstruction {
        core::View<core::PanoramicCamera> view;
        core::MixedGraph mg;
        core::MixedGraphPropertyTable props;
    };

}


class PanoRecProject : public Project {
public:
    explicit PanoRecProject(const QString & panoIm, QObject * parent) 
        : Project(parent), _panoImFileInfo(panoIm) {
    
        // initialize steps
        int stepLoad = _steps->addStep("Load Panorama", [this](){
            QImage im;
            im.load(_panoImFileInfo.absoluteFilePath());
            auto pim = vis::MakeCVMat(im);
            return data::PanoImage{ core::CreatePanoramicView(pim) };
        });



        // lines and vps
        int stepLinesVPs = _steps->addStep("Perspective Sampling, Extract Lines and Estimate Vanishing Points", [this](const data::PanoImage & im){
            auto & view = im.view;
            
            std::vector<core::PerspectiveCamera> cams;
            std::vector<std::vector<core::Classified<core::Line2>>> lines;
            std::vector<core::Vec3> vps;

            cams = core::CreateCubicFacedCameras(im.view.camera, im.view.image.rows, im.view.image.rows, im.view.image.rows * 0.4);
            lines.resize(cams.size());
            for (int i = 0; i < cams.size(); i++){
                auto pim = view.sampled(cams[i]).image;
                core::LineSegmentExtractor lineExtractor;
                lineExtractor.params().algorithm = core::LineSegmentExtractor::LSD;
                auto ls = lineExtractor(pim, 2, 300); // use pyramid
                lines[i].reserve(ls.size());
                for (auto & l : ls){
                    lines[i].push_back(core::ClassifyAs(l, -1));
                }
                //vis::Visualizer2D(pim) << ls << vis::manip2d::Show();

            }

            // estimate vp
            vps = core::EstimateVanishingPointsAndClassifyLines(cams, lines);
            // extract lines from segmentated region boundaries and classify them using estimated vps
            // make 3d lines
            std::vector<core::Line3> line3ds;
            for (int i = 0; i < cams.size(); i++){
                for (auto & l : lines[i]){
                    line3ds.emplace_back(normalize(cams[i].spatialDirection(l.component.first)),
                        normalize(cams[i].spatialDirection(l.component.second)));
                }
            }

            return data::LinesAndVPs{ std::move(cams), std::move(lines), std::move(line3ds), std::move(vps) };

        }, { stepLoad });


        // segmentation
        int stepSegmentation = _steps->addStep("Segmentation", [this](const data::PanoImage & im, const data::LinesAndVPs & linesVPs){
            
            auto & view = im.view;
            auto & line3ds = linesVPs.line3ds;

            core::Imagei segmentedImage;

            core::SegmentationExtractor segmenter;
            segmenter.params().algorithm = core::SegmentationExtractor::GraphCut;
            segmenter.params().c = 100.0;
            int segmentsNum = 0;
            std::tie(segmentedImage, segmentsNum) = segmenter(view.image, line3ds, view.camera, M_PI / 36.0);

            return data::Segmentation{ std::move(segmentedImage), segmentsNum };

        }, { stepLoad, stepLinesVPs });


        // reconstruction setup
        int stepReconstructionSetup = _steps->addStep("ReconstructionSetup", 
            [this](const data::PanoImage & im, const data::LinesAndVPs & linesVPs, const data::Segmentation & segs){
        
            core::MixedGraph mg;
            core::MixedGraphPropertyTable props;

            auto & view = im.view;
            auto & lines = linesVPs.lines;
            auto & cams = linesVPs.cams;
            auto & vps = linesVPs.vps;
            auto & segmentedImage = segs.segmentation;

            // append lines
            for (int i = 0; i < cams.size(); i++){
                core::AppendLines(mg, lines[i], cams[i], vps);
            }

            // append regions
            core::AppendRegions(mg, segmentedImage, view.camera, 0.01, 0.02, 3, 2);

            props = core::MakeMixedGraphPropertyTable(mg, vps);
            core::AttachPrincipleDirectionConstraints(mg, props, M_PI / 15.0);
            core::AttachWallConstriants(mg, props, M_PI / 30.0);

            return data::Reconstruction{ view, std::move(mg), std::move(props) };

        }, { stepLoad, stepLinesVPs, stepSegmentation });

        
        // reconstruction 1
        int stepReconstruction1 = _steps->addStep("Reconstruction 1",
            [this](const data::Reconstruction & lastRec){

            auto rec = lastRec;
            auto & mg = rec.mg;
            auto & props = rec.props;

            core::SolveVariablesUsingInversedDepths(mg, props);
            core::NormalizeVariables(mg, props);
            std::cout << "score = " << core::ComputeScore(mg, props) << std::endl;
            //core::Visualize(view, mg, props);
            core::LooseOrientationConstraintsOnComponents(mg, props, 0.2, 0.02, 0.1);

            core::SolveVariablesUsingInversedDepths(mg, props);
            core::NormalizeVariables(mg, props);
            std::cout << "score = " << core::ComputeScore(mg, props) << std::endl;
            //core::Visualize(view, mg, props);
            core::AttachFloorAndCeilingConstraints(mg, props, 0.1, 0.6);

            return rec;
        }, { stepReconstructionSetup });

        

    }



private:
    QFileInfo _panoImFileInfo;
    std::vector<core::PerspectiveCamera> _cams;
};



Project::Project(QObject * parent) 
: QObject(parent), _steps(std::make_unique<vis::Steps>()) {
}

Project::Project(const QString & projFile, QObject * parent) : QObject(parent) {
    loadFromDisk(projFile);
}

Project::~Project() {}



void Project::saveToDisk(const QString & filename) const {
    NOT_IMPLEMENTED_YET();
}

void Project::loadFromDisk(const QString & filename){
    _projectFileInfo = QFileInfo(filename);
    NOT_IMPLEMENTED_YET();
}




void Project::update() {

    NOT_IMPLEMENTED_YET();

    _steps->updateAll([this](int stepId){
        

        
    });

}


Project * Project::createProject(const QString & image, bool isPano, QObject * parent){
    if (isPano){
        Project * proj = new PanoRecProject(image, parent);
        return proj;
    }

    NOT_IMPLEMENTED_YET();
}