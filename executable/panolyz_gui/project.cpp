#include <QtGui>
#include <QtWidgets>

#include "../../src/core/basic_types.hpp"
#include "../../src/core/utilities.hpp"
#include "../../src/core/mixed_graph.hpp"
#include "../../src/vis/qt_glue.hpp"

#include "stepsdag.hpp"
#include "configuration.hpp"
#include "project.hpp"

using namespace panoramix;


class PanoRecProject : public Project {
public:
    explicit PanoRecProject(const QString & panoIm, QObject * parent) 
        : Project(parent), _panoImFileInfo(panoIm){
    
        // initialize steps
        int stepLoad = _steps->addStep("Panorama", [this](){
            QImage im;
            im.load(_panoImFileInfo.absoluteFilePath());
            im = im.convertToFormat(QImage::Format_RGB888);
            auto pim = vis::MakeCVMat(im);
            core::ResizeToHeight(pim, 900);
            bool b = core::MakePanorama(pim);
            Q_ASSERT(b);
            return PanoView{ core::CreatePanoramicView(pim) };
        });

        // lines and vps
        int stepLinesVPs = _steps->addStep(tr("Lines and Vanishing Points"), 
            [this](DataOfType<PanoView> & im){
            
            im.lockForRead();

            auto & view = im.content.view;
            
            std::vector<core::PerspectiveCamera> cams;
            std::vector<std::vector<core::Classified<core::Line2>>> lines;
            std::vector<core::Vec3> vps;

            cams = core::CreateCubicFacedCameras(view.camera, view.image.rows, view.image.rows, view.image.rows * 0.4);
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


            im.unlock();

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

            return LinesAndVPs{ std::move(cams), std::move(lines), std::move(line3ds), std::move(vps) };

        }, { stepLoad });


        // segmentation
        int stepSegmentation = _steps->addStep(tr("Segmentation"), 
            [this](DataOfType<PanoView> & im, DataOfType<LinesAndVPs> & linesVPs){
            
            im.lockForRead();
            linesVPs.lockForRead();

            auto & view = im.content.view;
            auto & line3ds = linesVPs.content.line3ds;

            core::Imagei segmentedImage;

            core::SegmentationExtractor segmenter;
            segmenter.params().algorithm = core::SegmentationExtractor::GraphCut;
            segmenter.params().c = 100.0;
            int segmentsNum = 0;

            std::tie(segmentedImage, segmentsNum) = segmenter(view.image, line3ds, view.camera, M_PI / 36.0);
            linesVPs.unlock();
            im.unlock();

            return Segmentation{ std::move(segmentedImage), segmentsNum };

        }, { stepLoad, stepLinesVPs });




        // reconstruction setup
        int stepReconstructionSetup = _steps->addStep(tr("Reconstruction Setup"), 
            [this](DataOfType<PanoView> & im, DataOfType<LinesAndVPs> & linesVPs, 
            DataOfType<Segmentation> & segs){
        
            core::MixedGraph mg;
            core::MixedGraphPropertyTable props;

            im.lockForRead();
            linesVPs.lockForRead();
            segs.lockForRead();

            auto & view = im.content.view;
            auto & lines = linesVPs.content.lines;
            auto & cams = linesVPs.content.cams;
            auto & vps = linesVPs.content.vps;
            auto segmentedImage = segs.content.segmentation;

            // append lines
            for (int i = 0; i < cams.size(); i++){
                core::AppendLines(mg, lines[i], cams[i], vps);
            }

            // append regions
            core::AppendRegions(mg, segmentedImage, view.camera, 0.01, 0.02, 3, 2);

            props = core::MakeMixedGraphPropertyTable(mg, vps);

            im.unlock();
            linesVPs.unlock();
            segs.unlock();

            core::AttachPrincipleDirectionConstraints(mg, props, M_PI / 15.0);
            core::AttachWallConstriants(mg, props, M_PI / 30.0);

            return ReconstructionSetup{ view, std::move(mg), std::move(props), std::move(segmentedImage) };

        }, { stepLoad, stepLinesVPs, stepSegmentation });






        
        // reconstruction
        int stepReconstruction = _steps->addStep("Reconstruction",
            [this](DataOfType<ReconstructionSetup> & lastRec){

            lastRec.lockForRead();
            Reconstruction rec = { lastRec.content.view, lastRec.content.mg, lastRec.content.props };
            lastRec.unlock();

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

            core::SolveVariablesUsingInversedDepths(mg, props);
            core::NormalizeVariables(mg, props);

            return rec;
        }, { stepReconstructionSetup });

        
        // add widgets and actions
        _widgets.clear();
        for (int i = 0; i < _steps->widgets().size(); i++){
            auto w = dynamic_cast<QWidget*>(_steps->widgets()[i]);
            if (w){
                w->setWindowTitle(_steps->stepNameAt(i));
                _widgets << w;
            } 
        }
    }



private:
    QFileInfo _panoImFileInfo;
    std::vector<core::PerspectiveCamera> _cams;
};


static int UnnamedProjectId = 1;
Project::Project(QObject * parent) 
: QObject(parent), _projectFileInfo(tr("Unnamed Project %1").arg(UnnamedProjectId ++))  {
    _steps = new StepsDAG(this);
    connect(_steps, SIGNAL(messageUpdated(QString)), this, SIGNAL(messageUpdated(QString)));
}

Project::Project(const QString & projFile, QObject * parent) : QObject(parent), _projectFileInfo(projFile) {
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




void Project::update(bool forceSourceStepUpdate) {

    _steps->updateAll(nullptr, forceSourceStepUpdate);

}


Project * Project::createProjectFromImage(const QString & image, bool isPano, QObject * parent){
    if (isPano){
        Project * proj = new PanoRecProject(image, parent);
        return proj;
    }

    NOT_IMPLEMENTED_YET();
}

Project * Project::loadProjectFromDisk(const QString & filename, QObject * parent){
    NOT_IMPLEMENTED_YET();
}