#include <QtGui>
#include <QtWidgets>

#include "../../src/core/basic_types.hpp"
#include "../../src/core/utilities.hpp"
#include "../../src/experimental/rl_graph.hpp"
#include "../../src/gui/qt_glue.hpp"

#include "stepsdag.hpp"
#include "configuration.hpp"
#include "project.hpp"

using namespace panoramix;
using PanoView = core::View<core::PanoramicCamera>;

namespace panoramix {
    namespace core{
        using namespace experimental;
    }
}

class PanoRecProject : public Project {
public:
    explicit PanoRecProject(const QString & panoIm, QObject * parent) 
        : Project(parent), _panoImFileInfo(panoIm){
    
        // initialize steps
        int stepLoad = _steps->addStep("Panorama", [this](){
            QImage im;
            im.load(_panoImFileInfo.absoluteFilePath());
            im = im.convertToFormat(QImage::Format_RGB888);
            auto pim = gui::MakeCVMat(im);
            core::ResizeToHeight(pim, 700);
            bool b = core::MakePanorama(pim);
            Q_ASSERT(b);
            return core::CreatePanoramicView(pim);
        });

        // lines and vps
        int stepLinesVPs = _steps->addStep(tr("Lines and Vanishing Points"), 
            [this](DataOfType<PanoView> & im){
            
            im.lockForRead();

            auto & view = im.content;
            
            std::vector<core::PerspectiveCamera> cams;
            std::vector<core::PerspectiveView> perspectiveViews;
            std::vector<std::vector<core::Classified<core::Line2>>> lines;
            std::vector<core::Vec3> vps;

            cams = core::CreateCubicFacedCameras(view.camera, view.image.rows, view.image.rows, view.image.rows * 0.4);
            perspectiveViews.resize(cams.size());
            lines.resize(cams.size());
            for (int i = 0; i < cams.size(); i++){
                perspectiveViews[i] = view.sampled(cams[i]);
                core::LineSegmentExtractor lineExtractor;
                lineExtractor.params().algorithm = core::LineSegmentExtractor::LSD;
                auto ls = lineExtractor(perspectiveViews[i].image, 2, 300); // use pyramid
                lines[i].reserve(ls.size());
                for (auto & l : ls){
                    lines[i].push_back(core::ClassifyAs(l, -1));
                }
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

            return LinesAndVPs{ 
                std::move(perspectiveViews), 
                std::move(lines), 
                std::move(line3ds), 
                std::move(vps) 
            };

        }, { stepLoad });


        // segmentation
        int stepSegmentation = _steps->addStep(tr("Segmentation"), 
            [this](DataOfType<PanoView> & im, DataOfType<LinesAndVPs> & linesVPs){
            
            im.lockForRead();
            linesVPs.lockForRead();

            auto & view = im.content;
            auto & line3ds = linesVPs.content.line3ds;

            core::Imagei segmentedImage;

            core::SegmentationExtractor segmenter;
            segmenter.params().algorithm = core::SegmentationExtractor::GraphCut;
            segmenter.params().useYUVColorSpace = false;
            segmenter.params().c = 100.0;
            int segmentsNum = 0;

            std::tie(segmentedImage, segmentsNum) = segmenter(view.image, line3ds, view.camera, M_PI / 36.0);
            linesVPs.unlock();
            im.unlock();

            return Segmentation{ std::move(segmentedImage), segmentsNum };

        }, { stepLoad, stepLinesVPs });




        // reconstruction setup
        setConf(tr("principle direction constraints angle"), M_PI / 30.0);
        setConf(tr("wall constraints angle"), M_PI / 60.0);
        int stepReconstructionSetup = _steps->addStep(tr("Reconstruction Setup"), 
            [this](DataOfType<PanoView> & im, DataOfType<LinesAndVPs> & linesVPs, 
            DataOfType<Segmentation> & segs) -> ReconstructionSetup<core::PanoramicCamera> {
        
            core::RLGraph mg;
            core::RLGraphPropertyTable props;

            im.lockForRead();
            linesVPs.lockForRead();
            segs.lockForRead();

            auto & view = im.content;
            auto & lines = linesVPs.content.lines;
            auto & perspectiveViews = linesVPs.content.perspectiveViews;
            auto & vps = linesVPs.content.vps;
            auto segmentedImage = segs.content.segmentation;

            // append lines
            for (int i = 0; i < perspectiveViews.size(); i++){
                core::AppendLines(mg, lines[i], perspectiveViews[i].camera, vps);
            }

            _confLock.lockForRead();

            // append regions
            core::AppendRegions(mg, segmentedImage, view.camera, 0.01, 0.02, 3, 2);

            props = core::MakeRLGraphPropertyTable(mg, vps);

            im.unlock();
            linesVPs.unlock();
            segs.unlock();

            core::AttachWallConstriants(mg, props, conf("wall constraints angle").value<double>());
            core::AttachPrincipleDirectionConstraints(mg, props, conf("principle direction constraints angle").value<double>());   
            AttachAnchorToCenterOfLargestRegion(mg, props, 1.0, 10.0);

            _confLock.unlock();

            return ReconstructionSetup<core::PanoramicCamera>{ 
                view, 
                    std::move(mg), std::move(props), 
                    std::move(segmentedImage) 
            };

        }, { stepLoad, stepLinesVPs, stepSegmentation });






        
        // reconstruction
        int stepReconstruction = _steps->addStep("Reconstruction",
            [this](DataOfType<ReconstructionSetup<core::PanoramicCamera>> & lastRec){

            lastRec.lockForRead();
            Reconstruction<core::PanoramicCamera> rec = { lastRec.content.view, lastRec.content.mg, lastRec.content.props };
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

using PerspView = panoramix::core::View<panoramix::core::PerspectiveCamera>;

class NormalRecProject : public Project {
public:
    explicit NormalRecProject(const QString & normalIm, QObject * parent) 
        : Project(parent), _imFileInfo(normalIm) {
        
        // initialize steps
        int stepLoad = _steps->addStep("Photo, lines and Vanishing Points", [this](){
            QImage im;
            im.load(_imFileInfo.absoluteFilePath());
            im = im.convertToFormat(QImage::Format_RGB888);
            auto cvim = gui::MakeCVMat(im);
            core::ResizeToHeight(cvim, 500);

            core::View<core::PerspectiveCamera> view;
            std::vector<core::Classified<core::Line2>> lines;
            std::vector<core::Classified<core::Line3>> line3ds;
            std::vector<core::Vec3> vps;

            double focal;
            core::VanishingPointsDetector::Params vpdParams(core::VanishingPointsDetector::TardifSimplified);
            view = core::CreatePerspectiveView(cvim, core::Point3(0, 0, 0), core::Point3(1, 0, 0), core::Point3(0, 0, -1),
                core::LineSegmentExtractor(), core::VanishingPointsDetector(vpdParams), &line3ds, &lines, &vps, &focal).unwrap();

            std::vector<core::Line3> pureLine3ds(line3ds.size());
            for (int i = 0; i < line3ds.size(); i++){
                pureLine3ds[i] = line3ds[i].component;
            }

            LinesAndVPs result;
            result.perspectiveViews = { std::move(view) };
            result.lines = { std::move(lines) };
            result.line3ds = std::move(pureLine3ds);
            result.vps = std::move(vps);
            return result;
            //return LinesAndVPs{ 
            //    std::vector<core::PerspectiveView>{ std::move(view) },
            //    std::vector<std::vector<core::Classified<core::Line2>>>{ std::move(lines) },
            //    std::move(pureLine3ds),
            //    std::move(vps) 
            //};
        });

        // segmentation
        int stepSegmentation = _steps->addStep(tr("Segmentation"),
            [this](DataOfType<LinesAndVPs> & linesVPs){

            linesVPs.lockForRead();

            assert(linesVPs.content.perspectiveViews.size() == 1);
            auto & view = linesVPs.content.perspectiveViews.front();
            std::vector<core::Line2> pureLines(linesVPs.content.lines.size());
            for (int i = 0; i < pureLines.size(); i++){
                pureLines[i] = linesVPs.content.lines.front()[i].component;
            }

            core::Imagei segmentedImage;

            core::SegmentationExtractor segmenter;
            segmenter.params().algorithm = core::SegmentationExtractor::GraphCut;
            segmenter.params().useYUVColorSpace = false;
            segmenter.params().c = 100.0;
            int segmentsNum = 0;

            std::tie(segmentedImage, segmentsNum) = segmenter(view.image, pureLines, view.image.cols / 100.0);
            linesVPs.unlock();

            return Segmentation{ std::move(segmentedImage), segmentsNum };

        }, { stepLoad });




        // reconstruction setup
        int stepReconstructionSetup = _steps->addStep(tr("Reconstruction Setup"),
            [this](DataOfType<LinesAndVPs> & linesVPs, DataOfType<Segmentation> & segs) -> ReconstructionSetup<core::PerspectiveCamera> {

            core::RLGraph mg;
            core::RLGraphPropertyTable props;

            linesVPs.lockForRead();
            segs.lockForRead();

            auto & lines = linesVPs.content.lines;
            auto & perspectiveViews = linesVPs.content.perspectiveViews;
            auto & vps = linesVPs.content.vps;
            auto segmentedImage = segs.content.segmentation;

            // append lines
            for (int i = 0; i < perspectiveViews.size(); i++)
                core::AppendLines(mg, lines[i], perspectiveViews[i].camera, vps);

            // append regions
            core::AppendRegions(mg, segmentedImage, perspectiveViews.front().camera, 0.001, 0.001, 3, 1);

            props = core::MakeRLGraphPropertyTable(mg, vps);

            linesVPs.unlock();
            segs.unlock();

            core::AttachPrincipleDirectionConstraints(mg, props, M_PI / 120.0);
            core::AttachWallConstriants(mg, props, M_PI / 100.0);
            AttachAnchorToCenterOfLargestRegion(mg, props, 1.0, 10.0);

            return ReconstructionSetup<core::PerspectiveCamera>{ 
                perspectiveViews.front(), 
                    std::move(mg), 
                    std::move(props), 
                    std::move(segmentedImage)
            };

        }, { stepLoad, stepSegmentation });







        // reconstruction
        int stepReconstruction = _steps->addStep("Reconstruction",
            [this](DataOfType<ReconstructionSetup<core::PerspectiveCamera>> & lastRec){

            lastRec.lockForRead();
            Reconstruction<core::PerspectiveCamera> rec = { 
                lastRec.content.view, lastRec.content.mg, lastRec.content.props
            };
            lastRec.unlock();

            auto & mg = rec.mg;
            auto & props = rec.props;

            core::SolveVariablesUsingInversedDepths(mg, props);
            core::NormalizeVariables(mg, props);
            std::cout << "score = " << core::ComputeScore(mg, props) << std::endl;
            //core::Visualize(view, mg, props);
            core::LooseOrientationConstraintsOnComponents(mg, props, 0.2, 0, 0.05);

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
    QFileInfo _imFileInfo;
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
    else {
        Project * proj = new NormalRecProject(image, parent);
        return proj;
    }

    NOT_IMPLEMENTED_YET();
}

Project * Project::loadProjectFromDisk(const QString & filename, QObject * parent){
    NOT_IMPLEMENTED_YET();
}