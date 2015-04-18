#include <QtGui>
#include <QtWidgets>

#include "../../src/core/basic_types.hpp"
#include "../../src/core/utilities.hpp"
#include "../../src/experimental/rl_graph.hpp"
#include "../../src/experimental/rl_graph_occlusion.hpp"
#include "../../src/experimental/tools.hpp"
#include "../../src/gui/qt_glue.hpp"

#include "stepsdag.hpp"
#include "configuration.hpp"
#include "project.hpp"

using namespace panoramix;
using namespace core;
using namespace experimental;

using PanoView = View<PanoramicCamera>;

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
            ResizeToHeight(pim, 700);
            bool b = MakePanorama(pim);
            Q_ASSERT(b);
            return CreatePanoramicView(pim);
        });

        // lines and vps
        int stepLinesVPs = _steps->addStep(tr("Lines and Vanishing Points"), 
            [this](DataOfType<PanoView> & im){
            
            im.lockForRead();

            auto & view = im.content;
            
            std::vector<PerspectiveCamera> cams;
            std::vector<PerspectiveView> perspectiveViews;
            std::vector<std::vector<Classified<Line2>>> lines;
            std::vector<Vec3> vps;

            cams = CreateCubicFacedCameras(view.camera, view.image.rows, view.image.rows, view.image.rows * 0.4);
            perspectiveViews.resize(cams.size());
            lines.resize(cams.size());
            for (int i = 0; i < cams.size(); i++){
                perspectiveViews[i] = view.sampled(cams[i]);
                LineSegmentExtractor lineExtractor;
                lineExtractor.params().algorithm = LineSegmentExtractor::LSD;
                auto ls = lineExtractor(perspectiveViews[i].image, 2, 300); // use pyramid
                lines[i].reserve(ls.size());
                for (auto & l : ls){
                    lines[i].push_back(ClassifyAs(l, -1));
                }
            }


            im.unlock();

            // estimate vp
            vps = EstimateVanishingPointsAndClassifyLines(cams, lines);
            // extract lines from segmentated region boundaries and classify them using estimated vps
            // make 3d lines
            std::vector<Line3> line3ds;
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

            Imagei segmentedImage;

            SegmentationExtractor segmenter;
            segmenter.params().algorithm = SegmentationExtractor::GraphCut;
            segmenter.params().useYUVColorSpace = false;
            segmenter.params().c = 30.0;
            segmenter.params().superpixelSizeSuggestion = 2000;
            int segmentsNum = 0;

            std::tie(segmentedImage, segmentsNum) = segmenter(view.image ,line3ds, view.camera, M_PI / 36.0);
            linesVPs.unlock();
            im.unlock();

            return Segmentation{ std::move(segmentedImage), segmentsNum };

        }, { stepLoad, stepLinesVPs });



        //// reconstruct lines
        //int stepReconstructLines = _steps->addStep("Reconstruct Lines",
        //    [this](DataOfType<PanoView> & im, DataOfType<LinesAndVPs> & linesVPs){

        //    im.lockForRead();
        //    linesVPs.lockForRead();

        //    auto & view = im.content;
        //    auto & lines = linesVPs.content.lines;
        //    auto & perspectiveViews = linesVPs.content.perspectiveViews;
        //    auto & vps = linesVPs.content.vps;

        //    // append lines
        //    RLGraph mg;
        //    for (int i = 0; i < perspectiveViews.size(); i++){
        //        AppendLines(mg, lines[i], perspectiveViews[i].camera, vps);
        //    }
        //    
        //    // CC
        //    auto ccids = MakeHandledTableForAllComponents<int>(mg);
        //    int ccnum = ConnectedComponents(mg, ccids);


        //    //props = MakeRLGraphPropertyTable(mg, vps);

        //    im.unlock();
        //    linesVPs.unlock();




        //}, { stepLoad, stepLinesVPs });



        // reconstruction setup
        setConf(tr("principle direction constraints angle"), M_PI / 20.0);
        setConf(tr("wall constraints angle"), M_PI / 20.0);
        int stepReconstructionSetup = _steps->addStep(tr("Reconstruction Setup"), 
            [this](DataOfType<PanoView> & im, DataOfType<LinesAndVPs> & linesVPs, 
            DataOfType<Segmentation> & segs) -> ReconstructionSetup<PanoramicCamera> {

            RLGraph mg;
            RLGraphControls controls;

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
                AppendLines(mg, lines[i], perspectiveViews[i].camera, vps); 
            }

            _confLock.lockForRead();

            // append regions
            auto segId2Rhs = AppendRegions(mg, segmentedImage, view.camera, 0.01, 0.02, 3, 2, true);
            controls = RLGraphControls(mg, vps);

            im.unlock();
            linesVPs.unlock();
            segs.unlock();

            AttachWallConstriants(mg, controls, conf("wall constraints angle").value<double>());
            AttachPrincipleDirectionConstraints(mg, controls, conf("principle direction constraints angle").value<double>());   

            _confLock.unlock();

            return ReconstructionSetup<PanoramicCamera>{ 
                view, 
                    std::move(mg), std::move(controls),
                    std::move(segmentedImage), std::move(segId2Rhs)
            };

        }, { stepLoad, stepLinesVPs, stepSegmentation });
        


        
        // reconstruction
        int stepReconstruction = _steps->addStep("Reconstruction",
            [this](DataOfType<ReconstructionSetup<PanoramicCamera>> & lastRec){

            lastRec.lockForRead();

            auto ccids = MakeHandledTableForAllComponents(lastRec.content.mg, -1);
            int ccnum = ConnectedComponents(lastRec.content.mg, ccids);
            std::vector<RLGraph> mgs = Decompose(lastRec.content.mg, ccids, ccnum);
            std::vector<RLGraphControls> cs = Decompose(lastRec.content.mg, lastRec.content.controls, ccids, ccnum);
            assert(cs.size() == mgs.size());

            Reconstruction<PanoramicCamera> rec;
            rec.view = lastRec.content.view;
            rec.gs.resize(mgs.size());

            lastRec.unlock();

            for (int i = 0; i < rec.gs.size(); i++){

                auto & mg = rec.gs[i].mg;
                auto & controls = rec.gs[i].controls;
                auto & vars = rec.gs[i].vars;

                mg = std::move(mgs[i]);
                controls = std::move(cs[i]);

                if (!AttachAnchorToCenterOfLargestRegionIfNoAnchorExists(mg, controls, 1.0, 1.0))
                    continue;

                ResetToSampledArmorAnchors(mg, controls, 0.1);
                vars = SolveVariablesBoundComponentAnchors(mg, controls, false, false, 1.0, 5.0, 100);
                NormalizeVariables(mg, controls, vars);
                std::cout << "score = " << Score(mg, controls, vars) << std::endl;

                LooseOrientationConstraintsOnComponents(mg, controls, vars, 0.2, 0.02, 0.1);
                if (!AttachAnchorToCenterOfLargestLineIfNoAnchorExists(mg, controls))
                    continue;

                
                vars = SolveVariablesBoundComponentAnchors(mg, controls, false, false, 1.0, 5.0, 100);
                NormalizeVariables(mg, controls, vars);

                AttachFloorAndCeilingConstraints(mg, controls, vars, 0.1, 0.6);

                if (!AttachAnchorToCenterOfLargestRegionIfNoAnchorExists(mg, controls))
                    continue;
                vars = SolveVariablesBoundComponentAnchors(mg, controls, false, false, 1.0, 5.0, 100);
                NormalizeVariables(mg, controls, vars);
            }

            return rec;
        }, { stepReconstructionSetup });



        // reconstruction refinement
        int stepReconstructionRefinement = _steps->addStep("Reconstruction Refinement",
            [this](DataOfType<Reconstruction<PanoramicCamera>> & rec){

            ReconstructionRefinement<PanoramicCamera> refinement;
            refinement.shapes.reserve(rec.content.gs.size());

            for (auto & g : rec.content.gs){
                auto & mg = g.mg;
                auto & controls = g.controls;
                auto & vars = g.vars;

                rec.lockForRead();
                refinement.view = rec.content.view;
                auto polygons = RegionPolygons(mg, controls, vars);
                int vertVPId = NearestDirectionId(controls.vanishingPoints);
                double medianDepth = MedianCenterDepth(mg, controls, vars);
                Vec3 vertDir = normalize(controls.vanishingPoints[vertVPId]);
                rec.unlock();

                auto range = EstimateEffectiveRangeAlongDirection(polygons, vertDir, medianDepth * 0.02, 0.7, -1e5, -1e5);

                std::vector<Chain3> chains;
                for (double x = range.first; x <= range.second; x += medianDepth * 0.02){
                    Plane3 cutplane(vertDir * x, vertDir);
                    auto loop = MakeSectionalPieces(polygons, cutplane);
                    if (loop.empty())
                        continue;
                    chains.push_back(MakeChain(loop));
                }

                LayeredShape3 shape;
                for (int i = 0; i < chains.size(); i++){
                    shape.layers.push_back(std::move(chains[i].points));
                }
                shape.normal = vertDir;
                refinement.shapes.push_back(std::move(shape));
            }

            return refinement;

        }, { stepReconstruction });


        
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
    std::vector<PerspectiveCamera> _cams;
};

using PerspView = View<PerspectiveCamera>;



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
            ResizeToHeight(cvim, 500);

            View<PerspectiveCamera> view;
            std::vector<Classified<Line2>> lines;
            std::vector<Classified<Line3>> line3ds;
            std::vector<Vec3> vps;

            double focal;
            VanishingPointsDetector::Params vpdParams(VanishingPointsDetector::TardifSimplified);
            view = CreatePerspectiveView(cvim, Point3(0, 0, 0), Point3(1, 0, 0), Point3(0, 0, -1),
                LineSegmentExtractor(), VanishingPointsDetector(vpdParams), &line3ds, &lines, &vps, &focal).unwrap();

            std::vector<Line3> pureLine3ds(line3ds.size());
            for (int i = 0; i < line3ds.size(); i++){
                pureLine3ds[i] = line3ds[i].component;
            }

            LinesAndVPs result;
            result.perspectiveViews = { std::move(view) };
            result.lines = { std::move(lines) };
            result.line3ds = std::move(pureLine3ds);
            result.vps = std::move(vps);
            return result;
        });

        // segmentation
        int stepSegmentation = _steps->addStep(tr("Segmentation"),
            [this](DataOfType<LinesAndVPs> & linesVPs){

            linesVPs.lockForRead();

            assert(linesVPs.content.perspectiveViews.size() == 1);
            auto & view = linesVPs.content.perspectiveViews.front();
            std::vector<Line2> pureLines(linesVPs.content.lines.size());
            for (int i = 0; i < pureLines.size(); i++){
                pureLines[i] = linesVPs.content.lines.front()[i].component;
            }

            Imagei segmentedImage;

            SegmentationExtractor segmenter;
            segmenter.params().algorithm = SegmentationExtractor::GraphCut;
            segmenter.params().useYUVColorSpace = false;
            segmenter.params().c = 100.0;
            int segmentsNum = 0;

            std::tie(segmentedImage, segmentsNum) = segmenter(view.image, pureLines, view.image.cols / 100.0);
            linesVPs.unlock();

            return Segmentation{ std::move(segmentedImage), segmentsNum };

        }, { stepLoad });




        // reconstruction setup
        int stepReconstructionSetup = _steps->addStep(tr("Reconstruction Setup"),
            [this](DataOfType<LinesAndVPs> & linesVPs, 
            DataOfType<Segmentation> & segs) -> ReconstructionSetup<PerspectiveCamera> {

            RLGraph mg;
            RLGraphControls controls;

            linesVPs.lockForRead();
            segs.lockForRead();

            auto & lines = linesVPs.content.lines;
            auto & perspectiveViews = linesVPs.content.perspectiveViews;
            auto & vps = linesVPs.content.vps;
            auto segmentedImage = segs.content.segmentation;

            // append lines
            for (int i = 0; i < perspectiveViews.size(); i++)
                AppendLines(mg, lines[i], perspectiveViews[i].camera, vps, 
                40.0 / 500.0, 
                100.0 / 500.0);

            _confLock.lockForRead();

            // append regions
            auto segId2Rhs = AppendRegions(mg, segmentedImage, perspectiveViews.front().camera, 0.01, 0.001, 1, 3, true);
            controls = RLGraphControls(mg, vps);

            linesVPs.unlock();
            segs.unlock();

            //AttachPrincipleDirectionConstraints(mg, controls, M_PI / 120.0);
            //AttachWallConstriants(mg, controls, M_PI / 100.0);

            _confLock.unlock();

            return ReconstructionSetup<PerspectiveCamera>{ 
                perspectiveViews.front(),
                    std::move(mg),
                    std::move(controls),
                    std::move(segmentedImage), std::move(segId2Rhs)
            };

        }, { stepLoad, stepSegmentation });





        // reconstruction
        int stepReconstruction = _steps->addStep("Reconstruction",
            [this](DataOfType<ReconstructionSetup<PerspectiveCamera>> & lastRec){

            auto & mg = lastRec.content.mg;
            auto & controls = lastRec.content.controls;

            auto ccids = MakeHandledTableForAllComponents(mg, -1);
            int ccnum = ConnectedComponents(mg, controls, ccids, [](const RLGraphConstraintControl & c){
                return c.used && c.weight > 0;
            });
            RLGraphOldToNew old2new;
            auto mgs = Decompose(mg, ccids, ccnum, &old2new);
            for (auto & o2n : old2new.container<RegionHandle>()){
                auto nrh = o2n.second.second;
                int ccid = o2n.second.first;
                assert(!nrh.invalid());
                assert(nrh.id < mgs[ccid].internalComponents<RegionData>().size());
            }
            auto cs = Decompose(mg, controls, ccids, ccnum);
            assert(mgs.size() == cs.size());

            Reconstruction<PerspectiveCamera> rec;
            rec.view = lastRec.content.view;
            rec.gs.resize(mgs.size());

            lastRec.unlock();

            for (int i = 0; i < rec.gs.size(); i++){

                auto & mg = rec.gs[i].mg;
                auto & controls = rec.gs[i].controls;
                auto & vars = rec.gs[i].vars;

                mg = std::move(mgs[i]);
                controls = std::move(cs[i]);

                if (!AttachAnchorToCenterOfLargestLineIfNoAnchorExists(mg, controls, 1.0, 1.0))
                    continue;

                vars = SolveVariablesBoundComponentAnchors(mg, controls, false, true, 1.0, 5.0, 100);
                NormalizeVariables(mg, controls, vars);
                std::cout << "score = " << Score(mg, controls, vars) << std::endl;
            }

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