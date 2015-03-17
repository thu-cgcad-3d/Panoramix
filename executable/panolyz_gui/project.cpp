#include <QtGui>
#include <QtWidgets>

#include "../../src/core/algorithms.hpp"
#include "../../src/core/basic_types.hpp"
#include "../../src/core/mixed_graph.hpp"
#include "../../src/vis/qt_glue.hpp"
#include "../../src/vis/visualizers.hpp"

#include "steps.hpp"
#include "project.hpp"

using namespace panoramix;






struct PanoView {
    core::View<core::PanoramicCamera> view;
};

QWidget * CreateBindingWidgetAndActions(DataOfType<PanoView> & pv,
    QList<QAction*> & actions, QWidget * parent) {

    class Widget : public QWidget {
    public:
        explicit Widget(DataOfType<PanoView> * d, QWidget * parent) : QWidget(parent), data(d) {
            setMinimumSize(300, 300);
        }
    protected:
        virtual void paintEvent(QPaintEvent * e) override {
            if (!data)
                return;
            QPainter painter(this);
            data->lockForRead();
            painter.drawImage(QPointF(), vis::MakeQImage(data->content.view.image));
            data->unlock();
        }
    private:
        DataOfType<PanoView>* data;
    };

    return new Widget(&pv, parent);
}

struct Segmentation {
    core::Imagei segmentation;
    int segmentsNum;
};


QWidget * CreateBindingWidgetAndActions(DataOfType<Segmentation> & segs,
    QList<QAction*> & actions, QWidget * parent) {

    class Widget : public QWidget {
    public:
        explicit Widget(DataOfType<Segmentation> * d, QWidget * parent) : QWidget(parent), data(d) {
            setMinimumSize(300, 300);
        }
    protected:
        virtual void paintEvent(QPaintEvent * e) override {
            if (!data)
                return;
            QPainter painter(this);
            data->lockForRead();
            // todo
            data->unlock();
        }
    private:
        DataOfType<Segmentation>* data;
        
    };

    return nullptr;
}



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

QWidget * CreateBindingWidgetAndActions(DataOfType<Reconstruction> & rec,
    QList<QAction*> & actions, QWidget * parent) {
    
    vis::ResourceStore::set("texture", rec.content.view.image);

    struct ComponentID {
        int handleID;
        bool isRegion;
    };

    vis::Visualizer viz("mixed graph optimizable");
    viz.renderOptions.bwColor = 1.0;
    viz.renderOptions.bwTexColor = 0.0;
    viz.installingOptions.discretizeOptions.colorTable = vis::ColorTableDescriptor::RGB;
    std::vector<std::pair<ComponentID, vis::Colored<vis::SpatialProjectedPolygon>>> spps;
    std::vector<vis::Colored<core::Line3>> lines;

    auto & mg = rec.content.mg;
    auto & props = rec.content.props;

    for (auto & c : mg.components<core::RegionData>()){
        if (!props[c.topo.hd].used)
            continue;
        auto uh = c.topo.hd;
        auto & region = c.data;
        vis::SpatialProjectedPolygon spp;
        // filter corners
        core::ForeachCompatibleWithLastElement(c.data.normalizedContours.front().begin(), c.data.normalizedContours.front().end(),
            std::back_inserter(spp.corners),
            [](const core::Vec3 & a, const core::Vec3 & b) -> bool {
            return core::AngleBetweenDirections(a, b) > M_PI / 1000.0;
        });
        if (spp.corners.size() < 3)
            continue;

        spp.projectionCenter = core::Point3(0, 0, 0);
        spp.plane = Instance(mg, props, uh);
        assert(!core::HasValue(spp.plane, core::IsInfOrNaN<double>));
        spps.emplace_back(ComponentID{ uh.id, true }, std::move(vis::ColorAs(spp, vis::ColorTag::Black)));
    }

    for (auto & c : mg.components<core::LineData>()){
        if (!props[c.topo.hd].used)
            continue;
        auto uh = c.topo.hd;
        auto & line = c.data;
        lines.push_back(vis::ColorAs(Instance(mg, props, uh), vis::ColorTag::Black));
    }

    viz.begin(spps).shaderSource(vis::OpenGLShaderSourceDescriptor::XPanorama).resource("texture").end();
    viz.installingOptions.discretizeOptions.color = vis::ColorTag::DarkGray;
    viz.installingOptions.lineWidth = 5.0;
    viz.add(lines);

    //std::vector<core::Line3> connectionLines;
    //for (auto & c : mg.constraints<core::RegionBoundaryData>()){
    //    if (!props[c.topo.hd].used)
    //        continue;
    //    auto & samples = c.data.normalizedSampledPoints;
    //    auto inst1 = Instance(mg, props, c.topo.component<0>());
    //    auto inst2 = Instance(mg, props, c.topo.component<1>());
    //    for (auto & ss : samples){
    //        for (auto & s : ss){
    //            double d1 = DepthAt(s, inst1);
    //            double d2 = DepthAt(s, inst2);
    //            connectionLines.emplace_back(normalize(s) * d1, normalize(s) * d2);
    //        }
    //    }
    //}
    //for (auto & c : mg.constraints<core::RegionLineConnectionData>()){
    //    if (!props[c.topo.hd].used)
    //        continue;
    //    auto inst1 = Instance(mg, props, c.topo.component<0>());
    //    auto inst2 = Instance(mg, props, c.topo.component<1>());
    //    for (auto & s : c.data.normalizedAnchors){
    //        double d1 = DepthAt(s, inst1);
    //        double d2 = DepthAt(s, inst2);
    //        connectionLines.emplace_back(normalize(s) * d1, normalize(s) * d2);
    //    }
    //}

    viz.installingOptions.discretizeOptions.color = vis::ColorTag::Black;
    viz.installingOptions.lineWidth = 1.0;
    //viz.add(connectionLines);

    viz.renderOptions.renderMode = vis::RenderModeFlag::Triangles | vis::RenderModeFlag::Lines;
    viz.renderOptions.backgroundColor = vis::ColorTag::White;
    viz.renderOptions.bwColor = 0.5;
    viz.renderOptions.bwTexColor = 0.5;
    viz.camera(core::PerspectiveCamera(1000, 800, 800, core::Point3(-1, 1, 1), core::Point3(0, 0, 0)));
    auto w = viz.createWidget(false, parent);
    for (auto a : w->actions()){
        actions << a;
    }
    return w;
}



class PanoRecProject : public Project {
public:
    explicit PanoRecProject(const QString & panoIm, QObject * parent) 
        : Project(parent), _panoImFileInfo(panoIm){
    
        // initialize steps
        int stepLoad = _steps->addStep("Load Panorama", [this](){
            QImage im;
            im.load(_panoImFileInfo.absoluteFilePath());
            im = im.convertToFormat(QImage::Format_RGB888);
            auto pim = vis::MakeCVMat(im);
            core::ResizeToMakeHeightUnder(pim, 900);
            return PanoView{ core::CreatePanoramicView(pim) };
        });

        // lines and vps
        int stepLinesVPs = _steps->addStep(tr("Perspective Sampling, Extract Lines and Estimate Vanishing Points"), 
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
        int stepReconstructionSetup = _steps->addStep("Reconstruction Setup", 
            [this](DataOfType<PanoView> & im, DataOfType<LinesAndVPs> & linesVPs, DataOfType<Segmentation> & segs){
        
            core::MixedGraph mg;
            core::MixedGraphPropertyTable props;

            im.lockForRead();
            linesVPs.lockForRead();
            segs.lockForRead();

            auto & view = im.content.view;
            auto & lines = linesVPs.content.lines;
            auto & cams = linesVPs.content.cams;
            auto & vps = linesVPs.content.vps;
            auto & segmentedImage = segs.content.segmentation;

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

            return Reconstruction{ view, std::move(mg), std::move(props) };

        }, { stepLoad, stepLinesVPs, stepSegmentation });

        
        // reconstruction 1
        int stepReconstruction1 = _steps->addStep("Reconstruction 1",
            [this](DataOfType<Reconstruction> & lastRec){

            lastRec.lockForRead();
            auto rec = lastRec.content;
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

            return rec;
        }, { stepReconstructionSetup });

        
        // add widgets and actions
        _widgets.clear();
        for (auto w : _steps->widgets()){
            if(w) _widgets << w;
        }
    }



private:
    QFileInfo _panoImFileInfo;
    std::vector<core::PerspectiveCamera> _cams;
};


static int UnnamedProjectId = 1;
Project::Project(QObject * parent) 
: QObject(parent), _steps(std::make_unique<Steps>()), _projectFileInfo(tr("Unnamed Project %1").arg(UnnamedProjectId ++))  {
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