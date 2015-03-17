#include <QtGui>
#include <QtWidgets>
#include <QtOpenGL>

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

WidgetLike * CreateBindingWidgetAndActions(DataOfType<PanoView> & pv,
    QList<QAction*> & actions, QWidget * parent) {

    class Widget : public QWidget, public WidgetLike {
    public:
        explicit Widget(DataOfType<PanoView> * d, QWidget * parent) : QWidget(parent), data(d) {
            setMinimumSize(300, 300);
        }

        virtual void refreshData() {}
        virtual void updatePainting() { update(); }
        virtual void showWidget() { show(); }
        virtual void hideWidget() { hide(); }

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


WidgetLike * CreateBindingWidgetAndActions(DataOfType<Segmentation> & segs,
    QList<QAction*> & actions, QWidget * parent) {

    class Widget : public QWidget, public WidgetLike {
    public:
        explicit Widget(DataOfType<Segmentation> * d, QWidget * parent) : QWidget(parent), data(d) {
            setMinimumSize(300, 300);
        }

        virtual void refreshData() {}
        virtual void updatePainting() { update(); }
        virtual void showWidget() { show(); }
        virtual void hideWidget() { hide(); }

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

WidgetLike * CreateBindingWidgetAndActions(DataOfType<Reconstruction> & rec,
    QList<QAction*> & actions, QWidget * parent) {
    
    class Widget : public QGLWidget, public WidgetLike {
    public:
        explicit Widget(DataOfType<Reconstruction> * d, QWidget * parent = nullptr)
            : QGLWidget(parent), data(d) {
            setMouseTracking(true);
            setAutoBufferSwap(false);
            grabKeyboard();
            setUpScene();
            _needsInitialization = true;
            setMinimumSize(300, 300);
        }

        virtual void refreshData() { setUpScene(); }
        virtual void updatePainting() { update(); }
        virtual void showWidget() { show(); }
        virtual void hideWidget() { hide(); }

    public:
        void setUpScene() {
            if (!data){
                scene.clear();
                return;
            }

            struct ComponentID {
                int handleID;
                bool isRegion;
            };


            vis::VisualObjectTree tree;
            auto activeOH = tree.addRoot(std::make_shared<vis::VisualObject>());

            data->lockForRead();
            vis::ResourceStore::set("texture", data->content.view.image);
            data->unlock();

            vis::Visualizer viz("mixed graph optimizable");
            viz.renderOptions.renderMode = vis::RenderModeFlag::Triangles;
            viz.renderOptions.bwColor = 1.0;
            viz.renderOptions.bwTexColor = 0.0;
            viz.installingOptions.discretizeOptions.colorTable = vis::ColorTableDescriptor::RGB;
            std::vector<std::pair<ComponentID, vis::Colored<vis::SpatialProjectedPolygon>>> spps;
            std::vector<vis::Colored<core::Line3>> lines;

            auto & mg = data->content.mg;
            auto & props = data->content.props;

            data->lockForRead();
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

            viz.installingOptions.discretizeOptions.color = vis::ColorTag::Black;
            viz.installingOptions.lineWidth = 1.0;
            //viz.add(connectionLines);

            viz.renderOptions.renderMode = vis::RenderModeFlag::Triangles | vis::RenderModeFlag::Lines;
            viz.renderOptions.backgroundColor = vis::ColorTag::White;
            viz.renderOptions.bwColor = 0.5;
            viz.renderOptions.bwTexColor = 0.5;
            viz.camera(core::PerspectiveCamera(1000, 800, 800, core::Point3(-1, 1, 1), core::Point3(0, 0, 0)));
            renderOptions = viz.renderOptions;

            data->unlock();

            scene.install(viz.tree());

            _needsInitializationLock.lockForWrite();
            _needsInitialization = true;
            _needsInitializationLock.unlock();
            //initializeGL();
        }

        void autoSetCamera() {
            auto sphere = scene.boundingBox().outerSphere();
            renderOptions.camera.resizeScreen(core::Size(width(), height()), false);
            renderOptions.camera.focusOn(sphere, true);
            update();
        }

    protected:
        void initializeGL() {
            if (!_needsInitialization)
                return;
            makeCurrent();
            glEnable(GL_MULTISAMPLE);
            GLint bufs;
            GLint samples;
            glGetIntegerv(GL_SAMPLE_BUFFERS, &bufs);
            glGetIntegerv(GL_SAMPLES, &samples);
            qDebug("Have %d buffers and %d samples", bufs, samples);
            qglClearColor(MakeQColor(renderOptions.backgroundColor));
            scene.initialize();

            _needsInitializationLock.lockForWrite();
            _needsInitialization = false;
            _needsInitializationLock.unlock();
        }

        void paintGL() {
            
            initializeGL();

            QPainter painter;
            painter.begin(this);

            painter.beginNativePainting();
            qglClearColor(MakeQColor(renderOptions.backgroundColor));

            glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);

            glFrontFace(GL_CCW); // face direction set to clockwise
            glEnable(GL_MULTISAMPLE);
            glEnable(GL_DEPTH_TEST);
            glEnable(GL_STENCIL_TEST);

            glEnable(GL_ALPHA_TEST);
            if (renderOptions.showInside){
                glEnable(GL_CULL_FACE);
            }
            else{
                glDisable(GL_CULL_FACE);
            }

            glEnable(GL_BLEND);
            glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);

            glEnable(GL_PROGRAM_POINT_SIZE);
            glEnable(GL_POINT_SPRITE);

            core::PerspectiveCamera & camera = renderOptions.camera;
            camera.resizeScreen(core::Size(width(), height()));

            scene.render(renderOptions);

            glDisable(GL_DEPTH_TEST);
            if (renderOptions.showInside){
                glDisable(GL_CULL_FACE);
            }

            painter.endNativePainting();
            swapBuffers();
        }

        void resizeGL(int w, int h) {
            core::PerspectiveCamera & camera = renderOptions.camera;
            camera.resizeScreen(core::Size(w, h));
            glViewport(0, 0, w, h);
        }

    protected:
        virtual void mousePressEvent(QMouseEvent * e) override {
            _lastPos = e->pos();
            if (e->buttons() & Qt::RightButton)
                setCursor(Qt::OpenHandCursor);
            else if (e->buttons() & Qt::MidButton)
                setCursor(Qt::SizeAllCursor);
            else if (e->buttons() & Qt::LeftButton){
                vis::VisualObjectHandle oh;
                vis::TriMesh::TriangleHandle t;
                std::tie(oh, t) = scene.pickOnScreen(renderOptions, core::Point2(e->pos().x(), e->pos().y()));
                if (oh.valid()){
                    int entityID = scene.tree().data(oh)->entityIDOfMeshTriangle(t);
                    if (e->modifiers() & Qt::ControlModifier){
                        scene.switchSelect(std::make_pair(oh, entityID));
                    }
                    else{
                        scene.clearSelection();
                        scene.select(std::make_pair(oh, entityID));
                    }
                    scene.tree().data(oh)
                        ->invokeCallbackFunction(vis::InteractionID::ClickLeftButton, scene.tree(), std::make_pair(oh, entityID));
                }
                else if (!(e->modifiers() & Qt::ControlModifier)){
                    scene.clearSelection();
                }
            }
            update();
        }

        virtual void mouseMoveEvent(QMouseEvent * e) override {
            QVector3D t(e->pos() - _lastPos);
            t.setX(-t.x());
            auto sphere = scene.boundingBox().outerSphere();
            if ((e->buttons() & Qt::RightButton) && !(e->modifiers() & Qt::ShiftModifier)) {
                core::Vec3 trans = t.x() * renderOptions.camera.rightward() + t.y() * renderOptions.camera.upward();
                trans *= 0.02;
                renderOptions.camera.moveEyeWithCenterFixed(trans, sphere, true, true);
                setCursor(Qt::ClosedHandCursor);
                update();
            }
            else if ((e->buttons() & Qt::MidButton) ||
                ((e->buttons() & Qt::RightButton) && (e->modifiers() & Qt::ShiftModifier))) {
                core::Vec3 trans = t.x() * renderOptions.camera.rightward() + t.y() * renderOptions.camera.upward();
                trans *= 0.02;
                renderOptions.camera.translate(trans, sphere, true);
                update();
            }
            _lastPos = e->pos();
        }

        virtual void wheelEvent(QWheelEvent * e) override {
            auto sphere = scene.boundingBox().outerSphere();
            double d = e->delta() / 10;
            double dist = core::Distance(renderOptions.camera.eye(), renderOptions.camera.center());
            core::Vec3 trans = d * dist / 1000.0 * renderOptions.camera.forward();
            renderOptions.camera.moveEyeWithCenterFixed(trans, sphere, false, true);
            update();
        }

        virtual void mouseReleaseEvent(QMouseEvent * e) override {
            unsetCursor();
        }

        virtual void keyPressEvent(QKeyEvent * e) override {
            if (e->key() == Qt::Key_Space){
                for (auto & n : scene.tree().nodes()){
                    if (n.exists){
                        for (int entityID : n.data->selectedEntities()){
                            n.data->invokeCallbackFunction(vis::InteractionID::PressSpace,
                                scene.tree(),
                                vis::VisualObjectEntityID{ n.topo.hd, entityID });
                        }
                    }
                }
            }
        }

    private:
        QPointF _lastPos;
        bool _needsInitialization;
        QReadWriteLock _needsInitializationLock;

    private:
        vis::RenderOptions renderOptions;
        //vis::VisualObjectInstallingOptions installingOptions;
        vis::VisualObjectScene scene;
        DataOfType<Reconstruction> * data;
    };

    return new Widget(&rec, parent);
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

            core::SolveVariablesUsingInversedDepths(mg, props);

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
    _steps = new Steps(this);
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