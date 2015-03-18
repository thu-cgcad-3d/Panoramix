#include <QtOpenGL>

#include "../../src/core/algorithms.hpp"
#include "../../src/vis/qt_glue.hpp"
#include "../../src/vis/visualizers.hpp"
#include "configuration.hpp"

using namespace panoramix;



template <class ContentT, class WidgetT>
class StepWidgetAdaptor : public StepWidgetInterface, public WidgetT {
    static_assert(std::is_base_of<QWidget, WidgetT>::value, "WidgetT must be derived from QWidget");
public:
    explicit StepWidgetAdaptor(DataOfType<ContentT> * d, QWidget * parent = nullptr) 
        : WidgetT(parent), _data(d) {
        setMinimumSize(500, 500);
    }
    virtual void refreshDataAsync() {}
    virtual void refreshData() {}
    virtual void updatePainting() { WidgetT::update(); }
    virtual void showWidget() { WidgetT::show(); }
    virtual void hideWidget() { WidgetT::hide(); }

protected:
    bool noData() const { return !_data; }
    DataOfType<ContentT> & data() const { return *_data; }
    void lockForRead() const { _data->lockForRead(); }
    void lockForWrite() const { _data->lockForWrite(); }
    void unlock() const { _data->unlock(); }
    ContentT & content() const { return _data->content; }

protected:
    DataOfType<ContentT> * _data;
};




StepWidgetInterface * CreateBindingWidgetAndActions(DataOfType<PanoView> & pv,
    QList<QAction*> & actions, QWidget * parent) {

    class Widget : public StepWidgetAdaptor<PanoView, QWidget> {
    public:
        explicit Widget(DataOfType<PanoView> * d, QWidget * parent)
            : StepWidgetAdaptor<PanoView, QWidget>(d, parent) {}
    protected:
        virtual void paintEvent(QPaintEvent * e) override {
            if (noData())
                return;
            QPainter painter(this);
            lockForRead();
            painter.drawImage(QPointF(), vis::MakeQImage(content().view.image));
            unlock();
        }
    };

    return new Widget(&pv, parent);
}






StepWidgetInterface * CreateBindingWidgetAndActions(DataOfType<Segmentation> & segs,
    QList<QAction*> & actions, QWidget * parent) {

    class Widget : public StepWidgetAdaptor<Segmentation, QWidget> {
    public:
        explicit Widget(DataOfType<Segmentation> * d, QWidget * parent)
            : StepWidgetAdaptor<Segmentation, QWidget>(d, parent){}
    protected:
        virtual void paintEvent(QPaintEvent * e) override {
            if (!data)
                return;
            QPainter painter(this);
            lockForRead();
            // todo
            unlock();
        }
    private:
        DataOfType<Segmentation>* data;
    };

    return new Widget(&segs, parent);
}







StepWidgetInterface * CreateBindingWidgetAndActions(DataOfType<LinesAndVPs> & segs,
    QList<QAction*> & actions, QWidget * parent){


    return nullptr;
}








StepWidgetInterface * CreateBindingWidgetAndActions(DataOfType<ReconstructionSetup> & segs,
    QList<QAction*> & actions, QWidget * parent){

    class Widget : public StepWidgetAdaptor<ReconstructionSetup, QWidget> {
    public:
        explicit Widget(DataOfType<ReconstructionSetup> * d, QWidget * parent) 
            : StepWidgetAdaptor<ReconstructionSetup, QWidget>(d, parent) {}

    protected:
        virtual void paintEvent(QPaintEvent * e) override {
            if (!data)
                return;
            QPainter painter(this);
            painter.setOpacity(1.0);
            if (_showIm){
                painter.drawImage(QPoint(), _im);
                painter.setOpacity(painter.opacity()*.5);
            }
            if (_showSeg)
                painter.drawImage(QPoint(), _segShow);
            painter.setOpacity(1.0);
            if (_showLabel){ painter.drawImage(QPoint(), _labelBg); }
        }

        void mousePressEvent(QMouseEvent *e) {
            if (_im.isNull())
                return;

            //static const QString setLabelTempl = "seglabel(%1)=%2";
            //static const QString getLabelTempl = "temp = seglabel(%1)";

            //int segIdx = segData_.pixel(e->pos().x(), e->pos().y());
            //matlab_ << getLabelTempl.arg(segIdx).toAscii();
            //int oldLabel = matlab_.get("temp").constReal()[0];

            ////if(labelLocks[oldLabel-1])
            ////	return;

            //redoCommands_.insert(commandIndex_, setLabelTempl.arg(segIdx).arg(currentLabel_));
            //undoCommands_.insert(commandIndex_, setLabelTempl.arg(segIdx).arg(oldLabel));
            //matlab_ << redoCommands_[commandIndex_].toAscii();
            //commandIndex_++;

            //refreshMasks();
            //update();

            //emit canUndoChanged(true);
            //emit canRedoChanged(false);
            //setWindowModified(true);
        }


    private:
        DataOfType<ReconstructionSetup>* data;
        bool _showIm, _showSeg, _showLabel;
        QImage _im;
        QImage _segShow;
        QImage _segData;
        QImage _labelBg;
    };


    return new Widget(&segs, parent);
}







StepWidgetInterface * CreateBindingWidgetAndActions(DataOfType<Reconstruction> & rec,
    QList<QAction*> & actions, QWidget * parent) {

    class Widget : public StepWidgetAdaptor<Reconstruction, QGLWidget> {
    public:
        explicit Widget(DataOfType<Reconstruction> * d, QWidget * parent = nullptr)
            : StepWidgetAdaptor<Reconstruction, QGLWidget>(d, parent){
            setMouseTracking(true);
            setAutoBufferSwap(false);
            grabKeyboard();
            //setUpScene();
            setMinimumSize(300, 300);
            _needsInitialization = true;
        }

        virtual void refreshDataAsync() override {
            if (noData()){
                scene.clear();
                return;
            }

            struct ComponentID {
                int handleID;
                bool isRegion;
            };

            vis::Visualizer viz("mixed graph optimizable");
            viz.installingOptions.discretizeOptions.colorTable = vis::ColorTableDescriptor::RGB;
            std::vector<std::pair<ComponentID, vis::Colored<vis::SpatialProjectedPolygon>>> spps;
            std::vector<vis::Colored<core::Line3>> lines;

            auto & mg = content().mg;
            auto & props = content().props;

            lockForRead();
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
            unlock();

            lockForRead();
            vis::ResourceStore::set("texture", content().view.image);
            unlock();

            viz.begin(spps).shaderSource(vis::OpenGLShaderSourceDescriptor::XPanorama).resource("texture").end();
            viz.installingOptions.discretizeOptions.color = vis::ColorTag::DarkGray;
            viz.installingOptions.lineWidth = 5.0;
            viz.add(lines);

            viz.installingOptions.discretizeOptions.color = vis::ColorTag::Black;
            viz.installingOptions.lineWidth = 1.0;
            //viz.add(connectionLines);

            viz.renderOptions.renderMode = vis::RenderModeFlag::Triangles /*| vis::RenderModeFlag::Lines*/;
            viz.renderOptions.backgroundColor = vis::ColorTag::White;
            viz.renderOptions.bwColor = 0.0;
            viz.renderOptions.bwTexColor = 1.0;
            viz.camera(core::PerspectiveCamera(1000, 800, 800, core::Point3(-1, 1, 1), core::Point3(0, 0, 0)));
            renderOptions = viz.renderOptions;

            scene.install(std::move(viz.tree()));

            _needsInitializationLock.lockForWrite();
            _needsInitialization = true;
            _needsInitializationLock.unlock();
        }

        virtual void refreshData() override { 
            initializeGL();
        }

    public:

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
        vis::VisualObjectScene scene;
        QReadWriteLock sceneLock;
    };

    return new Widget(&rec, parent);
}









