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



template <class ContentT>
class ImageViewer : public StepWidgetAdaptor<ContentT, QWidget> {
    using Base = StepWidgetAdaptor<ContentT, QWidget>;

public:
    explicit ImageViewer(DataOfType<ContentT> * d, QWidget * parent)
        : StepWidgetAdaptor<ContentT, QWidget>(d, parent), _scale(1.0) { setMouseTracking(true); }

protected:
    virtual void paintEvent(QPaintEvent * e) override {
        QPainter painter(this);
        painter.setBackground(Qt::white);
        painter.eraseRect(rect());

        if (noData())
            return;

        _imageLock.lockForRead();

        painter.resetTransform();
        painter.translate(rect().center() + _translate);
        painter.scale(_scale, _scale);
        painter.translate(-_image.rect().center());

        painter.fillRect(_image.rect().translated(QPoint(5, 5) / _scale), QColor(35, 30, 30, 100));
        painter.drawImage(_image.rect(), _image);
        _imageLock.unlock();

        painter.setPen(_strokePen);
        painter.drawPolyline(_stroke);
    }

    virtual void mousePressEvent(QMouseEvent * e) override {
        _lastPos = e->pos();
        if (e->buttons() & Qt::LeftButton || e->buttons() & Qt::MiddleButton){
            setCursor(Qt::SizeAllCursor);
        }
       /* if (e->buttons() & Qt::RightButton){
            auto p = positionOnImage(e->pos()).toPoint();
            if (_image.rect().contains(p))
                _image.setPixel(p, qRgb(0, 0, 0));
        }*/
        update();
        Base::mousePressEvent(e);
    }

    virtual void mouseMoveEvent(QMouseEvent * e) override {
        if (e->buttons() & Qt::LeftButton || e->buttons() & Qt::MiddleButton) {
            QPoint t(e->pos() - _lastPos);
            _translate += t;
            update();
        }
        _lastPos = e->pos();
        Base::mouseMoveEvent(e);
    }

    virtual void wheelEvent(QWheelEvent * e) override {
        double d = e->delta() / 500.0;
        if (_scale >= 50.0 && d >= 0)
            return;
        if (_scale <= 0.02 && d <= 0)
            return;
        _scale *= std::pow(2.0, d);
        update();
        Base::wheelEvent(e);
    }

    virtual void mouseReleaseEvent(QMouseEvent * e) override {
        unsetCursor();
        Base::mouseReleaseEvent(e);
    }

protected:
    inline QPointF positionOnImage(const QPointF & screenPos) const {
        // (x - imcenter) * scale + rect.center + translate
        return (screenPos - _translate - rect().center()) / _scale + _image.rect().center();
    }

protected:
    QReadWriteLock _imageLock;
    QImage _image;
    double _scale;
    QPoint _translate;
    QPoint _lastPos;
    QPolygonF _stroke;
    QPen _strokePen;
};

template <class ContentT, class RendererT>
inline StepWidgetInterface * CreateImageViewer(
    DataOfType<ContentT> * d, RendererT && r, QWidget * parent = nullptr) {

    using RendererType = std::decay_t<RendererT>;
    class Widget : public ImageViewer<ContentT> {
    public:
        explicit Widget(DataOfType<ContentT> * dd, RendererT && rr, QWidget * p) 
            : ImageViewer<ContentT>(dd, p), _renderer(std::forward<RendererT>(rr)) {}
        virtual void refreshDataAsync() {
            if (noData())
                return;
            QImage im = _renderer(data());
            _imageLock.lockForWrite();
            _image = im;
            _imageLock.unlock();
        }
    private:
        RendererType _renderer;
    };

    return new Widget(d, std::forward<RendererT>(r), parent);
}



StepWidgetInterface * CreateBindingWidgetAndActions(DataOfType<PanoView> & pv,
    QList<QAction*> & actions, QWidget * parent) {

    return CreateImageViewer(&pv, [](DataOfType<PanoView>& pv){
        pv.lockForRead();
        QImage im = vis::MakeQImage(pv.content.view.image);
        pv.unlock();
        return im;
    }, parent);
}






StepWidgetInterface * CreateBindingWidgetAndActions(DataOfType<Segmentation> & segs,
    QList<QAction*> & actions, QWidget * parent) {
    
    return CreateImageViewer(&segs, [](DataOfType<Segmentation> & segs){        
        segs.lockForRead();
        auto colorTable = vis::CreateRandomColorTableWithSize(segs.content.segmentsNum);
        core::Imageub3 im = colorTable(segs.content.segmentation);
        segs.unlock();
        return vis::MakeQImage(im);    
    }, parent);
}







StepWidgetInterface * CreateBindingWidgetAndActions(DataOfType<LinesAndVPs> & segs,
    QList<QAction*> & actions, QWidget * parent){
    return nullptr;
}








StepWidgetInterface * CreateBindingWidgetAndActions(DataOfType<ReconstructionSetup> & rec,
    QList<QAction*> & actions, QWidget * parent){

    class Widget : public ImageViewer<ReconstructionSetup> {

        void setNoPen() { _strokePen = QPen(Qt::transparent); _applyStroke = nullptr; }
        void setPen(const QColor & color, int strokeSize, 
            bool regionPropUsed, int regionPropOrientationClaz, int regionPropOrientationNotClaz) {
            _strokePen = QPen(color, strokeSize);
            _applyStroke = [=]() -> bool {
                bool changed = false;
                auto & recSetup = data();
                recSetup.lockForWrite();
                for (QPointF & p : _stroke){
                    int sz = strokeSize;
                    for (int dx = -sz; dx <= sz; dx++){
                        for (int dy = -sz; dy <= sz; dy++){
                            core::PixelLoc pp(p.x() + dx, p.y() + dy);
                            if (!core::Contains(recSetup.content.segmentation, pp))
                                continue;
                            int segId = recSetup.content.segmentation(pp);
                            core::RegionHandle rh(segId);
                            auto & prop = recSetup.content.props[rh];
                            if (prop.used == regionPropUsed && 
                                prop.orientationClaz == regionPropOrientationClaz && 
                                prop.orientationNotClaz == regionPropOrientationNotClaz)
                                continue;
                            changed = true;
                            prop.used = regionPropUsed;
                            prop.orientationClaz = regionPropOrientationClaz;
                            prop.orientationNotClaz = regionPropOrientationNotClaz;
                        }
                    }
                }               
                core::ResetVariables(recSetup.content.mg, recSetup.content.props);
                recSetup.unlock();
                return changed;
            };            
        }

    public:
        explicit Widget(DataOfType<ReconstructionSetup> * dd, QWidget * p) : ImageViewer<ReconstructionSetup>(dd, p) {            

            setContextMenuPolicy(Qt::ActionsContextMenu);
            QActionGroup * bas = new QActionGroup(this);
            
            QAction * defaultAction = nullptr;
            connect(defaultAction = bas->addAction(tr("No Brush")), &QAction::triggered, [this](){setNoPen(); });
            {
                QAction * sep = new QAction(bas);
                sep->setSeparator(true);
                bas->addAction(sep);
            }
            connect(bas->addAction(tr("Paint Void")), &QAction::triggered, [this](){ setPen(Qt::black, 3, false, -1, -1); });
            connect(bas->addAction(tr("Paint Free Plane")), &QAction::triggered, [this](){ setPen(Qt::gray, 3, true, -1, -1); });
            connect(bas->addAction(tr("Paint Toward Direction 1")), &QAction::triggered, [this](){ setPen(Qt::red, 3, true, 0, -1); });
            connect(bas->addAction(tr("Paint Toward Direction 2")), &QAction::triggered, [this](){ setPen(Qt::red, 3, true, 1, -1);  });
            connect(bas->addAction(tr("Paint Toward Direction 3")), &QAction::triggered, [this](){ setPen(Qt::red, 3, true, 2, -1); });
            connect(bas->addAction(tr("Paint Along Direction 1")), &QAction::triggered, [this](){ setPen(Qt::red, 3, true, -1, 0); });
            connect(bas->addAction(tr("Paint Along Direction 2")), &QAction::triggered, [this](){ setPen(Qt::red, 3, true, -1, 1); });
            connect(bas->addAction(tr("Paint Along Direction 3")), &QAction::triggered, [this](){ setPen(Qt::red, 3, true, -1, 2); });
            {
                QAction * sep = new QAction(bas);
                sep->setSeparator(true);
                bas->addAction(sep);
            }
            connect(bas->addAction(tr("Draw Occlusion")), &QAction::triggered, [this](){}); 
            connect(bas->addAction(tr("Remove Occlusion")), &QAction::triggered, [this](){});

            for (auto a : bas->actions()) a->setCheckable(true);

            bas->setExclusive(true);
            defaultAction->setChecked(true);

            addActions(bas->actions());
        }
        
    public:
        virtual void refreshDataAsync() {
            if (noData())
                return;

            auto & rec = data();
            rec.lockForRead();
            QImage transparentBg(rec.content.segmentation.cols, rec.content.segmentation.rows, QImage::Format::Format_RGB888);
            {
                transparentBg.fill(Qt::transparent);
                QPainter painter(&transparentBg);
                painter.fillRect(transparentBg.rect(), Qt::BrushStyle::HorPattern);
            }

            // render
            core::Imageub3 rendered = core::Imageub3::zeros(rec.content.segmentation.size());
            vis::ColorTable rgb = { vis::ColorTag::Red, vis::ColorTag::Green, vis::ColorTag::Blue };
            vis::ColorTable ymc = { vis::ColorTag::Yellow, vis::ColorTag::Magenta, vis::ColorTag::Cyan };
            double alpha = 0.3;
            for (auto it = rendered.begin(); it != rendered.end(); ++it){
                int regionId = rec.content.segmentation(it.pos());
                Q_ASSERT(regionId >= 0 && regionId < rec.content.mg.internalComponents<core::RegionData>().size());
                auto & prop = rec.content.props[core::RegionHandle(regionId)];
                if (!prop.used){
                    QRgb transPixel = transparentBg.pixel(it.pos().x, it.pos().y);
                    *it = core::Vec3b(qRed(transPixel), qGreen(transPixel), qBlue(transPixel));
                    continue;
                }
                vis::Color imColor = vis::ColorFromImage(rec.content.view.image, it.pos());
                if (prop.orientationClaz >= 0 && prop.orientationNotClaz == -1){
                    *it = imColor.blendWith(rgb[prop.orientationClaz], alpha);
                }
                else if (prop.orientationClaz == -1 && prop.orientationNotClaz >= 0){
                    *it = imColor.blendWith(ymc[prop.orientationNotClaz], alpha);
                }
                else{
                    *it = imColor;
                }
            }
            // render disconnected boundaries
            for (auto & b : rec.content.mg.constraints<core::RegionBoundaryData>()){
                if (rec.content.props[b.topo.hd].used)
                    continue;
                for (auto & e : b.data.normalizedEdges){
                    if (e.size() <= 1) continue;
                    for (int i = 1; i < e.size(); i++){
                        auto p1 = core::ToPixelLoc(rec.content.view.camera.screenProjection(e[i - 1]));
                        auto p2 = core::ToPixelLoc(rec.content.view.camera.screenProjection(e[i]));
                        cv::clipLine(cv::Rect(0, 0, rendered.cols, rendered.rows), p1, p2);
                        cv::line(rendered, p1, p2, vis::Color(vis::ColorTag::Black), 3);
                    }
                }
            }
            rec.unlock();
            _imageLock.lockForWrite();
            _image = vis::MakeQImage(rendered);
            _imageLock.unlock();
        }

    protected:
        virtual void mousePressEvent(QMouseEvent * e) override {
            if (_applyStroke && e->buttons() & Qt::LeftButton){
                _stroke << positionOnImage(e->pos());
                setCursor(Qt::CursorShape::OpenHandCursor);
                update();
            }
            else{
                ImageViewer<ReconstructionSetup>::mousePressEvent(e);
            }
        }

        virtual void mouseMoveEvent(QMouseEvent * e) override {
            if (_applyStroke && e->buttons() & Qt::LeftButton){
                auto pos = positionOnImage(e->pos());
                setCursor(Qt::CursorShape::ClosedHandCursor);
                if (_stroke.isEmpty() || (_stroke.last() - pos).manhattanLength() > 3){
                    _stroke << pos;
                    update();
                }
            }
            else{
                ImageViewer<ReconstructionSetup>::mouseMoveEvent(e);
            }
        }

        virtual void mouseReleaseEvent(QMouseEvent * e) override {
            if (_applyStroke && _applyStroke()){
                refreshDataAsync();
                refreshData();
                data().lockForWrite();
                data().setModified();
                data().unlock();
            }
            _stroke.clear();
            update();
            ImageViewer<ReconstructionSetup>::mouseReleaseEvent(e);
        }

    private:
        std::function<bool(void)> _applyStroke;
    };

    return new Widget(&rec, parent);
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
            setMinimumSize(600, 600);
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

            qDebug() << "loading texture";
            viz.begin(spps).shaderSource(vis::OpenGLShaderSourceDescriptor::XPanorama).resource("texture").end();
            qDebug() << "texture loaded";

            viz.installingOptions.discretizeOptions.color = vis::ColorTag::DarkGray;
            viz.installingOptions.lineWidth = 5.0;
            viz.add(lines);

            viz.installingOptions.discretizeOptions.color = vis::ColorTag::Black;
            viz.installingOptions.lineWidth = 1.0;

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









