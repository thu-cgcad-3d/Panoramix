#pragma once

#include <QtWidgets>

#include "../../src/gui/basic_types.hpp"
#include "../../src/gui/visualizers.hpp"
#include "stepsdag.hpp"

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





template <class ContentT, class SceneConstructorT>
class SceneViewer : public StepWidgetAdaptor<ContentT, QGLWidget> {
    using Base = StepWidgetAdaptor<ContentT, QGLWidget>;
public:
    explicit SceneViewer(DataOfType<ContentT> * d, SceneConstructorT && makeScene, QWidget * parent)
        : StepWidgetAdaptor<ContentT, QGLWidget>(d, parent), _makeScene(std::move(makeScene)) {
            setMouseTracking(true);
            setAutoBufferSwap(false);
            grabKeyboard();
            setMinimumSize(600, 600);
            _needsInitialization.component = true;
        }

    virtual void refreshDataAsync() override {
        if (noData()){
            _scene.lockForWrite();
            _scene.component.clear();
            _scene.unlock();
            return;
        }
        gui::Visualizer viz;
        _makeScene(data(), viz);
        _renderOptions = viz.renderOptions;

        _scene.lockForWrite();
        _scene.component.install(std::move(viz.tree()));
        _scene.unlock();

        _needsInitialization.lockForWrite();
        _needsInitialization.component = true;
        _needsInitialization.unlock();
    }

    virtual void refreshData() override {
        initializeGL();
    }

public:
    void autoSetCamera() {
        _scene.lockForRead();
        auto sphere = _scene.component.boundingBox().outerSphere();
        _scene.unlock();
        _renderOptions.camera.resizeScreen(core::Size(width(), height()), false);
        _renderOptions.camera.focusOn(sphere, true);
        update();
    }

protected:
    void initializeGL() {
        if (!_needsInitialization.component)
            return;
        makeCurrent();
        glEnable(GL_MULTISAMPLE);
        GLint bufs;
        GLint samples;
        glGetIntegerv(GL_SAMPLE_BUFFERS, &bufs);
        glGetIntegerv(GL_SAMPLES, &samples);
        qDebug("Have %d buffers and %d samples", bufs, samples);
        qglClearColor(MakeQColor(_renderOptions.backgroundColor));

        _scene.lockForWrite();
        _scene.component.initialize();
        _scene.unlock();

        _needsInitialization.lockForWrite();
        _needsInitialization.component = false;
        _needsInitialization.unlock();
    }

    void paintGL() {
        initializeGL();
        QPainter painter;
        painter.begin(this);
        painter.beginNativePainting();
        qglClearColor(MakeQColor(_renderOptions.backgroundColor));

        glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);

        core::PerspectiveCamera & camera = _renderOptions.camera;
        camera.resizeScreen(core::Size(width(), height()));

        _scene.lockForRead();
        _scene.component.render(_renderOptions);
        _scene.unlock();

        painter.endNativePainting();
        swapBuffers();
    }

    void resizeGL(int w, int h) {
        core::PerspectiveCamera & camera = _renderOptions.camera;
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
            gui::VisualObjectHandle oh;
            gui::TriMesh::TriangleHandle t;
            _scene.lockForWrite();
            std::tie(oh, t) = _scene.component.pickOnScreen(_renderOptions, core::Point2(e->pos().x(), e->pos().y()));
            if (oh.valid()){
                int entityID = _scene.component.tree().data(oh)->entityIDOfMeshTriangle(t);
                if (e->modifiers() & Qt::ControlModifier){
                    _scene.component.switchSelect(std::make_pair(oh, entityID));
                }
                else{
                    _scene.component.clearSelection();
                    _scene.component.select(std::make_pair(oh, entityID));
                }
                _scene.component.tree().data(oh)
                    ->invokeCallbackFunction(gui::InteractionID::ClickLeftButton, _scene.component.tree(), std::make_pair(oh, entityID));
            }
            else if (!(e->modifiers() & Qt::ControlModifier)){
                _scene.component.clearSelection();
            }
            _scene.unlock();
        }
        update();
    }

    virtual void mouseMoveEvent(QMouseEvent * e) override {
        QVector3D t(e->pos() - _lastPos);
        t.setX(-t.x());
        _scene.lockForRead();
        auto sphere = _scene.component.boundingBox().outerSphere();
        _scene.unlock();
        if ((e->buttons() & Qt::RightButton) && !(e->modifiers() & Qt::ShiftModifier)) {
            core::Vec3 trans = t.x() * _renderOptions.camera.rightward() + t.y() * _renderOptions.camera.upward();
            trans *= 0.02;
            _renderOptions.camera.moveEyeWithCenterFixed(trans, sphere, true, true);
            setCursor(Qt::ClosedHandCursor);
            update();
        }
        else if ((e->buttons() & Qt::MidButton) ||
            ((e->buttons() & Qt::RightButton) && (e->modifiers() & Qt::ShiftModifier))) {
            core::Vec3 trans = t.x() * _renderOptions.camera.rightward() + t.y() * _renderOptions.camera.upward();
            trans *= 0.02;
            _renderOptions.camera.translate(trans, sphere, true);
            update();
        }
        _lastPos = e->pos();
    }

    virtual void wheelEvent(QWheelEvent * e) override {
        _scene.lockForRead();
        auto sphere = _scene.component.boundingBox().outerSphere();
        _scene.unlock();
        double d = e->delta() / 10;
        double dist = core::Distance(_renderOptions.camera.eye(), _renderOptions.camera.center());
        core::Vec3 trans = d * dist / 1000.0 * _renderOptions.camera.forward();
        _renderOptions.camera.moveEyeWithCenterFixed(trans, sphere, false, true);
        update();
    }

    virtual void mouseReleaseEvent(QMouseEvent * e) override {
        unsetCursor();
    }

    virtual void keyPressEvent(QKeyEvent * e) override {
        if (e->key() == Qt::Key_Space){
            _scene.lockForRead();
            for (auto & n : _scene.component.tree().nodes()){
                if (n.exists){
                    for (int entityID : n.data->selectedEntities()){
                        n.data->invokeCallbackFunction(gui::InteractionID::PressSpace,
                            _scene.component.tree(),
                            gui::VisualObjectEntityID{ n.topo.hd, entityID });
                    }
                }
            }
            _scene.unlock();
        }
    }

private:
    QPointF _lastPos;
    LockableType<bool> _needsInitialization;

private:
    gui::RenderOptions _renderOptions;
    LockableType<gui::VisualObjectScene> _scene;
    SceneConstructorT _makeScene;
};

template <class ContentT, class SceneConstructorT>
inline StepWidgetInterface * CreateSceneViewer(
    DataOfType<ContentT> * d, SceneConstructorT && c, QWidget * parent = nullptr) {
    using SceneConstructorType = std::decay_t<SceneConstructorT>;
    return new SceneViewer<ContentT, SceneConstructorType>(d, std::forward<SceneConstructorT>(c), parent);
}











