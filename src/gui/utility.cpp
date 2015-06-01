#include <QtGui>
#include <QtOpenGL>
#include <QtWidgets>

#include "singleton.hpp"
#include "qt_glue.hpp"
#include "utility.hpp"
#include "scene.hpp"

namespace panoramix {
    namespace gui {

        int SelectFrom(const std::vector<std::string> & strs,
            const std::string & title, const std::string & text,
            int acceptId, int rejectId){
            Singleton::InitGui();
            QMessageBox mbox;
            std::vector<QAbstractButton*> buttons(strs.size(), nullptr);
            for (int i = 0; i < strs.size(); i++){
                buttons[i] = new QPushButton(QString::fromStdString(strs[i]));
                auto role = i == acceptId ? QMessageBox::ButtonRole::AcceptRole : 
                    (i == rejectId ? QMessageBox::ButtonRole::RejectRole : QMessageBox::ButtonRole::NoRole);
                mbox.addButton(buttons[i], role);
            }
            mbox.setWindowTitle(title.empty() ? QObject::tr("Make your decision") : QString::fromStdString(title));
            mbox.setText(text.empty() ? QObject::tr("Click on one of the buttons") : QString::fromStdString(text));
            mbox.exec();
            for (int i = 0; i < buttons.size(); i++){
                if (mbox.clickedButton() == buttons[i])
                    return i;
            }
            return -1;
        }

        core::Image PickAnImage(const std::string & dir, std::string * picked){
            Singleton::InitGui();
            auto filename = QFileDialog::getOpenFileName(nullptr, QObject::tr("Select an image file"), 
                QString::fromStdString(dir), 
                QObject::tr("Image Files (*.png;*.jpg;*.jpeg);;All Files (*.*)"));
            if (filename.isEmpty())
                return core::Image();
            if (picked){
                *picked = filename.toStdString();
            }
            return cv::imread(filename.toStdString());
        }

        std::vector<core::Image> PickImages(const std::string & dir){
            Singleton::InitGui();
            auto filenames = QFileDialog::getOpenFileNames(nullptr, QObject::tr("Select an image file"),
                QString::fromStdString(dir),
                QObject::tr("Image Files (*.png;*.jpg;*.jpeg);;All Files (*.*)"));
            std::vector<core::Image> ims;
            for (auto & filename : filenames){
                ims.push_back(cv::imread(filename.toStdString()));
            }
            return ims;
        }

        std::vector<core::Image> PickAllImagesFromAFolder(const std::string & dir){
            Singleton::InitGui();
            auto folder = QFileDialog::getExistingDirectory(nullptr, QObject::tr("Select a folder containing images"),
                QString::fromStdString(dir));
            NOT_IMPLEMENTED_YET();
        }

        Qt::PenStyle MakeQPenStyle(PenStyle ps) {
            return Qt::PenStyle(ps);
        }

        template <class BaseWidgetT>
        class PaintableWidget : public BaseWidgetT {
        public:
            explicit PaintableWidget(const std::vector<PenConfig> & pc, QWidget * parent = nullptr) 
                : BaseWidgetT(parent), _activePenId(-1), _penCursor(QPixmap(tr(":/icons/pencil_icon&16.png")), 0, 16), _pens(pc.size()) {

                setMouseTracking(true);

                // setup pen selection actions
                setContextMenuPolicy(Qt::ActionsContextMenu);
                QActionGroup * bas = new QActionGroup(this);

                QAction * defaultAction = nullptr;
                connect(defaultAction = bas->addAction(tr("No Brush")), &QAction::triggered, [this](){ _activePenId = -1; });
                {
                    QAction * sep = new QAction(bas);
                    sep->setSeparator(true);
                    bas->addAction(sep);
                }
                for (int i = 0; i < pc.size(); i++){
                    auto action = bas->addAction(QString::fromStdString(pc[i].name));
                    action->setToolTip(QString::fromStdString(pc[i].description));
                    action->setWhatsThis(QString::fromStdString(pc[i].description));
                    // draw icon
                    int sz = 16;
                    QImage image(sz, sz, QImage::Format::Format_ARGB32_Premultiplied);
                    image.fill(MakeQColor(gui::White));
                    _pens[i] = QPen(MakeQColor(pc[i].color), pc[i].thickness, MakeQPenStyle(pc[i].style));
                    QPainter painter(&image);
                    painter.setPen(_pens[i]);
                    painter.drawLine(QPointF(0, sz / 2), QPointF(sz, sz/2));
                    painter.end();
                    action->setIcon(QIcon(QPixmap::fromImage(image)));
                    connect(action, &QAction::triggered, [this, i] { _activePenId = i; });
                }

                for (auto a : bas->actions()) a->setCheckable(true);

                bas->setExclusive(true);
                defaultAction->setChecked(true);

                addActions(bas->actions());
            }

        protected:
            void mousePressEvent(QMouseEvent * e) override {
                _moved = false;
                if (e->buttons() & Qt::LeftButton){
                    if (_activePenId == -1)
                        return;
                    _stroke << e->pos();
                    setCursor(_penCursor);
                }
                else{
                    QWidget::mousePressEvent(e);
                }
            }

            void mouseMoveEvent(QMouseEvent * e) override {
                _moved = true;
                if (e->buttons() & Qt::LeftButton){
                    if (_activePenId == -1)
                        return;
                    auto pos = e->pos();
                    if (_stroke.isEmpty() || (_stroke.last() - pos).manhattanLength() > 3){
                        _stroke << pos;
                        update();
                    }
                }
                else{
                    QWidget::mouseMoveEvent(e);
                }
            }

            void mouseReleaseEvent(QMouseEvent * e) override {
                bool isClick = !_moved;
                _moved = false;
                unsetCursor();
                if (_activePenId == -1)
                    return;
                std::vector<core::Point2> points(_stroke.size());
                for (int i = 0; i < _stroke.size(); i++){
                    points[i][0] = _stroke[i].x();
                    points[i][1] = _stroke[i].y();
                }
                if (!points.empty()){
                    processStroke(points, _activePenId);
                    update();
                }
                _stroke.clear();
                update();
            }

            virtual void processStroke(const std::vector<core::Point2> & stroke, int penId) = 0;

        protected:
            QPolygonF _stroke;
            int _activePenId;
            std::vector<QPen> _pens;
            std::vector<PenConfig> _penConfigs;
            bool _moved;
            QCursor _penCursor;
        };



        void PaintWith(const std::function<core::Image()> & updater,
            const std::vector<PenConfig> & penConfigs,
            const std::function<bool(const std::vector<core::Point2> & polyline, int penId)> & callback){

            Singleton::InitGui();

            class Widget : public PaintableWidget<QWidget> {
                using BaseClass = PaintableWidget<QWidget>;
            public:
                explicit Widget(const std::function<core::Image()> & up,
                    const std::vector<PenConfig> & pc,
                    const std::function<bool(const std::vector<core::Point2> &, int)> & cb) 
                    : BaseClass(pc), _updater(up), _callback(cb) {
                    _scale = 1.0;
                    updateImage();
                }

                void updateImage() {
                    auto im = _updater();
                    _image = MakeQImage(im);
                }

            protected:
                void paintEvent(QPaintEvent * e) override {
                    QPainter painter(this);
                    painter.setBackground(Qt::white);
                    painter.eraseRect(rect());

                    painter.resetTransform();
                    painter.translate(rect().center() + _translate);
                    painter.scale(_scale, _scale);
                    painter.translate(-_image.rect().center());

                    painter.fillRect(_image.rect().translated(QPoint(5, 5) / _scale), QColor(35, 30, 30, 100));
                    painter.drawImage(_image.rect(), _image);
                    painter.resetTransform();

                    if (_activePenId >= 0){
                        auto pen = _pens[_activePenId];
                        painter.setPen(pen);                        
                    }
                    painter.setRenderHint(QPainter::Antialiasing);
                    painter.drawPolyline(_stroke);
                }

                void mousePressEvent(QMouseEvent * e) override {
                    if (e->buttons() & Qt::MidButton){
                        // move
                        _lastPos = e->pos();
                        setCursor(Qt::OpenHandCursor);
                        update();
                    }
                    else{
                        BaseClass::mousePressEvent(e);
                    }
                }

                void mouseMoveEvent(QMouseEvent * e) override {
                    if (e->buttons() & Qt::MidButton){
                        setCursor(Qt::ClosedHandCursor);
                        QPoint t(e->pos() - _lastPos);
                        _translate += t;
                        _lastPos = e->pos();
                        update();
                    }
                    else{
                        BaseClass::mouseMoveEvent(e);
                    }
                }

                void wheelEvent(QWheelEvent * e) override {
                    double d = e->delta() / 500.0;
                    if (_scale >= 50.0 && d >= 0)
                        return;
                    if (_scale <= 0.02 && d <= 0)
                        return;
                    _scale *= std::pow(2.0, d);
                    update();
                }

                void keyPressEvent(QKeyEvent * e) override {
                    if (e->key() == Qt::Key::Key_F5){
                        updateImage();
                        update();
                    }
                }

                inline QPointF positionOnImage(const QPointF & screenPos) const {
                    // (x - imcenter) * scale + rect.center + translate
                    return (screenPos - _translate - rect().center()) / _scale + _image.rect().center();
                }

                void processStroke(const std::vector<core::Point2> & stroke, int penId) {
                    if (_callback(stroke, penId)){
                        updateImage();
                        update();
                    }                    
                }

            private:
                QImage _image;
                double _scale;
                QPointF _translate;
                QPoint _lastPos;
                std::function<core::Image()> _updater;
                std::function<bool(const std::vector<core::Point2> & polyline, int penId)> _callback;
            };

            Widget w(updater, penConfigs, callback);
            w.show();

            Singleton::ContinueGui();
        }


        void VisualizeWithPanoramicOperation(const Scene & scene, const RenderOptions & options){

            using namespace panoramix::core;

            Singleton::InitGui();

            class Widget : public PaintableWidget<QGLWidget> {
                using BaseClass = PaintableWidget<QGLWidget>;
            public:
                explicit Widget(const Scene & scene, const RenderOptions & options,
                    QWidget * parent = nullptr)
                    : BaseClass({}, parent), _scene(scene), _options(options) {

                    setMouseTracking(true);
                    setAutoBufferSwap(false);
                    grabKeyboard();
                }

            protected:
                void initializeGL() {
                    makeCurrent();
                    glEnable(GL_MULTISAMPLE);
                    GLint bufs;
                    GLint samples;
                    glGetIntegerv(GL_SAMPLE_BUFFERS, &bufs);
                    glGetIntegerv(GL_SAMPLES, &samples);
                    qDebug("Have %d buffers and %d samples", bufs, samples);
                    qglClearColor(MakeQColor(_options.backgroundColor()));
                    _scene.initialize();
                }

                void resizeGL(int w, int h) {
                    core::PerspectiveCamera & camera = _options.camera();
                    camera.resizeScreen(core::Size(w, h));
                    glViewport(0, 0, w, h);
                }

                void paintEvent(QPaintEvent * e) override {
                    QPainter painter(this);
                    painter.beginNativePainting();
                    qglClearColor(MakeQColor(_options.backgroundColor()));
                    glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);

                    core::PerspectiveCamera & camera = _options.camera();
                    camera.resizeScreen(core::Size(width(), height()));

                    _scene.render(_options);

                    painter.endNativePainting();

                    if (_activePenId >= 0){
                        auto pen = _pens[_activePenId];
                        painter.setPen(pen);
                    }
                    painter.setRenderHint(QPainter::Antialiasing);
                    painter.drawPolyline(_stroke);
                    swapBuffers();
                }

                void mousePressEvent(QMouseEvent * e) {
                    if (e->buttons() & Qt::MidButton){
                        _lastPos = e->pos();
                        setCursor(Qt::OpenHandCursor);
                    }
                    else {
                        BaseClass::mousePressEvent(e);
                    }
                }

                void mouseMoveEvent(QMouseEvent * e) {
                    QVector3D t(e->pos() - _lastPos);
                    t.setX(-t.x());
                    if (e->buttons() & Qt::MidButton){
                        _options.camera().moveCenterWithEyeFixed(MakeCoreVec(t));
                        setCursor(Qt::ClosedHandCursor);
                        update();
                    }
                    else{
                        BaseClass::mouseMoveEvent(e);
                    }
                    _lastPos = e->pos();
                }

                void wheelEvent(QWheelEvent * e) {
                    _options.camera().setFocal(_options.camera().focal() * exp(e->delta() / 1000.0));
                    BaseClass::wheelEvent(e);
                    update();
                }

                void processStroke(const std::vector<core::Point2> & stroke, int penId) {
                }

            private:
                QPoint _lastPos;
                const Scene & _scene;
                RenderOptions _options;
            };

            Widget w(scene, options);
            w.show();

            Singleton::ContinueGui();

        }



        void PaintWithPanorama(const core::PanoramicView & view,
            const std::vector<PenConfig> & penConfigs,
            const std::function<bool(const std::vector<core::Point2> & polyline, int penId)> & callback){

            using namespace panoramix::core;

            Singleton::InitGui();

            class Widget : public PaintableWidget<QGLWidget> {
                using BaseClass = PaintableWidget<QGLWidget>;
            public:
                explicit Widget(const std::vector<PenConfig> & pc, const PanoramicView & v, 
                    const std::function<bool(const std::vector<core::Point2> & polyline, int penId)> & callback, 
                    QWidget * parent = nullptr)
                    : BaseClass(pc, parent), _view(v), _callback(callback) {

                    setMouseTracking(true);
                    setAutoBufferSwap(false);
                    grabKeyboard();


                    // build scene
                    SceneBuilder sb;
                    Sphere3 sp;
                    sp.center = Origin();
                    sp.radius = 1.0;
                    ResourceStore::set("tex", v.image);
                    sb.begin(sp).shaderSource(gui::OpenGLShaderSourceDescriptor::XPanorama).resource("tex").end();
                    _scene = sb.scene();

                    _options.camera() = PerspectiveCamera(800, 800);
                    _options.camera().setEye(Point3(0, 0, 0));
                    _options.camera().setCenter(v.camera.center());
                    _options.camera().setUp(- v.camera.up());
                }

            protected:
                void initializeGL() {
                    makeCurrent();
                    glEnable(GL_MULTISAMPLE);
                    GLint bufs;
                    GLint samples;
                    glGetIntegerv(GL_SAMPLE_BUFFERS, &bufs);
                    glGetIntegerv(GL_SAMPLES, &samples);
                    qDebug("Have %d buffers and %d samples", bufs, samples);
                    qglClearColor(MakeQColor(_options.backgroundColor()));
                    _scene.initialize();
                }

                void resizeGL(int w, int h) {
                    core::PerspectiveCamera & camera = _options.camera();
                    camera.resizeScreen(core::Size(w, h));
                    glViewport(0, 0, w, h);
                }

                void paintEvent(QPaintEvent * e) override {                   
                    QPainter painter;
                    painter.begin(this);
                    painter.beginNativePainting();
                    qglClearColor(MakeQColor(_options.backgroundColor()));
                    glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);

                    core::PerspectiveCamera & camera = _options.camera();
                    camera.resizeScreen(core::Size(width(), height()));

                    _scene.render(_options);

                    painter.endNativePainting();

                    if (_activePenId >= 0){
                        auto pen = _pens[_activePenId];
                        painter.setPen(pen);
                    }
                    painter.setRenderHint(QPainter::Antialiasing);
                    painter.drawPolyline(_stroke);
                    swapBuffers();
                }

                void mousePressEvent(QMouseEvent * e) {
                    if (e->buttons() & Qt::MidButton){
                        _lastPos = e->pos();
                        setCursor(Qt::OpenHandCursor);
                    }
                    else {
                        BaseClass::mousePressEvent(e);
                    }
                }

                void mouseMoveEvent(QMouseEvent * e) {
                    QVector3D t(e->pos() - _lastPos);
                    t.setX(-t.x());
                    if (e->buttons() & Qt::MidButton){
                        _options.camera().moveCenterWithEyeFixed(MakeCoreVec(t));
                        setCursor(Qt::ClosedHandCursor);
                        update();
                    }
                    else{
                        BaseClass::mouseMoveEvent(e);
                    }
                    _lastPos = e->pos();
                }

                void wheelEvent(QWheelEvent * e) {
                    _options.camera().setFocal(_options.camera().focal() * exp(e->delta() / 1000.0));
                    BaseClass::wheelEvent(e);
                    update();
                }

                void processStroke(const std::vector<core::Point2> & stroke, int penId) {
                    if (_callback(stroke, penId)){
                        update();
                    }
                }

            private:
                QPoint _lastPos;
                PanoramicView _view;
                Scene _scene;
                RenderOptions _options;
                std::function<bool(const std::vector<core::Point2> & polyline, int penId)> _callback;
            };

            Widget w(penConfigs, view, callback);
            w.show();

            Singleton::ContinueGui();
        }


    }
}
