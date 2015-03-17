#include <QtOpenGL>

#include "../core/utilities.hpp"
#include "qt_glue.hpp"
#include "visualizers.hpp"
#include "singleton.hpp"

#define GL_ALPHA_TEST 0x0BC0

namespace panoramix {
    namespace vis {

        using namespace core;


        // visualizer widget

        class VisualizerWidget : public QGLWidget {
        public:
            RenderOptions options;
            VisualObjectScene scene;

            VisualizerWidget(const Visualizer & v, QWidget * parent = nullptr) 
                : QGLWidget(parent), options(v.renderOptions), scene(v.tree()) {
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
                qglClearColor(MakeQColor(options.backgroundColor));   
                scene.initialize();
            }

            void paintGL() {
                QPainter painter;
                painter.begin(this);

                painter.beginNativePainting();
                qglClearColor(MakeQColor(options.backgroundColor));

                glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);

                glFrontFace(GL_CCW); // face direction set to clockwise
                glEnable(GL_MULTISAMPLE);
                glEnable(GL_DEPTH_TEST);
                glEnable(GL_STENCIL_TEST);

                glEnable(GL_ALPHA_TEST);
                if (options.showInside){
                    glEnable(GL_CULL_FACE);
                }
                else{
                    glDisable(GL_CULL_FACE);
                }

                glEnable(GL_BLEND);
                glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);

                glEnable(GL_PROGRAM_POINT_SIZE);
                glEnable(GL_POINT_SPRITE);

                core::PerspectiveCamera & camera = options.camera;
                camera.resizeScreen(core::Size(width(), height()));

                scene.render(options);

                glDisable(GL_DEPTH_TEST);
                if (options.showInside){
                    glDisable(GL_CULL_FACE);
                }

                painter.endNativePainting();
                swapBuffers();
            }

            void resizeGL(int w, int h) {
                core::PerspectiveCamera & camera = options.camera;
                camera.resizeScreen(core::Size(w, h));
                glViewport(0, 0, w, h);
            }

        public:
            void autoSetCamera() {
                auto sphere = scene.boundingBox().outerSphere();
                options.camera.resizeScreen(core::Size(width(), height()), false);
                options.camera.focusOn(sphere, true);
                update();
            }

        protected:
            virtual void mousePressEvent(QMouseEvent * e) override {
                _lastPos = e->pos();
                if (e->buttons() & Qt::RightButton)
                    setCursor(Qt::OpenHandCursor);
                else if (e->buttons() & Qt::MidButton)
                    setCursor(Qt::SizeAllCursor);
                else if (e->buttons() & Qt::LeftButton){
                    VisualObjectHandle oh;
                    TriMesh::TriangleHandle t;
                    std::tie(oh, t) = scene.pickOnScreen(options, core::Point2(e->pos().x(), e->pos().y()));                    
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
                            ->invokeCallbackFunction(InteractionID::ClickLeftButton, scene.tree(), std::make_pair(oh, entityID));
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
                    core::Vec3 trans = t.x() * options.camera.rightward() + t.y() * options.camera.upward();
                    trans *= 0.02;
                    options.camera.moveEyeWithCenterFixed(trans, sphere, true, true);
                    setCursor(Qt::ClosedHandCursor);
                    update();
                }
                else if ((e->buttons() & Qt::MidButton) || 
                    ((e->buttons() & Qt::RightButton) && (e->modifiers() & Qt::ShiftModifier))) {
                    core::Vec3 trans = t.x() * options.camera.rightward() + t.y() * options.camera.upward();
                    trans *= 0.02;
                    options.camera.translate(trans, sphere, true);
                    update();
                }
                _lastPos = e->pos();
            }

            virtual void wheelEvent(QWheelEvent * e) override {
                auto sphere = scene.boundingBox().outerSphere();
                double d = e->delta() / 10;
                double dist = core::Distance(options.camera.eye(), options.camera.center());
                core::Vec3 trans = d * dist / 1000.0 * options.camera.forward();
                options.camera.moveEyeWithCenterFixed(trans, sphere, false, true);
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
                                n.data->invokeCallbackFunction(InteractionID::PressSpace,
                                    scene.tree(),
                                    VisualObjectEntityID{ n.topo.hd, entityID });
                            }
                        }
                    }
                }
            }

        private:
            QPointF _lastPos;
        };

        
        template <class T, class = std::enable_if_t<std::is_floating_point<T>::value>>
        QWidget * MakeGuiAgent(core::Noted<T> & value, QWidget * parent = nullptr){
            QDoubleSpinBox * spinBox = new QDoubleSpinBox(parent);
            spinBox->setValue(value.component);
            spinBox->setSingleStep(0.01);
            spinBox->setDecimals(3);
            spinBox->setRange(0.0, 1.0);
            auto signal = static_cast<void (QDoubleSpinBox::*)(double)>(&QDoubleSpinBox::valueChanged);
            QObject::connect(spinBox, signal, [&value](double v){
                //std::cout << "value of " << value.note << " is set to " << v << "!" << std::endl;
                value.component = v;
            });
            return spinBox;
        }

        QWidget * MakeGuiAgent(core::Noted<bool> & value, QWidget * parent = nullptr){
            QCheckBox * checkBox = new QCheckBox(parent);
            checkBox->setCheckable(true);
            checkBox->setChecked(value.component);
            QObject::connect(checkBox, &QCheckBox::clicked, [&value, checkBox](){
                value.component = checkBox->isChecked();
            });
            return checkBox;
        }

        QWidget * MakeGuiAgent(core::Noted<Color> & value, QWidget * parent = nullptr){
            NOT_IMPLEMENTED_YET();
        }


        template <class ... Ts>
        void PopUpDialog(QWidget * parent, core::Noted<Ts> & ... values){
            QString names[] = { QString::fromStdString(values.note) ... };
            QWidget * agents[] = { MakeGuiAgent(values, nullptr) ... };
            QDialog dialog;
            QFormLayout * layout = new QFormLayout;
            for (int i = 0; i < sizeof...(Ts); i++){
                layout->addRow(names[i], agents[i]);
            }
            dialog.setLayout(layout);
            dialog.exec();
        }


        void PopUpGui(RenderOptions & options, QWidget * widget = nullptr){
            core::Noted<float> bwColor = core::NoteAs(options.bwColor, "Blend Weight of Color");
            core::Noted<float> bwTexColor = core::NoteAs(options.bwTexColor, "Blend Weight of Texture Color");
            core::Noted<bool> showInside = core::NoteAs(options.showInside, "Show Inside");
            core::Noted<bool> showPoints = core::NoteAs<bool>(options.renderMode & RenderModeFlag::Points, "Show Points");
            core::Noted<bool> showLines = core::NoteAs<bool>(options.renderMode & RenderModeFlag::Lines, "Show Lines");
            core::Noted<bool> showFaces = core::NoteAs<bool>(options.renderMode & RenderModeFlag::Triangles, "Show Faces");
            PopUpDialog(widget, bwColor, bwTexColor, showInside, showPoints, showLines, showFaces);
            options.bwColor = bwColor.component;
            options.bwTexColor = bwTexColor.component;
            options.showInside = showInside.component;
            options.renderMode = (showPoints.component ? RenderModeFlag::Points : RenderModeFlag::None)
                | (showLines.component ? RenderModeFlag::Lines : RenderModeFlag::None)
                | (showFaces.component ? RenderModeFlag::Triangles : RenderModeFlag::None);
        }



        QWidget * Visualizer::createWidget(bool autoSetCamera, QWidget * parent){
            VisualizerWidget * w = new VisualizerWidget(*this, parent);
            w->setMinimumSize(300, 300);
            w->resize(MakeQSize(renderOptions.camera.screenSize()));
            auto actionSettings = new QAction(QObject::tr("Settings"), nullptr);
            QObject::connect(actionSettings, &QAction::triggered, [w](){
                PopUpGui(w->options, w);
                w->update();
            });
            if (autoSetCamera) {
                w->autoSetCamera();
            }
            w->addAction(actionSettings);
            return w;
        }

        void Visualizer::show(bool doModal, bool autoSetCamera, Visualizer::CameraScalePolicy csp) {
            auto app = Singleton::InitGui();
            VisualizerWidget * w = new VisualizerWidget(*this);
            
            QMainWindow * mwin = new QMainWindow;
            mwin->setCentralWidget(w);
            mwin->setAttribute(Qt::WA_DeleteOnClose);
            mwin->resize(MakeQSize(renderOptions.camera.screenSize()));
            mwin->setWindowTitle(QString::fromStdString(renderOptions.winName));
            mwin->setWindowIcon(Singleton::DefaultIcon());
            mwin->setStyleSheet(Singleton::DefaultCSS());

            auto menuView = mwin->menuBar()->addMenu(QObject::tr("View"));
            auto actionSettings = menuView->addAction(QObject::tr("Settings"));
            QObject::connect(actionSettings, &QAction::triggered, [w](){
                PopUpGui(w->options, w);
                w->update();
            });
            auto menuAbout = mwin->menuBar()->addMenu(QObject::tr("About"));
            auto actionAbout = menuAbout->addAction(QObject::tr("About"));
            QObject::connect(actionAbout, &QAction::triggered, [mwin](){
                QMessageBox::about(mwin, QObject::tr("About this program"),
                    QObject::tr("Panoramix.Vis is the visulization module of project Panoramix developped by Yang Hao."));
            });
            mwin->statusBar()->show();


            auto palette = mwin->palette();
            palette.setColor(QPalette::Window, MakeQColor(renderOptions.backgroundColor));
            mwin->setPalette(palette);
            //qDebug() << mwin->styleSheet();
            if (autoSetCamera) {
                w->autoSetCamera();
            }
            mwin->show();
            if (doModal) {
                Singleton::ContinueGui(); // qApp->exec()
            }
        }

    }
}