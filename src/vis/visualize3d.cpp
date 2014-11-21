
#include <QtOpenGL>
#include <QtWidgets>
#include <QWidget>

#include <GL/GLU.h>

#include "../core/macros.hpp"
#include "../core/utilities.hpp"

#include "qt_glue.hpp"
#include "renderable_object_tree.hpp"
#include "singleton.hpp"
#include "visualize3D.hpp"

namespace panoramix {
    namespace vis {

        using namespace core;

        // visulizer private data
        struct GuiData {
            QMap<Visualizer3D *, QList<QWidget *>> widgetsTable;
        };

        static GuiData staticGuiData;

        // visualizer parameters
        Visualizer3D::Params::Params()
            :
            winName("Visualizer 3D"),
            backgroundColor(255, 255, 255),
            camera(700, 700, 700, core::Vec3(1, 1, 1) / 4, core::Vec3(0, 0, 0), core::Vec3(0, 0, -1)),
            renderMode(RenderModeFlag::All){
        }        

        Visualizer3D::Visualizer3D() : _root(std::make_shared<RenderableObject>())/*, _data(std::make_shared<PrivateData>())*/ {
            _activeObject = _root.get();
        }

        Visualizer3D::Visualizer3D(const Params & p, const DefaultRenderState & s)
            : params(p), defaultRenderState(s), _root(std::make_shared<RenderableObject>())/*, _data(std::make_shared<PrivateData>())*/ {
            _activeObject = _root.get();
        }

        Visualizer3D::~Visualizer3D() {
            staticGuiData.widgetsTable.remove(this);
        }

        void Visualizer3D::deactivateLast() {
            if (_activeObject != _root.get()) {
                _activeObject = _activeObject->parent();
            } else {
                qWarning() << "can not deactivate root!";
            }
        }

        // manipulators
        namespace manip3d {

            Manipulator<const std::string &> SetWindowName(const std::string & name) {
                return Manipulator<const std::string &>([](Visualizer3D & viz, const std::string & name) {
                    viz.params.winName = name;
                }, name);
            }

            Manipulator<Color> SetDefaultForegroundColor(Color color) {
                return Manipulator<Color>([](Visualizer3D & viz, Color c) {
                    viz.defaultRenderState.foregroundColor = c;
                }, color);
            }

            Manipulator<Color> SetBackgroundColor(Color color) {
                return Manipulator<Color>([](Visualizer3D & viz, Color c) {
                    viz.params.backgroundColor = c;
                }, color);
            }

            Manipulator<const PerspectiveCamera &> SetCamera(const PerspectiveCamera & camera) {
                return Manipulator<const PerspectiveCamera &>(
                    [](Visualizer3D & viz, const PerspectiveCamera & c) {
                    viz.params.camera = c;
                }, camera);
            }

            Manipulator<float> SetDefaultPointSize(float pointSize) {
                return Manipulator<float>(
                    [](Visualizer3D & viz, float t) {
                    viz.defaultRenderState.pointSize = t; },
                        pointSize);
            }

            Manipulator<float> SetDefaultLineWidth(float lineWidth) {
                return Manipulator<float>(
                    [](Visualizer3D & viz, float t) {
                    viz.defaultRenderState.lineWidth = t; },
                        lineWidth);
            }

            Manipulator<vis::ColorTable> SetDefaultColorTable(const vis::ColorTable & colorTable) {
                return Manipulator<ColorTable>(
                    [](Visualizer3D & viz, ColorTable d) {
                    viz.defaultRenderState.colorTable = d; },
                        colorTable);
            }

            Manipulator<RenderModeFlags> SetRenderMode(RenderModeFlags mode) {
                return Manipulator<RenderModeFlags>(
                    [](Visualizer3D & viz, RenderModeFlags d) {
                    viz.params.renderMode = d; },
                        mode);
            }


            // visualizer widget
            class Visualizer3DWidget : public QGLWidget {
            public:
                using Params = Visualizer3D::Params;

                Visualizer3DWidget(Visualizer3D & viz, QWidget * parent = nullptr)
                    : QGLWidget(parent), _params(viz.params), _renderableObjTree(viz.root()) {
                    setMouseTracking(true);
                    setAutoBufferSwap(false);
                    _boundingBox = BoundingBox(_renderableObjTree);
                }
                inline const Params & params() const { return _params; }
            protected:
                void initializeGL() {
                    makeCurrent();
                    
                    glEnable(GL_MULTISAMPLE);
                    GLint bufs;
                    GLint samples;
                    glGetIntegerv(GL_SAMPLE_BUFFERS, &bufs);
                    glGetIntegerv(GL_SAMPLES, &samples);
                    qDebug("Have %d buffers and %d samples", bufs, samples);

                    qglClearColor(MakeQColor(params().backgroundColor));
                    _renderableObjTree.initialize();
                }

                void paintGL() {
                    QPainter painter;
                    painter.begin(this);

                    painter.beginNativePainting();
                    qglClearColor(MakeQColor(params().backgroundColor));

                    glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);

                    glFrontFace(GL_CW); // face direction set to clockwise
                    //glCullFace(GL_FRONT); // specify whether front- or back-facing facets can be culled
                    //glEnable(GL_CULL_FACE);
                    glEnable(GL_MULTISAMPLE);
                    glEnable(GL_DEPTH_TEST);
                    glEnable(GL_STENCIL_TEST);

                    glEnable(GL_ALPHA_TEST);

                    glEnable(GL_BLEND);
                    glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);

                    core::PerspectiveCamera & camera = _params.camera;
                    camera.resizeScreen(core::SizeI(width(), height()));

                    _renderableObjTree.renderWithCamera(_params.renderMode, _params.camera);


                    glDisable(GL_DEPTH_TEST);
                    glDisable(GL_CULL_FACE);

                    painter.endNativePainting();
                    swapBuffers();
                }

                void resizeGL(int w, int h) {
                    core::PerspectiveCamera & camera = _params.camera;
                    camera.resizeScreen(core::Size(w, h));
                    glViewport(0, 0, w, h);
                }

            public:
                void autoSetCamera() {
                    // get bounding box for camera
                   /* auto minV = std::numeric_limits<double>::lowest();
                    auto maxV = std::numeric_limits<double>::max();
                    Point3 minC(maxV, maxV, maxV), maxC(minV, minV, minV);*/
                    auto sphere = _boundingBox.outerSphere();
                   /* sphere.center = Vec3(0, 0, 0);*/
                    _params.camera.focusOn(sphere, true);
                    update();
                }

            protected:
                virtual void mousePressEvent(QMouseEvent * e) override {
                    _lastPos = e->pos();
                    if (e->buttons() & Qt::RightButton)
                        setCursor(Qt::OpenHandCursor);
                    else if (e->buttons() & Qt::MidButton)
                        setCursor(Qt::SizeAllCursor);
                }

                virtual void mouseMoveEvent(QMouseEvent * e) override {
                    QVector3D t(e->pos() - _lastPos);
                    t.setX(-t.x());
                    auto sphere = _boundingBox.outerSphere();
                    //sphere.center = Vec3(0, 0, 0);
                    if (e->buttons() & Qt::RightButton) {
                        core::Vec3 trans = t.x() * _params.camera.rightward() + t.y() * _params.camera.upward();
                        trans *= 0.02;
                        _params.camera.moveEyeWithCenterFixed(trans, sphere, true, true);
                        setCursor(Qt::ClosedHandCursor);
                        update();
                    } else if (e->buttons() & Qt::MidButton) {
                        core::Vec3 trans = t.x() * _params.camera.rightward() + t.y() * _params.camera.upward();
                        trans *= 0.02;
                        _params.camera.translate(trans, sphere, true);
                        update();
                    }
                    _lastPos = e->pos();
                }

                virtual void wheelEvent(QWheelEvent * e) override {
                    auto sphere = _boundingBox.outerSphere();
                    /*sphere.center = Vec3(0, 0, 0);*/
                    double d = e->delta() / 10;
                    double dist = core::Distance(_params.camera.eye(), _params.camera.center());
                    core::Vec3 trans = d * dist/1000.0 * _params.camera.forward();
                    _params.camera.moveEyeWithCenterFixed(trans, sphere, false, true);
                    update();
                }

                virtual void mouseReleaseEvent(QMouseEvent * e) override {
                    unsetCursor();
                }

            private:
                Params _params;
                QPointF _lastPos;
                RenderableObjectTree _renderableObjTree;
                core::Box3 _boundingBox;
            };


            class Visualizer3DMainWindow : public QMainWindow {
            public:
                explicit Visualizer3DMainWindow(QWidget * parent = nullptr) : QMainWindow(parent) {
                    setupGui();
                }

                void setupGui() {
                    auto menuView = this->menuBar()->addMenu(tr("View"));
                    auto menuAbout = this->menuBar()->addMenu(tr("About"));
                    this->statusBar()->show();

                    auto actionAbout = menuAbout->addAction(tr("About"));
                    connect(actionAbout, &QAction::triggered, [this](){
                        QMessageBox::about(this, tr("About this program"), 
                            tr("Panoramix.Vis is the visulization module of project Panoramix developped by Yang Hao."));
                    });
                }
            };


            Manipulator<std::pair<bool, bool>> Show(bool doModel, bool autoSetCamera) {
                return Manipulator<std::pair<bool, bool>>(
                    [](Visualizer3D & viz, std::pair<bool, bool> doModelAndAutoSetCamera) {
                    bool doModal = doModelAndAutoSetCamera.first;
                    bool autoSetCamera = doModelAndAutoSetCamera.second;
                    auto app = Singleton::InitGui();
                    Visualizer3DWidget * w = new Visualizer3DWidget(viz);
                    Visualizer3DMainWindow * mwin = new Visualizer3DMainWindow();
                    mwin->setCentralWidget(w);
                    mwin->setAttribute(Qt::WA_DeleteOnClose);
                    staticGuiData.widgetsTable[&viz].append(mwin);
                    mwin->resize(MakeQSize(viz.params.camera.screenSize()));
                    mwin->setWindowTitle(QString::fromStdString(viz.params.winName));
                    mwin->setWindowIcon(Singleton::DefaultConfiguration().icon);
                    mwin->setStyleSheet(Singleton::DefaultConfiguration().css);
                    auto palette = mwin->palette();
                    palette.setColor(QPalette::Window, MakeQColor(viz.params.backgroundColor));
                    mwin->setPalette(palette);
                    //qDebug() << mwin->styleSheet();
                    if (autoSetCamera) {
                        w->autoSetCamera();
                    }
                    mwin->show();
                    if (doModal) {
                        Singleton::ContinueGui(); // qApp->exec()
                    }
                },  std::make_pair(doModel, autoSetCamera));
            }


            Manipulator<const Mat4 &> SetModelMatrix(const Mat4 & mat) {
                return Manipulator<const Mat4 &>(
                    [](Visualizer3D & viz, const Mat4 & m) {
                    viz.activeObject()->modelMatrix() = m;
                }, mat);
            }

            Manipulator<const Image &> SetTexture(const Image & tex) {
                return Manipulator<const Image &>(
                    [](Visualizer3D & viz, const Image & t) {
                    viz.activeObject()->setTexture(t);
                }, tex);
            }

        }

    }
}