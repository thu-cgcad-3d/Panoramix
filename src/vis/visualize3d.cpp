
#include <QtOpenGL>
#include <QtWidgets>
#include <QWidget>

#include <glut.h>
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
        struct Visualizer3D::PrivateData {
            ~PrivateData() {
                for (QWidget * w : widgets)
                    w->deleteLater();
            }
            QList<QWidget*> widgets;
        };


        // visualizer parameters
        Visualizer3D::Params::Params()
            :
            winName("Visualizer 3D"),
            backgroundColor(255, 255, 255),
            camera(700, 700, 200, core::Vec3(1, 1, 1) / 4, core::Vec3(0, 0, 0), core::Vec3(0, 0, -1)),
            renderMode(RenderModeFlag::All)
        {}        

        Visualizer3D::Visualizer3D() : _root(std::make_shared<RenderableObject>()), _data(std::make_shared<PrivateData>()) {
            _activeObject = _root.get();
        }

        Visualizer3D::Visualizer3D(const Params & p, const DefaultRenderState & s)
            : params(p), defaultRenderState(s), _root(std::make_shared<RenderableObject>()), _data(std::make_shared<PrivateData>()) {
            _activeObject = _root.get();
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

            Manipulator<const vis::ColorTable &> SetDefaultColorTable(const vis::ColorTable & colorTable) {
                return Manipulator<const ColorTable &>(
                    [](Visualizer3D & viz, const ColorTable & d) {
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

                Visualizer3DWidget(Visualizer3D & viz, QWidget * parent = 0)
                    : QGLWidget(parent), _params(viz.params), _renderableObjTree(viz.root()) {
                    setMouseTracking(true);
                    setAutoBufferSwap(false);
                    _boundingBox = BoundingBox(_renderableObjTree);
                }
                inline const Params & params() const { return _params; }
            protected:
                void initializeGL() {
                    makeCurrent();
                    qglClearColor(MakeQColor(params().backgroundColor));
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

            private:
                void moveCameraEyeWithCenterFixed(const QVector3D & t) {
                    core::PerspectiveCamera & camera = _params.camera;
                    QVector3D eye = MakeQVec(camera.eye());
                    QVector3D center = MakeQVec(camera.center());
                    QVector3D up = MakeQVec(camera.up());
                    QVector3D tt = t * (eye - center).length() * 0.002f;

                    QVector3D xv = QVector3D::crossProduct(center - eye, up).normalized();
                    QVector3D yv = QVector3D::crossProduct(xv, center - eye).normalized();
                    QVector3D xyTrans = xv * tt.x() + yv * tt.y();
                    double r = ((eye - center).length() - tt.z()) /
                        (eye + xyTrans - center).length();
                    eye = (eye + xyTrans - center) * r + center;
                    up = yv.normalized();
                    _params.camera.setEye(MakeCoreVec(eye), false);
                    _params.camera.setUp(MakeCoreVec(up), false);

                    auto meshCenter = MakeQVec(_boundingBox.center());
                    auto meshRadius = Line3(_boundingBox.minCorner, _boundingBox.maxCorner).length() / 2.0f;
                    auto nearPlane = (eye - meshCenter).length() - meshRadius;
                    nearPlane = nearPlane < 1e-3 ? 1e-3 : nearPlane;
                    auto farPlane = (eye - meshCenter).length() + meshRadius;
                    _params.camera.setNearAndFarPlanes(nearPlane, farPlane, true);
                }

                void moveCameraCenterAndCenter(const QVector3D & t) {
                    core::PerspectiveCamera & camera = _params.camera;
                    QVector3D eye = MakeQVec(camera.eye());
                    QVector3D center = MakeQVec(camera.center());
                    QVector3D up = MakeQVec(camera.up());
                    QVector3D tt = t * (eye - center).length() * 0.002;

                    QVector3D xv = QVector3D::crossProduct((center - eye), up).normalized();
                    QVector3D yv = QVector3D::crossProduct(xv, (center - eye)).normalized();
                    QVector3D zv = (center - eye).normalized();
                    QVector3D trans = xv * tt.x() + yv * tt.y() + zv * tt.z();
                    eye += trans;
                    center += trans;
                    _params.camera.setEye(MakeCoreVec(eye), false);
                    _params.camera.setCenter(MakeCoreVec(center), false);

                    auto meshCenter = MakeQVec(_boundingBox.center());
                    auto meshRadius = Line3(_boundingBox.minCorner, _boundingBox.maxCorner).length() / 2.0f;
                    auto nearPlane = (eye - meshCenter).length() - meshRadius;
                    nearPlane = nearPlane < 1e-3 ? 1e-3 : nearPlane;
                    auto farPlane = (eye - meshCenter).length() + meshRadius;
                    _params.camera.setNearAndFarPlanes(nearPlane, farPlane, true);
                }

            public:
                void autoSetCamera() {
                    auto & box = _boundingBox;
                    auto center = box.center();
                    auto radius = Line3(box.minCorner, box.maxCorner).length() / 2.0;
                    _params.camera.setCenter(center, false);
                    auto eyedirection = _params.camera.eye() - _params.camera.center();
                    eyedirection = eyedirection / core::norm(eyedirection) * radius * 0.8;
                    _params.camera.setEye(center + eyedirection, false);
                    _params.camera.setNearAndFarPlanes(radius / 2.0, radius * 2.0, true);
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
                    if (e->buttons() & Qt::RightButton) {
                        moveCameraEyeWithCenterFixed(t);
                        setCursor(Qt::ClosedHandCursor);
                        update();
                    } else if (e->buttons() & Qt::MidButton) {
                        moveCameraCenterAndCenter(t);
                        update();
                    }
                    _lastPos = e->pos();
                }

                virtual void wheelEvent(QWheelEvent * e) override {
                    moveCameraCenterAndCenter(QVector3D(0, 0, e->delta() / 10));
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


            Manipulator<std::pair<bool, bool>> Show(bool doModel, bool autoSetCamera) {
                return Manipulator<std::pair<bool, bool>>(
                    [](Visualizer3D & viz, std::pair<bool, bool> doModelAndAutoSetCamera) {
                    bool doModal = doModelAndAutoSetCamera.first;
                    bool autoSetCamera = doModelAndAutoSetCamera.second;
                    auto app = Singleton::InitGui();
                    Visualizer3DWidget * w = new Visualizer3DWidget(viz);
                    viz.data()->widgets.append(w);
                    w->resize(MakeQSize(viz.params.camera.screenSize()));
                    w->setWindowTitle(QString::fromStdString(viz.params.winName));
                    if (autoSetCamera) {
                        for (QWidget * w : viz.data()->widgets) {
                            auto v3dw = (Visualizer3DWidget*)w;
                            v3dw->autoSetCamera();
                        }
                    }
                    w->show();
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