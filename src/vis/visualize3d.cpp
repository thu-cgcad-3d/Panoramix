#include "visualize3D.hpp"

#include <QtOpenGL>
#include <QtWidgets>
#include <QWidget>

#include "opengl_object.hpp"

namespace panoramix {
    namespace vis {

        using namespace core;


        // visualizer parameters
        Visualizer3D::Params::Params()
            :
            winName("Image Feature Visualizer"),
            backgroundColor(200, 200, 2200),
            camera(600, 600, 300.0, core::Vec3(1, 1, 0), core::Vec3(0, 0, 0)),
            defaultColor(10, 10, 10),
            pointSize(10.0f),
            lineWidth(2.0f),
            colorTableDescriptor(ColorTableDescriptor::AllColors),
            renderMode(RenderModeFlag::All),
            modelMatrix(core::Mat4::eye())
        {}

        struct Visualizer3D::Entities {
            OpenGLMeshData mesh;
        };

        struct Visualizer3D::Widgets {
            QList<QWidget *> ws;
        };
        

        Visualizer3D::Visualizer3D(const Params & p) 
            : params(p), 
            _ents(std::make_shared<Entities>()),
            _widgets(std::make_shared<Widgets>()) {}

        Visualizer3D::~Visualizer3D() {
            qDeleteAll(_widgets->ws);
        }

        // manipulators
        namespace manip3d {

            Manipulator<const std::string &> SetWindowName(const std::string & name) {
                return Manipulator<const std::string &>(
                    [](Visualizer3D & viz, const std::string & name){
                    viz.params.winName = name; },
                        name);
            }

            Manipulator<const core::Color &> SetDefaultColor(const core::Color & color) {
                return Manipulator<const Color &>(
                    [](Visualizer3D & viz, const Color & c){
                    viz.params.defaultColor = c; },
                        color);
            }

            Manipulator<const core::Color &> SetBackgroundColor(const core::Color & color) {
                return Manipulator<const Color &>(
                    [](Visualizer3D & viz, const Color & c){
                    viz.params.backgroundColor = c; },
                        color);
            }

            Manipulator<const core::PerspectiveCamera &> SetCamera(const core::PerspectiveCamera & camera) {
                return Manipulator<const PerspectiveCamera &>(
                    [](Visualizer3D & viz, const PerspectiveCamera & c){
                    viz.params.camera = c; },
                        camera);
            }

            Manipulator<float> SetPointSize(float pointSize) {
                return Manipulator<float>(
                    [](Visualizer3D & viz, float t){
                    viz.params.pointSize = t; },
                        pointSize);
            }

            Manipulator<float> SetLineWidth(float lineWidth) {
                return Manipulator<float>(
                    [](Visualizer3D & viz, float t){
                    viz.params.lineWidth = t; },
                        lineWidth);
            }

            Manipulator<core::ColorTableDescriptor> SetColorTableDescriptor(core::ColorTableDescriptor descriptor) {
                return Manipulator<ColorTableDescriptor>(
                    [](Visualizer3D & viz, ColorTableDescriptor d){
                    viz.params.colorTableDescriptor = d; },
                        descriptor);
            }

            Manipulator<RenderModeFlags> SetRenderMode(RenderModeFlags mode) {
                return Manipulator<RenderModeFlags>(
                    [](Visualizer3D & viz, RenderModeFlags d){
                    viz.params.renderMode = d; },
                        mode);
            }

            Manipulator<const core::Mat4 &> SetModelMatrix(const core::Mat4 & mat) {
                return Manipulator<const core::Mat4 &>(
                    [](Visualizer3D & viz, const core::Mat4 & m){
                    viz.params.modelMatrix = m; },
                        mat);
            }

            void AutoSetCamera(Visualizer3D & viz) {
                auto box = viz.entities()->mesh.boundingBox();
                auto center = (box.first + box.second) / 2.0;
                auto radius = (box.second - box.first).length() / 2.0;
                viz.params.camera.setCenter(MakeCoreVec(center), false);
                auto eyedirection = viz.params.camera.eye() - viz.params.camera.center();
                eyedirection = eyedirection / core::norm(eyedirection) * radius * 1.5;
                viz.params.camera.setEye(MakeCoreVec(center) + eyedirection, false);
                viz.params.camera.setNearAndFarPlanes(radius / 2.0, radius * 4.0, true);
            }

            namespace {

                // visualizer widget
                class Visualizer3DWidget : public QGLWidget, protected QGLFunctions {
                public:
                    Visualizer3DWidget(Visualizer3D & viz, QWidget * parent = 0) 
                        : QGLWidget(parent), _viz(viz){
                        setAutoFillBackground(false);
                        setMouseTracking(true);

                        _meshBox = viz.entities()->mesh.boundingBox();
                    }
                protected:
                    void initializeGL() {
                        qglClearColor(MakeQColor(_viz.params.backgroundColor));
                        _trianglesObject = new OpenGLObject(this);
                        _trianglesObject->setUpShaders(OpenGLShaderSourceName::NormalTriangles);
                        _trianglesObject->setUpMesh(_viz.entities()->mesh);
                        _linesObject = new OpenGLObject(this);
                        _linesObject->setUpShaders(OpenGLShaderSourceName::NormalLines);
                        _linesObject->setUpMesh(_viz.entities()->mesh);
                        _pointsObject = new OpenGLObject(this);
                        _pointsObject->setUpShaders(OpenGLShaderSourceName::NormalPoints);
                        _pointsObject->setUpMesh(_viz.entities()->mesh);
                    }

                    void paintGL() {
                        glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);

                        glFrontFace(GL_CW); // face direction set to clockwise
                        glCullFace(GL_FRONT); // specify whether front- or back-facing facets can be culled
                        glEnable(GL_CULL_FACE);
                        glEnable(GL_DEPTH_TEST);
                        glEnable(GL_STENCIL_TEST);

                        glEnable(GL_ALPHA_TEST);

                        glEnable(GL_BLEND);
                        glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);
                        //glLineWidth(_viz.params.lineWidth);
                        core::PerspectiveCamera & camera = _viz.params.camera;
                        QMatrix4x4 modelMatrix = MakeQMatrix(_viz.params.modelMatrix);

                        if (_viz.params.renderMode & RenderModeFlag::Triangles){
                            _trianglesObject->render(RenderModeFlag::Triangles,
                                MakeQMatrix(camera.projectionMatrix()),
                                MakeQMatrix(camera.viewMatrix()),
                                modelMatrix);
                        }
                        if (_viz.params.renderMode & RenderModeFlag::Points){
                            _pointsObject->render(RenderModeFlag::Points,
                                MakeQMatrix(camera.projectionMatrix()),
                                MakeQMatrix(camera.viewMatrix()),
                                modelMatrix);
                        }
                        if (_viz.params.renderMode & RenderModeFlag::Lines){
                            _linesObject->render(RenderModeFlag::Lines,
                                MakeQMatrix(camera.projectionMatrix()),
                                MakeQMatrix(camera.viewMatrix()),
                                modelMatrix);
                        }
                        
                        glDisable(GL_DEPTH_TEST);
                        glDisable(GL_CULL_FACE);
                    }

                    void resizeGL(int w, int h) {
                        core::PerspectiveCamera & camera = _viz.params.camera;
                        camera.resizeScreen(core::Size(w, h));
                        glViewport(0, 0, w, h);
                    }

                private:
                    void moveCameraEyeWithCenterFixed(const QVector3D & t) {
                        core::PerspectiveCamera & camera = _viz.params.camera;
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
                        _viz.params.camera.setEye(MakeCoreVec(eye), false);
                        _viz.params.camera.setUp(MakeCoreVec(up), false);
                        
                        auto meshCenter = (_meshBox.first + _meshBox.second) / 2.0f;
                        auto meshRadius = (_meshBox.second - _meshBox.first).length() / 2.0f;
                        auto nearPlane = (eye - meshCenter).length() - meshRadius;
                        nearPlane = nearPlane < 1e-3 ? 1e-3 : nearPlane;
                        auto farPlane = (eye - meshCenter).length() + meshRadius;
                        _viz.params.camera.setNearAndFarPlanes(nearPlane, farPlane, true);
                    }

                    void moveCameraCenterAndCenter(const QVector3D & t) {
                        core::PerspectiveCamera & camera = _viz.params.camera;
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
                        _viz.params.camera.setEye(MakeCoreVec(eye), false);
                        _viz.params.camera.setCenter(MakeCoreVec(center), false);
                        
                        auto meshCenter = (_meshBox.first + _meshBox.second) / 2.0;
                        auto meshRadius = (_meshBox.second - _meshBox.first).length() / 2.0;
                        auto nearPlane = (eye - meshCenter).length() - meshRadius;
                        nearPlane = nearPlane < 1e-3 ? 1e-3 : nearPlane;
                        auto farPlane = (eye - meshCenter).length() + meshRadius;
                        _viz.params.camera.setNearAndFarPlanes(nearPlane, farPlane, true);
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
                        if (e->buttons() & Qt::RightButton){
                            moveCameraEyeWithCenterFixed(t);
                            setCursor(Qt::ClosedHandCursor);
                            update();
                        }
                        else if (e->buttons() & Qt::MidButton){
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
                    Visualizer3D & _viz;
                    QPointF _lastPos;
                    OpenGLObject * _linesObject;
                    OpenGLObject * _pointsObject;
                    OpenGLObject * _trianglesObject;
                    QPair<QVector3D, QVector3D> _meshBox;
                };
                

            }


            Manipulator<bool> Show(bool doModel) {
                return Manipulator<bool>(
                    [](Visualizer3D & viz, bool modal){                  

                    if (modal){
                        QDialog dialog;
                        dialog.setWindowTitle(QString::fromStdString(viz.params.winName));
                        QVBoxLayout * layout = new QVBoxLayout;
                        Visualizer3DWidget * w = new Visualizer3DWidget(viz, &dialog);
                        layout->addWidget(w);
                        layout->setContentsMargins(0, 0, 0, 0);
                        dialog.setLayout(layout);
                        dialog.resize(MakeQSize(viz.params.camera.screenSize()));
                        dialog.exec();
                    }
                    else {
                        QDialog * d = new QDialog;
                        QDialog & dialog = *d;
                        dialog.setWindowTitle(QString::fromStdString(viz.params.winName));
                        QVBoxLayout * layout = new QVBoxLayout;
                        Visualizer3DWidget * w = new Visualizer3DWidget(viz, &dialog);
                        layout->addWidget(w);
                        layout->setContentsMargins(0, 0, 0, 0);
                        dialog.setLayout(layout);
                        dialog.resize(MakeQSize(viz.params.camera.screenSize()));
                        viz.widgets()->ws.append(d);
                        d->show();
                    }

                },
                    doModel);
            }

        }


        Visualizer3D operator << (Visualizer3D viz, const core::Point3 & p) {
            OpenGLMeshData::Vertex v;
            v.position4 = MakeQVec(core::HPoint3(p, 1.0).toVector());
            v.color4 = MakeQVec(viz.params.defaultColor) / 255.0f;
            v.lineWidth1 = viz.params.lineWidth;
            v.pointSize1 = viz.params.pointSize;
            viz.entities()->mesh.addVertex(v);
            return viz;
        }


        Visualizer3D operator << (Visualizer3D viz, const core::Line3 & p) {
            auto & mesh = viz.entities()->mesh;
            core::Point3 ps[] = { p.first, p.second };
            OpenGLMeshData::Vertex vs[2];
            for (int i = 0; i < 2; i++){
                vs[i].position4 = MakeQVec(core::HPoint3(ps[i], 1.0).toVector());
                vs[i].color4 = MakeQVec(viz.params.defaultColor) / 255.0f;
                vs[i].lineWidth1 = viz.params.lineWidth;
                vs[i].pointSize1 = viz.params.pointSize;
            }
            viz.entities()->mesh.addIsolatedLine(vs[0], vs[1]);
            return viz;
        }

    }
}