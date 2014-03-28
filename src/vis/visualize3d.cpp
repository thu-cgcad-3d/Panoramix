#include "visualize3D.hpp"

#include <QtOpenGL>
#include <QtWidgets>
#include <QWidget>
//#include <gl/GL.h>
//#include <gl/GLU.h>

#include "opengl_object.hpp"
//#include "opengl_widgets.hpp"

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
            colorTableDescriptor(ColorTableDescriptor::AllColors)
        {}

        struct Visualizer3D::Entities {
            OpenGLMeshData mesh;
        };
        

        Visualizer3D::Visualizer3D(const Params & p) 
            : params(p), 
            _ents(std::make_shared<Entities>()) {}

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

            Manipulator<RenderModeFlags> SetRanderMode(RenderModeFlags mode) {
                return Manipulator<RenderModeFlags>(
                    [](Visualizer3D & viz, RenderModeFlags d){
                    viz.params.renderMode = d; },
                        mode);
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

                static const char *vsrc =
                    //"attribute lowp float typeFlag;\n" // points (=0) < 0.3; triangles (=1) > 0.7; lines (=0.5) otherwise
                    "attribute highp vec4 position;\n"
                    "attribute lowp float pointSize;\n"
                    "attribute highp vec3 normal;\n"
                    "attribute lowp vec4 color;\n"
                    "attribute lowp vec2 texCoord;\n"
                    "uniform highp mat4 viewMatrix;\n"
                    "uniform highp mat4 modelMatrix;\n"
                    "uniform highp mat4 projectionMatrix;\n"
                    "varying vec4 pixelColor;\n"
                    "void main(void)\n"
                    "{\n"
                    "    gl_Position = projectionMatrix * viewMatrix * modelMatrix * position;\n"
                    //"    if (typeFlag < 0.3) \n"
                    "       gl_PointSize = pointSize;\n"
                    "    pixelColor = color;\n"
                    "}\n";
                static const char *fsrc =
                    "varying lowp vec4 pixelColor;\n"
                    "void main(void)\n"
                    "{\n"
                    "    gl_FragColor = pixelColor;\n"
                    "}\n";


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
                        _object = new OpenGLObject(this);
                        _object->setUpShaders({ vsrc, fsrc });
                        _object->setUpMesh(_viz.entities()->mesh);
                    }

                    void paintGL() {
                        glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);

                        glFrontFace(GL_CW);
                        glCullFace(GL_FRONT);
                        glEnable(GL_CULL_FACE);
                        glEnable(GL_DEPTH_TEST);
                        glEnable(GL_STENCIL_TEST);
                        glEnable(GL_LINE_WIDTH);

                        glEnable(GL_BLEND);
                        glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);
                        glLineWidth(_viz.params.lineWidth);
                        core::PerspectiveCamera & camera = _viz.params.camera;
                        QMatrix4x4 modelMatrix;
                        _object->render(_viz.params.renderMode, 
                            MakeQMatrix(camera.projectionMatrix()), MakeQMatrix(camera.viewMatrix()), modelMatrix);
                        
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
                        //auto sc = core::norm(camera.eye() - camera.center());
                        QVector3D eye = MakeQVec(camera.eye());
                        QVector3D center = MakeQVec(camera.center());
                        QVector3D up = MakeQVec(camera.up());
                        QVector3D tt = t * (eye - center).length() * 0.002;

                        QVector3D xv = QVector3D::crossProduct(center - eye, up).normalized();
                        QVector3D yv = QVector3D::crossProduct(xv, center - eye).normalized();
                        QVector3D xyTrans = xv * tt.x() + yv * tt.y();
                        double r = ((eye - center).length() - tt.z()) /
                            (eye + xyTrans - center).length();
                        eye = (eye + xyTrans - center) * r + center;
                        up = yv.normalized();
                        _viz.params.camera.setEye(MakeCoreVec(eye), false);
                        _viz.params.camera.setUp(MakeCoreVec(up), false);
                        
                        auto meshCenter = (_meshBox.first + _meshBox.second) / 2.0;
                        auto meshRadius = (_meshBox.second - _meshBox.first).length() / 2.0;
                        auto nearPlane = (eye - meshCenter).length() - meshRadius;
                        nearPlane = nearPlane < 1e-3 ? 1e-3 : nearPlane;
                        auto farPlane = (eye - meshCenter).length() + meshRadius;
                        _viz.params.camera.setNearAndFarPlanes(nearPlane, farPlane, true);
                    }

                    void moveCameraCenterAndCenter(const QVector3D & t) {
                        core::PerspectiveCamera & camera = _viz.params.camera;
                        //float sc = std::max(camera.screenSize().width, camera.screenSize().height) * 0.001;
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
                    OpenGLObject * _object;
                    QPair<QVector3D, QVector3D> _meshBox;
                };
                

            }


            Manipulator<bool> Show(bool doModel) {
                return Manipulator<bool>(
                    [](Visualizer3D & viz, bool t){                  

                    QDialog dialog;
                    dialog.setWindowTitle(QString::fromStdString(viz.params.winName));
                    QVBoxLayout * layout = new QVBoxLayout;
                    Visualizer3DWidget * w = new Visualizer3DWidget(viz, &dialog);
                    layout->addWidget(w);
                    layout->setContentsMargins(0, 0, 0, 0);
                    dialog.setLayout(layout);
                    dialog.resize(MakeQSize(viz.params.camera.screenSize()));
                    dialog.exec();

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
            viz << p.first;
            mesh.vertices.back().lineWidth1 = viz.params.lineWidth;
            mesh.vertices.back().pointSize1 = viz.params.lineWidth;
            viz << p.second;
            mesh.vertices.back().lineWidth1 = viz.params.lineWidth;
            mesh.vertices.back().pointSize1 = viz.params.lineWidth;
            mesh.addLine(mesh.vertices.size() - 2, mesh.vertices.size() - 1);
            return viz;
        }

    }
}