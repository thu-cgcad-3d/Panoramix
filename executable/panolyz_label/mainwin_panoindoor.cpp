#include "../../src/gui/qttools.hpp"
#include "mainwin_panoindoor.hpp"


namespace panolyz {
    namespace PanoramaIndoor {

        using namespace pano::gui;

        Widget::Widget(QWidget * parent /*= nullptr*/) : QGLWidget(parent), _anno(nullptr), _state(Idle) {
            setMouseTracking(true);
            setAutoBufferSwap(false);
            grabKeyboard();

            _chain.closed = false;
            setContextMenuPolicy(Qt::ActionsContextMenu);
            QActionGroup * bas = new QActionGroup(this);

            QAction * defaultAction = nullptr;
            connect(defaultAction = bas->addAction(tr("Cancel")), &QAction::triggered, [this]() { 
                _state = Idle; 
                _chain.points.clear();
                rebuildStrokeScene();
                update();
            });            
            connect(defaultAction = bas->addAction(tr("Create Polygon")), &QAction::triggered, [this]() {
                _state = CreatingPolygon;
                _chain.points.clear();
                _chain.closed = false;
                rebuildStrokeScene();
                update();
            });
            connect(defaultAction = bas->addAction(tr("Create Line")), &QAction::triggered, [this]() {
                _state = CreatingLine;
                _chain.points.clear();
                _chain.closed = false;
                rebuildStrokeScene();
                update();
            });
        }

        Widget::~Widget() {

        }


        void Widget::setCurAnnotation(PanoIndoorAnnotation * anno) {
            _state = Idle;

            if (anno == nullptr) {
                _anno = anno;
                return;
            }

            assert(!anno.view.image.empty());
            auto im = anno->view.image;

            SceneBuilder sb;
            ResourceStore::set("tex", im);
            Sphere3 sp;
            sp.center = Origin();
            sp.radius = 10.0;
            sb.begin(sp).shaderSource(OpenGLShaderSourceDescriptor::XPanorama).resource("tex").end();
            _anno = anno;

            _options.panoramaAspectRatio(im.rows / float(im.cols));
            _options.panoramaHoriCenterRatio(0.5f);
            _options.camera(PerspectiveCamera(500, 500, Point2(250, 250), 200, Origin(), X(), -Z()));

            _imageScene = sb.scene();

            update();
        }


        void Widget::initializeGL() {
            if (!_anno) {
                return;
            }
            makeCurrent();
            glEnable(GL_MULTISAMPLE);
            GLint bufs;
            GLint samples;
            glGetIntegerv(GL_SAMPLE_BUFFERS, &bufs);
            glGetIntegerv(GL_SAMPLES, &samples);
            qDebug("Have %d buffers and %d samples", bufs, samples);
            qglClearColor(MakeQColor(_options.backgroundColor()));
            _imageScene.initialize();
        }

        void Widget::resizeGL(int w, int h) {
            core::PerspectiveCamera & camera = _options.camera();
            camera.resizeScreen(core::Size(w, h));
            glViewport(0, 0, w, h);
        }

        void Widget::paintEvent(QPaintEvent * e) {
            if (!_anno) {
                return;
            }

            QPainter painter(this);
            painter.beginNativePainting();
            qglClearColor(MakeQColor(_options.backgroundColor()));
            glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);

            core::PerspectiveCamera & camera = _options.camera();
            camera.resizeScreen(core::Size(width(), height()));

            _imageScene.render(_options);
            for (auto & s : _polygonLineScenes) {
                s.render(_options);
            }

            painter.endNativePainting();
            swapBuffers();
        }

        void Widget::mousePressEvent(QMouseEvent * e) {
            if (!_anno) {
                return;
            }
            if (e->buttons() & Qt::MidButton) {
                _lastPos = e->pos();
                setCursor(Qt::OpenHandCursor);
            } else {
                BaseClass::mousePressEvent(e);
            }
        }

        void Widget::mouseMoveEvent(QMouseEvent * e) {
            if (!_anno) {
                return;
            }
            QVector3D t(e->pos() - _lastPos);
            t.setX(-t.x());
            if (e->buttons() & Qt::MidButton) {
                _options.camera().moveCenterWithEyeFixed(MakeCoreVec(t));
                setCursor(Qt::ClosedHandCursor);
                update();
            } else {
                BaseClass::mouseMoveEvent(e);
            }
            _lastPos = e->pos();
        }

        void Widget::mouseReleaseEvent(QMouseEvent * e) {
            if (!_anno) {
                return;
            }
            unsetCursor();
            BaseClass::mouseReleaseEvent(e);
        }

        void Widget::wheelEvent(QWheelEvent * e) {
            if (!_anno) {
                return;
            }
            BaseClass::wheelEvent(e);
        }

        void Widget::keyPressEvent(QKeyEvent * e) {
            if (!_anno) {
                return;
            }
            BaseClass::keyPressEvent(e);
        }

        void Widget::rebuildPolygonLineScenes() {
            _polygonLineScenes.clear();
            
            // polygons
            assert(_anno->polygons.size() == _anno->polygonAlongVPIds.size() &&
                _anno->polygons.size() == _anno->polygonTowardVPIds.size() &&
                _anno->polygons.size() == _anno->polygonAreClutters.size());           
            
            for (int i = 0; i < _anno->polygons.size(); i++) {
                SceneBuilder sb;
                sb.installingOptions().defaultShaderSource = gui::OpenGLShaderSourceDescriptor::XTriangles;
                auto poly = _anno->polygons[i];
                for (auto & c : poly.corners) {
                    c = normalize(c);
                }
                Color color;

                ColorTable rgb = gui::RGB;

                if (_anno->polygonTowardVPIds[i] != -1) {
                    color = rgb[_anno->polygonTowardVPIds[i]];
                } else if (_anno->polygonAlongVPIds[i] != -1) {
                    if (_anno->polygonAlongVPIds[i] == _anno->vertVPId) {
                        color = gui::Gray;
                    } else {
                        color = gui::Yellow;
                    }
                } else if(_anno->polygonAreClutters[i]){
                    color = gui::Black;
                } else {
                    color = gui::White;
                }
                sb.add(ColorAs(poly, color));
                _polygonLineScenes.push_back(sb.scene());
            }

            // occlusion chains                   
            for (auto & occ : _anno->occlusions) {
                SceneBuilder sb;
                sb.lineWidth(5.0);
                sb.installingOptions().discretizeOptions.color = gui::Cyan;
                sb.installingOptions().defaultShaderSource = gui::OpenGLShaderSourceDescriptor::XLines;

                auto o = occ;
                for (auto & p : o.points) {
                    p = normalize(p) / 2.0;
                }
                sb.add(o);
                _polygonLineScenes.push_back(sb.scene());
            }
        }

        void Widget::rebuildStrokeScene() {
            SceneBuilder sb;
            sb.lineWidth(5.0);
            sb.installingOptions().discretizeOptions.color = gui::Orange;
            sb.installingOptions().defaultShaderSource = gui::OpenGLShaderSourceDescriptor::XLines;
            auto o = _chain;
            for (auto & p : o.points) {
                p = normalize(p) / 2.0;
            }
            sb.add(o);
            _strokeScene = sb.scene();
        }


        MainWin::MainWin() : QMainWindow(nullptr) {
            _w = new Widget;
            setCentralWidget(_w);
        }

        MainWin::~MainWin() {

        }

        void MainWin::selectFile(const QString & fname) {
            assert(!fname.isEmpty());
            QFileInfo finfo(fname);
            assert(finfo.exists());
            
            // get annotation filename
            auto annoFileName = finfo.fileName() + ".anno.cereal";
            QFileInfo annofinfo(annoFileName);

            // if exist, load it
            if (!annofinfo.exists() || !LoadFromDisk(annofinfo.absoluteFilePath().toStdString(), _anno)) {
                // can't load
                auto image = gui::MakeCVMat(QImage(fname));
                _anno.view = CreatePanoramicView(image);

                // collect lines in each view
                auto cams = CreateCubicFacedCameras(_anno.view.camera, image.rows, image.rows, image.rows * 0.4);
                std::vector<Line3> rawLine3s;
                for (int i = 0; i < cams.size(); i++) {
                    auto pim = _anno.view.sampled(cams[i]).image;
                    LineSegmentExtractor lineExtractor;
                    lineExtractor.params().algorithm = LineSegmentExtractor::LSD;
                    auto ls = lineExtractor(pim, 3, 300); // use pyramid
                    for (auto & l : ls) {
                        rawLine3s.emplace_back(normalize(cams[i].toSpace(l.first)),
                            normalize(cams[i].toSpace(l.second)));
                    }
                }
                rawLine3s = MergeLines(rawLine3s, DegreesToRadians(1)); 

                // estimate vp
                auto line3s = ClassifyEachAs(rawLine3s, -1);
                auto vps = EstimateVanishingPointsAndClassifyLines(line3s);
                int vertVPId = NearestDirectionId(vps, Vec3(0, 0, 1));

                _anno.lines = std::move(line3s);
                _anno.vps = std::move(vps);
                _anno.vertVPId = vertVPId;
            }  

            _w->setCurAnnotation(&_anno);
        }

    }
}
