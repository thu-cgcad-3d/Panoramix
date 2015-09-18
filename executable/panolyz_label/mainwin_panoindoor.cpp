#include "../../src/gui/qttools.hpp"
#include "mainwin_panoindoor.hpp"


namespace panolyz {
    namespace PanoramaIndoor {

        using namespace pano::gui;

        Widget::Widget(QWidget * parent /*= nullptr*/) : QGLWidget(parent), _anno(nullptr), _state(Idle) {
            _towardVPId = -1;
            _alongVPId = -1;
            _isClutter = false;

            setMouseTracking(true);
            setAutoBufferSwap(false);
            grabKeyboard();

            _chain.closed = false;
            setContextMenuPolicy(Qt::ActionsContextMenu);
            QActionGroup * bas = new QActionGroup(this);

            QAction * defaultAction = nullptr;
            connect(defaultAction = bas->addAction(tr("Cancel")), &QAction::triggered, [this]() { 
                _state = Idle; 
                clearStroke();
                update();
            });            
            connect(defaultAction = bas->addAction(tr("Create Polygon")), &QAction::triggered, [this]() {
                _state = CreatingPolygon;
                clearStroke();
                update();
            });
            connect(defaultAction = bas->addAction(tr("Create Occlusion")), &QAction::triggered, [this]() {
                _state = CreatingOcclusion;
                clearStroke();
                update();
            });

            {
                QAction * sep = new QAction(bas);
                sep->setSeparator(true);
                bas->addAction(sep);
            }
            
            connect(defaultAction = bas->addAction(tr("Polygon Brush: Toward Red VP")), &QAction::triggered, [this]() {
                _towardVPId = 0; _alongVPId = -1; _isClutter = false;
            });
            connect(defaultAction = bas->addAction(tr("Polygon Brush: Toward Green VP")), &QAction::triggered, [this]() {
                _towardVPId = 1; _alongVPId = -1; _isClutter = false;
            });
            connect(defaultAction = bas->addAction(tr("Polygon Brush: Toward Blue VP")), &QAction::triggered, [this]() {
                _towardVPId = 2; _alongVPId = -1; _isClutter = false;
            });

            connect(defaultAction = bas->addAction(tr("Polygon Brush: Vertical")), &QAction::triggered, [this]() {
                _towardVPId = -1; _alongVPId = _anno->vertVPId; _isClutter = false;
            });

            connect(defaultAction = bas->addAction(tr("Polygon Brush: Free")), &QAction::triggered, [this]() {
                _towardVPId = -1; _alongVPId = -1; _isClutter = false;
            });

            connect(defaultAction = bas->addAction(tr("Polygon Brush: Clutter")), &QAction::triggered, [this]() {
                _towardVPId = -1; _alongVPId = -1; _isClutter = true;
            });

            for (auto a : bas->actions()) a->setCheckable(true);
            bas->setExclusive(true);
            defaultAction->setChecked(true);
            addActions(bas->actions());
        }

        Widget::~Widget() {

        }


        void Widget::setCurAnnotation(PanoIndoorAnnotation * anno) {
            _state = Idle;

            assert(!anno->view.image.empty());
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
            _options.renderMode(gui::RenderModeFlag::All);

            _imageScene = sb.scene();

            clearStroke();
            rebuildPolygonLineScenes();
            update();
        }


        void Widget::paintEvent(QPaintEvent * e) {
            if (!_anno) {
                return;
            }
            _imageScene.initialize();
            for (auto & s : _polygonLineScenes) {
                s.initialize();
            }
            _strokeScene.initialize();

            QPainter painter(this);
            painter.beginNativePainting();

            glEnable(GL_MULTISAMPLE);
            GLint bufs;
            GLint samples;
            glGetIntegerv(GL_SAMPLE_BUFFERS, &bufs);
            glGetIntegerv(GL_SAMPLES, &samples);
            //qDebug("Have %d buffers and %d samples", bufs, samples);


            qglClearColor(MakeQColor(_options.backgroundColor()));
            glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);

            core::PerspectiveCamera & camera = _options.camera();
            camera.resizeScreen(core::Size(width(), height()));

            _imageScene.render(_options);
            for (auto & s : _polygonLineScenes) {
                s.render(_options);
            }
            _strokeScene.render(_options);

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
            } else if((e->buttons() & Qt::LeftButton) && _state != Idle){
                if (_state == CreatingOcclusion) {
                    _chain.points.push_back(normalize(_options.camera().toSpace(gui::MakeCoreVec(e->pos()))));
                } else if (_state == CreatingPolygon) {
                    _chain.points.push_back(normalize(_options.camera().toSpace(gui::MakeCoreVec(e->pos()))));
                }
                rebuildStrokeScene();
                update();
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

            update();
            BaseClass::wheelEvent(e);
        }

        void Widget::keyPressEvent(QKeyEvent * e) {
            if (!_anno) {
                return;
            }
            if (e->key() == Qt::Key_Space) {
                if (_state == CreatingOcclusion) {
                    acceptAsOcclusion();
                } else if (_state == CreatingPolygon) {
                    acceptAsPolygon(_towardVPId, _alongVPId, _isClutter);
                }
                clearStroke();
                rebuildPolygonLineScenes();
                update();
            } else if (e->key() == Qt::Key_Escape) {
                clearStroke();
                update();
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
                sb.installingOptions().lineWidth = 5.0;
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
            if (_chain.empty()) {
                _strokeScene = gui::Scene();
                return;
            }

            SceneBuilder sb;
            sb.installingOptions().lineWidth = 5.0;
            sb.installingOptions().discretizeOptions.color = gui::Orange;
            sb.installingOptions().defaultShaderSource = gui::OpenGLShaderSourceDescriptor::XLines;
            auto o = _chain;
            for (auto & p : o.points) {
                p = normalize(p) / 2.0;
            }
            sb.add(o);
            _strokeScene = sb.scene();
        }

        void Widget::clearStroke() {
            _chain.points.clear();
            _strokeScene = gui::Scene();
        }

        void Widget::acceptAsPolygon(int towardVPId, int alongVPId, bool clutter) {
            _anno->polygons.push_back(Polygon3(_chain));
            _anno->polygonTowardVPIds.push_back(towardVPId);
            _anno->polygonAlongVPIds.push_back(alongVPId);
            _anno->polygonAreClutters.push_back(clutter);
            _chain.points.clear();
            _strokeScene = gui::Scene();
            rebuildPolygonLineScenes();
        }

        void Widget::acceptAsOcclusion() {
            _anno->occlusions.push_back(std::move(_chain));
            _chain.points.clear();
            _strokeScene = gui::Scene();
            rebuildPolygonLineScenes();
        }


        MainWin::MainWin() : QMainWindow(nullptr) {
            _anno.reset();
            _w = new Widget(this);
            setCentralWidget(_w);

            QAction * defaultAction = nullptr;
            
            // menus...
            auto menuFile = menuBar()->addMenu(tr("File"));
            connect(defaultAction = menuFile->addAction(tr("Open Image (with Annotation if possible)")), &QAction::triggered, [this]() {
                auto fname = QFileDialog::getOpenFileName(this, tr("Select An Image File"), tr("H:\\DataSet"), 
                    tr("Image Files (*.png;*.jpg;*.jpeg);;All Files (*.*)"));
                if (fname.isEmpty())
                    return;
                selectFile(fname);
            });
            defaultAction->setShortcut(tr("Ctrl+o"));

            connect(defaultAction = menuFile->addAction(tr("Save Annotation")), &QAction::triggered, [this]() {
                if (_fname.isEmpty())
                    return;
                QFileInfo annofinfo(annoFileName());
                if (!SaveToDisk(annofinfo.absoluteFilePath().toStdString(), _anno)) {
                    QMessageBox::critical(this, tr("Saving Annotation Failed!"), tr("Can't Save To \"") + annofinfo.absoluteFilePath() + "\"!");
                } else {
                    qDebug() << "Success!";
                }
            });
            defaultAction->setShortcut(tr("Ctrl+s"));
        }

        MainWin::~MainWin() {

        }

        void MainWin::clear() {
            _fname.clear();
            _anno.reset();
        }

        void MainWin::selectFile(const QString & fname) {
            QFileInfo finfo(fname);
            assert(finfo.exists());

            _fname = fname;
            
            // get annotation filename
            QFileInfo annofinfo(annoFileName());

            // if exist, load it
            if (!annofinfo.exists() || !LoadFromDisk(annofinfo.absoluteFilePath().toStdString(), _anno)) {
                // can't load
                auto image = gui::MakeCVMat(QImage(fname));
                MakePanorama(image);
                ResizeToHeight(image, 700);

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
                OrderVanishingPoints(vps);
                int vertVPId = NearestDirectionId(vps, Vec3(0, 0, 1));

                _anno.lines = std::move(line3s);
                _anno.vps = std::move(vps);
                _anno.vertVPId = vertVPId;
            }  

            _w->setCurAnnotation(&_anno);
        }

        QString MainWin::annoFileName() const {
            QFileInfo finfo(_fname);
            if (!finfo.exists())
                return QString();

            auto annoFileName = finfo.absoluteFilePath() + ".anno.cereal";
            return annoFileName;
        }

    }
}
