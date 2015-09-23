#include "../gui/qttools.hpp"
#include "pi_graph_annotation_widgets.hpp"

namespace pano {
    namespace experimental {

        using namespace pano::gui;

        PIAnnotationWidget::PIAnnotationWidget(QWidget * parent /*= nullptr*/) : QGLWidget(parent), _anno(nullptr), _state(Idle) {
            _segControl.orientationClaz = -1;
            _segControl.orientationNotClaz = -1;
            _segControl.used = true;

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
                _segControl.orientationClaz = 0; _segControl.orientationNotClaz = -1; _segControl.used = true;
                _state = CreatingPolygon;
                clearStroke();
                update();
            });
            connect(defaultAction = bas->addAction(tr("Polygon Brush: Toward Green VP")), &QAction::triggered, [this]() {
                _segControl.orientationClaz = 1; _segControl.orientationNotClaz = -1; _segControl.used = true;
                _state = CreatingPolygon;
                clearStroke();
                update();
            });
            connect(defaultAction = bas->addAction(tr("Polygon Brush: Toward Blue VP")), &QAction::triggered, [this]() {
                _segControl.orientationClaz = 2; _segControl.orientationNotClaz = -1; _segControl.used = true;
                _state = CreatingPolygon;
                clearStroke();
                update();
            });

            connect(defaultAction = bas->addAction(tr("Polygon Brush: Vertical")), &QAction::triggered, [this]() {
                _segControl.orientationClaz = -1; _segControl.orientationNotClaz = _anno->vertVPId; _segControl.used = true;
                _state = CreatingPolygon;
                clearStroke();
                update();
            });

            connect(defaultAction = bas->addAction(tr("Polygon Brush: Free")), &QAction::triggered, [this]() {
                _segControl.orientationClaz = -1; _segControl.orientationNotClaz = -1; _segControl.used = true;
                _state = CreatingPolygon;
                clearStroke();
                update();
            });

            connect(defaultAction = bas->addAction(tr("Polygon Brush: Clutter")), &QAction::triggered, [this]() {
                _segControl.orientationClaz = -1; _segControl.orientationNotClaz = -1; _segControl.used = false;
                _state = CreatingPolygon;
                clearStroke();
                update();
            });

            {
                QAction * sep = new QAction(bas);
                sep->setSeparator(true);
                bas->addAction(sep);
            }

            for (auto a : bas->actions()) a->setCheckable(true);
            bas->setExclusive(true);
            defaultAction->setChecked(true);
            addActions(bas->actions());

            connect(defaultAction = new QAction(tr("Settings"), this), &QAction::triggered, [this]() {
                PopUpGui(_options, this);
                update();
            });
            addAction(defaultAction);
        }

        PIAnnotationWidget::~PIAnnotationWidget() {

        }

        void PIAnnotationWidget::setCurAnnotation(PIAnnotation * anno) {
            _state = Idle;
            _anno = anno;

            assert(!anno->view.image.empty());
            auto im = anno->view.image;

            _polygonsDeleted = std::vector<bool>(anno->polygons.size(), false);
            _occlusionsDeleted = std::vector<bool>(anno->occlusions.size(), false);

            SceneBuilder sb;
            ResourceStore::set("tex", im);
            Sphere3 sp;
            sp.center = Origin();
            sp.radius = 10.0;
            sb.begin(sp).shaderSource(OpenGLShaderSourceDescriptor::XPanorama).resource("tex").end();

            gui::ColorTable rgb = gui::RGBGreys;
            for (int i = 0; i < anno->vps.size(); i++) {
                sb.installingOptions().discretizeOptions.color = rgb[i].blendWith(gui::White, 0.4);
                sb.begin(normalize(anno->vps[i]) * 0.05).shaderSource(OpenGLShaderSourceDescriptor::XPoints).pointSize(40).end();
                sb.begin(-normalize(anno->vps[i]) * 0.05).shaderSource(OpenGLShaderSourceDescriptor::XPoints).pointSize(40).end();
            }

            _options.panoramaAspectRatio(im.rows / float(im.cols));
            _options.panoramaHoriCenterRatio(0.5f);
            _options.camera(PerspectiveCamera(500, 500, Point2(250, 250), 200, Origin(), X(), -Z()));
            _options.renderMode(gui::RenderModeFlag::All);

            _imageScene = sb.scene();

            clearStroke();
            rebuildLinesScene();
            rebuildPolygonScenes();
            rebuildOcclusionScenes();
            update();
        }


        void PIAnnotationWidget::paintEvent(QPaintEvent * e) {
            if (!_anno) {
                return;
            }
            _imageScene.initialize();
            for (auto & s : _polygonScenes) {
                s.initialize();
            }
            for (auto & s : _occlusionScenes) {
                s.initialize();
            }
            _linesScene.initialize();
            _strokeScene.initialize();

            QPainter painter(this);
            painter.beginNativePainting();

            glEnable(GL_ALPHA_TEST);
            glEnable(GL_MULTISAMPLE);
            GLint bufs;
            GLint samples;
            glGetIntegerv(GL_SAMPLE_BUFFERS, &bufs);
            glGetIntegerv(GL_SAMPLES, &samples);


            qglClearColor(MakeQColor(_options.backgroundColor()));
            glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);

            core::PerspectiveCamera & camera = _options.camera();
            camera.resizeScreen(core::Size(width(), height()));

            _imageScene.render(_options);
            for (auto & s : _polygonScenes) {
                s.render(_options);
            }
            for (auto & s : _occlusionScenes) {
                s.render(_options);
            }
            _linesScene.render(_options);
            _strokeScene.render(_options);

            painter.endNativePainting();
            swapBuffers();
        }

        void PIAnnotationWidget::mousePressEvent(QMouseEvent * e) {
            if (!_anno) {
                return;
            }
            if (e->buttons() & Qt::MidButton) {
                _lastPos = e->pos();
                setCursor(Qt::OpenHandCursor);
            } else if ((e->buttons() & Qt::LeftButton) && _state != Idle) {
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

        void PIAnnotationWidget::mouseMoveEvent(QMouseEvent * e) {
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

        void PIAnnotationWidget::mouseReleaseEvent(QMouseEvent * e) {
            if (!_anno) {
                return;
            }
            unsetCursor();
            BaseClass::mouseReleaseEvent(e);
        }

        void PIAnnotationWidget::wheelEvent(QWheelEvent * e) {
            if (!_anno) {
                return;
            }
            _options.camera().setFocal(_options.camera().focal() * exp(e->delta() / 1000.0));
            update();
            BaseClass::wheelEvent(e);
        }

        void PIAnnotationWidget::keyPressEvent(QKeyEvent * e) {
            if (!_anno) {
                return;
            }
            if (e->key() == Qt::Key_Space) {
                if (_state == CreatingOcclusion) {
                    acceptAsOcclusion();
                    rebuildOcclusionScenes();
                } else if (_state == CreatingPolygon) {
                    acceptAsPolygon(_segControl.orientationClaz, _segControl.orientationNotClaz, _segControl.used);
                    rebuildPolygonScenes();
                }
                clearStroke();
                update();
            } else if (e->key() == Qt::Key_Escape) {
                clearStroke();
                update();
            }
            BaseClass::keyPressEvent(e);
        }


        void PIAnnotationWidget::rebuildLinesScene() {
            SceneBuilder sb;
            sb.installingOptions().lineWidth = 3.0;
            sb.installingOptions().defaultShaderSource = gui::OpenGLShaderSourceDescriptor::XLines;
            ColorTable rgb = gui::RGB;
            rgb.exceptionalColor() = gui::Black;
            for (int i = 0; i < _anno->lines.size(); i++) {
                sb.add(ColorAs(normalize(_anno->lines[i].component) * 0.15, rgb[_anno->lines[i].claz]));
            }
            _linesScene = sb.scene();
        }

        void PIAnnotationWidget::rebuildPolygonScenes() {
            _polygonScenes.clear();
            // polygons       
            for (int i = 0; i < _anno->polygons.size(); i++) {
                SceneBuilder sb;
                sb.installingOptions().defaultShaderSource = gui::OpenGLShaderSourceDescriptor::XTriangles;
                auto poly = _anno->polygons[i];
                for (auto & c : poly.polygon.corners) {
                    c = normalize(c);
                }
                Color color;

                ColorTable rgb = gui::RGB;
                auto & p = _anno->polygons[i];
                if (p.control.orientationClaz != -1) {
                    color = rgb[p.control.orientationClaz];
                } else if (p.control.orientationNotClaz != -1) {
                    color = rgb[p.control.orientationNotClaz].blendWith(gui::White, 0.5);
                } else if (!p.control.used) {
                    color = gui::Black;
                } else {
                    color = gui::White;
                }
                sb.add(ColorAs(poly.polygon, color));
                _polygonScenes.push_back(sb.scene());
            }
        }

        void PIAnnotationWidget::rebuildOcclusionScenes() {
            _occlusionScenes.clear();
            // occlusion                   
            for (int i = 0; i < _anno->occlusions.size(); i++) {
                SceneBuilder sb;
                sb.installingOptions().lineWidth = 3.0;
                sb.installingOptions().pointSize = 10.0;
                sb.installingOptions().discretizeOptions.color = gui::Cyan;

                auto & chain = _anno->occlusions[i].chain;
                for (int i = 1; i < chain.size(); i++) {
                    auto p1 = normalize(chain[i - 1]) / 2.0;
                    auto p2 = normalize(chain[i]) / 2.0;
                    sb.begin(Line3(p1, p2)).shaderSource(gui::OpenGLShaderSourceDescriptor::XLines).end();
                    sb.begin(p2).shaderSource(gui::OpenGLShaderSourceDescriptor::XPoints).end();
                }
                _occlusionScenes.push_back(sb.scene());
            }
        }

        void PIAnnotationWidget::rebuildStrokeScene() {
            if (_chain.empty()) {
                _strokeScene = gui::Scene();
                return;
            }

            SceneBuilder sb;
            sb.installingOptions().lineWidth = 5.0;
            sb.installingOptions().pointSize = 7.0;
            sb.installingOptions().discretizeOptions.color = gui::Orange;

            auto & chain = _chain;
            for (int i = 1; i < chain.size(); i++) {
                auto p1 = normalize(chain[i - 1]) / 2.0;
                auto p2 = normalize(chain[i]) / 2.0;
                sb.begin(Line3(p1, p2)).shaderSource(gui::OpenGLShaderSourceDescriptor::XLines).end();
                sb.begin(p2).shaderSource(gui::OpenGLShaderSourceDescriptor::XPoints).end();
            }

            _strokeScene = sb.scene();
        }

        void PIAnnotationWidget::clearStroke() {
            _chain.clear();
            _strokeScene = gui::Scene();
        }

        void PIAnnotationWidget::acceptAsPolygon(int towardVPId, int alongVPId, bool used) {
            AnnotedPolygon poly;
            poly.control.orientationClaz = towardVPId;
            poly.control.orientationNotClaz = alongVPId;
            poly.control.used = used;
            poly.polygon = _chain;
            _anno->polygons.push_back(std::move(poly));
            _polygonsDeleted.push_back(false);
        }

        void PIAnnotationWidget::acceptAsOcclusion() {
            _anno->occlusions.push_back({ std::move(_chain) });
            _occlusionsDeleted.push_back(false);
        }



        //MainWin::MainWin() : QMainWindow(nullptr) {
        //    _anno.reset();
        //    _w = new PIAnnotationWidget(this);
        //    setCentralWidget(_w);

        //    QAction * defaultAction = nullptr;

        //    // menus...
        //    auto menuFile = menuBar()->addMenu(tr("File"));
        //    connect(defaultAction = menuFile->addAction(tr("Open Image (with Annotation if possible)")), &QAction::triggered, [this]() {
        //        auto fname = QFileDialog::getOpenFileName(this, tr("Select An Image File"), tr("H:\\DataSet"),
        //            tr("Image Files (*.png;*.jpg;*.jpeg);;All Files (*.*)"));
        //        if (fname.isEmpty())
        //            return;
        //        selectFile(fname);
        //    });
        //    defaultAction->setShortcut(tr("Ctrl+o"));

        //    connect(defaultAction = menuFile->addAction(tr("Save Annotation")), &QAction::triggered, [this]() {
        //        if (_fname.isEmpty())
        //            return;
        //        QFileInfo annofinfo(annoFileName());
        //        if (!SaveToDisk(annofinfo.absoluteFilePath().toStdString(), _anno)) {
        //            QMessageBox::critical(this, tr("Saving Annotation Failed!"), tr("Can't Save To \"") + annofinfo.absoluteFilePath() + "\"!");
        //        } else {
        //            qDebug() << "Success!";
        //        }
        //    });
        //    defaultAction->setShortcut(tr("Ctrl+s"));
        //}

        //MainWin::~MainWin() {

        //}

        //void MainWin::clear() {
        //    _fname.clear();
        //    _anno.reset();
        //}

        //void MainWin::selectFile(const QString & fname) {
        //    QFileInfo finfo(fname);
        //    assert(finfo.exists());

        //    _fname = fname;

        //    // get annotation filename
        //    QFileInfo annofinfo(annoFileName());

        //    // if exist, load it
        //    if (!annofinfo.exists() || !LoadFromDisk(annofinfo.absoluteFilePath().toStdString(), _anno)) {
        //        // can't load
        //        auto image = gui::MakeCVMat(QImage(fname));
        //        MakePanorama(image);
        //        ResizeToHeight(image, 700);

        //        _anno.view = CreatePanoramicView(image);

        //        // collect lines in each view
        //        auto cams = CreateCubicFacedCameras(_anno.view.camera, image.rows, image.rows, image.rows * 0.4);
        //        std::vector<Line3> rawLine3s;
        //        for (int i = 0; i < cams.size(); i++) {
        //            auto pim = _anno.view.sampled(cams[i]).image;
        //            LineSegmentExtractor lineExtractor;
        //            lineExtractor.params().algorithm = LineSegmentExtractor::LSD;
        //            auto ls = lineExtractor(pim, 3, 300); // use pyramid
        //            for (auto & l : ls) {
        //                rawLine3s.emplace_back(normalize(cams[i].toSpace(l.first)),
        //                    normalize(cams[i].toSpace(l.second)));
        //            }
        //        }
        //        rawLine3s = MergeLines(rawLine3s, DegreesToRadians(1));

        //        // estimate vp
        //        auto line3s = ClassifyEachAs(rawLine3s, -1);
        //        auto vps = EstimateVanishingPointsAndClassifyLines(line3s);
        //        OrderVanishingPoints(vps);
        //        int vertVPId = NearestDirectionId(vps, Vec3(0, 0, 1));

        //        _anno.lines = std::move(line3s);
        //        _anno.vps = std::move(vps);
        //        _anno.vertVPId = vertVPId;
        //    }

        //    _w->setCurAnnotation(&_anno);
        //}

        //QString MainWin::annoFileName() const {
        //    return QString::fromStdString(experimental::AnnotationFilePath(_fname.toStdString()));
        //}

    }
}
