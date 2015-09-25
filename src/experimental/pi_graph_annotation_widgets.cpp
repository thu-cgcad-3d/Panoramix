#include "../gui/qttools.hpp"
#include "pi_graph_annotation_widgets.hpp"

namespace pano {
    namespace experimental {

        using namespace pano::gui;

        //PIAnnotationWidget::PIAnnotationWidget(QWidget * parent /*= nullptr*/) : QGLWidget(parent), _anno(nullptr), _state(Idle) {
        //    _segControl.orientationClaz = -1;
        //    _segControl.orientationNotClaz = -1;
        //    _segControl.used = true;

        //    _showPolygons = _showLines = _showOcclusions = _showVPs = true;

        //    setMouseTracking(true);
        //    setAutoBufferSwap(false);
        //    grabKeyboard();

        //    _chain.closed = false;
        //    setContextMenuPolicy(Qt::ActionsContextMenu);
        //    QActionGroup * bas = new QActionGroup(this);

        //    QAction * defaultAction = nullptr;
        //    connect(defaultAction = bas->addAction(tr("Cancel")), &QAction::triggered, [this]() {
        //        _state = Idle;
        //        clearStroke();
        //        update();
        //    });

        //    connect(defaultAction = bas->addAction(tr("Create Lines")), &QAction::triggered, [this]() {
        //        _state = CreatingLine;
        //        clearStroke();
        //        update();
        //    });
        //    connect(defaultAction = bas->addAction(tr("Create Occlusion")), &QAction::triggered, [this]() {
        //        _state = CreatingOcclusion;
        //        clearStroke();
        //        update();
        //    });

        //    {
        //        QAction * sep = new QAction(bas);
        //        sep->setSeparator(true);
        //        bas->addAction(sep);
        //    }

        //    connect(defaultAction = bas->addAction(tr("Polygon Brush: Toward Red VP")), &QAction::triggered, [this]() {
        //        _segControl.orientationClaz = 0; _segControl.orientationNotClaz = -1; _segControl.used = true;
        //        _state = CreatingPolygon;
        //        clearStroke();
        //        update();
        //    });
        //    connect(defaultAction = bas->addAction(tr("Polygon Brush: Toward Green VP")), &QAction::triggered, [this]() {
        //        _segControl.orientationClaz = 1; _segControl.orientationNotClaz = -1; _segControl.used = true;
        //        _state = CreatingPolygon;
        //        clearStroke();
        //        update();
        //    });
        //    connect(defaultAction = bas->addAction(tr("Polygon Brush: Toward Blue VP")), &QAction::triggered, [this]() {
        //        _segControl.orientationClaz = 2; _segControl.orientationNotClaz = -1; _segControl.used = true;
        //        _state = CreatingPolygon;
        //        clearStroke();
        //        update();
        //    });

        //    connect(defaultAction = bas->addAction(tr("Polygon Brush: Vertical")), &QAction::triggered, [this]() {
        //        _segControl.orientationClaz = -1; _segControl.orientationNotClaz = _anno->vertVPId; _segControl.used = true;
        //        _state = CreatingPolygon;
        //        clearStroke();
        //        update();
        //    });

        //    connect(defaultAction = bas->addAction(tr("Polygon Brush: Free")), &QAction::triggered, [this]() {
        //        _segControl.orientationClaz = -1; _segControl.orientationNotClaz = -1; _segControl.used = true;
        //        _state = CreatingPolygon;
        //        clearStroke();
        //        update();
        //    });

        //    connect(defaultAction = bas->addAction(tr("Polygon Brush: Clutter")), &QAction::triggered, [this]() {
        //        _segControl.orientationClaz = -1; _segControl.orientationNotClaz = -1; _segControl.used = false;
        //        _state = CreatingPolygon;
        //        clearStroke();
        //        update();
        //    });

        //    {
        //        QAction * sep = new QAction(bas);
        //        sep->setSeparator(true);
        //        bas->addAction(sep);
        //    }

        //    for (auto a : bas->actions()) a->setCheckable(true);
        //    bas->setExclusive(true);
        //    defaultAction->setChecked(true);
        //    addActions(bas->actions());


        //    connect(defaultAction = new QAction(tr("Show Lines"), this), &QAction::triggered, [this](bool checked) {
        //        _showLines = checked;
        //        update();
        //    });
        //    defaultAction->setCheckable(true);
        //    defaultAction->setChecked(_showLines);
        //    addAction(defaultAction);

        //    connect(defaultAction = new QAction(tr("Show Polygons"), this), &QAction::triggered, [this](bool checked) {
        //        _showPolygons = checked;
        //        update();
        //    });
        //    defaultAction->setCheckable(true);
        //    defaultAction->setChecked(_showPolygons);
        //    addAction(defaultAction);

        //    connect(defaultAction = new QAction(tr("Show Occlusions"), this), &QAction::triggered, [this](bool checked) {
        //        _showOcclusions = checked;
        //        update();
        //    });
        //    defaultAction->setCheckable(true);
        //    defaultAction->setChecked(_showOcclusions);
        //    addAction(defaultAction);

        //    connect(defaultAction = new QAction(tr("Show VPs"), this), &QAction::triggered, [this](bool checked) {
        //        _showVPs = checked;
        //        update();
        //    });
        //    defaultAction->setCheckable(true);
        //    defaultAction->setChecked(_showVPs);
        //    addAction(defaultAction);

        //    {
        //        QAction * sep = new QAction(bas);
        //        sep->setSeparator(true);
        //        bas->addAction(sep);
        //    }

        //    connect(defaultAction = new QAction(tr("Settings"), this), &QAction::triggered, [this]() {
        //        PopUpGui(_options, this);
        //        update();
        //    });
        //    addAction(defaultAction);
        //}

        //PIAnnotationWidget::~PIAnnotationWidget() {

        //}

        //void PIAnnotationWidget::setCurAnnotation(PIAnnotation * anno) {
        //    _state = Idle;
        //    _anno = anno;

        //    assert(!anno->view.image.empty());
        //    auto im = anno->view.image;

        //    _polygonsDeleted = std::vector<bool>(anno->polygons.size(), false);
        //    _occlusionsDeleted = std::vector<bool>(anno->occlusions.size(), false);

        //    SceneBuilder sb;
        //    ResourceStore::set("tex", im);
        //    Sphere3 sp;
        //    sp.center = Origin();
        //    sp.radius = 10.0;
        //    sb.begin(sp).shaderSource(OpenGLShaderSourceDescriptor::XPanorama).resource("tex").end();

        //    gui::ColorTable rgb = gui::RGBGreys;
        //    for (int i = 0; i < anno->vps.size(); i++) {
        //        sb.installingOptions().discretizeOptions.color = rgb[i].blendWith(gui::White, 0.4);
        //        sb.begin(normalize(anno->vps[i]) * 0.05).shaderSource(OpenGLShaderSourceDescriptor::XPoints).pointSize(40).end();
        //        sb.begin(-normalize(anno->vps[i]) * 0.05).shaderSource(OpenGLShaderSourceDescriptor::XPoints).pointSize(40).end();
        //    }

        //    _options.panoramaAspectRatio(im.rows / float(im.cols));
        //    _options.panoramaHoriCenterRatio(0.5f);
        //    _options.camera(PerspectiveCamera(500, 500, Point2(250, 250), 200, Origin(), X(), -Z()));
        //    _options.renderMode(gui::RenderModeFlag::All);

        //    _imageScene = sb.scene();

        //    clearStroke();
        //    rebuildLinesScene();
        //    rebuildPolygonScenes();
        //    rebuildOcclusionScenes();
        //    update();
        //}


        //void PIAnnotationWidget::paintEvent(QPaintEvent * e) {
        //    if (!_anno) {
        //        return;
        //    }
        //    _imageScene.initialize();
        //    for (auto & s : _polygonScenes) {
        //        s.initialize();
        //    }
        //    for (auto & s : _occlusionScenes) {
        //        s.initialize();
        //    }
        //    _linesScene.initialize();
        //    _strokeScene.initialize();

        //    QPainter painter(this);
        //    painter.beginNativePainting();

        //    glEnable(GL_ALPHA_TEST);
        //    glEnable(GL_MULTISAMPLE);
        //    GLint bufs;
        //    GLint samples;
        //    glGetIntegerv(GL_SAMPLE_BUFFERS, &bufs);
        //    glGetIntegerv(GL_SAMPLES, &samples);


        //    qglClearColor(MakeQColor(_options.backgroundColor()));
        //    glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);

        //    core::PerspectiveCamera & camera = _options.camera();
        //    camera.resizeScreen(core::Size(width(), height()));

        //    _imageScene.render(_options);
        //    if (_showPolygons) {
        //        for (auto & s : _polygonScenes) {
        //            s.render(_options);
        //        }
        //    }
        //    if (_showOcclusions) {
        //        for (auto & s : _occlusionScenes) {
        //            s.render(_options);
        //        }
        //    }
        //    if (_showLines) {
        //        _linesScene.render(_options);
        //    }            
        //    _strokeScene.render(_options);

        //    painter.endNativePainting();
        //    swapBuffers();
        //}

        //void PIAnnotationWidget::mousePressEvent(QMouseEvent * e) {
        //    if (!_anno) {
        //        return;
        //    }
        //    if (e->buttons() & Qt::MidButton) {
        //        _lastPos = e->pos();
        //        setCursor(Qt::OpenHandCursor);
        //    } else if ((e->buttons() & Qt::LeftButton) && _state != Idle) {
        //        if (_state == CreatingOcclusion) {
        //            _chain.points.push_back(normalize(_options.camera().toSpace(gui::MakeCoreVec(e->pos()))));
        //        } else if (_state == CreatingPolygon) {
        //            _chain.points.push_back(normalize(_options.camera().toSpace(gui::MakeCoreVec(e->pos()))));
        //        } else if (_state == CreatingLine) {
        //            _chain.points.push_back(normalize(_options.camera().toSpace(gui::MakeCoreVec(e->pos()))));
        //        }
        //        rebuildStrokeScene();
        //        update();
        //    } else {
        //        BaseClass::mousePressEvent(e);
        //    }
        //}

        //void PIAnnotationWidget::mouseMoveEvent(QMouseEvent * e) {
        //    if (!_anno) {
        //        return;
        //    }
        //    QVector3D t(e->pos() - _lastPos);
        //    t.setX(-t.x());
        //    if (e->buttons() & Qt::MidButton) {
        //        _options.camera().moveCenterWithEyeFixed(MakeCoreVec(t));
        //        setCursor(Qt::ClosedHandCursor);
        //        update();
        //    } else {
        //        BaseClass::mouseMoveEvent(e);
        //    }
        //    _lastPos = e->pos();
        //}

        //void PIAnnotationWidget::mouseReleaseEvent(QMouseEvent * e) {
        //    if (!_anno) {
        //        return;
        //    }
        //    unsetCursor();
        //    BaseClass::mouseReleaseEvent(e);
        //}

        //void PIAnnotationWidget::wheelEvent(QWheelEvent * e) {
        //    if (!_anno) {
        //        return;
        //    }
        //    _options.camera().setFocal(_options.camera().focal() * exp(e->delta() / 1000.0));
        //    update();
        //    BaseClass::wheelEvent(e);
        //}

        //void PIAnnotationWidget::keyPressEvent(QKeyEvent * e) {
        //    if (!_anno) {
        //        return;
        //    }
        //    if (e->key() == Qt::Key_Space) {
        //        if (_state == CreatingOcclusion) {
        //            acceptAsOcclusion();
        //            rebuildOcclusionScenes();
        //        } else if (_state == CreatingPolygon) {
        //            acceptAsPolygon(_segControl.orientationClaz, _segControl.orientationNotClaz, _segControl.used);
        //            rebuildPolygonScenes();
        //        } else if (_state == CreatingLine) {
        //            acceptAsLines();
        //            rebuildLinesScene();
        //        }
        //        clearStroke();
        //        update();
        //    } else if (e->key() == Qt::Key_Escape) {
        //        clearStroke();
        //        update();
        //    }
        //    BaseClass::keyPressEvent(e);
        //}


        //void PIAnnotationWidget::rebuildLinesScene() {
        //    SceneBuilder sb;
        //    sb.installingOptions().lineWidth = 3.0;
        //    sb.installingOptions().defaultShaderSource = gui::OpenGLShaderSourceDescriptor::XLines;
        //    ColorTable rgb = gui::RGB;
        //    rgb.exceptionalColor() = gui::Black;
        //    for (int i = 0; i < _anno->lines.size(); i++) {
        //        sb.add(ColorAs(normalize(_anno->lines[i].component) * 0.15, rgb[_anno->lines[i].claz]));
        //    }
        //    _linesScene = sb.scene();
        //}

        //void PIAnnotationWidget::rebuildPolygonScenes() {
        //    _polygonScenes.clear();
        //    // polygons       
        //    for (int i = 0; i < _anno->polygons.size(); i++) {
        //        SceneBuilder sb;
        //        sb.installingOptions().defaultShaderSource = gui::OpenGLShaderSourceDescriptor::XTriangles;
        //        auto poly = _anno->polygons[i];
        //        for (auto & c : poly.polygon.corners) {
        //            c = normalize(c);
        //        }
        //        Color color;

        //        ColorTable rgb = gui::RGB;
        //        auto & p = _anno->polygons[i];
        //        if (p.control.orientationClaz != -1) {
        //            color = rgb[p.control.orientationClaz];
        //        } else if (p.control.orientationNotClaz != -1) {
        //            color = rgb[p.control.orientationNotClaz].blendWith(gui::White, 0.5);
        //        } else if (!p.control.used) {
        //            color = gui::Black;
        //        } else {
        //            color = gui::White;
        //        }
        //        sb.add(ColorAs(poly.polygon, color));
        //        _polygonScenes.push_back(sb.scene());
        //    }
        //}

        //void PIAnnotationWidget::rebuildOcclusionScenes() {
        //    _occlusionScenes.clear();
        //    // occlusion                   
        //    for (int i = 0; i < _anno->occlusions.size(); i++) {
        //        SceneBuilder sb;
        //        sb.installingOptions().lineWidth = 3.0;
        //        sb.installingOptions().pointSize = 10.0;
        //        sb.installingOptions().discretizeOptions.color = gui::Cyan;

        //        auto & chain = _anno->occlusions[i].chain;
        //        for (int i = 1; i < chain.size(); i++) {
        //            auto p1 = normalize(chain[i - 1]) / 2.0;
        //            auto p2 = normalize(chain[i]) / 2.0;
        //            sb.begin(Line3(p1, p2)).shaderSource(gui::OpenGLShaderSourceDescriptor::XLines).end();
        //            sb.begin(p2).shaderSource(gui::OpenGLShaderSourceDescriptor::XPoints).end();
        //        }
        //        _occlusionScenes.push_back(sb.scene());
        //    }
        //}

        //void PIAnnotationWidget::rebuildStrokeScene() {
        //    if (_chain.empty()) {
        //        _strokeScene = gui::Scene();
        //        return;
        //    }

        //    SceneBuilder sb;
        //    sb.installingOptions().lineWidth = 5.0;
        //    sb.installingOptions().pointSize = 7.0;
        //    sb.installingOptions().discretizeOptions.color = gui::Orange;

        //    auto & chain = _chain;
        //    for (int i = 1; i < chain.size(); i++) {
        //        auto p1 = normalize(chain[i - 1]) / 2.0;
        //        auto p2 = normalize(chain[i]) / 2.0;
        //        sb.begin(Line3(p1, p2)).shaderSource(gui::OpenGLShaderSourceDescriptor::XLines).end();
        //        sb.begin(p2).shaderSource(gui::OpenGLShaderSourceDescriptor::XPoints).end();
        //    }

        //    _strokeScene = sb.scene();
        //}

        //void PIAnnotationWidget::clearStroke() {
        //    _chain.clear();
        //    _strokeScene = gui::Scene();
        //}

        //void PIAnnotationWidget::acceptAsPolygon(int towardVPId, int alongVPId, bool used) {
        //    if (_chain.empty()) {
        //        return;
        //    }
        //    AnnotedPolygon poly;
        //    poly.control.orientationClaz = towardVPId;
        //    poly.control.orientationNotClaz = alongVPId;
        //    poly.control.used = used;
        //    poly.polygon = _chain;
        //    for (auto & c : poly.polygon.corners) {
        //        c = normalize(c) * (poly.control.used ? 0.5 : 0.4);
        //    }
        //    _anno->polygons.push_back(std::move(poly));
        //    _polygonsDeleted.push_back(false);
        //}

        //void PIAnnotationWidget::acceptAsOcclusion() {
        //    if (_chain.empty()) {
        //        return;
        //    }
        //    _anno->occlusions.push_back({ std::move(_chain) });
        //    _occlusionsDeleted.push_back(false);
        //}

        //void PIAnnotationWidget::acceptAsLines() {
        //    if (_chain.empty()) {
        //        return;
        //    }
        //    for (int i = 1; i < _chain.size(); i++) {
        //        _anno->lines.push_back(ClassifyAs(_chain.edge(i), -1));
        //    }
        //}




        static const double visualDepthImage = 100.0;
        static const double visualDepthFace = 9.0;
        static const double visualDepthBorder = 1.0;
        static const double visualDepthClutter = 0.7;
        static const double visualDepthVP = 0.6;

        static const double visualBorderWidth = 5;
        static const double visualCornerRadius = 10.0;

        PILayoutAnnotationWidget::PILayoutAnnotationWidget(QWidget * parent /*= nullptr*/) 
            : QGLWidget(parent), _anno(nullptr), _state(Idle) {

            _showLayouts = _showVPs = true;
            _cornerClicked = _borderClicked = _faceClicked = -1;

            setMouseTracking(true);
            setAutoBufferSwap(false);
            grabKeyboard();

            _stroke.closed = false;
            setContextMenuPolicy(Qt::ActionsContextMenu);


            QAction * defaultAction = nullptr;           

            // set mode
            {
                QActionGroup * bas = new QActionGroup(this);
                connect(defaultAction = bas->addAction(tr("Pick Mode")), &QAction::triggered, [this]() {
                    _state = Idle;
                    clearStroke();
                    update();
                });
                connect(defaultAction = bas->addAction(tr("Draw Border")), &QAction::triggered, [this]() {
                    _state = DrawingBorder;
                    clearStroke();
                    update();
                });
                connect(defaultAction = bas->addAction(tr("Draw Clutter")), &QAction::triggered, [this]() {
                    _state = DrawingClutter;
                    clearStroke();
                    update();
                });
                for (auto a : bas->actions()) a->setCheckable(true);
                bas->setExclusive(true);
                addActions(bas->actions());
            }


            {
                QAction * sep = new QAction(this);
                sep->setSeparator(true);
                addAction(sep);
            }

            // set face orientations
            connect(defaultAction = new QAction(tr("Make Face Toward Red VP"), this), &QAction::triggered, [this]() {
                if (_faceClicked != -1) {
                    _anno->face2control[_faceClicked] = { 0, -1, true };
                    rebuildLayoutScene();
                }                
                update();
            });
            addAction(defaultAction);
            connect(defaultAction = new QAction(tr("Make Face Toward Green VP"), this), &QAction::triggered, [this]() {
                if (_faceClicked != -1) {
                    _anno->face2control[_faceClicked] = { 1, -1, true };
                    rebuildLayoutScene();
                }
                update();
            });
            addAction(defaultAction);
            connect(defaultAction = new QAction(tr("Make Face Toward Blue VP"), this), &QAction::triggered, [this]() {
                if (_faceClicked != -1) {
                    _anno->face2control[_faceClicked] = { 2, -1, true };
                    rebuildLayoutScene();
                }
                update();
            });
            addAction(defaultAction);
            connect(defaultAction = new QAction(tr("Make Face Along Red VP"), this), &QAction::triggered, [this]() {
                if (_faceClicked != -1) {
                    _anno->face2control[_faceClicked] = { -1, 0, true };
                    rebuildLayoutScene();
                }
                update();
            });
            addAction(defaultAction);
            connect(defaultAction = new QAction(tr("Make Face Along Green VP"), this), &QAction::triggered, [this]() {
                if (_faceClicked != -1) {
                    _anno->face2control[_faceClicked] = { -1, 1, true };
                    rebuildLayoutScene();
                }
                update();
            });
            addAction(defaultAction);
            connect(defaultAction = new QAction(tr("Make Face Along Blue VP"), this), &QAction::triggered, [this]() {
                if (_faceClicked != -1) {
                    _anno->face2control[_faceClicked] = { -1, 2, true };
                    rebuildLayoutScene();
                }
                update();
            });
            addAction(defaultAction);

            {
                QAction * sep = new QAction(this);
                sep->setSeparator(true);
                addAction(sep);
            }

            // set border occlusion
            connect(defaultAction = new QAction(tr("Make Border Disconnected"), this), &QAction::triggered, [this]() {
                if (_borderClicked != -1) {
                    _anno->border2connected[_borderClicked] = false;
                    rebuildLayoutScene();
                }
                update();
            });
            addAction(defaultAction);
            connect(defaultAction = new QAction(tr("Make Border Connected"), this), &QAction::triggered, [this]() {
                if (_borderClicked != -1) {
                    _anno->border2connected[_borderClicked] = true;
                    rebuildLayoutScene();
                }
                update();
            });
            addAction(defaultAction);


            // visibility
            connect(defaultAction = new QAction(tr("Show Layouts"), this), &QAction::triggered, [this](bool checked) {
                _showLayouts = checked;
                update();
            });
            defaultAction->setCheckable(true);
            defaultAction->setChecked(_showLayouts);
            connect(defaultAction = new QAction(tr("Show VPs"), this), &QAction::triggered, [this](bool checked) {
                _showVPs = checked;
                update();
            });
            defaultAction->setCheckable(true);
            defaultAction->setChecked(_showVPs);
            addAction(defaultAction);

            {
                QAction * sep = new QAction(this);
                sep->setSeparator(true);
                addAction(sep);
            }

            // settings
            connect(defaultAction = new QAction(tr("Settings"), this), &QAction::triggered, [this]() {
                PopUpGui(_options, this);
                update();
            });
            addAction(defaultAction);
            connect(defaultAction = new QAction(tr("Rebuild Faces"), this), &QAction::triggered, [this]() {
                _anno->generateFacesWithBorders();
                rebuildLayoutScene();
                clearStroke();
                update();
            });
            addAction(defaultAction);

        }

        PILayoutAnnotationWidget::~PILayoutAnnotationWidget() {

        }

        void PILayoutAnnotationWidget::setCurAnnotation(PILayoutAnnotation * anno) {
            _state = Idle;
            _anno = anno;
            _lastHitCornerId = -1;

            _cornerClicked = _borderClicked = _faceClicked = -1;

            assert(!anno->view.image.empty());
            auto im = anno->view.image;

            // build image scene
            SceneBuilder sb;
            ResourceStore::set("tex", anno->rectifiedImage);
            Sphere3 sp;
            sp.center = Origin();
            sp.radius = visualDepthImage;
            sb.begin(sp).shaderSource(OpenGLShaderSourceDescriptor::XPanorama).resource("tex").end();
            _imageScene = sb.scene();


            // build vps scene
            sb.clear();
            gui::ColorTable rgb = gui::RGBGreys;
            for (int i = 0; i < anno->vps.size(); i++) {
                sb.installingOptions().discretizeOptions.color = rgb[i].blendWith(gui::White, 0.4);
                sb.begin(normalize(anno->vps[i]) * visualDepthVP).shaderSource(OpenGLShaderSourceDescriptor::XPoints).pointSize(40).end();
                sb.begin(-normalize(anno->vps[i]) * visualDepthVP).shaderSource(OpenGLShaderSourceDescriptor::XPoints).pointSize(40).end();
            }
            _vpsScene = sb.scene();

            _options.panoramaAspectRatio(im.rows / float(im.cols));
            _options.panoramaHoriCenterRatio(0.5f);
            _options.camera(PerspectiveCamera(500, 500, Point2(250, 250), 200, Origin(), X(), -Z()));
            _options.renderMode(gui::RenderModeFlag::All);
            
            rebuildLayoutScene();
            rebuildCluttersScene();

            update();
        }

        void PILayoutAnnotationWidget::paintEvent(QPaintEvent * e) {
            if (!_anno) {
                return;
            }
            _imageScene.initialize();
            _layoutScene.initialize();
            _cluttersScene.initialize();
            _vpsScene.initialize();
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
            if (_showLayouts) {
                _layoutScene.render(_options);
            }
            if (_showClutters) {
                _cluttersScene.render(_options);
            }
            if (_showVPs) {
                _vpsScene.render(_options);
            }
            _strokeScene.render(_options);

            painter.endNativePainting();
            swapBuffers();
        }


        void PILayoutAnnotationWidget::mousePressEvent(QMouseEvent * e) {
            if (!_anno) {
                return;
            }
            if (e->buttons() & Qt::MidButton) {
                _lastPos = e->pos();
                setCursor(Qt::OpenHandCursor);
            } else if ((e->buttons() & Qt::LeftButton) && _showLayouts) {
                Vec3 dir = normalize(_options.camera().toSpace(MakeCoreVec(e->pos())));
                Ray3 ray(_options.camera().eye(), dir);
                
                _cornerClicked = _borderClicked = _faceClicked = -1;

                // select layout components and invoke callback functions
                auto ents = _layoutScene.pickOnScreen(_options, core::Point2(e->pos().x(), e->pos().y()));
                _layoutScene.invokeCallbackFunctions(InteractionID::ClickLeftButton, ents, false);
                
                if (_state == Idle) {

                    // select something
                    if (e->modifiers() & Qt::ControlModifier) {
                        for (auto & ent : ents) {
                            _layoutScene.switchSelect(ent);
                        }
                    } else {
                        _layoutScene.clearSelection();
                        for (auto & ent : ents) {
                            _layoutScene.select(ent);
                        }
                    }

                }else if (_state == DrawingBorder) {
                    auto & corners = _anno->corners;               
                    if (_cornerClicked != -1) { // yes, click on a corner
                        if (_lastHitCornerId != -1) {
                            // connect a border!
                            _anno->addBorder(_lastHitCornerId, _cornerClicked);
                            _lastHitCornerId = _cornerClicked;
                        } else {
                            _lastHitCornerId = _cornerClicked;
                        }
                    } else if (_borderClicked != -1) { // yes, click on a border
                        // make a new corner here
                        // break this border into two
                        _anno->corners.push_back(dir);
                        int newCorner = _anno->corners.size() - 1;
                        _anno->splitBorderBy(_borderClicked, newCorner);
                        if (_lastHitCornerId != -1) {
                            // connect a border!
                            _anno->addBorder(_lastHitCornerId, newCorner);
                            _lastHitCornerId = newCorner;
                        } else {
                            _lastHitCornerId = newCorner;
                        }
                    } else { // a new corner not on any border
                        _anno->corners.push_back(dir);
                        int newCorner = _anno->corners.size() - 1;
                        if (_lastHitCornerId != -1) { 
                            // connect a border!
                            _anno->addBorder(_lastHitCornerId, newCorner);
                            _lastHitCornerId = newCorner;
                        } else {
                            _lastHitCornerId = newCorner;
                        }
                    }
                    rebuildLayoutScene();
                } else if (_state == DrawingClutter) {

                    _stroke.points.push_back(normalize(_options.camera().toSpace(gui::MakeCoreVec(e->pos()))));
                
                }
                rebuildStrokeScene();
                update();
            } else {
                QGLWidget::mousePressEvent(e);
            }
        }

        void PILayoutAnnotationWidget::mouseMoveEvent(QMouseEvent * e) {
            if (!_anno) {
                return;
            }
            QVector3D t(e->pos() - _lastPos);
            t.setX(-t.x());
            if (e->buttons() & Qt::MidButton) {
                _options.camera().moveCenterWithEyeFixed(MakeCoreVec(t));
                setCursor(Qt::ClosedHandCursor);
                update();
            } else if (e->buttons() & Qt::LeftButton) {
                if (_cornerClicked != -1) {
                    _anno->corners[_cornerClicked] = normalize(_options.camera().toSpace(gui::MakeCoreVec(e->pos())));
                    rebuildLayoutScene();
                    update();
                }            
            } else {
                QGLWidget::mouseMoveEvent(e);
            }
            _lastPos = e->pos();
        }

        void PILayoutAnnotationWidget::mouseReleaseEvent(QMouseEvent * e) {
            if (!_anno) {
                return;
            }
            unsetCursor();
            QGLWidget::mouseReleaseEvent(e);
        }

        void PILayoutAnnotationWidget::wheelEvent(QWheelEvent * e) {
            if (!_anno) {
                return;
            }
            _options.camera().setFocal(_options.camera().focal() * exp(e->delta() / 1000.0));
            update();
            QGLWidget::wheelEvent(e);
        }

        void PILayoutAnnotationWidget::keyPressEvent(QKeyEvent * e) {
            if (!_anno) {
                return;
            }
            if (e->key() == Qt::Key_Space) {
                if (_state == DrawingClutter) {
                    acceptClutter();
                    rebuildCluttersScene();
                }
                clearStroke();
                update();
            } else if (e->key() == Qt::Key_Escape) {
                _state = Idle;
                clearStroke();
                update();
            }
            QGLWidget::keyPressEvent(e);
        }

        void PILayoutAnnotationWidget::rebuildLayoutScene() {
            SceneBuilder sb;
            sb.installingOptions().lineWidth = 10.0;
            sb.installingOptions().pointSize = 30.0;
            sb.installingOptions().discretizeOptions.color = gui::Cyan;
            
            auto & corners = _anno->corners;
            
            // corners
            _cornerPoints.resize(corners.size());
            for (int i = 0; i < corners.size(); i++) {
                _cornerPoints[i].component.component = normalize(corners[i]) * visualDepthBorder;
                _cornerPoints[i].component.color = gui::White;
                _cornerPoints[i].decoration = i;
            }
            sb.begin(_cornerPoints, [this](gui::InteractionID iid, const core::Decorated<Colored<Point3>, int> & point) {
                std::cout << "corner " << point.decoration << " is clicked" << std::endl;
                _cornerClicked = point.decoration;
                _borderClicked = -1;
                _faceClicked = -1;
            }).shaderSource(gui::OpenGLShaderSourceDescriptor::XPoints).end();

            // borders
            _borderLines.clear();
            _borderLines.resize(_anno->border2corners.size());
            for (int i = 0; i < _anno->border2corners.size(); i++) {
                _borderLines[i].decoration = i;
                auto & cs = _anno->border2corners[i];
                _borderLines[i].component.component.first = normalize(corners[cs.first]) * visualDepthBorder;
                _borderLines[i].component.component.second = normalize(corners[cs.second]) * visualDepthBorder;
                _borderLines[i].component.color = _anno->border2connected[i] ? gui::White : gui::Gray;
            }
            sb.begin(_borderLines, [this](gui::InteractionID iid, const core::Decorated<Colored<Line3>, int> & line) {
                std::cout << "border " << line.decoration << " is clicked" << std::endl;
                _cornerClicked = -1;
                _borderClicked = line.decoration;
                _faceClicked = -1;
            }).shaderSource(gui::OpenGLShaderSourceDescriptor::XLines).end();

            // faces
            _facePolygons.clear();
            // rebuild faces
            _facePolygons.resize(_anno->face2control.size());
            const static gui::ColorTable rgbtable = gui::RGB;
            for (int i = 0; i < _anno->face2control.size(); i++) {
                _facePolygons[i].decoration = i;
                const auto & cs = _anno->face2corners[i];
                auto & control = _anno->face2control[i];
                _facePolygons[i].component.color = control.dof() == 3 
                    ? gui::White 
                    : (control.dof() == 2 
                        ? rgbtable[control.orientationNotClaz].blendWith(gui::Blue, 0.2) 
                        : rgbtable[control.orientationClaz]);
                _facePolygons[i].component.component.normal = Origin();
                for (int corner : cs) {
                    _facePolygons[i].component.component.corners.push_back(normalize(corners[corner]) * visualDepthFace);
                    _facePolygons[i].component.component.normal += _facePolygons[i].component.component.corners.back();
                }
                _facePolygons[i].component.component.normal /= norm(_facePolygons[i].component.component.normal);
            }
            sb.begin(_facePolygons, [this](gui::InteractionID iid, const core::Decorated<Colored<Polygon3>, int> & polygon) {
                std::cout << "face " << polygon.decoration << " is clicked" << std::endl;
                _cornerClicked = -1;
                _borderClicked = -1;
                _faceClicked = polygon.decoration;
            }).shaderSource(gui::OpenGLShaderSourceDescriptor::XPanorama).end();

            _layoutScene = sb.scene();
        }


        void PILayoutAnnotationWidget::rebuildCluttersScene() {

        }

        void PILayoutAnnotationWidget::rebuildStrokeScene() {

        }



        void PILayoutAnnotationWidget::acceptClutter() {

        }


        void PILayoutAnnotationWidget::clearStroke() {
            _lastHitCornerId = -1;
            _stroke.clear();
            _strokeScene = gui::Scene();
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
