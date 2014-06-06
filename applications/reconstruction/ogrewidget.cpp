

#include "ogrewidget.hpp"

using Ogre::Camera;
using Ogre::SceneManager;
using Ogre::ResourceGroupManager;
using Ogre::Light;
using Ogre::String;
using Ogre::NameValuePairList;
using Ogre::Entity;
using Ogre::SceneNode;
using Ogre::Viewport;
using Ogre::Real;

Ogre::ColourValue QgreColourFromQColor(QColor c) {
    return Ogre::ColourValue(c.redF(), c.greenF(), c.blueF(), c.alphaF());
}


OgreWidget::OgreWidget(QWidget * parent) : QGLWidget(parent) {
    // Create an instance of the OGRE Root Class
    _root = new Ogre::Root;

    // Configures the application
    if (!_root->restoreConfig()) {
        _root->showConfigDialog();
        _root->saveConfig();
    }
    _root->initialise(false);
}

OgreWidget::~OgreWidget() {
    _root->shutdown();
    delete _root;
    destroy();
}

void OgreWidget::initializeGL() {

    {
        //== Creating and Acquiring Ogre Window ==//
        // Get the parameters of the window QT created
        Ogre::String winHandle;
        winHandle += Ogre::StringConverter::toString((unsigned long)(this->parentWidget()->winId()));

        Ogre::NameValuePairList params;

        // code for Windows and Linux
        params["parentWindowHandle"] = winHandle;
        _window = _root->createRenderWindow("QOgreWidget_RenderWindow",
            this->width(),
            this->height(),
            false,
            &params);

        _window->setActive(true);
        WId ogreWinId = 0x0;
        _window->getCustomAttribute("WINDOW", &ogreWinId);
        assert(ogreWinId);

        // bug fix, extract geometry
        QRect geo = this->frameGeometry();

        // create new window
        this->create(ogreWinId);

        // set geometrie infos to new window
        this->setGeometry(geo);

        setAutoBufferSwap(false);
        setAttribute(Qt::WA_PaintOnScreen, true);
        setAttribute(Qt::WA_NoBackground);
    }

    //== Ogre Initialization ==//
    //Ogre::SceneType scene_manager_type = Ogre::ST_EXTERIOR_CLOSE;
    Ogre::SceneType scene_manager_type = Ogre::ST_GENERIC;

    _sceneMgr = _root->createSceneManager(scene_manager_type);
    _sceneMgr->setAmbientLight(QgreColourFromQColor(Qt::gray));

    // Create a new _camera
    _camera = _sceneMgr->createCamera("Camera");
    _camera->setPosition(Ogre::Vector3(30, 30, 30));
    _camera->lookAt(Ogre::Vector3(0, 0, 0));
    _camera->setNearClipDistance(5);

    // Add our model to our resources and index it
    ResourceGroupManager::getSingleton().addResourceLocation("E:/Tools/ogre/v1-9/media/packs/SinBad.zip", "Zip");
    ResourceGroupManager::getSingleton().addResourceLocation("E:/Tools/ogre/v1-9/media/packs/dragon.zip", "Zip");
    ResourceGroupManager::getSingleton().addResourceLocation("E:/Tools/ogre/v1-9/media/models/", "FileSystem");
    ResourceGroupManager::getSingleton().initialiseAllResourceGroups();

    //    
    {
        Light* light1 = _sceneMgr->createLight("Light1");
        light1->setType(Ogre::Light::LT_POINT);
        // Set Light Color
        light1->setDiffuseColour(QgreColourFromQColor(Qt::white));
        // Set Light Reflective Color
        light1->setSpecularColour(QgreColourFromQColor(Qt::green));
        // Set Light (Range, Brightness, Fade Speed, Rapid Fade Speed)
        light1->setAttenuation(1000, 0.5, 0.045, 0.0);

        //
        Entity* lightEnt = _sceneMgr->createEntity("LightEntity", "sphere.mesh");
        SceneNode* lightNode = _sceneMgr->createSceneNode("LightNode");
        lightNode->attachObject(lightEnt);
        lightNode->attachObject(light1);
        lightNode->setScale(0.01f, 0.01f, 0.01f);
        lightNode->setPosition(0, 4, 10);
        _sceneMgr->getRootSceneNode()->addChild(lightNode);
    }
    {
        Light* light2 = _sceneMgr->createLight("Light2");
        light2->setType(Ogre::Light::LT_POINT);
        // Set Light Color
        light2->setDiffuseColour(1.0f, 0, 1.0f);
        // Set Light Reflective Color
        light2->setSpecularColour(QgreColourFromQColor(Qt::yellow));
        // Set Light (Range, Brightness, Fade Speed, Rapid Fade Speed)
        light2->setAttenuation(100, 0.1, 0.045, 0.0);

        Entity* lightEnt = _sceneMgr->createEntity("LightEntity2", "sphere.mesh");
        SceneNode* lightNode = _sceneMgr->createSceneNode("LightNode2");
        lightNode->attachObject(lightEnt);
        lightNode->attachObject(light2);
        lightNode->setScale(0.01f, 0.01f, 0.01f);
        lightNode->setPosition(-10, 10, 0);
        _sceneMgr->getRootSceneNode()->addChild(lightNode);
    }


    // Using the _camera create a _viewport and set the background color to black
    _viewport = _window->addViewport(_camera);
    _viewport->setBackgroundColour(QgreColourFromQColor(Qt::gray));

    // Use the _viewport to set the aspect ratio of the _camera
    _camera->setAspectRatio(Real(_viewport->getActualWidth()) /
        Real(_viewport->getActualHeight()));

    // Create an instance of our model and add it to the scene
    {
        Entity* ent = _sceneMgr->createEntity("SinBad.mesh");
        SceneNode* entNode = _sceneMgr->createSceneNode("Character");
        entNode->attachObject(ent);
        _sceneMgr->getRootSceneNode()->addChild(entNode);
        entNode->setPosition(0, 0, 0);
        _focusedNode = entNode;
    }
    {
        Entity* ent = _sceneMgr->createEntity("dragon.mesh");

        auto r = ent->getBoundingRadius();
        auto c = ent->getBoundingBox().getCenter();

        SceneNode* entNode = _sceneMgr->createSceneNode("Dragon");
        entNode->attachObject(ent);
        entNode->scale(0.05, 0.05, 0.05);
        entNode->rotate(Ogre::Vector3(1, -1, 1), Ogre::Radian(M_PI_2));

        _sceneMgr->getRootSceneNode()->addChild(entNode);
        entNode->setPosition(0, 0, 0);
    }

    // Create a plane
    {
       
    }

    assert(_focusedNode);
    _camera->lookAt(_focusedNode->getPosition());
}


void OgreWidget::resizeGL(int width, int height) {
    assert(_window);
    if (_camera && _viewport) {
        _camera->setAutoAspectRatio(true);
    }
    _window->reposition(this->pos().x(),
        this->pos().y());
    _window->resize(width, height);
    paintGL();
}

void OgreWidget::paintGL() {
    // Be sure to call "OgreWidget->repaint();" to call paintGL
    assert(_window);
    _root->renderOneFrame();
}

void OgreWidget::mousePressEvent(QMouseEvent * e) {
    _lastPos = e->pos();
    if (e->buttons() & Qt::RightButton)
        setCursor(Qt::OpenHandCursor);
    else if (e->buttons() & Qt::MidButton)
        setCursor(Qt::SizeAllCursor);
}

void OgreWidget::mouseMoveEvent(QMouseEvent * e) {
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

void OgreWidget::wheelEvent(QWheelEvent * e) {
    moveCameraCenterAndCenter(QVector3D(0, 0, e->delta() / 10));
    update();
}

void OgreWidget::mouseReleaseEvent(QMouseEvent *) {
    unsetCursor();
}

void OgreWidget::moveCameraEyeWithCenterFixed(const QVector3D & t) {
    auto & camera = *_camera;
    auto eye = camera.getPosition();
    auto center = _focusedNode->getPosition();
    auto up = camera.getUp();
    auto tt = t * (eye - center).length() * 0.002f;

    auto xv = (center - eye).crossProduct(up).normalisedCopy();
    auto yv = xv.crossProduct(center - eye).normalisedCopy();
    auto xyTrans = xv * tt.x() + yv * tt.y();
    double r = ((eye - center).length() - tt.z()) /
        (eye + xyTrans - center).length();
    eye = (eye + xyTrans - center) * r + center;
    up = yv.normalisedCopy();

    camera.setPosition(eye);
    camera.lookAt(center);
}

void OgreWidget::moveCameraCenterAndCenter(const QVector3D & t) {
    auto & camera = *_camera;
    auto eye = camera.getPosition();
    auto center = _focusedNode->getPosition();
    auto up = camera.getUp();
    auto tt = t * (eye - center).length() * 0.002f;

    auto xv = (center - eye).crossProduct(up).normalisedCopy();
    auto yv = xv.crossProduct(center - eye).normalisedCopy();
    auto zv = (center - eye).normalisedCopy();
    auto trans = xv * tt.x() + yv * tt.y() + zv * tt.z();
    eye += trans;
    center += trans;
    camera.setPosition(eye);
    camera.lookAt(center);
}


