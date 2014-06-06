#include <QtGui>
#include <QtOpenGL>

#include <Ogre.h>
#include <OgreFrameListener.h>

#include <iostream>

#include "gtest/gtest.h"


class OgreWidget : public QGLWidget {
public:
    OgreWidget(QWidget * parent = nullptr) : QGLWidget(parent) {
        // Create an instance of the OGRE Root Class
        _root = new Ogre::Root;

        // Configures the application
        //if (!root->restoreConfig())
        _root->showConfigDialog();
        _root->saveConfig();
        _root->initialise(false);
    }

    virtual ~OgreWidget() {
        _root->shutdown();
        delete _root;
        destroy();
    }

    virtual void initializeGL() {

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
        _sceneMgr->setAmbientLight(Ogre::ColourValue(0, 0, 0));

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
            light1->setDiffuseColour(1.0f, 0, 1.0f);
            // Set Light Reflective Color
            light1->setSpecularColour(1.0f, 0.0f, 0.0f);
            // Set Light (Range, Brightness, Fade Speed, Rapid Fade Speed)
            light1->setAttenuation(10, 0.5, 0.045, 0.0);

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
            light2->setSpecularColour(1.0f, 0.0f, 0.0f);
            // Set Light (Range, Brightness, Fade Speed, Rapid Fade Speed)
            light2->setAttenuation(100, 0.05, 0.045, 0.0);

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
        _viewport->setBackgroundColour(Ogre::ColourValue(0.0, 0.0, 0.0));

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
    }


    virtual void resizeGL(int width, int height) {
        assert(_window);
        _window->reposition(this->pos().x(),
            this->pos().y());
        _window->resize(width, height);
        paintGL();
    }

    virtual void paintGL() {
        // Be sure to call "OgreWidget->repaint();" to call paintGL
        //swapBuffers();
        assert(_window);
        _root->renderOneFrame();
    }

    Ogre::Root * _root;
    Ogre::RenderWindow *_window;
    Ogre::Camera *_camera;
    Ogre::Viewport *_viewport;
    Ogre::SceneManager *_sceneMgr;
};



TEST(Ogre, Basic) {

    using namespace Ogre;

    //
    class MyFrameListener : public FrameListener {
    public:
        bool frameStarted(const FrameEvent &evt) { return true; }
        bool frameEnded(const FrameEvent &evt) { return true; }
        bool frameRenderingQueued(const FrameEvent &evt) { return true; }
    };

    // Create an instance of the OGRE Root Class
    Root* root = new Root;

    // Configures the application
    //if (!root->restoreConfig())
    root->showConfigDialog();
    root->saveConfig();

    // Create a render window
    RenderWindow* window = root->initialise(true, "Tutorial 1");
    //window->setFullscreen(false, 400, 400);

    // Create a new scene manager.
    SceneManager* sceneManager = root->createSceneManager(ST_GENERIC);
    sceneManager->setAmbientLight(Ogre::ColourValue(0.0, 0.0, 0.0));

    // Create a new camera
    Camera* camera = sceneManager->createCamera("Camera");
    camera->setPosition(Ogre::Vector3(30, 30, 30));
    camera->lookAt(Ogre::Vector3(0, 0, 0));
    camera->setNearClipDistance(5);

    // Add our model to our resources and index it
    ResourceGroupManager::getSingleton().addResourceLocation("E:/Tools/ogre/v1-9/media/packs/SinBad.zip", "Zip");
    ResourceGroupManager::getSingleton().addResourceLocation("E:/Tools/ogre/v1-9/media/packs/dragon.zip", "Zip");
    ResourceGroupManager::getSingleton().addResourceLocation("E:/Tools/ogre/v1-9/media/models/", "FileSystem");
    ResourceGroupManager::getSingleton().initialiseAllResourceGroups();

    //    
    {
        Light* light1 = sceneManager->createLight("Light1");
        light1->setType(Ogre::Light::LT_POINT);
        // Set Light Color
        light1->setDiffuseColour(1.0f, 0, 1.0f);
        // Set Light Reflective Color
        light1->setSpecularColour(1.0f, 0.0f, 0.0f);
        // Set Light (Range, Brightness, Fade Speed, Rapid Fade Speed)
        light1->setAttenuation(10, 0.5, 0.045, 0.0);

        //
        Entity* lightEnt = sceneManager->createEntity("LightEntity", "sphere.mesh");
        SceneNode* lightNode = sceneManager->createSceneNode("LightNode");
        lightNode->attachObject(lightEnt);
        lightNode->attachObject(light1);
        lightNode->setScale(0.01f, 0.01f, 0.01f);
        lightNode->setPosition(0, 4, 10);
        sceneManager->getRootSceneNode()->addChild(lightNode);
    }
    {
        Light* light2 = sceneManager->createLight("Light2");
        light2->setType(Ogre::Light::LT_POINT);
        // Set Light Color
        light2->setDiffuseColour(1.0f, 0, 1.0f);
        // Set Light Reflective Color
        light2->setSpecularColour(1.0f, 0.0f, 0.0f);
        // Set Light (Range, Brightness, Fade Speed, Rapid Fade Speed)
        light2->setAttenuation(100, 0.05, 0.045, 0.0);

        Entity* lightEnt = sceneManager->createEntity("LightEntity2", "sphere.mesh");
        SceneNode* lightNode = sceneManager->createSceneNode("LightNode2");
        lightNode->attachObject(lightEnt);
        lightNode->attachObject(light2);
        lightNode->setScale(0.01f, 0.01f, 0.01f);
        lightNode->setPosition(-10, 10, 0);
        sceneManager->getRootSceneNode()->addChild(lightNode);
    }


    // Using the camera create a viewport and set the background color to black
    Viewport* viewport = window->addViewport(camera);
    viewport->setBackgroundColour(Ogre::ColourValue(0.0, 0.0, 0.0));

    // Use the viewport to set the aspect ratio of the camera
    camera->setAspectRatio(Real(viewport->getActualWidth()) /
        Real(viewport->getActualHeight()));

    // Create an instance of our model and add it to the scene
    {
        Entity* ent = sceneManager->createEntity("SinBad.mesh");
        SceneNode* entNode = sceneManager->createSceneNode("Character");
        entNode->attachObject(ent);
        sceneManager->getRootSceneNode()->addChild(entNode);
        entNode->setPosition(0, 0, 0);
    }
    {
        Entity* ent = sceneManager->createEntity("dragon.mesh");

        auto r = ent->getBoundingRadius();
        auto c = ent->getBoundingBox().getCenter();

        SceneNode* entNode = sceneManager->createSceneNode("Dragon");
        entNode->attachObject(ent);
        entNode->scale(0.05, 0.05, 0.05);
        entNode->rotate(Ogre::Vector3(1, -1, 1), Ogre::Radian(M_PI_2));

        sceneManager->getRootSceneNode()->addChild(entNode);
        entNode->setPosition(0, 0, 0);
    }

    // Create an instance of the MyFrameListener Class and add it to the root object
    MyFrameListener* myListener = new MyFrameListener();
    root->addFrameListener(myListener);

    // Tell root to start rendering
    root->startRendering();

    // Cleanup
    delete myListener;
    delete root;

}


TEST(Ogre, Qt) {
    char * name = "app";
    int i = 1;
    QApplication app(i, &name);

    QWidget window;

    window.resize(800, 600);
    window.setWindowTitle("Simple example");

    OgreWidget* ogreWidget = new OgreWidget;

    QVBoxLayout *layout = new QVBoxLayout;
    layout->setMargin(0);
    layout->addWidget(ogreWidget);

    window.setLayout(layout);
    window.show();

    app.exec();
}


int main(int argc, char * argv[], char * envp[]) {
    testing::InitGoogleTest(&argc, argv);
    testing::FLAGS_gtest_filter = "*Qt";
    return RUN_ALL_TESTS();
}