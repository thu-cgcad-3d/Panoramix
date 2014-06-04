#include <Ogre.h>
#include <OgreFrameListener.h>

//#include <QGLWidget>
//#include <QtGui>
//#include <QtWidgets>
//#include <QApplication>


#include <iostream>




//namespace {
//
//    class OgreWidget : public QGLWidget
//    {
//        //Q_OBJECT;
//
//    public:
//        OgreWidget(QWidget *parent = 0) :
//            QGLWidget(parent),
//            mOgreWindow(NULL), mViewport(nullptr)
//        {
//            init();
//        }
//
//        virtual ~OgreWidget()
//        {
//            mOgreRoot->shutdown();
//            delete mOgreRoot;
//            destroy();
//        }
//
//    protected:
//        virtual void initializeGL() {
//            //== Creating and Acquiring Ogre Window ==//
//            // Get the parameters of the window QT created
//            Ogre::String winHandle;
//            winHandle += Ogre::StringConverter::toString((unsigned long)(this->parentWidget()->winId()));
//
//            Ogre::NameValuePairList params;
//
//            // code for Windows and Linux
//            params["parentWindowHandle"] = winHandle;
//            mOgreWindow = mOgreRoot->createRenderWindow("QOgreWidget_RenderWindow",
//                this->width(),
//                this->height(),
//                false,
//                &params);
//
//            mOgreWindow->setActive(true);
//            WId ogreWinId = 0x0;
//            mOgreWindow->getCustomAttribute("WINDOW", &ogreWinId);
//            assert(ogreWinId);
//
//            // bug fix, extract geometry
//            QRect geo = this->frameGeometry();
//
//            // create new window
//            this->create(ogreWinId);
//
//            // set geometrie infos to new window
//            this->setGeometry(geo);
//
//            setAutoBufferSwap(false);
//            setAttribute(Qt::WA_PaintOnScreen, true);
//            setAttribute(Qt::WA_NoBackground);
//
//            //== Ogre Initialization ==//
//            Ogre::SceneType scene_manager_type = Ogre::ST_EXTERIOR_CLOSE;
//
//            mSceneMgr = mOgreRoot->createSceneManager(scene_manager_type);
//            mSceneMgr->setAmbientLight(Ogre::ColourValue(1, 1, 1));
//
//            mCamera = mSceneMgr->createCamera("QOgreWidget_Cam");
//            mCamera->setPosition(Ogre::Vector3(0, 10, 0));
//            mCamera->lookAt(Ogre::Vector3(0, 0, 0));
//            mCamera->setNearClipDistance(20.0);
//
//
//
//            // setup resource
//            Ogre::ConfigFile cf;
//            cf.load("resources.cfg");
//            // Go through all sections & settings in the file
//            Ogre::ConfigFile::SectionIterator seci = cf.getSectionIterator();
//
//            Ogre::String secName, typeName, archName;
//            while (seci.hasMoreElements())
//            {
//                secName = seci.peekNextKey();
//                Ogre::ConfigFile::SettingsMultiMap *settings = seci.getNext();
//                Ogre::ConfigFile::SettingsMultiMap::iterator i;
//                for (i = settings->begin(); i != settings->end(); ++i)
//                {
//                    typeName = i->first;
//                    archName = i->second;
//                    Ogre::ResourceGroupManager::getSingleton().addResourceLocation(
//                        archName, typeName, secName);
//                }
//            }
//
//            // Set default mipmap level (note: some APIs ignore this)
//            Ogre::TextureManager::getSingleton().setDefaultNumMipmaps(5);
//            // initialise all resource groups
//            Ogre::ResourceGroupManager::getSingleton().initialiseAllResourceGroups();
//
//            auto sg = Ogre::ResourceGroupManager::getSingleton().getResourceGroups();
//            for (auto s : sg){
//                std::cout << s << std::endl;
//            }
//
//            {
//                Ogre::Entity* ogreHead = mSceneMgr->createEntity("Head", "ogrehead.mesh", "Popular");
//                Ogre::SceneNode* headNode = mSceneMgr->getRootSceneNode()->createChildSceneNode("HeadNode");
//                headNode->attachObject(ogreHead);
//
//                headNode->yaw(Ogre::Degree(-90));
//
//                Ogre::Entity* ogreHead2 = mSceneMgr->createEntity("Head2", "ogrehead.mesh", "Popular");
//                Ogre::SceneNode* headNode2 = mSceneMgr->getRootSceneNode()->createChildSceneNode("HeadNode2", Ogre::Vector3(100, 0, 0));
//                headNode2->attachObject(ogreHead2);
//
//                headNode2->pitch(Ogre::Degree(-90));
//                headNode2->scale(1, 2, 1);
//
//                Ogre::Entity* ogreHead3 = mSceneMgr->createEntity("Head3", "ogrehead.mesh", "Popular");
//                Ogre::SceneNode* headNode3 = mSceneMgr->getRootSceneNode()->createChildSceneNode("HeadNode3", Ogre::Vector3(200, 0, 0));
//                headNode3->attachObject(ogreHead3);
//
//                headNode3->roll(Ogre::Degree(-90));
//            }
//
//            mViewport = mOgreWindow->addViewport(mCamera);
//            mViewport->setBackgroundColour(Ogre::ColourValue(0.8, 0.8, 1));
//
//        }
//
//        virtual void resizeGL(int width, int height){
//            assert(mOgreWindow);
//            if (mViewport)
//                mViewport->setDimensions(0, 0, width, height);
//            mOgreWindow->reposition(this->pos().x(),
//                this->pos().y());
//            mOgreWindow->resize(width, height);
//            paintGL();
//        }
//
//        virtual void paintGL() {
//            // Be sure to call "OgreWidget->repaint();" to call paintGL
//            //swapBuffers();
//            assert(mOgreWindow);
//            mOgreRoot->renderOneFrame();
//        }
//
//        void init(){
//            // create the main ogre object
//            mOgreRoot = new Ogre::Root;
//
//            if (!mOgreRoot->restoreConfig())
//                mOgreRoot->showConfigDialog();
//            mOgreRoot->saveConfig();
//
//            // setup a renderer
//            const Ogre::RenderSystemList & renderers = mOgreRoot->getAvailableRenderers();
//            assert(!renderers.empty()); // we need at least one renderer to do anything useful
//
//            Ogre::RenderSystem *renderSystem;
//            renderSystem = chooseRenderer(renderers);
//
//            assert(renderSystem); // user might pass back a null renderer, which would be bad!
//
//            mOgreRoot->setRenderSystem(renderSystem);
//
//            QString dimensions = QString("%1x%2")
//                .arg(this->width())
//                .arg(this->height());
//
//            //Ogre::String s = dimensions.toStdString();
//            renderSystem->setConfigOption("Video Mode", "640x480");
//            renderSystem->setConfigOption("VSync", "true");
//
//            // initialize without creating window
//            mOgreRoot->getRenderSystem()->setConfigOption("Full Screen", "No");
//            mOgreRoot->saveConfig();
//            mOgreRoot->initialise(false); // don't create a window
//        }
//
//        virtual Ogre::RenderSystem* chooseRenderer(const Ogre::RenderSystemList & renderers) {
//            return renderers.front();
//        }
//
//        Ogre::Root * mOgreRoot;
//        Ogre::RenderWindow *mOgreWindow;
//        Ogre::Camera *mCamera;
//        Ogre::Viewport *mViewport;
//        Ogre::SceneManager *mSceneMgr;
//    };
//
//    void OGRETest(){
//
//        char * name = "app";
//        int i = 1;
//        QApplication app(i, &name);
//
//        QWidget window;
//
//        window.resize(800, 600);
//        window.setWindowTitle("Simple example");
//
//        OgreWidget* ogreWidget = new OgreWidget;
//
//        QVBoxLayout *layout = new QVBoxLayout;
//        layout->setMargin(0);
//        layout->addWidget(ogreWidget);
//
//        window.setLayout(layout);
//        window.show();
//
//        app.exec();
//
//    }
//
//}



using namespace Ogre;

//
class MyFrameListener : public FrameListener {
public:
    bool frameStarted(const FrameEvent &evt);
    bool frameEnded(const FrameEvent &evt);
    bool frameRenderingQueued(const FrameEvent &evt);
};

//
bool MyFrameListener::frameStarted(const FrameEvent &evt) {
    std::cout << "Frame Started" << std::endl;
    return true;
}

//
bool MyFrameListener::frameEnded(const FrameEvent &evt) {
    std::cout << "Frame Ended" << std::endl;
    return true;
}

//
bool MyFrameListener::frameRenderingQueued(const FrameEvent &evt) {
    std::cout << "Frame Queued" << std::endl;
    return true;
}

int main(void)
{
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
    camera->setPosition(Ogre::Vector3(0, 0, 15));
    camera->lookAt(Ogre::Vector3(0, 0, 0));
    camera->setNearClipDistance(5);

    // Add our model to our resources and index it
    ResourceGroupManager::getSingleton().addResourceLocation("E:/Tools/ogre/v1-9/media/packs/Sinbad.zip", "Zip");
    ResourceGroupManager::getSingleton().addResourceLocation("E:/Tools/ogre/v1-9/media/models/", "FileSystem");
    ResourceGroupManager::getSingleton().initialiseAllResourceGroups();

    //
    Light* light1 = sceneManager->createLight("Light1");
    light1->setType(Ogre::Light::LT_POINT);
    // Set Light Color
    light1->setDiffuseColour(1.0f, 1.0f, 1.0f);
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


    // Using the camera create a viewport and set the background color to black
    Viewport* viewport = window->addViewport(camera);
    viewport->setBackgroundColour(Ogre::ColourValue(0.0, 0.0, 0.0));

    // Use the viewport to set the aspect ratio of the camera
    camera->setAspectRatio(Real(viewport->getActualWidth()) /
        Real(viewport->getActualHeight()));

    // Create an instance of our model and add it to the scene
    Entity* ent = sceneManager->createEntity("Sinbad.mesh");
    SceneNode* entNode = sceneManager->createSceneNode("Character");
    entNode->attachObject(ent);
    sceneManager->getRootSceneNode()->addChild(entNode);
    entNode->setPosition(0, 0, 0);

    // Create an instance of the MyFrameListener Class and add it to the root object
    MyFrameListener* myListener = new MyFrameListener();
    root->addFrameListener(myListener);

    // Tell root to start rendering
    root->startRendering();

    // Cleanup
    delete myListener;
    delete root;

    return 0;
}