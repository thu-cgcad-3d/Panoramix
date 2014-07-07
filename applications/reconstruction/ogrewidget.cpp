
#include "../../src/core/mesh_maker.hpp"
#include "../../src/rec/reconstruction_engine_visualize.hpp"

#include "ogrewidget.hpp"

using namespace panoramix;

using Ogre::Camera;
using Ogre::SceneManager;
using Ogre::ResourceGroupManager;
using Ogre::MeshManager;
using Ogre::MeshPtr;
using Ogre::MaterialPtr;
using Ogre::MaterialManager;
using Ogre::Light;
using Ogre::String;
using Ogre::NameValuePairList;
using Ogre::Entity;
using Ogre::SceneNode;
using Ogre::Viewport;
using Ogre::Real;
using Ogre::SubMesh;
using Ogre::RenderSystem;
using Ogre::RGBA;
using Ogre::ColourValue;
using Ogre::VertexData;
using Ogre::VertexDeclaration;
using Ogre::VertexElement;
using Ogre::Root;


Ogre::ColourValue ToQgreColour(QColor c) {
    return Ogre::ColourValue(c.redF(), c.greenF(), c.blueF(), c.alphaF());
}


OgreWidget::OgreWidget(QWidget * parent) : QGLWidget(parent) {
    // Create an instance of the OGRE Root Class
    _root = new Ogre::Root;

    // Configures the application
    //if (!_root->restoreConfig()) {
        _root->showConfigDialog();
        _root->saveConfig();
    //}
    _root->initialise(false);
}

OgreWidget::~OgreWidget() {
    _root->shutdown();
    delete _root;
    destroy();
}


void OgreWidget::setupPanorama(const QString & filename) {
    cv::Mat panorama = cv::imread(filename.toStdString());
    cv::resize(panorama, panorama, cv::Size(2000, 1000));
    core::PanoramicCamera originCam(panorama.cols / M_PI / 2.0);

    std::vector<core::PerspectiveCamera> cams;
    core::Mesh<core::Vec3> cameraStand;
    core::MakeQuadFacedSphere(cameraStand, 6, 12);
    for (auto & v : cameraStand.vertices()) {
        core::Vec3 direction = v.data;
        if (core::AngleBetweenDirections(direction, core::Vec3(0, 0, 1)) <= 0.1 ||
            core::AngleBetweenDirections(direction, core::Vec3(0, 0, -1)) <= 0.1) {
            //cams.emplace_back(700, 700, originCam.focal(), core::Vec3(0, 0, 0), direction, core::Vec3(0, 1, 0));
            continue;
        } else {
            cams.emplace_back(700, 700, originCam.focal(), core::Vec3(0, 0, 0), direction, core::Vec3(0, 0, -1));
        }
    }


    /// insert into views net
    rec::ReconstructionEngine::Params params;
    params.mjWeightT = 2.0;
    params.intersectionConstraintLineDistanceAngleThreshold = 0.05;
    params.incidenceConstraintLineDistanceAngleThreshold = 0.2;
    params.mergeLineDistanceAngleThreshold = 0.05;
    rec::ReconstructionEngine net(params);

    net.insertPanorama(panorama, cams, originCam);

    #pragma omp parallel for
    for (int i = 0; i < net.views().internalElements<0>().size(); i++) {
        auto viewHandle = rec::ReconstructionEngine::ViewHandle(i);
        net.computeFeatures(viewHandle);
        net.buildRegionNet(viewHandle);
    }

    std::cout << "calibrating camera and classifying lines ...";
    net.estimateVanishingPointsAndClassifyLines();
    net.rectifySpatialLines();
        
}

void OgreWidget::createCube(const Ogre::String & name) {
    /// Create the mesh via the MeshManager
    Ogre::MeshPtr msh = MeshManager::getSingleton().createManual(name, ResourceGroupManager::DEFAULT_RESOURCE_GROUP_NAME);

    /// Create one submesh
    SubMesh* sub = msh->createSubMesh();

    const float sqrt13 = 0.577350269f; /* sqrt(1/3) */

    /// Define the vertices (8 vertices, each consisting of 2 groups of 3 floats
    const size_t nVertices = 8;
    const size_t vbufCount = 3 * 2 * nVertices;
    float vertices[vbufCount] = {
        -100.0, 100.0, -100.0,        //0 position
        -sqrt13, sqrt13, -sqrt13,     //0 normal
        100.0, 100.0, -100.0,         //1 position
        sqrt13, sqrt13, -sqrt13,      //1 normal
        100.0, -100.0, -100.0,        //2 position
        sqrt13, -sqrt13, -sqrt13,     //2 normal
        -100.0, -100.0, -100.0,       //3 position
        -sqrt13, -sqrt13, -sqrt13,    //3 normal
        -100.0, 100.0, 100.0,         //4 position
        -sqrt13, sqrt13, sqrt13,      //4 normal
        100.0, 100.0, 100.0,          //5 position
        sqrt13, sqrt13, sqrt13,       //5 normal
        100.0, -100.0, 100.0,         //6 position
        sqrt13, -sqrt13, sqrt13,      //6 normal
        -100.0, -100.0, 100.0,        //7 position
        -sqrt13, -sqrt13, sqrt13,     //7 normal
    };

    RenderSystem* rs = Root::getSingleton().getRenderSystem();
    RGBA colours[nVertices];
    RGBA *pColour = colours;
    // Use render system to convert colour value since colour packing varies
    rs->convertColourValue(ColourValue(1.0, 0.0, 0.0), pColour++); //0 colour
    rs->convertColourValue(ColourValue(1.0, 1.0, 0.0), pColour++); //1 colour
    rs->convertColourValue(ColourValue(0.0, 1.0, 0.0), pColour++); //2 colour
    rs->convertColourValue(ColourValue(0.0, 0.0, 0.0), pColour++); //3 colour
    rs->convertColourValue(ColourValue(1.0, 0.0, 1.0), pColour++); //4 colour
    rs->convertColourValue(ColourValue(1.0, 1.0, 1.0), pColour++); //5 colour
    rs->convertColourValue(ColourValue(0.0, 1.0, 1.0), pColour++); //6 colour
    rs->convertColourValue(ColourValue(0.0, 0.0, 1.0), pColour++); //7 colour


    /// Define 12 triangles (two triangles per cube face)
    /// The values in this table refer to vertices in the above table
    const size_t ibufCount = 36;
    unsigned short faces[ibufCount] = {
        0, 2, 3,
        0, 1, 2,
        1, 6, 2,
        1, 5, 6,
        4, 6, 5,
        4, 7, 6,
        0, 7, 4,
        0, 3, 7,
        0, 5, 1,
        0, 4, 5,
        2, 7, 3,
        2, 6, 7
    };

    /// Create vertex data structure for 8 vertices shared between submeshes
    msh->sharedVertexData = new VertexData();
    msh->sharedVertexData->vertexCount = nVertices;

    /// Create declaration (memory format) of vertex data
    VertexDeclaration* decl = msh->sharedVertexData->vertexDeclaration;
    size_t offset = 0;




    // 1st buffer
    decl->addElement(0, offset, Ogre::VET_FLOAT3, Ogre::VES_POSITION);
    offset += VertexElement::getTypeSize(Ogre::VET_FLOAT3);
    decl->addElement(0, offset, Ogre::VET_FLOAT3, Ogre::VES_NORMAL);
    offset += VertexElement::getTypeSize(Ogre::VET_FLOAT3);
    /// Allocate vertex buffer of the requested number of vertices (vertexCount) 
    /// and bytes per vertex (offset)
    Ogre::HardwareVertexBufferSharedPtr vbuf =
        Ogre::HardwareBufferManager::getSingleton().createVertexBuffer(
        offset, msh->sharedVertexData->vertexCount, Ogre::HardwareBuffer::HBU_STATIC_WRITE_ONLY);
    /// Upload the vertex data to the card
    vbuf->writeData(0, vbuf->getSizeInBytes(), vertices, true);

    /// Set vertex buffer binding so buffer 0 is bound to our vertex buffer
    Ogre::VertexBufferBinding* bind = msh->sharedVertexData->vertexBufferBinding;
    bind->setBinding(0, vbuf);




    // 2nd buffer
    offset = 0;
    decl->addElement(1, offset, Ogre::VET_COLOUR, Ogre::VES_DIFFUSE);
    offset += VertexElement::getTypeSize(Ogre::VET_COLOUR);
    /// Allocate vertex buffer of the requested number of vertices (vertexCount) 
    /// and bytes per vertex (offset)
    vbuf = Ogre::HardwareBufferManager::getSingleton().createVertexBuffer(
        offset, msh->sharedVertexData->vertexCount, Ogre::HardwareBuffer::HBU_STATIC_WRITE_ONLY);
    /// Upload the vertex data to the card
    vbuf->writeData(0, vbuf->getSizeInBytes(), colours, true);

    /// Set vertex buffer binding so buffer 1 is bound to our colour buffer
    bind->setBinding(1, vbuf);



    // index buffer
    /// Allocate index buffer of the requested number of vertices (ibufCount) 
    Ogre::HardwareIndexBufferSharedPtr ibuf = Ogre::HardwareBufferManager::getSingleton().
        createIndexBuffer(
        Ogre::HardwareIndexBuffer::IT_16BIT,
        ibufCount,
        Ogre::HardwareBuffer::HBU_STATIC_WRITE_ONLY);

    /// Upload the index data to the card
    ibuf->writeData(0, ibuf->getSizeInBytes(), faces, true);

    /// Set parameters of the submesh
    sub->useSharedVertices = true;
    sub->indexData->indexBuffer = ibuf;
    sub->indexData->indexCount = ibufCount;
    sub->indexData->indexStart = 0;

    /// Set bounding information (for culling)
    msh->_setBounds(Ogre::AxisAlignedBox(-100, -100, -100, 100, 100, 100));
    msh->_setBoundingSphereRadius(Ogre::Math::Sqrt(3 * 100 * 100));

    /// Notify -Mesh object that it has been loaded
    msh->load();

}

void OgreWidget::createSphere(const Ogre::String & name, const float r, const int nRings, const int nSegments) {
    using namespace Ogre;

    MeshPtr pSphere = MeshManager::getSingleton().createManual(name, ResourceGroupManager::DEFAULT_RESOURCE_GROUP_NAME);
    SubMesh *pSphereVertex = pSphere->createSubMesh();

    pSphere->sharedVertexData = new VertexData();
    VertexData* vertexData = pSphere->sharedVertexData;

    // define the vertex format
    VertexDeclaration* vertexDecl = vertexData->vertexDeclaration;
    size_t currOffset = 0;
    // positions
    vertexDecl->addElement(0, currOffset, VET_FLOAT3, VES_POSITION);
    currOffset += VertexElement::getTypeSize(VET_FLOAT3);
    // normals
    vertexDecl->addElement(0, currOffset, VET_FLOAT3, VES_NORMAL);
    currOffset += VertexElement::getTypeSize(VET_FLOAT3);
    // two dimensional texture coordinates
    vertexDecl->addElement(0, currOffset, VET_FLOAT2, VES_TEXTURE_COORDINATES, 0);
    currOffset += VertexElement::getTypeSize(VET_FLOAT2);

    // allocate the vertex buffer
    vertexData->vertexCount = (nRings + 1) * (nSegments + 1);
    HardwareVertexBufferSharedPtr vBuf = HardwareBufferManager::getSingleton().createVertexBuffer(vertexDecl->getVertexSize(0), 
        vertexData->vertexCount, HardwareBuffer::HBU_STATIC_WRITE_ONLY, false);
    VertexBufferBinding* binding = vertexData->vertexBufferBinding;
    binding->setBinding(0, vBuf);
    float* pVertex = static_cast<float*>(vBuf->lock(HardwareBuffer::HBL_DISCARD));


    // allocate index buffer
    pSphereVertex->indexData->indexCount = 6 * nRings * (nSegments + 1);
    pSphereVertex->indexData->indexBuffer = HardwareBufferManager::getSingleton().createIndexBuffer(HardwareIndexBuffer::IT_16BIT, pSphereVertex->indexData->indexCount, HardwareBuffer::HBU_STATIC_WRITE_ONLY, false);
    HardwareIndexBufferSharedPtr iBuf = pSphereVertex->indexData->indexBuffer;
    unsigned short* pIndices = static_cast<unsigned short*>(iBuf->lock(HardwareBuffer::HBL_DISCARD));

    float fDeltaRingAngle = (Math::PI / nRings);
    float fDeltaSegAngle = (2 * Math::PI / nSegments);
    unsigned short wVerticeIndex = 0;

    // Generate the group of rings for the sphere
    for (int ring = 0; ring <= nRings; ring++) {
        float r0 = r * sinf(ring * fDeltaRingAngle);
        float y0 = r * cosf(ring * fDeltaRingAngle);

        // Generate the group of segments for the current ring
        for (int seg = 0; seg <= nSegments; seg++) {
            float x0 = r0 * sinf(seg * fDeltaSegAngle);
            float z0 = r0 * cosf(seg * fDeltaSegAngle);

            // Add one vertex to the strip which makes up the sphere
            *pVertex++ = x0;
            *pVertex++ = y0;
            *pVertex++ = z0;

            Vector3 vNormal = Vector3(x0, y0, z0).normalisedCopy();
            *pVertex++ = vNormal.x;
            *pVertex++ = vNormal.y;
            *pVertex++ = vNormal.z;

            *pVertex++ = (float)seg / (float)nSegments;
            *pVertex++ = (float)ring / (float)nRings;

            if (ring != nRings) {
                // each vertex (except the last) has six indices pointing to it
                *pIndices++ = wVerticeIndex + nSegments + 1;
                *pIndices++ = wVerticeIndex;
                *pIndices++ = wVerticeIndex + nSegments;
                *pIndices++ = wVerticeIndex + nSegments + 1;
                *pIndices++ = wVerticeIndex + 1;
                *pIndices++ = wVerticeIndex;
                wVerticeIndex++;
            }
        }; // end for seg
    } // end for ring

    // Unlock
    vBuf->unlock();
    iBuf->unlock();
    // Generate face list
    pSphereVertex->useSharedVertices = true;

    // the original code was missing this line:
    pSphere->_setBounds(AxisAlignedBox(Vector3(-r, -r, -r), Vector3(r, r, r)), false);
    pSphere->_setBoundingSphereRadius(r);
    // this line makes clear the mesh is loaded (avoids memory leaks)
    pSphere->load();
}

void OgreWidget::createCamera() {
    // Create a new _camera
    _camera = _sceneMgr->createCamera("Camera");
    _camera->setPosition(Ogre::Vector3(200, 200, 200));
    _camera->lookAt(Ogre::Vector3(0, 0, 0));
    _camera->setNearClipDistance(5);
}

void OgreWidget::createViewports() {
    // Using the _camera create a _viewport and set the background color to black
    _viewport = _window->addViewport(_camera);
    _viewport->setBackgroundColour(ToQgreColour(Qt::black));

    // Use the _viewport to set the aspect ratio of the _camera
    _camera->setAutoAspectRatio(true);
}

void OgreWidget::prepareResources() {
    ResourceGroupManager::getSingleton().addResourceLocation(OGRE_MEDIA_FOLDER "/packs/SinBad.zip", "Zip");
    ResourceGroupManager::getSingleton().addResourceLocation(OGRE_MEDIA_FOLDER "/packs/dragon.zip", "Zip");
    ResourceGroupManager::getSingleton().addResourceLocation(OGRE_MEDIA_FOLDER "/models/", "FileSystem");
    //ResourceGroupManager::getSingleton().addResourceLocation(OGRE_MEDIA_FOLDER "/materials/scripts/", "FileSystem");
   // ResourceGroupManager::getSingleton().addResourceLocation(OGRE_MEDIA_FOLDER "/materials/programs/Cg/", "FileSystem", ResourceGroupManager::DEFAULT_RESOURCE_GROUP_NAME, true);
    //ResourceGroupManager::getSingleton().addResourceLocation(OGRE_MEDIA_FOLDER "/materials/programs/glsl/", "FileSystem", ResourceGroupManager::DEFAULT_RESOURCE_GROUP_NAME, true);
    ResourceGroupManager::getSingleton().addResourceLocation(OGRE_MEDIA_FOLDER "/materials/textures/", "FileSystem");
    ResourceGroupManager::getSingleton().addResourceLocation(PROJECT_DATA_DIR_STR "/", "FileSystem");

    ResourceGroupManager::getSingleton().initialiseAllResourceGroups();
}

void OgreWidget::prepareMeshes() {
    createCube("Cube");
    createSphere("Sphere", 100, 50, 50);
}


void OgreWidget::prepareMaterials() {

}

void OgreWidget::createScene() {

    _sceneMgr->setAmbientLight(Ogre::ColourValue(0, 0, 0));
    _sceneMgr->setShadowTechnique(Ogre::SHADOWTYPE_STENCIL_ADDITIVE);

    // entities
    //Ogre::Entity* entCharacter = _sceneMgr->createEntity("Hero", "Ninja.mesh");
    //entCharacter->setCastShadows(true);
    //auto characterNode = _sceneMgr->getRootSceneNode()->createChildSceneNode();
    //characterNode->attachObject(entCharacter);

    Ogre::Plane plane(Ogre::Vector3::UNIT_Y, 0);
    Ogre::MeshManager::getSingleton().createPlane("ground", Ogre::ResourceGroupManager::DEFAULT_RESOURCE_GROUP_NAME,
        plane, 1500, 1500, 20, 20, true, 1, 5, 5, Ogre::Vector3::UNIT_Z);
    Ogre::Entity* entGround = _sceneMgr->createEntity("GroundEntity", "ground");
    _sceneMgr->getRootSceneNode()->createChildSceneNode()->attachObject(entGround);

    entGround->setMaterialName("Examples/Rockwall");
    //entGround->setMaterialName("Panorama/Test");
    entGround->setCastShadows(false);

    

    // custom cube
    Entity* entCube = _sceneMgr->createEntity("cc", "Cube");
    entCube->setMaterialName("Panorama/Panorama");

    SceneNode* cubeNode = _sceneMgr->getRootSceneNode()->createChildSceneNode();
    cubeNode->setPosition(0, 100, 0);
    cubeNode->attachObject(entCube);
    _focusedNode = cubeNode;


    // lights
    Ogre::Light* pointLight = _sceneMgr->createLight("pointLight");
    pointLight->setType(Ogre::Light::LT_POINT);
    pointLight->setPosition(Ogre::Vector3(0, 350, 250));

    pointLight->setDiffuseColour(ToQgreColour(Qt::gray));
    pointLight->setSpecularColour(ToQgreColour(Qt::gray));

    Ogre::Light* directionalLight = _sceneMgr->createLight("directionalLight");
    directionalLight->setType(Ogre::Light::LT_DIRECTIONAL);
    directionalLight->setDiffuseColour(ToQgreColour(Qt::darkYellow));
    directionalLight->setSpecularColour(ToQgreColour(Qt::yellow));

    directionalLight->setDirection(Ogre::Vector3(0, -1, 1));

    Ogre::Light* spotLight = _sceneMgr->createLight("spotLight");
    spotLight->setType(Ogre::Light::LT_SPOTLIGHT);
    spotLight->setDiffuseColour(ToQgreColour(Qt::cyan));
    spotLight->setSpecularColour(ToQgreColour(Qt::cyan));
    spotLight->setAttenuation(1000, 0.2, 1, 0);

    spotLight->setDirection(-1, -1, 0);
    spotLight->setPosition(Ogre::Vector3(300, 300, 0));

    spotLight->setSpotlightRange(Ogre::Degree(35), Ogre::Degree(50));


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
        params["FSAA"] = "16";
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
    _sceneMgr->setAmbientLight(ToQgreColour(Qt::gray));
    _sceneMgr->setShadowTechnique(Ogre::SHADOWDETAILTYPE_ADDITIVE);

    createCamera();
    createViewports();
    prepareResources();
    prepareMaterials();
    prepareMeshes();
    createScene();   

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
}


