#pragma once

#include <Ogre.h>
#include <OgreFrameListener.h>

#include <QtGui>
#include <QtOpenGL>

#include "../../src/rec/reconstruction_engine.hpp"

class OgreWidget : public QGLWidget {
public:
    OgreWidget(QWidget * parent = nullptr);
    virtual ~OgreWidget();

    void setupPanorama(const QString & filename);

private:
    static void createCube(const Ogre::String & name);
    static void createSphere(const Ogre::String & name, const float r, const int nRings = 16, const int nSegments = 16);

    void createCamera();
    void createViewports();
    void prepareResources();
    void prepareMeshes();
    void prepareMaterials();
    void createScene();

protected:
    virtual void initializeGL();
    virtual void resizeGL(int width, int height);
    virtual void paintGL();

protected:
    void mousePressEvent(QMouseEvent *) override;
    void mouseMoveEvent(QMouseEvent *) override;
    void mouseReleaseEvent(QMouseEvent *) override;
    void wheelEvent(QWheelEvent *) override;

private:
    void moveCameraEyeWithCenterFixed(const QVector3D & t);
    void moveCameraCenterAndCenter(const QVector3D & t);

    Ogre::Root * _root;
    Ogre::RenderWindow *_window;
    Ogre::Camera *_camera;
    Ogre::SceneNode *_focusedNode;
    Ogre::Viewport *_viewport;
    Ogre::SceneManager *_sceneMgr;
    QPointF _lastPos;

    panoramix::rec::ReconstructionEngine _viewsNet;
};