#pragma once

#include <QtGui>
#include <QtOpenGL>

#include <Ogre.h>
#include <OgreFrameListener.h>

class OgreWidget : public QGLWidget {
public:
    OgreWidget(QWidget * parent = nullptr);
    virtual ~OgreWidget();

private:
    void createCamera();
    void createViewports();
    void prepareResources();
    void prepareMeshes();
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
};