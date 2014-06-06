#pragma once

#include <QtGui>
#include <QtOpenGL>

#include <Ogre.h>
#include <OgreFrameListener.h>

class OgreWidget : public QGLWidget {
public:
    OgreWidget(QWidget * parent = nullptr);
    virtual ~OgreWidget();

    virtual void initializeGL();
    virtual void resizeGL(int width, int height);
    virtual void paintGL();

    Ogre::Root * _root;
    Ogre::RenderWindow *_window;
    Ogre::Camera *_camera;
    Ogre::Viewport *_viewport;
    Ogre::SceneManager *_sceneMgr;
};