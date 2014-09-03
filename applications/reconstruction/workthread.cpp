
#include "../../src/core/feature.hpp"
#include "../../src/core/mesh_maker.hpp"
#include "../../src/vis/qt_glue.hpp"

#include "../../src/rec/reconstruction_engine.hpp"

#include "workthread.hpp"

using namespace panoramix;

WorkThread::WorkThread(QObject *parent)
    : QThread(parent){
    
    core::Mesh<core::Vec3> cameraStand;
    core::MakeQuadFacedSphere(cameraStand, 4, 7);
    for (auto & v : cameraStand.vertices()){
        _camDirections.append(vis::MakeQVec(v.data));
    }

}

WorkThread::~WorkThread(){

}

void WorkThread::run(){
    QThread::msleep(10000);
}

void WorkThread::reconstructPanorama(const QImage & pan){




}