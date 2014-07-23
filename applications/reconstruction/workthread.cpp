
#include "../../src/core/mesh_maker.hpp"
#include "../../src/vis/qt_glue.hpp"
#include "../../src/rec/reconstruction_engine_visualize.hpp"

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

    cv::Mat panorama = panoramix::vis::MakeCVMat(pan);
    cv::resize(panorama, panorama, cv::Size(2000, 1000));
    core::PanoramicCamera originCam(panorama.cols / M_PI / 2.0);

    std::vector<core::PerspectiveCamera> cams;
    for (auto & v : _camDirections) {
        core::Vec3 direction = vis::MakeCoreVec(v);
        cams.emplace_back(700, 700, originCam.focal(), core::Vec3(0, 0, 0), direction, core::Vec3(0, 0, -1));
    }

    /// insert into views net
    rec::ReconstructionEngine::Params params;
    rec::ReconstructionEngine net(params);

    net.insertPanorama(panorama, cams, originCam);

#pragma omp parallel for
    for (int i = 0; i < net.views().internalElements<0>().size(); i++) {
        auto viewHandle = rec::ReconstructionEngine::ViewHandle(i);
        net.computeFeatures(viewHandle);
    }

    std::cout << "calibrating camera and classifying lines ...";
    net.estimateVanishingPointsAndClassifyLines();
    net.reconstructLinesAndFaces();

    // extract reconstructed regions
    using namespace core;
    using namespace rec;

    std::vector<std::vector<Point3>> regions;
    regions.reserve(net.constraints().internalElements<0>().size());
    for (auto & vd : net.constraints().elements<0>()){
        if (vd.data.type == ReconstructionEngine::ComponentData::Type::Region){
            auto theta = deriv::MakeCoreVec(vd.data.asRegion.thetaExpr.expression().result());
            auto & cam = net.views().data(vd.data.asRegion.viewHandle).camera;
            auto & regionData =
                net.views().data(vd.data.asRegion.viewHandle).regionNet->regions().data(vd.data.asRegion.regionHandle);
            auto & outCountours = regionData.contours.back();
            std::vector<Point3> spatialRegion(outCountours.size());
            for (int i = 0; i < outCountours.size(); i++){
                Point3 dir = cam.spatialDirection(outCountours[i]);
                dir /= norm(dir);
                dir *= (theta.dot(theta) / theta.dot(dir));
                spatialRegion[i] = dir;
            }
            regions.push_back(std::move(spatialRegion));
        }
    }



}