#include "utilities.hpp"
#include "cameras.hpp"

namespace panoramix {
    namespace core {

        PerspectiveCamera::PerspectiveCamera(int w, int h, double focal, const Vec3 & eye,
            const Vec3 & center, const Vec3 & up, double near, double far)
            : _screenW(w), _screenH(h), _focal(focal), _eye(eye), _center(center), _up(up), _near(near), _far(far) {
            updateMatrices();
        }

        void PerspectiveCamera::updateMatrices() {
            _viewMatrix = MakeMat4LookAt(_eye, _center, _up);

            double verticalViewAngle = atan(_screenH / 2.0 / _focal) * 2;
            double aspect = double(_screenW) / double(_screenH);
            _projectionMatrix = MakeMat4Perspective(verticalViewAngle, aspect, _near, _far);

            _viewProjectionMatrix = _projectionMatrix * _viewMatrix;
            //_viewProjectionMatrixInv = _viewProjectionMatrix.inv(cv::DECOMP_LU);
        }

        Vec2 PerspectiveCamera::screenProjection(const Vec3 & p3) const {
            Vec4 p4(p3(0), p3(1), p3(2), 1);
            Vec4 position = _viewProjectionMatrix * p4;
            double xratio = position(0) / position(3) / 2;
            double yratio = position(1) / position(3) / 2;
            double x = (xratio + 0.5) * _screenW;
            double y = _screenH - (yratio + 0.5) * _screenH;
            return Vec2(x, y);
        }

        bool PerspectiveCamera::isVisibleOnScreen(const Vec3 & p3d) const {
            Vec4 p4(p3d(0), p3d(1), p3d(2), 1);
            Vec4 position = _viewProjectionMatrix * p4;
            return position(3) > 0 && position(2) > 0;
        }

        HPoint2 PerspectiveCamera::screenProjectionInHPoint(const Vec3 & p3) const {
            Vec4 p4(p3(0), p3(1), p3(2), 1);
            Vec4 position = _viewProjectionMatrix * p4;
            double xratio = position(0) / 2;
            double yratio = position(1) / 2;
            double zratio = position(3);

            double x = (xratio + 0.5 * zratio) * _screenW;
            double y = _screenH * zratio - (yratio + 0.5 * zratio) * _screenH;
            return HPoint2({ x, y }, zratio);
        }

        Vec3 PerspectiveCamera::spatialDirection(const Vec2 & p2d) const {
            double xratio = (p2d(0) / _screenW - 0.5) * 2;
            double yratio = ((_screenH - p2d(1)) / _screenH - 0.5) * 2;
            Vec4 position(xratio, yratio, 1, 1);
            Vec4 realPosition;
            bool solvable = cv::solve(_viewProjectionMatrix, position, realPosition);// _viewProjectionMatrixInv * position;
            assert(solvable);
            return Vec3(realPosition(0) / realPosition(3),
                realPosition(1) / realPosition(3),
                realPosition(2) / realPosition(3));
        }

        void PerspectiveCamera::resizeScreen(const Size & sz, bool updateMat) {
            if (_screenH == sz.height && _screenW == sz.width)
                return;
            _screenH = sz.height;
            _screenW = sz.width;
            if (updateMat)
                updateMatrices();
        }

        void PerspectiveCamera::setFocal(double f, bool updateMat) {
            if (f == _focal)
                return;
            _focal = f;
            if (updateMat)
                updateMatrices();
        }

        void PerspectiveCamera::setEye(const Vec3 & e, bool updateMat) {
            if (_eye == e)
                return;
            _eye = e;
            if (updateMat)
                updateMatrices();
        }

        void PerspectiveCamera::setCenter(const Vec3 & c, bool updateMat) {
            if (_center == c)
                return;
            _center = c;
            if (updateMat)
                updateMatrices();
        }

        void PerspectiveCamera::setUp(const Vec3 & up, bool updateMat) {
            if (_up == up)
                return;
            _up = up;
            if (updateMat)
                updateMatrices();
        }

        void PerspectiveCamera::setNearAndFarPlanes(double near, double far, bool updateMat) {
            if (_near == near && _far == far)
                return;
            _near = near;
            _far = far;
            if (updateMat)
                updateMatrices();
        }

        inline void AdjustNearAndFar(double & n, double & f, const Sphere3 & target, const Vec3 & eye){
            n = BoundBetween(norm(target.center - eye) - target.radius, 1e-2, 1e8);
            f = BoundBetween(norm(target.center - eye) + target.radius, 1e2, 1e8);
        }

        void PerspectiveCamera::focusOn(const Sphere3 & target, bool updateMat) {
            _center = target.center;
            auto eyedirection = _eye - _center;
            eyedirection = eyedirection / core::norm(eyedirection) * target.radius * 0.8;
            _eye = _center + eyedirection;
            AdjustNearAndFar(_near, _far, target, _eye);
            if (updateMat)
                updateMatrices();
        }

        void PerspectiveCamera::translate(const Vec3 & t, const Sphere3 & target, bool updateMat){
            _eye += t;
            _center += t;
            AdjustNearAndFar(_near, _far, target, _eye);
            if (updateMat)
                updateMatrices();
        }

        void PerspectiveCamera::moveEyeWithCenterFixed(const Vec3 & t, const Sphere3 & target, bool distanceFixed, bool updateMat){
            double dist = norm(_eye - _center);
            _eye += t;
            if (distanceFixed){
                _eye = normalize(_eye - _center) * dist + _center;
            }
            AdjustNearAndFar(_near, _far, target, _eye);
            if (updateMat)
                updateMatrices();
        }





        PanoramicCamera::PanoramicCamera(double focal, const Vec3 & eye,
            const Vec3 & center, const Vec3 & up)
            : _focal(focal), _eye(eye), _center(center), _up(up) {
            _xaxis = (_center - _eye); _xaxis /= core::norm(_xaxis);
            _yaxis = _up.cross(_xaxis); _yaxis /= core::norm(_yaxis);
            _zaxis = _xaxis.cross(_yaxis);
        }

        Vec2 PanoramicCamera::screenProjection(const Vec3 & p3) const {
            double xx = p3.dot(_xaxis);
            double yy = p3.dot(_yaxis);
            double zz = p3.dot(_zaxis);
            GeoCoord pg = core::Vec3(xx, yy, zz);
            auto sz = screenSize();
            double x = (pg.longitude + M_PI) / 2.0 / M_PI * sz.width;
            double y = (pg.latitude + M_PI_2) / M_PI * sz.height;
            return Vec2(x, y);
        }

        Vec3 PanoramicCamera::spatialDirection(const Vec2 & p2d) const {
            auto sz = screenSize();
            double longi = p2d(0) / double(sz.width) * 2 * M_PI - M_PI;
            double lati = p2d(1) / double(sz.height) * M_PI - M_PI_2;
            Vec3 dd = (GeoCoord(longi, lati).toVector());
            return dd(0) * _xaxis + dd(1) * _yaxis + dd(2) * _zaxis;
        }




        PartialPanoramicCamera::PartialPanoramicCamera(int w, int h, double focal, const Vec3 & eye,
            const Vec3 & center, const Vec3 & up)
            : _screenW(w), _screenH(h), _focal(focal), _eye(eye), _center(center), _up(up){
            _xaxis = (_center - _eye); _xaxis /= core::norm(_xaxis);
            _yaxis = _up.cross(_xaxis); _yaxis /= core::norm(_yaxis);
            _zaxis = _xaxis.cross(_yaxis);
        }

        PartialPanoramicCamera::PartialPanoramicCamera(const PanoramicCamera & panoCam, 
            int w, int h) 
            : _screenW(w), _screenH(h), _focal(panoCam.focal()), _eye(panoCam.eye()), _center(panoCam.center()), _up(panoCam.up()){
            _xaxis = (_center - _eye); _xaxis /= core::norm(_xaxis);
            _yaxis = _up.cross(_xaxis); _yaxis /= core::norm(_yaxis);
            _zaxis = _xaxis.cross(_yaxis);
        }

        Vec2 PartialPanoramicCamera::screenProjection(const Vec3 & p3) const {
            double xx = (p3 - _eye).dot(_xaxis);
            double yy = (p3 - _eye).dot(_yaxis);
            double zz = (p3 - _eye).dot(_zaxis);
            GeoCoord pg = core::Vec3(xx, yy, zz);

            double halfLongitudeAngleBound = _screenW / 2.0 / _focal;
            double halfLatitudeAngleBound = _screenH / 2.0 / _focal;

            double x = (pg.longitude + halfLongitudeAngleBound) * _focal;
            double y = (pg.latitude + halfLatitudeAngleBound) * _focal;

            return Vec2(x, y);
        }

        bool PartialPanoramicCamera::isVisibleOnScreen(const Vec3 & p3) const{
            double xx = (p3 - _eye).dot(_xaxis);
            double yy = (p3 - _eye).dot(_yaxis);
            double zz = (p3 - _eye).dot(_zaxis);
            GeoCoord pg = core::Vec3(xx, yy, zz);

            double halfLongitudeAngleBound = _screenW / 2.0 / _focal;
            double halfLatitudeAngleBound = _screenH / 2.0 / _focal;

            double x = (pg.longitude + halfLongitudeAngleBound) * _focal;
            double y = (pg.latitude + halfLatitudeAngleBound) * _focal;

            return IsBetween(x, 0, _screenW) && IsBetween(y, 0, _screenH);
        }

        Vec3 PartialPanoramicCamera::spatialDirection(const Vec2 & p2d) const {
            double halfLongitudeAngleBound = _screenW / 2.0 / _focal;
            double halfLatitudeAngleBound = _screenH / 2.0 / _focal;

            double longi = p2d(0) / _focal - halfLongitudeAngleBound;
            double lati = p2d(1) / _focal - halfLatitudeAngleBound;

            Vec3 dd = (GeoCoord(longi, lati).toVector());
            return dd(0) * _xaxis + dd(1) * _yaxis + dd(2) * _zaxis;
        }










        std::vector<PerspectiveCamera> CreateHorizontalPerspectiveCameras(const PanoramicCamera & panoCam, int num,
            int width, int height, double focal){
            std::vector<PerspectiveCamera> cams(num);
            Vec3 x, y;
            std::tie(x, y) = ProposeXYDirectionsFromZDirection(-panoCam.up());
            for (int i = 0; i < num; i++){
                double angle = M_PI * 2.0 / num * i;
                Vec3 direction = sin(angle) * x + cos(angle) * y;
                cams[i] = PerspectiveCamera(width, height, focal, panoCam.eye(), panoCam.eye() + direction, -panoCam.up(), 0.01, 1e4);
            }
            return cams;
        }

        std::vector<PerspectiveCamera> CreateHorizontalPerspectiveCameras(const PanoramicCamera & panoCam,
            const std::vector<Vec3> & dirs,
            int width, int height, double focal, double angleThreshold){
            std::vector<PerspectiveCamera> cams;
            for (int i = 0; i < dirs.size(); i++){
                if (AngleBetweenUndirectedVectors(dirs[i], panoCam.up()) < angleThreshold){
                    continue;
                }
                cams.emplace_back(width, height, focal, panoCam.eye(), panoCam.eye() + dirs[i], -panoCam.up(), 0.01, 1e4);
                cams.emplace_back(width, height, focal, panoCam.eye(), panoCam.eye() - dirs[i], -panoCam.up(), 0.01, 1e4);
            }
            return cams;
        }

        std::vector<PerspectiveCamera> CreateCubicFacedCameras(const PanoramicCamera & panoCam,
            int width/* = 500*/, int height/* = 500*/, double focal/* = 250.0*/){
            Vec3 x, y;
            std::tie(x, y) = ProposeXYDirectionsFromZDirection(-panoCam.up());
            std::vector<core::PerspectiveCamera> cams = {
                core::PerspectiveCamera(width, height, focal, panoCam.eye(), panoCam.eye() + Vec3{ 1, 0, 0 },  -panoCam.up()),
                core::PerspectiveCamera(width, height, focal, panoCam.eye(), panoCam.eye() + Vec3{ 0, 1, 0 },  -panoCam.up()),
                core::PerspectiveCamera(width, height, focal, panoCam.eye(), panoCam.eye() + Vec3{ -1, 0, 0 }, -panoCam.up()),
                core::PerspectiveCamera(width, height, focal, panoCam.eye(), panoCam.eye() + Vec3{ 0, -1, 0 }, -panoCam.up()),
                core::PerspectiveCamera(width, height, focal, panoCam.eye(), panoCam.eye() + Vec3{ 0, 0, 1 }, x),
                core::PerspectiveCamera(width, height, focal, panoCam.eye(), panoCam.eye() + Vec3{ 0, 0, -1 }, x)
            };
            return cams;
        }




        View<PanoramicCamera> CreatePanoramicView(const Image & panorama,
            const Point3 & eye, const Point3 & center, const Vec3 & up) {
            assert(abs(panorama.cols - panorama.rows * 2) < panorama.rows / 10.0f);
            return View<PanoramicCamera>{
                panorama,
                    PanoramicCamera(panorama.cols / M_PI / 2.0, eye, center, up)
            };
        }

        View<PerspectiveCamera> CreatePerspectiveView(const Image & perspectiveImage,
            const Point3 & eye, const Point3 & center, const Vec3 & up,
            const LineSegmentExtractor & lse, const VanishingPointsDetector & vpd,
            std::vector<Classified<Line3>> * line3sPtr,
            std::vector<Classified<Line2>> * line2sPtr,
            std::vector<Vec3> * vpsPtr,
            double * focalPtr){

            auto lines = lse(perspectiveImage, 2);
            std::vector<HPoint2> vps;
            double focal;
            std::vector<int> lineClasses;
            std::tie(vps, focal, lineClasses) = vpd(lines, perspectiveImage.size()).unwrap();
            assert(vps.size() >= 3);

            // only reserve first 3 vps
            vps.erase(vps.begin() + 3, vps.end());
            for (auto & c : lineClasses){
                if (c >= 3) c = -1;
            }

            View<PerspectiveCamera> view;
            view.image = perspectiveImage;
            view.camera = PerspectiveCamera(perspectiveImage.cols, perspectiveImage.rows, focal, eye, center, up);
                
            if (line3sPtr){
                std::vector<Classified<Line3>> line3s(lines.size());
                for (int i = 0; i < lines.size(); i++){
                    line3s[i].component.first = view.camera.spatialDirection(lines[i].first);
                    line3s[i].component.second = view.camera.spatialDirection(lines[i].second);
                    line3s[i].claz = lineClasses[i];
                }
                *line3sPtr = std::move(line3s);
            }
            if (line2sPtr){
                auto line2s = ClassifyEachAs(lines, -1);
                for (int i = 0; i < lines.size(); i++)
                    line2s[i].claz = lineClasses[i];
                *line2sPtr = std::move(line2s);
            }
            if (vpsPtr){
                std::vector<Vec3> vp3s(vps.size());
                for (int i = 0; i < vps.size(); i++){
                    vp3s[i] = normalize(view.camera.spatialDirection(vps[i].value()));
                }
                *vpsPtr = std::move(vp3s);
            }
            if (focalPtr){
                *focalPtr = focal;
            }

            return view;

        }



    }
}