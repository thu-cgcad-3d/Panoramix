#pragma once

#include "basic_types.hpp"
#include "feature.hpp"

namespace pano {
namespace core {

// perspective camera
class PerspectiveCamera {
public:
  PerspectiveCamera();
  PerspectiveCamera(int w, int h);
  PerspectiveCamera(int w, int h, const Point2 &pp, const Vec2 &focalxy,
                    const Point3 &eye = Point3(0, 0, 0),
                    const Point3 &center = Point3(1, 0, 0),
                    const Vec3 &up = Vec3(0, 0, -1), double nearPlane = 0.01,
                    double farPlane = 1e4);
  PerspectiveCamera(int w, int h, const Point2 &pp, double focal,
                    const Point3 &eye = Point3(0, 0, 0),
                    const Point3 &center = Point3(1, 0, 0),
                    const Vec3 &up = Vec3(0, 0, -1), double nearPlane = 0.01,
                    double farPlane = 1e4);

  inline Sizei screenSize() const {
    return Sizei(static_cast<int>(_screenW), static_cast<int>(_screenH));
  }
  inline double screenWidth() const { return _screenW; }
  inline double screenHeight() const { return _screenH; }
  inline const Point2 &principlePoint() const { return _principlePoint; }
  inline double aspect() const { return double(_screenW) / double(_screenH); }
  inline double focalX() const { return _focalxy[0]; }
  inline double focalY() const { return _focalxy[1]; }
  inline double focal() const {
    assert(_focalxy[0] == _focalxy[1]);
    return _focalxy[0];
  }
  inline const Point3 &eye() const { return _eye; }
  inline const Point3 &center() const { return _center; }
  inline const Vec3 &up() const { return _up; }
  inline double nearPlane() const { return _near; }
  inline double farPlane() const { return _far; }
  Point2 toScreen(const Point3 &p3d) const;
  bool isVisibleOnScreen(const Point3 &p3d) const;
  HPoint2 toScreenInHPoint(const Point3 &p3d) const;
  Point3 toSpace(const Point2 &p2d) const;
  inline Point3 toSpace(const Pixel &p) const {
    return toSpace(Point2(p.x, p.y));
  }
  inline Vec3 direction(const Point2 &p2d) const { return toSpace(p2d) - _eye; }
  inline Vec3 direction(const Pixel &p) const { return direction(p); }

  inline const Mat4 &viewMatrix() const { return _viewMatrix; }
  inline const Mat4 &projectionMatrix() const { return _projectionMatrix; }
  inline const Mat4 &viewProjectionMatrix() const {
    return _viewProjectionMatrix;
  }

  // operations
  void resizeScreen(const Size &sz, bool updateMat = true);
  void setPrinciplePoint(const Point2 &pp, bool updateMat = true);
  void setFocalX(double fx, bool updateMat = true);
  void setFocalY(double fy, bool updateMat = true);
  void setFocal(double f, bool updateMat = true);
  void setEye(const Point3 &e, bool updateMat = true);
  void setCenter(const Point3 &c, bool updateMat = true);
  void setUp(const Vec3 &up, bool updateMat = true);
  void setNearAndFarPlanes(double nearPlane, double farPlane,
                           bool updateMat = true);

  // advanced functions
  inline Vec3 forward() const { return normalize(_center - _eye); }
  inline Vec3 backward() const { return normalize(_eye - _center); }
  inline Vec3 upward() const { return normalize(_up); }
  inline Vec3 downward() const { return -upward(); }
  inline Vec3 leftward() const { return normalize(_up.cross(forward())); }
  inline Vec3 rightward() const { return -leftward(); }

  void focusOn(const Sphere3 &target, bool updateMat = true);
  void translate(const Vec3 &t, const Sphere3 &target, bool updateMat = true);
  void moveEyeWithCenterFixed(const Vec3 &t, const Sphere3 &target,
                              bool distanceFixed = false,
                              bool updateMat = true);
  void moveCenterWithEyeFixed(const Vec3 &t, bool updateMat = true);

  // compare
  bool operator==(const PerspectiveCamera &cam) const;

protected:
  void updateMatrices();

protected:
  double _screenW, _screenH;
  Point2 _principlePoint;
  Vec2 _focalxy;
  double _near, _far;
  Vec3 _eye, _center, _up;
  Mat4 _viewMatrix, _projectionMatrix, _viewProjectionMatrix;

  template <class Archive> inline void serialize(Archive &ar) {
    ar(_screenW, _screenH);
    ar(_principlePoint);
    ar(_focalxy, _near, _far);
    ar(_eye, _center, _up);
    ar(_viewMatrix, _projectionMatrix, _viewProjectionMatrix);
  }
  friend class cereal::access;
};

// panoramic camera
class PanoramicCamera {
public:
  explicit PanoramicCamera(double focal = 250,
                           const Point3 &eye = Vec3(0, 0, 0),
                           const Point3 &center = Vec3(1, 0, 0),
                           const Vec3 &up = Vec3(0, 0, 1));

  inline Sizei screenSize() const {
    return Sizei(static_cast<int>(_focal * 2 * M_PI),
                 static_cast<int>(_focal * M_PI));
  }
  inline double focal() const { return _focal; }
  inline const Point3 &eye() const { return _eye; }
  inline const Point3 &center() const { return _center; }
  inline const Vec3 &up() const { return _up; }
  Point2 toScreen(const Point3 &p3d) const;
  bool isVisibleOnScreen(const Point3 &p3d) const { return true; }
  inline HPoint2 toScreenInHPoint(const Point3 &p3d) const {
    return HPoint2(toScreen(p3d), 1.0);
  }
  Point3 toSpace(const Point2 &p2d) const;
  inline Point3 toSpace(const Pixel &p) const {
    return toSpace(Vec2(p.x, p.y));
  }
  inline Point3 toSpace(const HPoint2 &p) const { return toSpace(p.value()); }
  Vec3 direction(const Point2 &p2d) const;
  inline Vec3 direction(const Pixel &p) const { return direction(p); }

private:
  double _focal;
  Vec3 _eye, _center, _up;
  Vec3 _xaxis, _yaxis, _zaxis;

  template <class Archive> inline void serialize(Archive &ar) {
    ar(_focal);
    ar(_eye, _center, _up);
    ar(_xaxis, _yaxis, _zaxis);
  }
  friend class cereal::access;
};

// partial panoramic camera
class PartialPanoramicCamera {
public:
  explicit PartialPanoramicCamera(int w = 500, int h = 500, double focal = 250,
                                  const Point3 &eye = Vec3(0, 0, 0),
                                  const Point3 &center = Vec3(1, 0, 0),
                                  const Vec3 &up = Vec3(0, 0, -1));
  explicit PartialPanoramicCamera(const PanoramicCamera &panoCam, int w = 500,
                                  int h = 500);

  inline Sizei screenSize() const {
    return Size(static_cast<int>(_screenW), static_cast<int>(_screenH));
  }
  inline double focal() const { return _focal; }
  inline const Point3 &eye() const { return _eye; }
  inline const Point3 &center() const { return _center; }
  inline const Vec3 &up() const { return _up; }
  Point2 toScreen(const Point3 &p3d) const;
  bool isVisibleOnScreen(const Point3 &p3d) const;
  inline HPoint2 toScreenInHPoint(const Point3 &p3d) const {
    return HPoint2(toScreen(p3d), 1.0);
  }
  Point3 toSpace(const Point2 &p2d) const;
  inline Point3 toSpace(const Pixel &p) const {
    return toSpace(Vec2(p.x, p.y));
  }
  Vec3 direction(const Point2 &p2d) const;
  inline Vec3 direction(const Pixel &p) const { return direction(p); }
  PanoramicCamera toPanoramic() const {
    return PanoramicCamera(_focal, _eye, _center, _up);
  }

private:
  double _screenW, _screenH;
  double _focal;
  Vec3 _eye, _center, _up;
  Vec3 _xaxis, _yaxis, _zaxis;

  template <class Archive> inline void serialize(Archive &ar) {
    ar(_screenW, _screenH);
    ar(_focal);
    ar(_eye, _center, _up);
    ar(_xaxis, _yaxis, _zaxis);
  }
  friend class cereal::access;
};


namespace {
template <class T> struct IsCameraImpl {
  template <class TT>
  static auto test(int) -> decltype(
      std::declval<TT>().eye(),
      std::declval<TT>().toScreen(std::declval<core::Point3>()),
      std::declval<TT>().toScreenInHPoint(std::declval<core::Point3>()),
      std::declval<TT>().isVisibleOnScreen(std::declval<core::Point3>()),
      std::declval<TT>().toSpace(std::declval<core::Point2>()),
      std::declval<TT>().toSpace(std::declval<core::Pixel>()),
      std::declval<TT>().screenSize(), std::true_type()) {
    return std::true_type();
  }
  template <class> static std::false_type test(...) {
    return std::false_type();
  }
  static const bool value =
      std::is_same<decltype(test<T>(0)), std::true_type>::value;
};
}

// judge whether T is a camera type
template <class T>
struct IsCamera : std::integral_constant<bool, IsCameraImpl<T>::value> {};

// sample image from image using camera conversion
template <class OutCameraT, class InCameraT> class CameraSampler {
  static_assert(IsCamera<OutCameraT>::value && IsCamera<InCameraT>::value,
                "OutCameraT and InCameraT should both be cameras!");

public:
  template <class OCamT, class ICamT>
  CameraSampler(OCamT &&outCam, ICamT &&inCam)
      : _outCam(std::forward<OCamT>(outCam)),
        _inCam(std::forward<ICamT>(inCam)) {
    assert(outCam.eye() == inCam.eye());
    auto outCamSize = _outCam.screenSize();
    _mapx = cv::Mat::zeros(outCamSize, CV_32FC1);
    _mapy = cv::Mat::zeros(outCamSize, CV_32FC1);
    for (int j = 0; j < outCamSize.height; j++) {
      for (int i = 0; i < outCamSize.width; i++) {
        Vec2 screenp(i, j);
        Vec3 p3 = _outCam.toSpace(screenp);
        if (!_inCam.isVisibleOnScreen(p3)) {
          _mapx.at<float>(j, i) = -1;
          _mapy.at<float>(j, i) = -1;
          continue;
        }
        Vec2 screenpOnInCam = _inCam.toScreen(p3);
        _mapx.at<float>(j, i) = static_cast<float>(screenpOnInCam(0));
        _mapy.at<float>(j, i) = static_cast<float>(screenpOnInCam(1));
      }
    }
  }

  Image operator()(const Image &inputIm, int borderMode = cv::BORDER_REPLICATE,
                   const cv::Scalar &borderValue = cv::Scalar(0, 0, 0,
                                                              0)) const {
    Image outputIm;
    cv::remap(inputIm, outputIm, _mapx, _mapy,
              inputIm.channels() <= 4 ? cv::INTER_LINEAR : cv::INTER_NEAREST,
              borderMode, borderValue);
    return outputIm;
  }

  template <class T>
  ImageOf<T> operator()(const ImageOf<T> &inputIm,
                        int borderMode = cv::BORDER_REPLICATE,
                        const T &borderValue = T()) const {
    ImageOf<T> outputIm;
    cv::remap(inputIm, outputIm, _mapx, _mapy, cv::INTER_NEAREST, borderMode,
              borderValue);
    return outputIm;
  }

  template <class T, int N, class = std::enable_if_t<(N <= 4)>>
  ImageOf<Vec<T, N>>
  operator()(const ImageOf<Vec<T, N>> &inputIm,
             int borderMode = cv::BORDER_REPLICATE,
             const Vec<T, N> &borderValue = Vec<T, N>()) const {
    Image outputIm;
    cv::Scalar bv;
    for (int i = 0; i < N; i++) {
      bv[i] = borderValue[i];
    }
    cv::remap(inputIm, outputIm, _mapx, _mapy, cv::INTER_NEAREST, borderMode,
              bv);
    return outputIm;
  }

  template <class T, int N, class = std::enable_if_t<(N > 4)>, class = void>
  ImageOf<Vec<T, N>>
  operator()(const ImageOf<Vec<T, N>> &inputIm,
             int borderMode = cv::BORDER_REPLICATE,
             const Vec<T, N> &borderValue = Vec<T, N>()) const {
    std::vector<Image> channels;
    cv::split(inputIm, channels);
    for (int i = 0; i < N; i++) {
      auto &c = channels[i];
      cv::remap(c, c, _mapx, _mapy, cv::INTER_NEAREST, borderMode,
                borderValue[i]);
    }
    ImageOf<Vec<T, N>> result;
    cv::merge(channels, result);
    return result;
  }

private:
  OutCameraT _outCam;
  InCameraT _inCam;
  cv::Mat _mapx, _mapy;
};

template <class OutCameraT, class InCameraT>
inline CameraSampler<std::decay_t<OutCameraT>, std::decay_t<InCameraT>>
MakeCameraSampler(OutCameraT &&outCam, InCameraT &&inCam) {
  return CameraSampler<std::decay_t<OutCameraT>, std::decay_t<InCameraT>>(
      std::forward<OutCameraT>(outCam), std::forward<InCameraT>(inCam));
}

template <class CameraT, class = std::enable_if_t<IsCamera<CameraT>::value>>
struct DimensionConvertor {
  const CameraT &cam;
  DimensionConvertor(const CameraT &c) : cam(c) {}

  Point3 toSpace(const Point2 &p) const { return cam.toSpace(p); }
  Point2 toScreen(const Point3 &p) const { return cam.toScreen(p); }

  Line3 toSpace(const Line2 &line2) const {
    return Line3(cam.toSpace(line2.first), cam.toSpace(line2.second));
  }
  Line2 toScreen(const Line3 &line3) const {
    return Line2(cam.toScreen(line3.first), cam.toScreen(line3.second));
  }

  Ray3 toSpace(const Ray2 &ray) const {
    return Line3(cam.toSpace(ray.anchor),
                 cam.toSpace(ray.anchor + ray.direction))
        .ray();
  }
  Ray2 toScreen(const Ray3 &ray) const {
    return toScreen(Line3(ray.anchor, ray.anchor + ray.direction)).ray();
  }

  Chain3 toSpace(const Chain2 &chain) const {
    Chain3 c;
    c.closed = chain.closed;
    c.points.resize(chain.points.size());
    for (int i = 0; i < c.points.size(); i++) {
      c.points[i] = toSpace(chain.points[i]);
    }
    return c;
  }
  Chain2 toScreen(const Chain3 &chain) const {
    Chain2 c;
    c.closed = chain.closed;
    c.points.resize(chain.points.size());
    for (int i = 0; i < c.points.size(); i++) {
      c.points[i] = toScreen(chain.points[i]);
    }
    return c;
  }

  template <class T, class = std::enable_if_t<IsDecorated<T>::value>>
  struct DecoratedToSpace {
    using type = decltype(std::declval<T>().converted(
        std::declval<DimensionConvertor<CameraT>>().toSpace(
            std::declval<T>().component)));
  };
  template <class T, class = std::enable_if_t<IsDecorated<T>::value>>
  struct DecoratedToScreen {
    using type = decltype(std::declval<T>().converted(
        std::declval<DimensionConvertor<CameraT>>().toScreen(
            std::declval<T>().component)));
  };

  template <class T> struct GeneralToSpace {
    using type =
        std::decay_t<decltype(std::declval<DimensionConvertor<CameraT>>()
                                  .toSpace(std::declval<T>()))>;
  };
  template <class T> struct GeneralToScreen {
    using type =
        std::decay_t<decltype(std::declval<DimensionConvertor<CameraT>>()
                                  .toScreen(std::declval<T>()))>;
  };

  template <class T, class = std::enable_if_t<IsDecorated<T>::value>>
  typename DecoratedToSpace<T>::type toSpace(const T &c) const {
    return c.converted(toSpace(c.component));
  }
  template <class T, class = std::enable_if_t<IsDecorated<T>::value>>
  typename DecoratedToScreen<T>::type toScreen(const T &c) const {
    return c.converted(toScreen(c.component));
  }

  template <class T>
  std::vector<typename GeneralToSpace<T>::type>
  toSpace(const std::vector<T> &v) const {
    std::vector<typename GeneralToSpace<T>::type> results(v.size());
    for (int i = 0; i < v.size(); i++)
      results[i] = toSpace(v[i]);
    return results;
  }
  template <class T>
  std::vector<typename GeneralToScreen<T>::type>
  toScreen(const std::vector<T> &v) const {
    std::vector<typename GeneralToScreen<T>::type> results(v.size());
    for (int i = 0; i < v.size(); i++)
      results[i] = toScreen(v[i]);
    return results;
  }
};

template <class CameraT>
DimensionConvertor<CameraT> AsDimensionConvertor(const CameraT &cam) {
  return DimensionConvertor<CameraT>(cam);
}

// create horizontal cameras
std::vector<PerspectiveCamera>
CreateHorizontalPerspectiveCameras(const PanoramicCamera &panoCam, int num = 16,
                                   int width = 500, int height = 500,
                                   double focal = 250.0);

std::vector<PerspectiveCamera> CreateHorizontalPerspectiveCameras(
    const PanoramicCamera &panoCam, const std::vector<Vec3> &dirs,
    int width = 500, int height = 500, double focal = 250.0,
    double angleThreshold = 0.1);

std::vector<PerspectiveCamera>
CreatePanoContextCameras(const PanoramicCamera &panoCam, int width = 500,
                         int height = 500, double focal = 250.0);

// create cameras toward cubic faces
std::vector<PerspectiveCamera>
CreateCubicFacedCameras(const PanoramicCamera &panoCam, int width = 500,
                        int height = 500, double focal = 250.0);

// view class
template <class CameraT, class ImageT = Image,
          class = std::enable_if_t<IsCamera<CameraT>::value>>
struct View {
  ImageT image;
  CameraT camera;

  View() {}
  View(const ImageT &im, const CameraT &cam) : image(im), camera(cam) {}
  View(const View<CameraT, Image> &v) : image(v.image), camera(v.camera) {}
  template <class T>
  View(const View<CameraT, ImageOf<T>> &v) : image(v.image), camera(v.camera) {}

  template <
      class AnotherCameraT,
      class = std::enable_if_t<IsCamera<std::decay_t<AnotherCameraT>>::value>>
  inline View<std::decay_t<AnotherCameraT>, ImageT>
  sampled(AnotherCameraT &&cam) const {
    View<std::decay_t<AnotherCameraT>, ImageT> v;
    v.image = MakeCameraSampler(cam, camera)(image);
    v.camera = std::forward<AnotherCameraT>(cam);
    return v;
  }

  template <class Archiver> void serialize(Archiver &ar) { ar(image, camera); }
};

template <class CameraT, class ImageT>
View<CameraT, ImageT> MakeView(const ImageT &im, const CameraT &c) {
  assert(im.size() == c.screenSize());
  return View<CameraT, ImageT>(im, c);
}

template <class OutCameraT, class InCameraT, class T,
          class = std::enable_if_t<IsCamera<std::decay_t<InCameraT>>::value &&
                                   IsCamera<std::decay_t<OutCameraT>>::value>>
View<OutCameraT, ImageOf<T>>
Combine(const OutCameraT &camera,
        const std::vector<View<InCameraT, ImageOf<T>>> &views) {
  if (views.empty()) {
    return View<OutCameraT, ImageOf<T>>();
  }
  ImageOf<T> converted = ImageOf<T>::zeros(camera.screenSize());
  static const int channels = cv::DataType<T>::channels;
  Imagef counts(camera.screenSize(), 0.0f);
  for (int i = 0; i < views.size(); i++) {
    auto sampler = MakeCameraSampler(camera, views[i].camera);
    auto piece = sampler(views[i].image, cv::BORDER_CONSTANT);
    converted += piece;
    counts += sampler(Imagef::ones(views[i].image.size()), cv::BORDER_CONSTANT);
  }
  View<OutCameraT, ImageOf<T>> v;
  v.camera = camera;
  v.image = ImageOf<T>::zeros(camera.screenSize());
  for (auto it = converted.begin(); it != converted.end(); ++it) {
    float count = counts(it.pos());
    v.image(it.pos()) = *it / std::max(count, 1.0f);
  }
  return v;
}

template <class OutCameraT, class InCameraT, class T, class W,
          class = std::enable_if_t<IsCamera<std::decay_t<InCameraT>>::value &&
                                   IsCamera<std::decay_t<OutCameraT>>::value>>
View<OutCameraT, ImageOf<T>>
Combine(const OutCameraT &camera,
        const std::vector<Weighted<View<InCameraT, ImageOf<T>>, W>> &views) {
  if (views.empty()) {
    return View<OutCameraT, ImageOf<T>>();
  }
  static const int channels = cv::DataType<T>::channels;
  using channel_type = typename cv::DataType<T>::channel_type;
  ImageOf<T> converted = ImageOf<T>::zeros(camera.screenSize());
  Imagef counts(camera.screenSize(), 0.0f);
  for (int i = 0; i < views.size(); i++) {
    auto sampler = MakeCameraSampler(camera, views[i].component.camera);
    auto piece = sampler(views[i].component.image, cv::BORDER_CONSTANT);
    for (auto it = converted.begin(); it != converted.end(); ++it) {
      *it += (piece(it.pos()) * views[i].weight());
    }
    counts +=
        sampler(Imagef(views[i].component.image.size(), views[i].weight()),
                cv::BORDER_CONSTANT);
  }
  View<OutCameraT, ImageOf<T>> v;
  v.camera = camera;
  v.image = ImageOf<T>::zeros(camera.screenSize());
  for (auto it = converted.begin(); it != converted.end(); ++it) {
    float count = counts(it.pos());
    v.image(it.pos()) = *it / std::max(count, 1.0f);
  }
  return v;
}

template <class CameraT, class InCameraT, class T>
T MeanInMask(const View<InCameraT, ImageOf<T>> &source,
             const View<CameraT, Imageub> &maskView) {
  ImageOf<T> converted = source.sampled(maskView.camera).image;
  int votes = 0;
  T featureSum;
  for (auto it = maskView.image.begin(); it != maskView.image.end(); ++it) {
    if (!*it) {
      continue;
    }
    featureSum += converted(it.pos());
    votes += 1;
  }
  return featureSum * (1.0 / std::max(votes, 1));
}

using PerspectiveView = View<PerspectiveCamera>;
using PanoramicView = View<PanoramicCamera>;
using PartialPanoramicView = View<PartialPanoramicCamera>;

PanoramicCamera CreatePanoramicCamera(const Image &panorama,
                                      const Point3 &eye = Point3(0, 0, 0),
                                      const Point3 &center = Point3(1, 0, 0),
                                      const Vec3 &up = Vec3(0, 0, 1));

Failable<PerspectiveCamera> CreatePerspectiveCamera(
    const Image &perspectiveImage, const Point3 &eye = Point3(0, 0, 0),
    const Point3 &center = Point3(1, 0, 0), const Vec3 &up = Vec3(0, 0, -1),
    const LineSegmentExtractor &lse = LineSegmentExtractor(),
    const VanishingPointsDetector &vpd = VanishingPointsDetector(),
    std::vector<Classified<Line3>> *line3s = nullptr,
    std::vector<Classified<Line2>> *line2s = nullptr,
    std::vector<Vec3> *vps = nullptr, double *focal = nullptr);

PerspectiveCamera CreatePerspectiveCamera(
    const Image &perspectiveImage, const std::vector<HPoint2> &vps,
    const Point3 &eye = Point3(0, 0, 0), const Point3 &center = Point3(1, 0, 0),
    const Vec3 &up = Vec3(0, 0, -1));

// create panoramic view
inline PanoramicView CreatePanoramicView(const Image &panorama,
                                         const Point3 &eye = Point3(0, 0, 0),
                                         const Point3 &center = Point3(1, 0, 0),
                                         const Vec3 &up = Vec3(0, 0, 1)) {
  return PanoramicView(panorama,
                       CreatePanoramicCamera(panorama, eye, center, up));
}
template <class T>
inline View<PanoramicCamera, ImageOf<T>> CreatePanoramicView(
    const ImageOf<T> &panorama, const Point3 &eye = Point3(0, 0, 0),
    const Point3 &center = Point3(1, 0, 0), const Vec3 &up = Vec3(0, 0, 1)) {
  return View<PanoramicCamera, ImageOf<T>>(
      panorama, CreatePanoramicCamera(panorama, eye, center, up));
}

// create perspective view
inline Failable<PerspectiveView> CreatePerspectiveView(
    const Image &perspectiveImage, const Point3 &eye = Point3(0, 0, 0),
    const Point3 &center = Point3(1, 0, 0), const Vec3 &up = Vec3(0, 0, -1),
    const LineSegmentExtractor &lse = LineSegmentExtractor(),
    const VanishingPointsDetector &vpd = VanishingPointsDetector(),
    std::vector<Classified<Line3>> *line3s = nullptr,
    std::vector<Classified<Line2>> *line2s = nullptr,
    std::vector<Vec3> *vps = nullptr, double *focal = nullptr) {
  auto cam = CreatePerspectiveCamera(perspectiveImage, eye, center, up, lse,
                                     vpd, line3s, line2s, vps, focal);
  if (cam.failed())
    return nullptr;
  return PerspectiveView(perspectiveImage, cam.unwrap());
}
template <class T>
inline Failable<View<PerspectiveCamera, ImageOf<T>>> CreatePerspectiveView(
    const ImageOf<T> &perspectiveImage, const Point3 &eye = Point3(0, 0, 0),
    const Point3 &center = Point3(1, 0, 0), const Vec3 &up = Vec3(0, 0, -1),
    const LineSegmentExtractor &lse = LineSegmentExtractor(),
    const VanishingPointsDetector &vpd = VanishingPointsDetector(),
    std::vector<Classified<Line3>> *line3s = nullptr,
    std::vector<Classified<Line2>> *line2s = nullptr,
    std::vector<Vec3> *vps = nullptr, double *focal = nullptr) {
  auto cam = CreatePerspectiveCamera(perspectiveImage, eye, center, up, lse,
                                     vpd, line3s, line2s, vps, focal);
  if (cam.failed())
    return nullptr;
  return View<PerspectiveCamera, ImageOf<T>>(perspectiveImage, cam.unwrap());
}

inline PerspectiveView CreatePerspectiveView(
    const Image &perspectiveImage, const std::vector<HPoint2> &vps,
    const Point3 &eye = Point3(0, 0, 0), const Point3 &center = Point3(1, 0, 0),
    const Vec3 &up = Vec3(0, 0, -1)) {
  return PerspectiveView(
      perspectiveImage,
      CreatePerspectiveCamera(perspectiveImage, vps, eye, center, up));
}
template <class T>
inline View<PerspectiveCamera, ImageOf<T>> CreatePerspectiveView(
    const ImageOf<T> &perspectiveImage, const std::vector<HPoint2> &vps,
    const Point3 &eye = Point3(0, 0, 0), const Point3 &center = Point3(1, 0, 0),
    const Vec3 &up = Vec3(0, 0, -1)) {
  return View<PerspectiveCamera, ImageOf<T>>(
      perspectiveImage,
      CreatePerspectiveCamera(perspectiveImage, vps, eye, center, up));
}

template <class CameraT, class = std::enable_if_t<IsCamera<CameraT>::value>>
Failable<View<CameraT>> CreateView(const Image &image,
                                   const Point3 &eye = Point3(0, 0, 0),
                                   const Point3 &center = Point3(1, 0, 0),
                                   const Vec3 &up = Vec3(0, 0, -1));

template <class CameraT, class T,
          class = std::enable_if_t<IsCamera<CameraT>::value>>
Failable<View<CameraT, ImageOf<T>>>
CreateView(const ImageOf<T> &image, const Point3 &eye = Point3(0, 0, 0),
           const Point3 &center = Point3(1, 0, 0),
           const Vec3 &up = Vec3(0, 0, -1)) {
  auto v = CreateView<CameraT>((const Image &)image, eye, center, up);
  if (v.failed())
    return nullptr;
  return View<CameraT, ImageOf<T>>(image, v.unwrap.camera);
}

template <>
inline Failable<View<PanoramicCamera>>
CreateView(const Image &image, const Point3 &eye, const Point3 &center,
           const Vec3 &up) {
  return CreatePanoramicView(image, eye, center, up);
}

template <>
inline Failable<View<PerspectiveCamera>>
CreateView(const Image &image, const Point3 &eye, const Point3 &center,
           const Vec3 &up) {
  return CreatePerspectiveView(image, eye, center, up);
}
}
}
