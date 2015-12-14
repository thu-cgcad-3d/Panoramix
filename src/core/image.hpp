#pragma once

#include "geometry.hpp"

namespace pano {
namespace core {

using Image = cv::Mat;

template <class T> using ImageOf = cv::Mat_<T>;
using Imageb = ImageOf<bool>;
using Image3b = ImageOf<Vec<bool, 3>>;
using Imageub = ImageOf<uint8_t>;
using Image3ub = ImageOf<Vec<uint8_t, 3>>;
using Imagei = ImageOf<int>;
using Image3i = ImageOf<Vec<int, 3>>;
using Imagef = ImageOf<float>;
using Image3f = ImageOf<Vec<float, 3>>;
using Imaged = ImageOf<double>;
using Image3d = ImageOf<Vec<double, 3>>;
using Image5d = ImageOf<Vec<double, 5>>;
using Image6d = ImageOf<Vec<double, 6>>;
using Image7d = ImageOf<Vec<double, 7>>;

template <> struct IsNotContainerByHand<Image> : yes {};
template <class T> struct IsNotContainerByHand<ImageOf<T>> : yes {};

namespace {
template <class To, class From>
inline ImageOf<To> VecCastPrivate(const ImageOf<From> &im, std::false_type) {
  ImageOf<To> cim(im.size());
  for (auto it = im.begin(); it != im.end(); ++it) {
    cim(it.pos()) = static_cast<To>(*it);
  }
  return cim;
}
template <class To>
inline ImageOf<To> VecCastPrivate(const ImageOf<To> &v, std::true_type) {
  return v.clone();
}

template <class To, class From, int N>
inline ImageOf<Vec<To, N>> VecCastPrivate(const ImageOf<Vec<From, N>> &im,
                                          std::false_type) {
  ImageOf<Vec<To, N>> cim(im.size());
  for (auto it = im.begin(); it != im.end(); ++it) {
    cim(it.pos()) = ecast<To>(*it);
  }
  return cim;
}
template <class To, int N>
inline ImageOf<Vec<To, N>> VecCastPrivate(const ImageOf<Vec<To, N>> &v,
                                          std::true_type) {
  return v.clone();
}
}

template <class To, class From>
inline ImageOf<To> ecast(const ImageOf<From> &v) {
  return VecCastPrivate<To>(
      v, std::integral_constant<bool, std::is_same<To, From>::value>());
}

template <class To, class From, int N>
inline ImageOf<Vec<To, N>> ecast(const ImageOf<Vec<From, N>> &v) {
  return VecCastPrivate<To>(
      v, std::integral_constant<bool, std::is_same<To, From>::value>());
}

using Pixel = cv::Point;
template <class T> inline Vec<T, 2> ecast(const Pixel &p) {
  return Vec<T, 2>(static_cast<T>(p.x), static_cast<T>(p.y));
}

template <class T = Vec<uint8_t, 3>>
inline ImageOf<T> ImageRead(const std::string &filename) {
  return cv::imread(filename);
}
template <class T>
inline bool ImageWrite(const std::string &filename, const ImageOf<T> &im) {
  return cv::imwrite(filename, im);
}

template <class T = int> inline T Area(const Image &im) {
  return im.cols * im.rows;
}
template <class T = double> inline Point<T, 2> Center(const Image &im) {
  return Point<T, 2>(im.cols / 2.0, im.rows / 2.0);
}

void ClipToSquare(Image &im);
Imageb ClipToDisk(Image &im);
Imageb Rotate(Image &im, double angle);

template <class T> inline void ReverseCols(ImageOf<T> &im) {
  for (int i = 0; i < im.cols / 2; i++) {
    for (int k = 0; k < im.rows; k++) {
      std::swap(im(k, i), im(k, im.cols - 1 - i));
    }
  }
}

template <class T> inline void ReverseRows(ImageOf<T> &im) {
  for (int i = 0; i < im.rows / 2; i++) {
    for (int k = 0; k < im.cols; k++) {
      std::swap(im(i, k), im(im.rows - 1 - i, k));
    }
  }
}

void ResizeToWidth(Image &im, int width);
void ResizeToHeight(Image &im, int height);
void ResizeToMakeWidthUnder(Image &im, int widthUpperBound);
void ResizeToMakeHeightUnder(Image &im, int heightUpperBound);
bool MayBeAPanorama(const Image &im);
bool MakePanorama(Image &im, int horiCenter = -1, bool *extendedOnTop = nullptr,
                  bool *extendedOnBottom = nullptr);

std::pair<Pixel, Pixel> MinMaxLocOfImage(const Image &im);
std::pair<double, double> MinMaxValOfImage(const Image &im);
template <class T> inline T Mean(const ImageOf<T> &im, const Imageub &mask) {
  T sum = T();
  int count = 0;
  for (auto it = im.begin(); it != im.end(); ++it) {
    if (!mask.empty() && !mask(it.pos())) {
      continue;
    }
    sum += (*it);
    count++;
  }
  return sum / count;
}

template <class T> inline Pixel ToPixel(const Point<T, 2> &p) {
  return Pixel(static_cast<int>(p[0]), static_cast<int>(p[1]));
}

Pixel PixelFromGeoCoord(const GeoCoord &p, int longidiv, int latidiv);
GeoCoord GeoCoordFromPixel(const Pixel &pixel, int longidiv, int latidiv);

// dense mat
template <class T> using DenseMat = cv::Mat_<T>;
using DenseMati = DenseMat<int>;
using DenseMatd = DenseMat<double>;

// sparse mat
template <class T> using SparseMat = cv::SparseMat_<T>;
using SparseMatd = SparseMat<double>;
template <class T> struct SparseMatElement {
  using ValueType = T;
  int row, col;
  T value;
  inline SparseMatElement() : row(-1), col(-1) {}
  inline SparseMatElement(int r, int c, T v) : row(r), col(c), value(v) {}
};
using SparseMatElementd = SparseMatElement<double>;
template <class SparseMatElementIteratorT,
          class T = typename std::iterator_traits<
              SparseMatElementIteratorT>::value_type::ValueType>
inline SparseMat<T> MakeSparseMatFromElements(int row, int col,
                                              SparseMatElementIteratorT &&begin,
                                              SparseMatElementIteratorT &&end) {
  int dims[] = {row, col};
  SparseMat<T> mat(2, dims);
  while (begin != end) {
    mat.ref(begin->row, begin->col) = begin->value;
    ++begin;
  }
  return mat;
}
}
}
