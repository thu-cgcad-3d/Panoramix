#include "basic_types.hpp"

#include "../core/utility.hpp"

namespace pano {
namespace gui {

using namespace core;

Color::Color(ColorTag tag) {
  switch (tag) {
  case Transparent:
    _rgba = Vec4i(0, 0, 0, 0);
    break;

  case White:
    _rgba = Vec4i(255, 255, 255, 255);
    break;
  case Black:
    _rgba = Vec4i(0, 0, 0, 255);
    break;

  case DimGray:
    _rgba = Vec4i(105, 105, 105, 255);
    break;
  case Gray:
    _rgba = Vec4i(128, 128, 128, 255);
    break;
  case DarkGray:
    _rgba = Vec4i(169, 169, 169, 255);
    break;
  case Silver:
    _rgba = Vec4i(192, 192, 192, 255);
    break;
  case LightGray:
    _rgba = Vec4i(211, 211, 211, 255);
    break;

  case Red:
    _rgba = Vec4i(255, 0, 0, 255);
    break;
  case Green:
    _rgba = Vec4i(0, 255, 0, 255);
    break;
  case Blue:
    _rgba = Vec4i(0, 0, 255, 255);
    break;

  case Yellow:
    _rgba = Vec4i(255, 255, 0, 255);
    break;
  case Magenta:
    _rgba = Vec4i(255, 0, 255, 255);
    break;
  case Cyan:
    _rgba = Vec4i(0, 255, 255, 255);
    break;
  case Orange:
    _rgba = Vec4i(255, 165, 0, 255);
    break;
  default:
    _rgba = Vec4i(255, 255, 255, 255);
  }
}

Color::Color(const std::uint8_t *data, int cvType) {
  const int16_t *datai16 = (const int16_t *)(data);
  const int32_t *datai32 = (const int32_t *)(data);
  const float *dataf32 = (const float *)(data);
  const double *dataf64 = (const double *)(data);

  switch (cvType) {
  case CV_8UC1:
    _rgba = Vec4i(data[0], data[0], data[0], 255);
    break;
  case CV_8UC3:
    _rgba = Vec4i(data[0], data[1], data[2], 255);
    break;
  case CV_8UC4:
    _rgba = Vec4i(data[0], data[1], data[2], data[3]);
    break;

  case CV_16SC1:
    _rgba = Vec4i(datai16[0], datai16[0], datai16[0], 255);
    break;
  case CV_16SC3:
    _rgba = Vec4i(datai16[0], datai16[1], datai16[2], 255);
    break;
  case CV_16SC4:
    _rgba = Vec4i(datai16[0], datai16[1], datai16[2], datai16[3]);
    break;

  case CV_32SC1:
    _rgba = Vec4i(datai32[0], datai32[0], datai32[0], 255);
    break;
  case CV_32SC3:
    _rgba = Vec4i(datai32[0], datai32[1], datai32[2], 255);
    break;
  case CV_32SC4:
    _rgba = Vec4i(datai32[0], datai32[1], datai32[2], datai32[3]);
    break;

  case CV_32FC1:
    _rgba = Vec4i(dataf32[0], dataf32[0], dataf32[0], 1) * 255;
    break;
  case CV_32FC3:
    _rgba = Vec4i(dataf32[0], dataf32[1], dataf32[2], 1) * 255;
    break;
  case CV_32FC4:
    _rgba = Vec4i(dataf32[0], dataf32[1], dataf32[2], dataf32[3]) * 255;
    break;

  case CV_64FC1:
    _rgba = Vec4i(dataf64[0], dataf64[0], dataf64[0], 1) * 255;
    break;
  case CV_64FC3:
    _rgba = Vec4i(dataf64[0], dataf64[1], dataf64[2], 1) * 255;
    break;
  case CV_64FC4:
    _rgba = Vec4i(dataf64[0], dataf64[1], dataf64[2], dataf64[3]) * 255;
    break;
  default:
    std::cerr << "cannot convert this cv type to vis::Color [cvType = "
              << cvType << "]!" << std::endl;
    break;
  }
}

std::ostream &operator<<(std::ostream &os, ColorTag ct) {
  switch (ct) {
  case ColorTag::Transparent:
    os << "Transparent";
    break;

  case ColorTag::White:
    os << "White";
    break;
  case ColorTag::Black:
    os << "Black";
    break;

  case ColorTag::DimGray:
    os << "DimGray";
    break;
  case ColorTag::Gray:
    os << "Gray";
    break;
  case ColorTag::DarkGray:
    os << "DarkGray";
    break;
  case ColorTag::Silver:
    os << "Silver";
    break;
  case ColorTag::LightGray:
    os << "LightGray";
    break;

  case ColorTag::Red:
    os << "Red";
    break;
  case ColorTag::Green:
    os << "Green";
    break;
  case ColorTag::Blue:
    os << "Blue";
    break;

  case ColorTag::Yellow:
    os << "Yellow";
    break;
  case ColorTag::Magenta:
    os << "Magenta";
    break;
  case ColorTag::Cyan:
    os << "Cyan";
    break;
  case ColorTag::Orange:
    os << "Orange";
    break;
  default:
    os << "Unknown Color";
    break;
  }

  return os;
}

Color ColorFromHSV(double h, double s, double v, double a) {
  int i;
  double f, p, q, t;
  if (s == 0.0) {
    // achromatic (grey)
    return Color(v, v, v, a);
  }
  h *= 6.0; // sector 0 to 5
  i = floor(h);
  f = h - i; // factorial part of h
  p = v * (1 - s);
  q = v * (1 - s * f);
  t = v * (1 - s * (1 - f));
  switch (i) {
  case 0:
    return Color(v, t, p, a);
  case 1:
    return Color(q, v, p, a);
  case 2:
    return Color(p, v, t, a);
  case 3:
    return Color(p, q, v, a);
  case 4:
    return Color(t, p, v, a);
  default:
    return Color(v, p, q, a);
  }
}

Color RandomColor() { return Color(rand() % 255, rand() % 255, rand() % 255); }

namespace {

using ColorTableInit = std::pair<std::vector<Color>, Color>;
static const ColorTableInit allColorTable = {{White, Black,

                                              DimGray, Gray, DarkGray, Silver,
                                              LightGray,

                                              Red, Green, Blue,

                                              Yellow, Magenta, Cyan, Orange},
                                             Transparent};

static const ColorTableInit allColorExcludingWhiteTable = {
    {// ColorTag::White,
     Black,

     DimGray, Gray, DarkGray, Silver, LightGray,

     Red, Green, Blue,

     Yellow, Magenta, Cyan, Orange},
    Transparent};

static const ColorTableInit allColorExcludingBlackTable = {
    {White,
     // Black,

     DimGray, Gray, DarkGray, Silver, LightGray,

     Red, Green, Blue,

     Yellow, Magenta, Cyan, Orange},
    Transparent};

static const ColorTableInit allColorExcludingWhiteAndBlackTable = {
    {// White,
     // Black,

     DimGray, Gray, DarkGray, Silver, LightGray,

     Red, Green, Blue,

     Yellow, Magenta, Cyan, Orange},
    Transparent};

static const ColorTableInit RGBColorTable = {{Red, Green, Blue}, White};

static const ColorTableInit RGBGreysColorTable = {{Red, Green, Blue,

                                                   DimGray, Gray, DarkGray,
                                                   Silver, LightGray},
                                                  White};
}

ColorTable::ColorTable(ColorTableDescriptor descriptor) {
  const ColorTableInit *data = nullptr;
  switch (descriptor) {
  case ColorTableDescriptor::RGB:
    data = &RGBColorTable;
    break;
  case ColorTableDescriptor::AllColorsExcludingBlack:
    data = &allColorExcludingBlackTable;
    break;
  case ColorTableDescriptor::AllColorsExcludingWhite:
    data = &allColorExcludingWhiteTable;
    break;
  case ColorTableDescriptor::AllColorsExcludingWhiteAndBlack:
    data = &allColorExcludingWhiteAndBlackTable;
    break;
  case ColorTableDescriptor::AllColors:
    data = &allColorTable;
    break;
  default:
    data = &RGBGreysColorTable;
  }
  _colors = data->first;
  _exceptionalColor = data->second;
}

// ColorTable::ColorTable(std::initializer_list<ColorTag> ctags, ColorTag
// exceptColor) {
//    _colors.reserve(ctags.size());
//    for (auto ct : ctags) {
//        _colors.push_back(ct);
//    }
//    _exceptionalColor = exceptColor;
//}

core::Image3ub ColorTable::operator()(const core::Imagei &indexIm) const {
  core::Image3ub im(indexIm.size());
  for (auto i = indexIm.begin(); i != indexIm.end(); ++i) {
    core::Vec3ub color = (*this)[*i];
    im(i.pos()) = color;
  }
  return im;
}

ColorTable &ColorTable::randomize() {
  std::random_shuffle(_colors.begin(), _colors.end());
  return *this;
}

ColorTable &ColorTable::appendRandomizedColors(size_t sz) {
  int dimSplit = std::max(int(sqrt(sz)), 3);
  std::vector<Color> colors;
  colors.reserve(dimSplit * dimSplit * dimSplit - dimSplit);
  for (int i = 0; i < dimSplit; i++) {
    for (int j = 0; j < dimSplit; j++) {
      for (int k = 0; k < dimSplit; k++) {
        if (i == j && j == k)
          continue;
        colors.push_back(
            Color(1.0 * i / dimSplit, 1.0 * j / dimSplit, 1.0 * k / dimSplit));
      }
    }
  }
  assert(colors.size() > sz);
  std::random_shuffle(colors.begin(), colors.end());
  _colors.insert(_colors.end(), colors.begin(), colors.begin() + sz);
  return *this;
}

ColorTable &ColorTable::appendRandomizedGreyColors(size_t size) {
  core::Vec3 full(255, 255, 255);
  std::vector<Color> colors(size);
  for (int i = 0; i < size; i++) {
    colors[i] = Color(double(i) * full / double(size));
  }
  std::random_shuffle(colors.begin(), colors.end());
  _colors.insert(_colors.end(), colors.begin(), colors.end());
  return *this;
}

ColorTable CreateGreyColorTableWithSize(int sz, const Color &exceptColor) {
  core::Vec3 full(1, 1, 1);
  std::vector<Color> colors(sz);
  for (int i = 0; i < sz; i++) {
    colors[i] = Color(double(i) * full / double(sz));
  }
  return ColorTable(colors, exceptColor);
}

ColorTable CreateRandomColorTableWithSize(int sz, const Color &exceptColor) {
  int dimSplit = std::max(int(sqrt(sz)), 3);
  std::vector<Color> colors;
  colors.reserve(dimSplit * dimSplit * dimSplit - dimSplit);
  for (int i = 0; i < dimSplit; i++) {
    for (int j = 0; j < dimSplit; j++) {
      for (int k = 0; k < dimSplit; k++) {
        if (i == j && j == k)
          continue;
        colors.push_back(
            Color(1.0 * i / dimSplit, 1.0 * j / dimSplit, 1.0 * k / dimSplit));
      }
    }
  }
  assert(colors.size() > sz);
  std::random_shuffle(colors.begin(), colors.end());
  return ColorTable(colors.begin(), colors.begin() + sz, exceptColor);
}

namespace {
double interpolate(double val, double y0, double x0, double y1, double x1) {
  return (val - x0) * (y1 - y0) / (x1 - x0) + y0;
}
double base(double val) {
  if (val <= -0.75)
    return 0;
  else if (val <= -0.25)
    return interpolate(val, 0.0, -0.75, 1.0, -0.25);
  else if (val <= 0.25)
    return 1.0;
  else if (val <= 0.75)
    return interpolate(val, 1.0, 0.25, 0.0, 0.75);
  else
    return 0.0;
}
double red(double gray) { return base(gray - 0.5); }
double green(double gray) { return base(gray); }
double blue(double gray) { return base(gray + 0.5); }
}

ColorTable CreateJetColorTableWithSize(int sz, const Color &exceptColor) {
  std::vector<Color> colors(sz);
  for (int i = 0; i < sz; i++) {
    double val = i / double(sz) * 2.0 - 1.0;
    colors[i] = Color(red(val), green(val), blue(val));
  }
  return ColorTable(std::move(colors), exceptColor);
}

// opengl shader source
const std::pair<std::string, std::string> &
PredefinedShaderSource(OpenGLShaderSourceDescriptor name) {

  static const std::pair<std::string, std::string> defaultPointsShaderSource = {
      "#version 120\n"
      "attribute highp vec4 position;\n"
      "attribute highp vec3 normal;\n"
      "attribute lowp vec4 color;\n"
      "attribute lowp vec2 texCoord;\n"
      "uniform highp mat4 matrix;\n"
      "uniform float pointSize;\n"
      "varying vec4 pixelColor;\n"
      "void main(void)\n"
      "{\n"
      "    gl_Position = matrix * position;\n"
      "    gl_PointSize = pointSize;\n"
      "    pixelColor = color;\n"
      "}\n",

      "#version 120\n"
      "varying lowp vec4 pixelColor;\n"
      "void main(void)\n"
      "{\n"
      "   gl_FragColor = pixelColor;\n"
      "   float distance = length(gl_PointCoord - vec2(0.5));\n"
      "   if(distance > 0.4 && distance <= 0.5)\n"
      "       gl_FragColor.a = 1.0 - (distance - 0.4) * 0.1;\n"
      "   else if(distance > 0.5)\n"
      "       discard;\n"
      "}\n"};

  static const std::pair<std::string, std::string> defaultLinesShaderSource = {
      "#version 120\n"
      "attribute lowp vec4 position;\n"
      "attribute lowp vec3 normal;\n"
      "attribute lowp vec4 color;\n"
      "attribute lowp vec2 texCoord;\n"
      "uniform highp mat4 matrix;\n"
      "uniform float pointSize;\n"
      "varying vec4 pixelColor;\n"
      "void main(void)\n"
      "{\n"
      "    gl_Position = matrix * position;\n"
      "    pixelColor = color;\n"
      "}\n",

      "#version 120\n"
      "varying lowp vec4 pixelColor;\n"
      "void main(void)\n"
      "{\n"
      "    gl_FragColor = pixelColor;\n"
      "}\n"};

  static const std::pair<std::string, std::string>
      defaultTrianglesShaderSource = {
          "#version 120\n"
          "attribute highp vec4 position;\n"
          "attribute highp vec3 normal;\n"
          "attribute lowp vec4 color;\n"
          "attribute lowp vec2 texCoord;\n"
          "uniform highp mat4 matrix;\n"
          "uniform float pointSize;\n"
          "varying vec4 pixelColor;\n"
          "void main(void)\n"
          "{\n"
          "    gl_Position = matrix * position;\n"
          "    highp vec4 transformedNormal = viewMatrix * modelMatrix * "
          "vec4(normal, 1.0);\n"
          "    highp vec3 transformedNormal3 = transformedNormal.xyz / "
          "transformedNormal.w;\n"
          "    pixelColor = abs(dot(transformedNormal3 / "
          "length(transformedNormal), vec3(1.0, 0.0, 0.0))) * vec4(1.0, 1.0, "
          "1.0, 1.0);\n"
          "}\n",

          "#version 120\n"
          "varying lowp vec4 pixelColor;\n"
          "void main(void)\n"
          "{\n"
          "    gl_FragColor = pixelColor;\n"
          "}\n"};

  static const std::pair<std::string, std::string> panoramaShaderSource = {
      "#version 120\n"
      "attribute highp vec3 position;\n"
      "attribute highp vec3 normal;\n"
      "attribute highp vec4 color;\n"
      "uniform highp mat4 matrix;\n"
      "varying highp vec3 pixelPosition;\n"
      "varying highp vec3 pixelNormal;\n"
      "varying highp vec4 pixelColor;\n"
      "void main(void)\n"
      "{\n"
      "    pixelPosition = position.xyz;\n"
      "    pixelNormal = normal;\n"
      "    pixelColor = color;\n"
      "    gl_Position = matrix * vec4(position, 1.0);\n"
      "}\n",

      // 3.14159265358979323846264338327950288
      "uniform sampler2D tex;\n"
      "uniform highp vec3 panoramaCenter;\n"
      "varying highp vec3 pixelPosition;\n"
      "varying highp vec3 pixelNormal;\n"
      "varying highp vec4 pixelColor;\n"
      "void main(void)\n"
      "{\n"
      "    highp vec3 direction = pixelPosition - panoramaCenter;\n"
      "    highp float longi = atan(direction.y, direction.x);\n"
      "    highp float lati = asin(direction.z / length(direction));\n"
      "    highp vec2 texCoord = vec2(longi / 3.1415926535897932 / 2.0 + 0.5, "
      "- lati / 3.1415926535897932 + 0.5);\n"
      "    gl_FragColor = texture2D(tex, texCoord) * 1.0 + pixelColor * 0.0;\n"
      "    gl_FragColor.a = 0.7;\n"
      "}\n"};

  static const std::pair<std::string, std::string> xPointsShaderSource = {
      "#version 130\n"

      "attribute highp vec4 position;\n"
      "attribute highp vec3 normal;\n"
      "attribute lowp vec4 color;\n"
      "attribute lowp vec2 texCoord;\n"
      "attribute uint isSelected;\n"

      "uniform highp mat4 viewMatrix;\n"
      "uniform highp mat4 modelMatrix;\n"
      "uniform highp mat4 projectionMatrix;\n"

      "varying lowp vec4 pixelColor;\n"
      "varying lowp vec2 pixelTexCoord;\n"
      "varying lowp float pixelSelection;\n"

      "void main(void)\n"
      "{\n"
      "    gl_Position = projectionMatrix * viewMatrix * modelMatrix * "
      "position;\n"
      "    pixelColor = color;\n"
      "    pixelTexCoord = texCoord;\n"
      "    pixelSelection = isSelected == 0u ? 0.0 : 1.0;\n"
      "}\n",

      "#version 130\n"

      "uniform sampler2D tex;\n"
      "uniform lowp vec4 globalColor;\n"

      "uniform lowp float bwColor;\n"
      "uniform lowp float bwTexColor;\n"

      "varying lowp vec4 pixelColor;\n"
      "varying lowp vec2 pixelTexCoord;\n"
      "varying lowp float pixelSelection;\n"

      "void main(void)\n"
      "{\n"
      "   lowp vec4 texColor = texture2D(tex, pixelTexCoord);\n"
      //"   gl_FragColor = (pixelColor * bwColor + texColor * bwTexColor)"
      //        "       / (bwColor + bwTexColor);\n"
      "   gl_FragColor = (pixelColor * 1.0 + texColor * 0.0);\n"
      "   gl_FragColor.a = 1.0 - pixelSelection * 0.5;\n"
      "   float distance = length(gl_PointCoord - vec2(0.5));\n"
      "   if(distance > 0.4 && distance <= 0.5)\n"
      "       gl_FragColor.a = (1.0 - (distance - 0.4) * 0.1) * "
      "gl_FragColor.a;\n"
      "   else if(distance > 0.5)\n"
      "       discard;\n"
      "}\n"};

  static const std::pair<std::string, std::string> xLinesShaderSource = {
      "#version 130\n"

      "attribute highp vec4 position;\n"
      "attribute highp vec3 normal;\n"
      "attribute lowp vec4 color;\n"
      "attribute lowp vec2 texCoord;\n"
      "attribute uint isSelected;\n"

      "uniform highp mat4 viewMatrix;\n"
      "uniform highp mat4 modelMatrix;\n"
      "uniform highp mat4 projectionMatrix;\n"

      "varying lowp vec4 pixelColor;\n"
      "varying lowp vec2 pixelTexCoord;\n"
      "varying lowp float pixelSelection;\n"

      "void main(void)\n"
      "{\n"
      "    gl_Position = projectionMatrix * viewMatrix * modelMatrix * "
      "position;\n"
      "    pixelColor = color;\n"
      "    pixelTexCoord = texCoord;\n"
      "    pixelSelection = isSelected == 0u ? 0.0 : 1.0;\n"
      "}\n",

      "#version 130\n"

      "uniform sampler2D tex;\n"
      "uniform lowp vec4 globalColor;\n"

      "uniform lowp float bwColor;\n"
      "uniform lowp float bwTexColor;\n"

      "varying lowp vec4 pixelColor;\n"
      "varying lowp vec2 pixelTexCoord;\n"
      "varying lowp float pixelSelection;\n"

      "void main(void)\n"
      "{\n"
      "   lowp vec4 texColor = texture2D(tex, pixelTexCoord);\n"
      "   gl_FragColor = (pixelColor * 1.0 + texColor * 0.0)"
      "     ;\n"
      "   gl_FragColor.a = 1.0 - pixelSelection * 0.5;\n"
      "}\n"};

  static const std::pair<std::string, std::string> xTrianglesShaderSource = {
      "#version 130\n"

      "attribute highp vec4 position;\n"
      "attribute highp vec3 normal;\n"
      "attribute lowp vec4 color;\n"
      "attribute lowp vec2 texCoord;\n"
      "attribute uint isSelected;\n"

      "uniform highp mat4 viewMatrix;\n"
      "uniform highp mat4 modelMatrix;\n"
      "uniform highp mat4 projectionMatrix;\n"

      "varying lowp vec4 pixelColor;\n"
      "varying lowp vec2 pixelTexCoord;\n"
      "varying lowp float pixelSelection;\n"

      "void main(void)\n"
      "{\n"
      "    gl_Position = projectionMatrix * viewMatrix * modelMatrix * "
      "position;\n"
      "    pixelColor = color;\n"
      "    pixelTexCoord = texCoord;\n"
      "    pixelSelection = isSelected == 0u ? 0.0 : 1.0;\n"
      "}\n",

      "#version 130\n"

      "uniform sampler2D tex;\n"
      "uniform lowp vec4 globalColor;\n"

      "uniform lowp float bwColor;\n"
      "uniform lowp float bwTexColor;\n"

      "varying lowp vec4 pixelColor;\n"
      "varying lowp vec2 pixelTexCoord;\n"
      "varying lowp float pixelSelection;\n"

      "void main(void)\n"
      "{\n"
      "    lowp vec4 texColor = texture2D(tex, pixelTexCoord);\n"
      "    gl_FragColor = (pixelColor * bwColor + texColor * bwTexColor)"
      "       / (bwColor + bwTexColor);\n"
      "    gl_FragColor.a = 1.0 - pixelSelection * 0.5;\n"
      "}\n"};

  static const std::pair<std::string, std::string> xPanoramaShaderSource = {
      "#version 130\n"

      "attribute highp vec3 position;\n"
      "attribute highp vec3 normal;\n"
      "attribute highp vec4 color;\n"
      "attribute lowp vec2 texCoord;\n"
      "attribute uint isSelected;\n"

      "uniform highp mat4 viewMatrix;\n"
      "uniform highp mat4 modelMatrix;\n"
      "uniform highp mat4 projectionMatrix;\n"

      "varying highp vec3 pixelPosition;\n"
      "varying highp vec3 pixelNormal;\n"
      "varying highp vec4 pixelColor;\n"
      "varying lowp float pixelSelection;\n"

      "void main(void)\n"
      "{\n"
      "    pixelPosition = position.xyz;\n"
      "    pixelNormal = normal;\n"
      "    pixelColor = color;\n"
      "    gl_Position = projectionMatrix * viewMatrix * modelMatrix * "
      "vec4(position, 1.0);\n"
      "    pixelSelection = isSelected == 0u ? 0.0 : 1.0;\n"
      "}\n",

      // 3.14159265358979323846264338327950288
      "#version 130\n"

      "uniform sampler2D tex;\n"

      "uniform lowp float bwColor;\n"
      "uniform lowp float bwTexColor;\n"
      "uniform highp vec3 panoramaCenter;\n"
      "uniform highp float panoramaHoriCenterRatio;\n"
      "uniform highp float panoramaAspectRatio;\n" // height/width

      "varying highp vec3 pixelPosition;\n"
      "varying highp vec3 pixelNormal;\n"
      "varying highp vec4 pixelColor;\n"
      "varying lowp float pixelSelection;\n"

      "void main(void)\n"
      "{\n"
      "    highp vec3 direction = pixelPosition - panoramaCenter;\n"
      "    highp float longi = atan(direction.y, direction.x);\n"
      "    highp float lati = asin(direction.z / length(direction));\n"
      "    highp float globalx = longi / 3.1415926535897932 / 2.0 + 0.5;\n"
      "    highp float globaly = - lati / 3.1415926535897932 + 0.5;\n"
      //"    highp vec2 texCoord = vec2(longi / 3.1415926535897932 / 2.0 + 0.5,
      //- lati / 3.1415926535897932 + 0.5);\n"
      "    highp float texx = globalx;\n"
      "    highp float texy =  (globaly / 2.0 - (0.25 - "
      "panoramaHoriCenterRatio * panoramaAspectRatio)) / panoramaAspectRatio;\n"
      "    if(texy < 0.0 || texy > 1.0){ discard; } \n"
      "    highp vec2 texCoord = vec2(texx, texy);\n"
      "    lowp vec4 texColor = texture2D(tex, texCoord);\n"
      "    gl_FragColor = (pixelColor * bwColor + texColor * bwTexColor)"
      "       / (bwColor + bwTexColor);\n"
      "    gl_FragColor.a = 1.0 - pixelSelection * 0.5;\n"
      "}\n"};

  switch (name) {
  case OpenGLShaderSourceDescriptor::DefaultPoints:
    return defaultPointsShaderSource;
  case OpenGLShaderSourceDescriptor::DefaultLines:
    return defaultLinesShaderSource;
  case OpenGLShaderSourceDescriptor::DefaultTriangles:
    return defaultTrianglesShaderSource;
  case OpenGLShaderSourceDescriptor::Panorama:
    return panoramaShaderSource;

  case OpenGLShaderSourceDescriptor::XPoints:
    return xPointsShaderSource;
  case OpenGLShaderSourceDescriptor::XLines:
    return xLinesShaderSource;
  case OpenGLShaderSourceDescriptor::XTriangles:
    return xTrianglesShaderSource;
  case OpenGLShaderSourceDescriptor::XPanorama:
    return xPanoramaShaderSource;
  default:
    return defaultTrianglesShaderSource;
  }
}

OpenGLShaderSource::OpenGLShaderSource(OpenGLShaderSourceDescriptor d) {
  std::tie(_vshaderSrc, _fshaderSrc) = PredefinedShaderSource(d);
}
}

namespace core {
Box3 BoundingBox(const gui::SpatialProjectedPolygon &spp) {
  std::vector<Vec3> cs(spp.corners.size());
  for (int i = 0; i < spp.corners.size(); i++) {
    cs[i] =
        IntersectionOfLineAndPlane(
            Ray3(spp.projectionCenter, spp.corners[i] - spp.projectionCenter),
            spp.plane)
            .position;
  }
  return BoundingBoxOfContainer(cs);
}
}
}