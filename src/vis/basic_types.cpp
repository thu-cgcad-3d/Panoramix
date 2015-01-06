#include "basic_types.hpp"

#include "../core/utilities.hpp"
#include "../core/algorithms.hpp"

namespace panoramix {
    namespace vis {

        using namespace core;

        const std::vector<ColorTag> & AllColorTags() {
            static const std::vector<ColorTag> _allColorTags = {
                ColorTag::Transparent,
                ColorTag::White,
                ColorTag::Gray,
                ColorTag::Red,
                ColorTag::Green,
                ColorTag::Blue,
                ColorTag::Yellow,
                ColorTag::Magenta,
                ColorTag::Cyan,
                ColorTag::Orange
            };
            return _allColorTags;
        }

        std::ostream & operator << (std::ostream & os, ColorTag ct) {
            switch (ct){
            case ColorTag::Transparent: os << "Transparent"; break;

            case ColorTag::White: os << "White"; break;
            case ColorTag::Black: os << "Black"; break;

            case ColorTag::DimGray: os << "DimGray"; break;
            case ColorTag::Gray: os << "Gray"; break;
            case ColorTag::DarkGray: os << "DarkGray"; break;
            case ColorTag::Silver: os << "Silver"; break;
            case ColorTag::LightGray: os << "LightGray"; break;

            case ColorTag::Red: os << "Red"; break;
            case ColorTag::Green: os << "Green"; break;
            case ColorTag::Blue: os << "Blue"; break;

            case ColorTag::Yellow: os << "Yellow"; break;
            case ColorTag::Magenta: os << "Magenta"; break;
            case ColorTag::Cyan: os << "Cyan"; break;
            case ColorTag::Orange: os << "Orange"; break;
            default:
                os << "Unknown Color"; break;
            }

            return os;
        }

        Color ColorFromRGB(double r, double g, double b, double a) {
            return Color(b, g, r, a);
        }

        Color ColorFromHSV(double h, double s, double v, double a) {
            int i;
            double f, p, q, t;
            if (s == 0) {
                // achromatic (grey)
                return ColorFromRGB(v * 255, v * 255, v * 255, a);
            }
            h *= 360.0;
            h /= 60;			// sector 0 to 5
            i = floor(h);
            f = h - i;			// factorial part of h
            p = v * (1 - s);
            q = v * (1 - s * f);
            t = v * (1 - s * (1 - f));
            switch (i) {
            case 0: return ColorFromRGB(v*255, t*255, p*255, a);
            case 1: return ColorFromRGB(q*255, v*255, p*255, a);
            case 2: return ColorFromRGB(p*255, v*255, t*255, a);
            case 3: return ColorFromRGB(p*255, q*255, v*255, a);
            case 4: return ColorFromRGB(t*255, p*255, v*255, a);
            default:return ColorFromRGB(v*255, p*255, q*255, a);
            }
        }

        Color ColorFromTag(ColorTag t) {
            switch (t){
            case ColorTag::Transparent: return ColorFromRGB(0, 0, 0, 0);

            case ColorTag::White: return ColorFromRGB(255, 255, 255);
            case ColorTag::Black: return ColorFromRGB(0, 0, 0);

            case ColorTag::DimGray: return ColorFromRGB(105, 105, 105);
            case ColorTag::Gray: return ColorFromRGB(128, 128, 128);
            case ColorTag::DarkGray: return ColorFromRGB(169, 169, 169);
            case ColorTag::Silver: return ColorFromRGB(192, 192, 192);
            case ColorTag::LightGray: return ColorFromRGB(211, 211, 211);

            case ColorTag::Red: return ColorFromRGB(255, 0, 0);
            case ColorTag::Green: return ColorFromRGB(0, 255, 0);
            case ColorTag::Blue: return ColorFromRGB(0, 0, 255);

            case ColorTag::Yellow: return ColorFromRGB(255, 255, 0);
            case ColorTag::Magenta: return ColorFromRGB(255, 0, 255);
            case ColorTag::Cyan: return ColorFromRGB(0, 255, 255);
            case ColorTag::Orange: return ColorFromRGB(255, 165, 0);
            default:
                return Color(255, 255, 255);
            }
        }

        Color RandomColor() {
            return Color(rand() % 255, rand() % 255, rand() % 255);
        }

        core::Vec3b ToVec3b(const Color & c) {
            return core::Vec3b(static_cast<uint8_t>(c[0]), static_cast<uint8_t>(c[1]), static_cast<uint8_t>(c[2]));
        }


        ColorTable::ColorTable(ColorTableDescriptor descriptor) {
            const auto & predefined = PredefinedColorTable(descriptor);
            _colors = predefined._colors;
            _exceptionalColor = predefined._exceptionalColor;
        }
        ColorTable::ColorTable(std::initializer_list<ColorTag> ctags, ColorTag exceptColor) {
            _colors.reserve(ctags.size());
            for (auto ct : ctags) {
                _colors.push_back(ColorFromTag(ct));
            }
            _exceptionalColor = ColorFromTag(exceptColor);
        }

        void ColorTable::randomize() { std::random_shuffle(_colors.begin(), _colors.end()); }

        const ColorTable & PredefinedColorTable(ColorTableDescriptor descriptor) {
           
            static const ColorTable allColorTable = {
                {
                    ColorTag::White,
                    ColorTag::Black,

                    ColorTag::DimGray,
                    ColorTag::Gray,
                    ColorTag::DarkGray,
                    ColorTag::Silver,
                    ColorTag::LightGray,

                    ColorTag::Red,
                    ColorTag::Green,
                    ColorTag::Blue,

                    ColorTag::Yellow,
                    ColorTag::Magenta,
                    ColorTag::Cyan,
                    ColorTag::Orange
                }, 
                ColorTag::Transparent
            };

            static const ColorTable allColorExcludingWhiteTable = {
                {
                    //ColorTag::White,
                    ColorTag::Black,

                    ColorTag::DimGray,
                    ColorTag::Gray,
                    ColorTag::DarkGray,
                    ColorTag::Silver,
                    ColorTag::LightGray,

                    ColorTag::Red,
                    ColorTag::Green,
                    ColorTag::Blue,

                    ColorTag::Yellow,
                    ColorTag::Magenta,
                    ColorTag::Cyan,
                    ColorTag::Orange
                },
                ColorTag::Transparent
            };

            static const ColorTable allColorExcludingBlackTable = {
                {
                    ColorTag::White,
                    //ColorTag::Black,

                    ColorTag::DimGray,
                    ColorTag::Gray,
                    ColorTag::DarkGray,
                    ColorTag::Silver,
                    ColorTag::LightGray,

                    ColorTag::Red,
                    ColorTag::Green,
                    ColorTag::Blue,

                    ColorTag::Yellow,
                    ColorTag::Magenta,
                    ColorTag::Cyan,
                    ColorTag::Orange
                },
                ColorTag::Transparent
            };

            static const ColorTable allColorExcludingWhiteAndBlackTable = {
                {
                    //ColorTag::White,
                    //ColorTag::Black,

                    ColorTag::DimGray,
                    ColorTag::Gray,
                    ColorTag::DarkGray,
                    ColorTag::Silver,
                    ColorTag::LightGray,

                    ColorTag::Red,
                    ColorTag::Green,
                    ColorTag::Blue,

                    ColorTag::Yellow,
                    ColorTag::Magenta,
                    ColorTag::Cyan,
                    ColorTag::Orange
                },
                ColorTag::Transparent
            };

            static const ColorTable RGBColorTable = {
                {
                    ColorTag::Red,
                    ColorTag::Green,
                    ColorTag::Blue
                },
                ColorTag::White
            };

            switch (descriptor) {
            case ColorTableDescriptor::RGB: return RGBColorTable;
            case ColorTableDescriptor::AllColorsExcludingBlack: return allColorExcludingBlackTable;
            case ColorTableDescriptor::AllColorsExcludingWhite: return allColorExcludingWhiteTable;
            case ColorTableDescriptor::AllColorsExcludingWhiteAndBlack: return allColorExcludingWhiteAndBlackTable;
            case ColorTableDescriptor::AllColors: return allColorTable;
            default: return allColorTable;
            }
            
        }

        ColorTable CreateGreyColorTableWithSize(int sz) {
            auto exeptColor = ColorFromTag(ColorTag::Blue);
            Color full(255, 255, 255);
            std::vector<Color> colors(sz);
            for (int i = 0; i < sz; i++){
                colors[i] = double(i) * full / double(sz);
            }
            return ColorTable(colors, exeptColor);
        }

        ColorTable CreateRandomColorTableWithSize(int sz, const Color & exceptColor) {
            int dimSplit = std::max(int(sqrt(sz)), 3);
            std::vector<Color> colors;
            colors.reserve(dimSplit * dimSplit * dimSplit - dimSplit);
            for (int i = 0; i < dimSplit; i++){
                for (int j = 0; j < dimSplit; j++){
                    for (int k = 0; k < dimSplit; k++){
                        if (i == j && j == k)
                            continue;
                        colors.push_back(Color(255.0 * i / dimSplit, 255.0 * j / dimSplit, 255.0 * k / dimSplit));
                    }
                }
            }
            assert(colors.size() > sz);
            std::random_shuffle(colors.begin(), colors.end());
            return ColorTable(colors.begin(), colors.begin() + sz, exceptColor);
        }


        OpenGLShaderSource::OpenGLShaderSource(OpenGLShaderSourceDescriptor d) {
            auto & ss = PredefinedShaderSource(d);
            _vshaderSrc = ss._vshaderSrc;
            _fshaderSrc = ss._fshaderSrc;
        }


        // opengl shader source 
        const OpenGLShaderSource & PredefinedShaderSource(OpenGLShaderSourceDescriptor name) {

            static const OpenGLShaderSource defaultPointsShaderSource = {
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
                "}\n"
            };

            static const OpenGLShaderSource defaultLinesShaderSource = {
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
                "}\n"
            };

            static const OpenGLShaderSource defaultTrianglesShaderSource = {
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
                "    highp vec4 transformedNormal = viewMatrix * modelMatrix * vec4(normal, 1.0);\n"
                "    highp vec3 transformedNormal3 = transformedNormal.xyz / transformedNormal.w;\n"
                "    pixelColor = abs(dot(transformedNormal3 / length(transformedNormal), vec3(1.0, 0.0, 0.0))) * vec4(1.0, 1.0, 1.0, 1.0);\n"
                "}\n",

                "#version 120\n"
                "varying lowp vec4 pixelColor;\n"
                "void main(void)\n"
                "{\n"
                "    gl_FragColor = pixelColor;\n"
                "}\n"
            };

            static const OpenGLShaderSource panoramaShaderSource = {
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
                "}\n"
                ,

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
                "    highp vec2 texCoord = vec2(longi / 3.1415926535897932 / 2.0 + 0.5, - lati / 3.1415926535897932 + 0.5);\n"
                "    gl_FragColor = texture2D(tex, texCoord) * 1.0 + pixelColor * 0.0;\n"
                "    gl_FragColor.a = 0.7;\n"
                "}\n"
            };








            static const OpenGLShaderSource xPointsShaderSource = {
                "#version 130\n"

                "attribute highp vec4 position;\n"
                "attribute highp vec3 normal;\n"
                "attribute lowp vec4 color;\n"
                "attribute lowp vec2 texCoord;\n"

                "uniform highp mat4 viewMatrix;\n"
                "uniform highp mat4 modelMatrix;\n"
                "uniform highp mat4 projectionMatrix;\n"

                "varying lowp vec4 pixelColor;\n"
                "varying lowp vec2 pixelTexCoord;\n"

                "void main(void)\n"
                "{\n"
                "    gl_Position = projectionMatrix * viewMatrix * modelMatrix * position;\n"
                "    pixelColor = color;\n"
                "    pixelTexCoord = texCoord;\n"
                "}\n",



                "#version 130\n"

                "uniform sampler2D tex;\n"
                "uniform lowp vec4 globalColor;\n"

                "uniform lowp float bwColor;\n"
                "uniform lowp float bwTexColor;\n"

                "varying lowp vec4 pixelColor;\n"
                "varying lowp vec2 pixelTexCoord;\n"

                "void main(void)\n"
                "{\n"
                "   lowp vec4 texColor = texture2D(tex, pixelTexCoord);\n"
                "   gl_FragColor = (pixelColor * bwColor + texColor * bwTexColor)" 
                        "       / (bwColor + bwTexColor);\n"
                "   float distance = length(gl_PointCoord - vec2(0.5));\n"
                "   if(distance > 0.4 && distance <= 0.5)\n"
                "       gl_FragColor.a = 1.0 - (distance - 0.4) * 0.1;\n"
                "   else if(distance > 0.5)\n"
                "       discard;\n"
                "}\n"
            };

            static const OpenGLShaderSource xLinesShaderSource = {
                "#version 130\n"

                "attribute highp vec4 position;\n"
                "attribute highp vec3 normal;\n"
                "attribute lowp vec4 color;\n"
                "attribute lowp vec2 texCoord;\n"

                "uniform highp mat4 viewMatrix;\n"
                "uniform highp mat4 modelMatrix;\n"
                "uniform highp mat4 projectionMatrix;\n"

                "varying lowp vec4 pixelColor;\n"
                "varying lowp vec2 pixelTexCoord;\n"

                "void main(void)\n"
                "{\n"
                "    gl_Position = projectionMatrix * viewMatrix * modelMatrix * position;\n"
                "    pixelColor = color;\n"
                "    pixelTexCoord = texCoord;\n"
                "}\n",



                "#version 130\n"

                "uniform sampler2D tex;\n"
                "uniform lowp vec4 globalColor;\n"

                "uniform lowp float bwColor;\n"
                "uniform lowp float bwTexColor;\n"

                "varying lowp vec4 pixelColor;\n"
                "varying lowp vec2 pixelTexCoord;\n"

                "void main(void)\n"
                "{\n"
                "   lowp vec4 texColor = texture2D(tex, pixelTexCoord);\n"
                "   gl_FragColor = (pixelColor * bwColor + texColor * bwTexColor)"
                "       / (bwColor + bwTexColor);\n"
                "}\n"
            };

            static const OpenGLShaderSource xTrianglesShaderSource = {
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
                "    gl_Position = projectionMatrix * viewMatrix * modelMatrix * position;\n"
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
                "}\n"
            };

            static const OpenGLShaderSource xPanoramaShaderSource = {
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
                "    gl_Position = projectionMatrix * viewMatrix * modelMatrix * vec4(position, 1.0);\n"
                "    pixelSelection = isSelected == 0u ? 0.0 : 1.0;\n"
                "}\n"
                ,

                // 3.14159265358979323846264338327950288
                "#version 130\n"

                "uniform sampler2D tex;\n"

                "uniform lowp float bwColor;\n"
                "uniform lowp float bwTexColor;\n"
                "uniform highp vec3 panoramaCenter;\n"

                "varying highp vec3 pixelPosition;\n"
                "varying highp vec3 pixelNormal;\n"
                "varying highp vec4 pixelColor;\n"
                "varying lowp float pixelSelection;\n"

                "void main(void)\n"
                "{\n"
                "    highp vec3 direction = pixelPosition - panoramaCenter;\n"
                "    highp float longi = atan(direction.y, direction.x);\n"
                "    highp float lati = asin(direction.z / length(direction));\n"
                "    highp vec2 texCoord = vec2(longi / 3.1415926535897932 / 2.0 + 0.5, - lati / 3.1415926535897932 + 0.5);\n"
                "    lowp vec4 texColor = texture2D(tex, texCoord);\n"
                "    gl_FragColor = (pixelColor * bwColor + texColor * bwTexColor)"
                "       / (bwColor + bwTexColor);\n"
                "    gl_FragColor.a = 1.0 - pixelSelection * 0.5;\n"
                "}\n"
            };



            switch (name) {
            case OpenGLShaderSourceDescriptor::DefaultPoints: return defaultPointsShaderSource;
            case OpenGLShaderSourceDescriptor::DefaultLines: return defaultLinesShaderSource;
            case OpenGLShaderSourceDescriptor::DefaultTriangles: return defaultTrianglesShaderSource;
            case OpenGLShaderSourceDescriptor::Panorama: return panoramaShaderSource;

            case OpenGLShaderSourceDescriptor::XPoints: return xPointsShaderSource;
            case OpenGLShaderSourceDescriptor::XLines: return xLinesShaderSource;
            case OpenGLShaderSourceDescriptor::XTriangles: return xTrianglesShaderSource;
            case OpenGLShaderSourceDescriptor::XPanorama: return xPanoramaShaderSource;
            default:
                return defaultTrianglesShaderSource;
            }
        }




       

        namespace {

            inline Vec3 ToVec3Affine(const Vec4 & v4) {
                return Vec3(v4[0], v4[1], v4[2]) / v4[3];
            }

            inline Vec2 ToVec2(const Vec3 & v3) {
                return Vec2(v3[0], v3[1]);
            }

        }



        // tri mesh implementation
        TriMesh::Vertex::Vertex()
            : position(0, 0, 0, 1), normal(0, 0, 0), color(0, 0, 0, 1), texCoord(0, 0), entityIndex(-1), isSelected(false) {
        }


        TriMesh::VertHandle TriMesh::addVertex(const TriMesh::Vertex & v) {
            vertices.push_back(v);
            iPoints.push_back(static_cast<TriMesh::VertHandle>(vertices.size() - 1));
            return iPoints.back();
        }

        TriMesh::LineHandle TriMesh::addLine(TriMesh::VertHandle v1, TriMesh::VertHandle v2) {
            iLines.push_back(v1);
            iLines.push_back(v2);
            return iLines.size() / 2;
        }

        TriMesh::LineHandle TriMesh::addIsolatedLine(const Vertex & v1, const Vertex & v2) {
            vertices.push_back(v1);
            iLines.push_back(vertices.size() - 1);
            vertices.push_back(v2);
            iLines.push_back(vertices.size() - 1);
            return iLines.size() / 2;
        }

        size_t TriMesh::numberOfLines() const {
            return iLines.size() / 2;
        }

        void TriMesh::fetchLineVerts(TriMesh::LineHandle l, TriMesh::VertHandle & v1, TriMesh::VertHandle & v2) const {
            v1 = iLines[l * 2];
            v2 = iLines[l * 2 + 1];
        }

        TriMesh::TriangleHandle TriMesh::addTriangle(TriMesh::VertHandle v1, TriMesh::VertHandle v2, TriMesh::VertHandle v3) {
            iTriangles.push_back(v1);
            iTriangles.push_back(v2);
            iTriangles.push_back(v3);
            return iTriangles.size() / 3;
        }

        TriMesh::TriangleHandle TriMesh::addIsolatedTriangle(const Vertex & v1, const Vertex & v2, const Vertex & v3) {
            vertices.push_back(v1);
            iTriangles.push_back(vertices.size() - 1);
            vertices.push_back(v2);
            iTriangles.push_back(vertices.size() - 1);
            vertices.push_back(v3);
            iTriangles.push_back(vertices.size() - 1);
            return iTriangles.size() / 3;
        }

        size_t TriMesh::numberOfTriangles() const {
            return iTriangles.size() / 3;
        }

        void TriMesh::fetchTriangleVerts(TriMesh::TriangleHandle t, TriMesh::VertHandle & v1, TriMesh::VertHandle & v2, TriMesh::VertHandle & v3) const {
            v1 = iTriangles[t * 3];
            v2 = iTriangles[t * 3 + 1];
            v3 = iTriangles[t * 3 + 2];
        }

        void TriMesh::addQuad(TriMesh::VertHandle v1, TriMesh::VertHandle v2, TriMesh::VertHandle v3, TriMesh::VertHandle v4) {
            addTriangle(v1, v2, v3);
            addTriangle(v1, v3, v4);
        }

        void TriMesh::addPolygon(const std::vector<TriMesh::VertHandle> & vhs) {
            assert(vhs.size() >= 3);
            // get normal direction
            Vec3 normal = normalize(
                (ToVec3Affine(vertices[vhs[1]].position) - ToVec3Affine(vertices[vhs[0]].position)).cross(
                (ToVec3Affine(vertices[vhs[2]].position) - ToVec3Affine(vertices[vhs[1]].position)))
                );

            TriangulatePolygon(vhs.begin(), vhs.end(), [this, &normal](VertHandle vh) {
                Vec3 v = ToVec3Affine(vertices[vh].position);
                return ToVec2(v - v.dot(normal) * normal);
            }, [this](VertHandle a, VertHandle b, VertHandle c){
                addTriangle(a, b, c);
            });
        }

        void TriMesh::clear() {
            vertices.clear();
            iPoints.clear();
            iLines.clear();
            iTriangles.clear();
        }

        Box3 TriMesh::boundingBox() const {
            if (vertices.empty())
                return Box3();
            Box3 box(ToVec3Affine(vertices.front().position), ToVec3Affine(vertices.front().position));
            for (auto & v : vertices) {
                auto p = ToVec3Affine(v.position);
                box = box | BoundingBox(p);
            }
            return box;
        }


      


    }

    namespace core {

        Box3 BoundingBox(const vis::SpatialProjectedPolygon & spp){
            std::vector<Vec3> cs(spp.corners.size());
            for (int i = 0; i < spp.corners.size(); i++){
                cs[i] = IntersectionOfLineAndPlane(InfiniteLine3(spp.projectionCenter, spp.corners[i] - spp.projectionCenter), 
                    spp.plane).position;
            }
            return BoundingBoxOfContainer(cs);
        }

    }

}