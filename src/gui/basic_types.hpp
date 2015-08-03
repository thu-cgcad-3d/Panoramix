#pragma once


#include "../core/basic_types.hpp"
 
namespace pano {
    namespace gui {

        // color
        enum ColorTag {
            Transparent,

            White,
            Black,

            DimGray,
            Gray,
            DarkGray,
            Silver,
            LightGray,

            Red,
            Green,
            Blue,

            Yellow,
            Magenta,
            Cyan,
            Orange
        };

        class Color {
        public:
            inline Color() : _rgba(255, 255, 255, 255){}
            inline Color(int r, int g, int b, int a = 255) 
                : _rgba(r, g, b, a) {}
            inline Color(double r, double g, double b, double a = 1.0)
                : _rgba(static_cast<int>(r * 255), static_cast<int>(g * 255),
                static_cast<int>(b * 255), static_cast<int>(a * 255)) {}

            // from vec4
            template <class T, class = std::enable_if_t<std::is_integral<T>::value>>
            inline Color(const core::Vec<T, 4> & v) : _rgba(v) {}

            template <class T, class = std::enable_if_t<std::is_floating_point<T>::value>, class = void>
            inline Color(const core::Vec<T, 4> & v) : _rgba(v * 255) {}

            // from vec3
            template <class T, class = std::enable_if_t<std::is_integral<T>::value>>
            inline Color(const core::Vec<T, 3> & v, T a = 255) 
                : _rgba(static_cast<int>(v[0]), static_cast<int>(v[1]), 
                static_cast<int>(v[2]), a) {}

            template <class T, class = std::enable_if_t<std::is_floating_point<T>::value>, class = void>
            inline Color(const core::Vec<T, 3> & v, T a = 1.0) 
                : _rgba(static_cast<int>(v[0] * 255), static_cast<int>(v[1] * 255),
                static_cast<int>(v[2] * 255), a * 255) {}

            // from tag
            Color(ColorTag tag);
            // from raw data
            Color(const std::uint8_t * data, int cvType);
            
        public:
            inline int red() const { return _rgba[0]; }
            inline int green() const { return _rgba[1]; }
            inline int blue() const { return _rgba[2]; }
            inline int alpha() const { return _rgba[3]; }

            inline float redf() const { return _rgba[0] / 255.0f; }
            inline float greenf() const { return _rgba[1] / 255.0f; }
            inline float bluef() const { return _rgba[2] / 255.0f; }
            inline float alphaf() const { return _rgba[3] / 255.0f; }

            inline bool isTransparent() const { return _rgba[3] == 0; }

            // to cv::Scalar (bgra)
            inline operator cv::Scalar() const { return cv::Scalar(_rgba[2], _rgba[1], _rgba[0], _rgba[3]); }
            
            // to vec4            
            template <class T, class = std::enable_if_t<std::is_integral<T>::value>>
            inline operator core::Vec<T, 4>() const { return _rgba; }

            template <class T, class = std::enable_if_t<std::is_floating_point<T>::value>, class = void>
            inline operator core::Vec<T, 4>() const { return _rgba / 255.0; }

            // to vec3
            template <class T, class = std::enable_if_t<std::is_integral<T>::value>>
            inline operator core::Vec<T, 3>() const { return core::Vec<T, 3>(_rgba[0], _rgba[1], _rgba[2]); }

            template <class T, class = std::enable_if_t<std::is_floating_point<T>::value>, class = void>
            inline operator core::Vec<T, 3>() const { 
                return core::Vec<T, 3>(_rgba[0] / 255.0, _rgba[1] / 255.0, _rgba[2] / 255.0); 
            }

            inline bool operator == (const Color & color) const { return _rgba == color._rgba; }
            
            template <class T, class = std::enable_if_t<std::is_arithmetic<T>::value>>
            inline Color & operator *= (T s) { _rgba *= s; return *this; }
            template <class T, class = std::enable_if_t<std::is_arithmetic<T>::value>>
            inline Color & operator /= (T s) { _rgba /= s; return *this; }

            inline Color blendWith(const Color & c, double alpha) const { return _rgba * (1.0 - alpha) + c._rgba * alpha; }

            template <class Archive> inline void serialize(Archive & ar) { ar(_rgba); }

        private:
            core::Vec4i _rgba;
        };

        template <class T, class = std::enable_if_t<std::is_arithmetic<T>::value>>
        inline Color operator * (const Color & c, T d) { Color cc = c; cc *= d; return cc; }
        template <class T, class = std::enable_if_t<std::is_arithmetic<T>::value>>
        inline Color operator * (T d, const Color & c) { Color cc = c; cc *= d; return cc; }
        template <class T, class = std::enable_if_t<std::is_arithmetic<T>::value>>
        inline Color operator / (const Color & c, T d) { Color cc = c; cc /= d; return cc; }


        std::ostream & operator << (std::ostream & os, ColorTag ct);
        Color ColorFromHSV(double h, double s, double v, double a = 1.0);
        Color RandomColor();
        inline Color ColorFromImage(const core::Image & im, core::Pixel p) { return Color(im.ptr(p.y, p.x), im.type()); }



        // line style
        // keep sync with Qt::PenStyle
        enum PenStyle {
            NoPen,
            SolidLine,	//1	A plain line.
            DashLine,	//2	Dashes separated by a few pixels.
            DotLine,	    //3	Dots separated by a few pixels.
            DashDotLine,	//4	Alternate dots and dashes.
            DashDotDotLine,	//5	One dash, two dots, one dash, two dots.
            CustomDashLine
        };

        struct PenConfig {
            std::string name;
            std::string description;
            double thickness;
            Color color;
            PenStyle style;
        };


        // color table
        enum ColorTableDescriptor {
            RGB,
            AllColors,
            AllColorsExcludingWhite,
            AllColorsExcludingBlack,
            AllColorsExcludingWhiteAndBlack,
            RGBGreys
        };

        class ColorTable {
        public:
            inline ColorTable() 
                : _exceptionalColor(ColorTag::Transparent) {}
            inline ColorTable(const std::vector<Color> & ctable, const Color & exceptColor = ColorTag::Transparent)
                : _colors(ctable), _exceptionalColor(exceptColor) {}
            inline ColorTable(std::vector<Color> && ctable, const Color & exceptColor = ColorTag::Transparent)
                : _colors(std::move(ctable)), _exceptionalColor(exceptColor) {}

            ColorTable(ColorTableDescriptor descriptor);

            inline ColorTable(std::initializer_list<Color> c, const Color & exceptColor = ColorTag::Transparent) 
                : _colors(c), _exceptionalColor(exceptColor) {}

            template <class ColorIteratorT, class = std::enable_if_t<std::is_same<std::iterator_traits<ColorIteratorT>::value_type, Color>::value>>
            inline ColorTable(ColorIteratorT begin, ColorIteratorT end, const Color & exceptColor = ColorTag::Transparent)
                : _colors(begin, end), _exceptionalColor(exceptColor) {}


            inline ColorTable(ColorTable && ctable) 
                : _colors(std::move(ctable._colors)), _exceptionalColor(std::move(ctable._exceptionalColor)) {}
            ColorTable & operator = (ColorTable && ctable) { 
                _colors = std::move(ctable._colors); 
                _exceptionalColor = std::move(ctable._exceptionalColor); 
                return *this;
            }


        public:
            const std::vector<Color> & colors() const { return _colors; }
            size_t size() const { return _colors.size(); }
            const Color & exceptionalColor() const { return _exceptionalColor; }
            Color & exceptoinalColor() { return _exceptionalColor; }
            const Color & operator[](int claz) const { return claz < 0 ? _exceptionalColor : _colors[claz]; }
            Color & operator[](int claz) { return claz < 0 ? _exceptionalColor : _colors[claz]; }
            bool empty() const { return _colors.empty(); }

            core::Image3ub operator()(const core::Imagei & indexIm) const;

            const Color & roundedAt(int claz) const { return claz < 0 ? _exceptionalColor : _colors[claz % _colors.size()]; }

            ColorTable & randomize();
            ColorTable & appendRandomizedColors(size_t size);
            ColorTable & appendRandomizedGreyColors(size_t size);

            template <class Archive> inline void serialize(Archive & ar) { ar(_colors, _exceptionalColor); }

        private:
            std::vector<Color> _colors;
            Color _exceptionalColor;
        };

        ColorTable CreateGreyColorTableWithSize(int sz);
        ColorTable CreateRandomColorTableWithSize(int sz, const Color & exceptColor = ColorTag::Transparent);


        // render mode
        using RenderModeFlags = int8_t;
        enum RenderModeFlag : RenderModeFlags {
            None        = 0,
            Points      = 1,
            Lines       = 1 << 1,
            Triangles   = 1 << 2,
            All         = (1 << 3) - 1
        };
        inline RenderModeFlags operator | (RenderModeFlag f1, RenderModeFlag f2) {
            return static_cast<RenderModeFlags>(f1) | static_cast<RenderModeFlags>(f2);
        }



        // opengl shader source
        enum class OpenGLShaderSourceDescriptor {
            DefaultPoints,
            DefaultLines,
            DefaultTriangles,
            Panorama,

            XPoints,
            XLines,
            XTriangles,
            XPanorama
        };

        class OpenGLShaderSource {
        public:
            OpenGLShaderSource(OpenGLShaderSourceDescriptor d = OpenGLShaderSourceDescriptor::XTriangles);

            template <class StringT1, class StringT2>
            OpenGLShaderSource(StringT1 && vs, StringT2 && fs) : _vshaderSrc(std::forward<StringT1>(vs)), _fshaderSrc(std::forward<StringT2>(fs)) {}

            const std::string & vertexShaderSource() const { return _vshaderSrc; }
            const std::string & fragmentShaderSource() const { return _fshaderSrc; }

            template <class Archive> inline void serialize(Archive & ar) { ar(_vshaderSrc, _fshaderSrc); }

        private:
            std::string _vshaderSrc, _fshaderSrc;
        };




        // spatial projected polygon for panorama reconstruction
        struct SpatialProjectedPolygon {
            std::vector<core::Vec3> corners;
            core::Point3 projectionCenter;
            core::Plane3 plane;
        };


        template <class T>
        struct Colored {
            T component;
            Color color;

            template <class Archive> 
            inline void serialize(Archive & ar) {
                ar(component, color);
            }
        };

        template <class T>
        inline Colored<std::decay_t<T>> ColorAs(T && comp, const Color & color){
            return Colored<std::decay_t<T>>{std::forward<T>(comp), color};
        }
    }


    namespace core {

        Box3 BoundingBox(const gui::SpatialProjectedPolygon & spp);

    }

}

#define DECL_PROPERTY(claz, type, name) \
    private: type _##name;  \
    public: inline const type & name() const { return _##name; }\
    public: inline type & name() { return _##name; }\
    public: inline claz & name(const type & v) { _##name = v; return *this; }

 
