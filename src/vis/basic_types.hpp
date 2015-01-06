#ifndef PANORAMIX_VIS_BASIC_TYPES_HPP
#define PANORAMIX_VIS_BASIC_TYPES_HPP

#include "../core/basic_types.hpp"
 
namespace panoramix {
    namespace vis {

        // color
        using Color = cv::Scalar;
        enum class ColorTag {
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

        const std::vector<ColorTag> & AllColorTags();
        std::ostream & operator << (std::ostream & os, ColorTag ct);
        Color ColorFromRGB(double r, double g, double b, double a = 255.0);
        Color ColorFromHSV(double h, double s, double v, double a = 255.0);
        Color ColorFromTag(ColorTag t);
        Color RandomColor();
        core::Vec3b ToVec3b(const Color & c);
        

        // line style
        enum class PenStyle {
            NoPen,
            SolidLine,	//1	A plain line.
            DashLine,	//2	Dashes separated by a few pixels.
            DotLine,	    //3	Dots separated by a few pixels.
            DashDotLine,	//4	Alternate dots and dashes.
            DashDotDotLine,	//5	One dash, two dots, one dash, two dots.
            CustomDashLine
        };

        // color table
        enum class ColorTableDescriptor {
            RGB,
            AllColors,
            AllColorsExcludingWhite,
            AllColorsExcludingBlack,
            AllColorsExcludingWhiteAndBlack
        };

        class ColorTable {
        public:
            inline ColorTable() 
                : _exceptionalColor(ColorFromTag(ColorTag::Transparent)) {}
            inline ColorTable(const std::vector<Color> & ctable, const Color & exceptColor = ColorFromTag(ColorTag::Transparent))
                : _colors(ctable), _exceptionalColor(exceptColor) {}

            ColorTable(ColorTableDescriptor descriptor);

            inline ColorTable(std::initializer_list<Color> c, const Color & exceptColor) 
                : _colors(c), _exceptionalColor(exceptColor) {}

            ColorTable(std::initializer_list<ColorTag> ctags, ColorTag exceptColor);

            template <class ColorIteratorT, class = std::enable_if_t<std::is_same<std::iterator_traits<ColorIteratorT>::value_type, Color>::value>>
            inline ColorTable(ColorIteratorT begin, ColorIteratorT end, const Color & exceptColor = ColorFromTag(ColorTag::Transparent))
                : _colors(begin, end), _exceptionalColor(exceptColor) {}

        public:
            const std::vector<Color> & colors() const { return _colors; }
            size_t size() const { return _colors.size(); }
            const Color & exceptionalColor() const { return _exceptionalColor; }
            Color & exceptoinalColor() { return _exceptionalColor; }
            const Color & operator[](int claz) const { return claz < 0 ? _exceptionalColor : _colors[claz]; }
            Color & operator[](int claz) { return claz < 0 ? _exceptionalColor : _colors[claz]; }
            bool empty() const { return _colors.empty(); }

            const Color & roundedAt(int claz) const { return claz < 0 ? _exceptionalColor : _colors[claz % _colors.size()]; }

            void randomize();

        private:
            std::vector<Color> _colors;
            Color _exceptionalColor;
        };

        const ColorTable & PredefinedColorTable(ColorTableDescriptor descriptor);
        ColorTable CreateGreyColorTableWithSize(int sz);
        ColorTable CreateRandomColorTableWithSize(int sz, const Color & exceptColor = ColorFromTag(ColorTag::Transparent));


        // render mode
        using RenderModeFlags = int8_t;
        enum RenderModeFlag : RenderModeFlags {
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
            OpenGLShaderSource(OpenGLShaderSourceDescriptor d = OpenGLShaderSourceDescriptor::DefaultTriangles);
            inline OpenGLShaderSource(const char * vs, const char * fs) : _vshaderSrc(vs), _fshaderSrc(fs) {}
            inline OpenGLShaderSource(const std::string & vs, const std::string & fs) :_vshaderSrc(vs), _fshaderSrc(fs) {}
            const std::string & vertexShaderSource() const { return _vshaderSrc; }
            const std::string & fragmentShaderSource() const { return _fshaderSrc; }

        private:
            std::string _vshaderSrc, _fshaderSrc;
        };

        const OpenGLShaderSource & PredefinedShaderSource(OpenGLShaderSourceDescriptor descriptor);




        // spatial projected polygon for panorama reconstruction
        struct SpatialProjectedPolygon {
            std::vector<core::Vec3> corners;
            core::Point3 projectionCenter;
            core::Plane3 plane;
        };


        // triangular mesh
        struct TriMesh {

            struct Vertex {
                Vertex();
                core::Vec4f position;
                core::Vec3f normal;
                core::Vec4f color; // the intrinsic color
                core::Vec2f texCoord;
                int entityIndex; // index in container
                uint8_t isSelected;
            };

            using VertHandle = uint32_t;
            using LineHandle = uint32_t;
            using TriangleHandle = uint32_t;

            std::vector<Vertex> vertices;
            std::vector<VertHandle> iPoints;
            std::vector<LineHandle> iLines;
            std::vector<TriangleHandle> iTriangles;


            VertHandle addVertex(const Vertex & v);

            LineHandle addLine(VertHandle v1, VertHandle v2);
            LineHandle addIsolatedLine(const Vertex & v1, const Vertex & v2);
            size_t numberOfLines() const;
            void fetchLineVerts(LineHandle l, VertHandle & v1, VertHandle & v2) const;

            TriangleHandle addTriangle(VertHandle v1, VertHandle v2, VertHandle v3);
            TriangleHandle addIsolatedTriangle(const Vertex & v1, const Vertex & v2, const Vertex & v3);
            size_t numberOfTriangles() const;
            void fetchTriangleVerts(TriangleHandle t, VertHandle & v1, VertHandle & v2, VertHandle & v3) const;

            void addQuad(VertHandle v1, VertHandle v2, VertHandle v3, VertHandle v4);
            void addPolygon(const std::vector<VertHandle> & vhs);

            void clear();

            core::Box3 boundingBox() const;     

        };


        template <class T>
        struct Colored {
            T component;
            Color color;
        };

        template <class T>
        inline Colored<T> ColorAs(const T & comp, const Color & color){
            return Colored<T>{comp, color};
        }
    }


    namespace core {

        Box3 BoundingBox(const vis::SpatialProjectedPolygon & spp);

        inline Box3 BoundingBox(const vis::TriMesh & m) {
            return m.boundingBox();
        }

    }
}
 
#endif