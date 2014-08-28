#ifndef PANORAMIX_VIS_BASIC_TYPES_HPP
#define PANORAMIX_VIS_BASIC_TYPES_HPP

#include <cstdint>
#include <opencv2/opencv.hpp>
 
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
        Color ColorFromTag(ColorTag t);
        

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


        enum class ColorTableDescriptor : int8_t {
            RGB,
            WRGB,
            AllColors
        };
        const std::vector<Color> &
            PredefinedColorTable(ColorTableDescriptor descriptor = ColorTableDescriptor::AllColors);



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

    }
}
 
#endif