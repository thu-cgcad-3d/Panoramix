#include "basic_types.hpp"

namespace panoramix {
    namespace vis {

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

        namespace {
            inline Color rgb(int R, int G, int B) {
                return Color(B, G, R);
            }
            inline Color rgba(int R, int G, int B, int A) {
                return Color(B, G, R, A);
            }
        }

        Color ColorFromRGB(double r, double g, double b, double a) {
            return Color(b, g, r, a);
        }

        Color ColorFromTag(ColorTag t) {
            switch (t){
            case ColorTag::Transparent: return rgba(0, 0, 0, 0);

            case ColorTag::White: return rgb(255, 255, 255);
            case ColorTag::Black: return rgb(0, 0, 0);

            case ColorTag::DimGray: return rgb(105, 105, 105);
            case ColorTag::Gray: return rgb(128, 128, 128);
            case ColorTag::DarkGray: return rgb(169, 169, 169);
            case ColorTag::Silver: return rgb(192, 192, 192);
            case ColorTag::LightGray: return rgb(211, 211, 211);

            case ColorTag::Red: return rgb(255, 0, 0);
            case ColorTag::Green: return rgb(0, 255, 0);
            case ColorTag::Blue: return rgb(0, 0, 255);

            case ColorTag::Yellow: return rgb(255, 255, 0);
            case ColorTag::Magenta: return rgb(255, 0, 255);
            case ColorTag::Cyan: return rgb(0, 255, 255);
            case ColorTag::Orange: return rgb(255, 165, 0);
            default:
                return Color(255, 255, 255);
            }
        }

        Color RandomColor() {
            return Color(rand() % 255, rand() % 255, rand() % 255);
        }

        const std::vector<Color> & PredefinedColorTable(ColorTableDescriptor descriptor) {
            static const std::vector<Color> allColorTable = {
                ColorFromTag(ColorTag::White),
                ColorFromTag(ColorTag::Black),
                               
                ColorFromTag(ColorTag::DimGray),
                ColorFromTag(ColorTag::Gray),
                ColorFromTag(ColorTag::DarkGray),
                ColorFromTag(ColorTag::Silver),
                ColorFromTag(ColorTag::LightGray),

                ColorFromTag(ColorTag::Red),
                ColorFromTag(ColorTag::Green),
                ColorFromTag(ColorTag::Blue),

                ColorFromTag(ColorTag::Yellow),
                ColorFromTag(ColorTag::Magenta),
                ColorFromTag(ColorTag::Cyan),
                ColorFromTag(ColorTag::Orange)
            };

            static const std::vector<Color> allColorExcludingWhiteTable = {
                //ColorFromTag(ColorTag::White),
                ColorFromTag(ColorTag::Black),

                ColorFromTag(ColorTag::DimGray),
                ColorFromTag(ColorTag::Gray),
                ColorFromTag(ColorTag::DarkGray),
                ColorFromTag(ColorTag::Silver),
                ColorFromTag(ColorTag::LightGray),

                ColorFromTag(ColorTag::Red),
                ColorFromTag(ColorTag::Green),
                ColorFromTag(ColorTag::Blue),

                ColorFromTag(ColorTag::Yellow),
                ColorFromTag(ColorTag::Magenta),
                ColorFromTag(ColorTag::Cyan),
                ColorFromTag(ColorTag::Orange)
            };

            static const std::vector<Color> allColorExcludingBlackTable = {
                ColorFromTag(ColorTag::White),
                //ColorFromTag(ColorTag::Black),

                ColorFromTag(ColorTag::DimGray),
                ColorFromTag(ColorTag::Gray),
                ColorFromTag(ColorTag::DarkGray),
                ColorFromTag(ColorTag::Silver),
                ColorFromTag(ColorTag::LightGray),

                ColorFromTag(ColorTag::Red),
                ColorFromTag(ColorTag::Green),
                ColorFromTag(ColorTag::Blue),

                ColorFromTag(ColorTag::Yellow),
                ColorFromTag(ColorTag::Magenta),
                ColorFromTag(ColorTag::Cyan),
                ColorFromTag(ColorTag::Orange)
            };

            static const std::vector<Color> allColorIncludingTransparentTable = {
                ColorFromTag(ColorTag::Transparent),

                ColorFromTag(ColorTag::White),
                ColorFromTag(ColorTag::Black),

                ColorFromTag(ColorTag::DimGray),
                ColorFromTag(ColorTag::Gray),
                ColorFromTag(ColorTag::DarkGray),
                ColorFromTag(ColorTag::Silver),
                ColorFromTag(ColorTag::LightGray),

                ColorFromTag(ColorTag::Red),
                ColorFromTag(ColorTag::Green),
                ColorFromTag(ColorTag::Blue),

                ColorFromTag(ColorTag::Yellow),
                ColorFromTag(ColorTag::Magenta),
                ColorFromTag(ColorTag::Cyan),
                ColorFromTag(ColorTag::Orange)
            };


            static const std::vector<Color> WRGBColorTable = {
                ColorFromTag(ColorTag::White),
                ColorFromTag(ColorTag::Red),
                ColorFromTag(ColorTag::Green),
                ColorFromTag(ColorTag::Blue),
            };
            static const std::vector<Color> RGBColorTable = {
                ColorFromTag(ColorTag::Red),
                ColorFromTag(ColorTag::Green),
                ColorFromTag(ColorTag::Blue),
            };
            
            switch (descriptor){
            case ColorTableDescriptor::WRGB: return WRGBColorTable;
            case ColorTableDescriptor::RGB: return RGBColorTable;
            case ColorTableDescriptor::AllColorsExcludingBlack: return allColorExcludingBlackTable;
            case ColorTableDescriptor::AllColorsExcludingWhite: return allColorExcludingWhiteTable;
            case ColorTableDescriptor::AllColors: return allColorTable;
            case ColorTableDescriptor::AllColorsIncludingTransparent: return allColorIncludingTransparentTable;
            default: return allColorTable;
            }
        }

    }
}