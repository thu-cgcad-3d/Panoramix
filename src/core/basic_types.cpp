#include "basic_types.hpp"

namespace panoramix {
    namespace core {

        namespace {
            inline Color rgb(int R, int G, int B) {
                return Color(G, B, R);
            }
            inline Color rgba(int R, int G, int B, int A) {
                return Color(G, B, R, A);
            }
        }

        Color ColorFromTag(ColorTag t) {
            switch (t){
            case ColorTag::Transparent: return rgba(0, 0, 0, 0);
            case ColorTag::White: return rgb(255, 255, 255);
            case ColorTag::Black: return rgb(0, 0, 0);
            case ColorTag::Gray: return rgb(128, 128, 128);
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

    }
}