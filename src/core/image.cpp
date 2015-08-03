#include "macros.hpp"
#include "image.hpp"

namespace pano {
    namespace core {

        void ClipToSquare(Image & image) {
            int mindim = std::min(image.cols, image.rows);
            image = image(cv::Range(image.rows / 2 - mindim / 2, image.rows / 2 + mindim / 2),
                cv::Range(image.cols / 2 - mindim / 2, image.cols / 2 + mindim / 2));
        }

        Imageb ClipToDisk(Image & im) {
            NOT_IMPLEMENTED_YET();
        }

        Imageb Rotate(Image & im, double angle) {
            int len = std::max(im.cols, im.rows);
            cv::Mat r = cv::getRotationMatrix2D(Center<float>(im), angle / M_PI * 180, 1.0);
            Imageb mask(im.size(), true);
            cv::warpAffine(im, im, r, cv::Size(len, len));
            cv::warpAffine(mask, mask, r, cv::Size(len, len));
            return mask;
        }

        void ResizeToWidth(Image & im, int width) {
            if (im.cols == width)
                return;
            cv::resize(im, im, cv::Size(width, static_cast<int>(im.rows * width / im.cols)));
        }

        void ResizeToHeight(Image & im, int height) {
            if (im.rows == height)
                return;
            cv::resize(im, im, cv::Size(static_cast<int>(im.cols * height / im.rows), height));
        }

        void ResizeToMakeWidthUnder(Image & im, int widthUpperBound) {
            if (im.cols <= widthUpperBound)
                return;
            cv::resize(im, im, cv::Size(widthUpperBound, static_cast<int>(im.rows * widthUpperBound / im.cols)));
        }

        void ResizeToMakeHeightUnder(Image & im, int heightUpperBound) {
            if (im.rows <= heightUpperBound)
                return;
            cv::resize(im, im, cv::Size(static_cast<int>(im.cols * heightUpperBound / im.rows), heightUpperBound));
        }

        bool MayBeAPanorama(const Image & im) {
            if (abs(im.cols - im.rows * 2) > im.rows / 10.0f)
                return false;
            // check boundary pixels
            /* for (int x = 0; x < im.cols; x++){
            const uchar* p1 = im.ptr(0, x);
            const uchar* p2 = im.ptr(im.rows - 1, x);
            for (int k = 0; k < im.elemSize(); k++){
            NOT_IMPLEMENTED_YET();
            }
            }*/
            //NOT_IMPLEMENTED_YET();
            THERE_ARE_BUGS_HERE("check continuity on borders");
            return true;
        }

        bool MakePanorama(Image & im) {
            if (im.cols < im.rows * 2)
                return false;
            if (im.cols == im.rows * 2)
                return true;
            Image pim = Image::zeros(im.cols / 2, im.cols, im.type());
            im.copyTo(pim(cv::Rect(0, (pim.rows - im.rows) / 2, pim.cols, im.rows)));
            im = pim;
            return true;
        }


        std::pair<Pixel, Pixel> MinMaxLocOfImage(const Image & im) {
            Pixel minLoc, maxLoc;
            double minVal, maxVal;
            cv::minMaxLoc(im, &minVal, &maxVal, &minLoc, &maxLoc);
            return std::make_pair(minLoc, maxLoc);
        }

        std::pair<double, double> MinMaxValOfImage(const Image & im) {
            Pixel minLoc, maxLoc;
            double minVal, maxVal;
            cv::minMaxLoc(im, &minVal, &maxVal, &minLoc, &maxLoc);
            return std::make_pair(minVal, maxVal);
        }



        Pixel PixelFromGeoCoord(const GeoCoord & p, int longidiv, int latidiv) {
            int longtid = static_cast<int>((p.longitude + M_PI) * longidiv / M_PI / 2);
            int latid = static_cast<int>((p.latitude + M_PI_2) * latidiv / M_PI);
            longtid = (longtid % longidiv + longidiv) % longidiv;
            latid = (latid % latidiv + latidiv) % latidiv;
            return Pixel(longtid, latid);
        }

        GeoCoord GeoCoordFromPixel(const cv::Point & pixel, int longidiv, int latidiv) {
            return GeoCoord{ pixel.x * M_PI * 2 / longidiv - M_PI, pixel.y * M_PI / latidiv - M_PI_2 };
        }

    }
}