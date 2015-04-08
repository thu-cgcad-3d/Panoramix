#include "algorithms.hpp"
#include "basic_types.hpp"

namespace panoramix {
    namespace core {

        void ResizeToWidth(Image & im, int width){
            if (im.cols == width)
                return;
            cv::resize(im, im, cv::Size(width, static_cast<int>(im.rows * width / im.cols)));
        }

        void ResizeToHeight(Image & im, int height){
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

        bool MayBeAPanorama(const Image & im){
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

        bool MakePanorama(Image & im){
            if (im.cols < im.rows * 2)
                return false;
            if (im.cols == im.rows * 2)
                return true;
            Image pim = Image::zeros(im.cols / 2, im.cols, im.type());
            im.copyTo(pim(cv::Rect(0, (pim.rows - im.rows) / 2, pim.cols, im.rows)));
            im = pim;
            return true;
        }


        std::pair<PixelLoc, PixelLoc> MinMaxLocOfImage(const Image & im) {
            PixelLoc minLoc, maxLoc;
            double minVal, maxVal;
            cv::minMaxLoc(im, &minVal, &maxVal, &minLoc, &maxLoc);
            return std::make_pair(minLoc, maxLoc);
        }
        
        std::pair<double, double> MinMaxValOfImage(const Image & im) {
            PixelLoc minLoc, maxLoc;
            double minVal, maxVal;
            cv::minMaxLoc(im, &minVal, &maxVal, &minLoc, &maxLoc);
            return std::make_pair(minVal, maxVal);
        }



        PixelLoc PixelLocFromGeoCoord(const GeoCoord & p, int longidiv, int latidiv) {
            int longtid = static_cast<int>((p.longitude + M_PI) * longidiv / M_PI / 2);
            int latid = static_cast<int>((p.latitude + M_PI_2) * latidiv / M_PI);
            longtid = (longtid % longidiv + longidiv) % longidiv;
            latid = (latid % latidiv + latidiv) % latidiv;
            return PixelLoc(longtid, latid);
        }

        GeoCoord GeoCoordFromPixelLoc(const cv::Point & pixel, int longidiv, int latidiv) {
            return GeoCoord{ pixel.x * M_PI * 2 / longidiv - M_PI, pixel.y * M_PI / latidiv - M_PI_2 };
        }

        template <class T>
        T AreaTemplated(const Polygon<T, 3> & polygon){
            T area = 0.0;
            Vec<T, 3> x, y;
            std::tie(x, y) = ProposeXYDirectionsFromZDirection(polygon.normal);
            TriangulatePolygon(polygon.corners.begin(), polygon.corners.end(), [&](const Point<T, 3> & v) {
                auto proj = v - v.dot(polygon.normal) * polygon.normal;
                return Vec<T, 2>(proj.dot(x), proj.dot(y));
            }, [&area](const Point<T, 3> & a, const Point<T, 3> & b, const Point<T, 3> & c){
                double aa = AreaOfTriangle<T, 3>(a, b, c);
                area += IsInfOrNaN(aa) ? 0.0 : aa;
            });
            return area;
        }

        double Area(const Polygon3 & polygon) { return AreaTemplated(polygon); }
        float Area(const Polygon3f & polygon) { return AreaTemplated(polygon); }

    }
}