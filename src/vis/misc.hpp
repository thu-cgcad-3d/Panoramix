#ifndef PANORAMIX_VIS_MISC_HPP
#define PANORAMIX_VIS_MISC_HPP

#include <QtCore>
#include <QtGui>

#include "../core/basic_types.hpp"

namespace panoramix {
    namespace vis {

        // QVectorXD
        template <class T>
        inline QVector2D MakeQVec(const core::Vec<T, 2> & v) {
            return QVector2D(static_cast<float>(v[0]), static_cast<float>(v[1]));
        }

        template <class T>
        inline QVector3D MakeQVec(const core::Vec<T, 3> & v) {
            return QVector3D(static_cast<float>(v[0]), static_cast<float>(v[1]), static_cast<float>(v[2]));
        }

        template <class T>
        inline QVector4D MakeQVec(const core::Vec<T, 4> & v) {
            return QVector4D(static_cast<float>(v[0]), static_cast<float>(v[1]), static_cast<float>(v[2]), static_cast<float>(v[3]));
        }



        QImage MakeQImage(const cv::Mat & im);
        inline QPixmap MakeQPixmap(const cv::Mat & im) {
            return QPixmap::fromImage(MakeQImage(im));
        }

        cv::Mat MakeCVMat(const QImage & im, bool inCloneImageData = true);
        inline cv::Mat MakeCVMat(const QPixmap & im, bool inCloneImageData = true) {
            return MakeCVMat(im.toImage(), inCloneImageData);
        }


    }
}
 
#endif