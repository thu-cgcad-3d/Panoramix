#ifndef PANORAMIX_VIS_MISC_HPP
#define PANORAMIX_VIS_MISC_HPP

#include <QtCore>
#include <QtGui>

#include "../core/basic_types.hpp"

namespace panoramix {
    namespace vis {

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



        // color
        inline QRgb MakeQRgb(const core::Color & c) { 
            return qRgba(static_cast<int>(c[2]), static_cast<int>(c[1]), 
                static_cast<int>(c[0]), static_cast<int>(c[3])); 
        }

        inline QColor MakeQColor(const core::Color & c) {
            return QColor(MakeQRgb(c));
        }





        // vector
        template <class T>
        inline QVector2D MakeQVec(const core::Vec<T, 2> & v) {
            return QVector2D(static_cast<float>(v[0]), static_cast<float>(v[1]));
        }

        template <class T>
        inline QVector3D MakeQVec(const core::Vec<T, 3> & v) {
            return QVector3D(static_cast<float>(v[0]), static_cast<float>(v[1]), 
                static_cast<float>(v[2]));
        }
        
        inline core::Vec3 MakeCoreVec(const QVector3D & v) {
            return core::Vec3(v.x(), v.y(), v.z());
        }

        template <class T>
        inline QVector4D MakeQVec(const core::Vec<T, 4> & v) {
            return QVector4D(static_cast<float>(v[0]), static_cast<float>(v[1]), 
                static_cast<float>(v[2]), static_cast<float>(v[3]));
        }

        inline QVector4D MakeQVec(const core::Color & v) {
            return QVector4D(static_cast<float>(v[0]), static_cast<float>(v[1]),
                static_cast<float>(v[2]), static_cast<float>(v[3] == 0.0 ? 255.0f : v[3]));
        }


        // matrix
        template <class T, int M, int N>
        QGenericMatrix<N, M, T> MakeQMatrix(const core::Mat<T, M, N> & m) {
            QGenericMatrix<N, M, T> mat;
            for (int i = 0; i < M; i++){
                for (int j = 0; j < N; j++){
                    mat(i, j) = m(i, j);
                }
            }
            return mat;
        }
        QMatrix4x4 MakeQMatrix(const core::Mat<float, 4, 4> & m);
        QMatrix4x4 MakeQMatrix(const core::Mat<double, 4, 4> & m);


        // point
        template <class T>
        inline QPointF MakeQPointF(const core::Point<T, 2> & p) {
            return QPointF(static_cast<float>(p[0]), static_cast<float>(p[1]));
        }
        inline QPoint MakeQPoint(const core::PixelLoc & p) {
            return QPoint(p.x, p.y);
        }



        // size
        inline QSizeF MakeQSizeF(const core::Size & sz) {
            return QSizeF(sz.width, sz.height);
        }
        inline QSize MakeQSize(const core::SizeI & sz) {
            return QSize(sz.width, sz.height);
        }



        // image
        QImage MakeQImage(const core::Image & im);
        inline QPixmap MakeQPixmap(const core::Image & im) {
            return QPixmap::fromImage(MakeQImage(im));
        }

        core::Image MakeCVMat(const QImage & im, bool inCloneImageData = true);
        inline core::Image MakeCVMat(const QPixmap & im, bool inCloneImageData = true) {
            return MakeCVMat(im.toImage(), inCloneImageData);
        }


    }
}
 
#endif