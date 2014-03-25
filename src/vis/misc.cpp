#include "misc.hpp"

#include <QtDebug>

namespace panoramix {
    namespace vis {

        namespace {
            inline QVector<QRgb> MakeColorTable() {
                QVector<QRgb>  sColorTable;
                for (int i = 0; i < 256; ++i)
                    sColorTable.push_back(qRgb(i, i, i));
                return sColorTable;
            }
        }

        QImage MakeQImage(const cv::Mat & inMat) {
            switch (inMat.type()) {
                // 8-bit, 4 channel
            case CV_8UC4:
                return QImage(inMat.data, inMat.cols, inMat.rows, inMat.step, QImage::Format_RGB32);
                // 8-bit, 3 channel
            case CV_8UC3:
                return QImage(inMat.data, inMat.cols, inMat.rows, inMat.step, QImage::Format_RGB888).rgbSwapped();
                // 8-bit, 1 channel
            case CV_8UC1:{
                             static const QVector<QRgb> sColorTable = MakeColorTable();
                             QImage image(inMat.data, inMat.cols, inMat.rows, inMat.step, QImage::Format_Indexed8);
                             image.setColorTable(sColorTable);
                             return image;
            }
            default:
                qWarning() << "MakeQImage() - cv::Mat image type not handled in switch:" << inMat.type();
                break;
            }
            return QImage();
        }

        cv::Mat MakeCVMat(const QImage & inImage, bool inCloneImageData) {
            switch (inImage.format()){
                // 8-bit, 4 channel
            case QImage::Format_RGB32:{
                                          cv::Mat mat(inImage.height(), inImage.width(), CV_8UC4, 
                                              const_cast<uchar*>(inImage.bits()), inImage.bytesPerLine());
                                          return (inCloneImageData ? mat.clone() : mat);
            }
                // 8-bit, 3 channel
            case QImage::Format_RGB888:{
                                           if (!inCloneImageData)
                                               qWarning() << "MakeCVMat() - Conversion requires cloning since we use a temporary QImage";
                                           QImage swapped = inImage.rgbSwapped();
                                           return cv::Mat(swapped.height(), swapped.width(), CV_8UC3, 
                                               const_cast<uchar*>(swapped.bits()), swapped.bytesPerLine()).clone();
            }
                // 8-bit, 1 channel
            case QImage::Format_Indexed8:{
                                             cv::Mat mat(inImage.height(), inImage.width(), CV_8UC1, 
                                                 const_cast<uchar*>(inImage.bits()), inImage.bytesPerLine());
                                             return (inCloneImageData ? mat.clone() : mat);
            }
            default:
                qWarning() << "AMakeCVMat() - QImage format not handled in switch:" << inImage.format();
                break;
            }
            return cv::Mat();
        }

    }
}