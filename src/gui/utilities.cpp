#include <QtGui>
#include <QtWidgets>

#include "singleton.hpp"
#include "qt_glue.hpp"
#include "utitilies.hpp"

namespace panoramix {
    namespace gui {

        core::Image PickAnImage(const std::string & dir){
            Singleton::InitGui();
            auto filename = QFileDialog::getOpenFileName(nullptr, QObject::tr("Select an image file"), 
                QString::fromStdString(dir), 
                QObject::tr("Image Files (*.png;*.jpg;*.jpeg);;All Files (*.*)"));
            if (filename.isEmpty())
                return core::Image();
            return cv::imread(filename.toStdString());
        }

        std::vector<core::Image> PickImages(const std::string & dir){
            Singleton::InitGui();
            auto filenames = QFileDialog::getOpenFileNames(nullptr, QObject::tr("Select an image file"),
                QString::fromStdString(dir),
                QObject::tr("Image Files (*.png;*.jpg;*.jpeg);;All Files (*.*)"));
            std::vector<core::Image> ims;
            for (auto & filename : filenames){
                ims.push_back(cv::imread(filename.toStdString()));
            }
            return ims;
        }

        std::vector<core::Image> PickAllImagesFromAFolder(const std::string & dir){
            Singleton::InitGui();
            auto folder = QFileDialog::getExistingDirectory(nullptr, QObject::tr("Select a folder containing images"),
                QString::fromStdString(dir));
            NOT_IMPLEMENTED_YET();
        }

    }
}
