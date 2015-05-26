#include <QtGui>
#include <QtWidgets>

#include "singleton.hpp"
#include "qt_glue.hpp"
#include "utility.hpp"

namespace panoramix {
    namespace gui {

        int SelectFrom(const std::vector<std::string> & strs,
            const std::string & title, const std::string & text,
            int acceptId, int rejectId){
            Singleton::InitGui();
            QMessageBox mbox;
            std::vector<QAbstractButton*> buttons(strs.size(), nullptr);
            for (int i = 0; i < strs.size(); i++){
                buttons[i] = new QPushButton(QString::fromStdString(strs[i]));
                auto role = i == acceptId ? QMessageBox::ButtonRole::AcceptRole : 
                    (i == rejectId ? QMessageBox::ButtonRole::RejectRole : QMessageBox::ButtonRole::NoRole);
                mbox.addButton(buttons[i], role);
            }
            mbox.setWindowTitle(title.empty() ? QObject::tr("Make your decision") : QString::fromStdString(title));
            mbox.setText(text.empty() ? QObject::tr("Click on one of the buttons") : QString::fromStdString(text));
            mbox.exec();
            for (int i = 0; i < buttons.size(); i++){
                if (mbox.clickedButton() == buttons[i])
                    return i;
            }
            return -1;
        }

        core::Image PickAnImage(const std::string & dir, std::string * picked){
            Singleton::InitGui();
            auto filename = QFileDialog::getOpenFileName(nullptr, QObject::tr("Select an image file"), 
                QString::fromStdString(dir), 
                QObject::tr("Image Files (*.png;*.jpg;*.jpeg);;All Files (*.*)"));
            if (filename.isEmpty())
                return core::Image();
            if (picked){
                *picked = filename.toStdString();
            }
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

        Qt::PenStyle MakeQPenStyle(PenStyle ps) {
            return Qt::PenStyle(ps);
        }

        void PaintWith(const std::function<core::Image()> & updater,
            const std::vector<PenConfig> & penConfigs,
            const std::function<bool(const std::vector<core::Point2> & polyline, int penId)> & callback){

            Singleton::InitGui();

            class Widget : public QWidget {
            public:
                explicit Widget(const std::vector<PenConfig> & pc) : QWidget(), penConfigs(pc), pens(pc.size()) {
                    setMouseTracking(true);

                    // setup pen selection actions
                    setContextMenuPolicy(Qt::ActionsContextMenu);
                    QActionGroup * bas = new QActionGroup(this);

                    QAction * defaultAction = nullptr;
                    connect(defaultAction = bas->addAction(tr("No Brush")), &QAction::triggered, [this](){ activePenId = -1; });
                    {
                        QAction * sep = new QAction(bas);
                        sep->setSeparator(true);
                        bas->addAction(sep);
                    }
                    for (int i = 0; i < pc.size(); i++){
                        auto action = bas->addAction(QString::fromStdString(pc[i].name));
                        // draw icon
                        QImage image(64, 64, QImage::Format::Format_ARGB32_Premultiplied);
                        image.fill(MakeQColor(gui::White));
                        pens[i] =  QPen(MakeQColor(pc[i].color), pc[i].thickness, MakeQPenStyle(pc[i].style));
                        QPainter painter(&image);
                        painter.setPen(pens[i]);
                        painter.drawLine(QPointF(0, 32), QPointF(64, 32));
                        painter.end();
                        action->setIcon(QIcon(QPixmap::fromImage(image)));
                        connect(action, &QAction::triggered, [this, i] { activePenId = i; });
                    }

                    for (auto a : bas->actions()) a->setCheckable(true);

                    bas->setExclusive(true);
                    defaultAction->setChecked(true);

                    addActions(bas->actions());    
                    
                }

            protected:
                void paintEvent(QPaintEvent * e) override {

                }

                void mousePressEvent(QMouseEvent * e) override {
                    if (e->buttons() & Qt::LeftButton){
                        if (activePenId == -1)
                            return;
                        lastPos = e->pos();

                    }
                    else if (e->buttons() & Qt::MidButton){
                        lastPos = e->pos();

                    }
                    else{
                        QWidget::mousePressEvent(e);
                    }
                }

                void mouseMoveEvent(QMouseEvent * e) override {
                    if (e->buttons() & Qt::LeftButton){

                    }
                    else if (e->buttons() & Qt::MidButton){

                    }
                    else{
                        QWidget::mouseMoveEvent(e);
                    }
                }

                void mouseReleaseEvent(QMouseEvent * e) override {

                }

                void wheelEvent(QWheelEvent * e) override {

                }

                void keyPressEvent(QKeyEvent * e) override {
                    
                }

            public:
                QImage image;
                double scale;
                QPointF translate;
                QPolygonF stroke;
                QPointF lastPos;
                int activePenId;
                std::vector<QPen> pens;
                std::vector<PenConfig> penConfigs;
                std::function<core::Image()> updater;
                std::function<bool(const std::vector<core::Point2> & polyline, int penId)> callback;
            };

            Widget w(penConfigs);
            w.updater = updater;
            w.callback = callback;
            w.show();

            Singleton::ContinueGui();
        }

    }
}
