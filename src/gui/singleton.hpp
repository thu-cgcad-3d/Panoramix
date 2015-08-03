#pragma once


class QIcon;
class QString;
class QApplication;

namespace pano {
    namespace gui {

        struct Singleton {

            static const QIcon & DefaultIcon();
            static const QString & DefaultCSS();

            static QApplication* InitGui(int argc, char ** argv);
            static QApplication* InitGui();

            static int ContinueGui();
        };

    }
}
 