#ifndef PANORAMIX_GUI_SINGLETON_HPP
#define PANORAMIX_GUI_SINGLETON_HPP

class QIcon;
class QString;
class QApplication;

namespace panoramix {
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
 
#endif