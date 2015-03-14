#ifndef PANORAMIX_VIS_SINGLETON_HPP
#define PANORAMIX_VIS_SINGLETON_HPP

class QIcon;
class QString;
class QApplication;

namespace panoramix {
    namespace vis {

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