#ifndef PANORAMIX_VIS_SINGLETON_HPP
#define PANORAMIX_VIS_SINGLETON_HPP

class QIcon;
class QApplication;

namespace panoramix {
    namespace vis {

        struct Singleton {

            struct Configuration {
                QIcon icon;
                QString css;
            };

            static const Configuration & DefaultConfiguration();

            static QApplication* InitGui(int argc, char ** argv);
            static QApplication* InitGui();

            static void ContinueGui();
        };

    }
}
 
#endif