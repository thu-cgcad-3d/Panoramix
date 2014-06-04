#ifndef PANORAMIX_VIS_SINGLETON_HPP
#define PANORAMIX_VIS_SINGLETON_HPP

class QApplication;

namespace panoramix {
    namespace vis {

        struct Singleton {
            static QApplication* InitGui(int argc, char ** argv);
            static QApplication* InitGui();

            static void ContinueGui();
        };

    }
}
 
#endif