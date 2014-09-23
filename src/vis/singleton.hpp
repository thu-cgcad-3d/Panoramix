#ifndef PANORAMIX_VIS_SINGLETON_HPP
#define PANORAMIX_VIS_SINGLETON_HPP

class QGuiApplication;

namespace panoramix {
    namespace vis {

        struct Singleton {
            static QGuiApplication* InitGui(int argc, char ** argv);
            static QGuiApplication* InitGui();

            static void ContinueGui();
        };

    }
}
 
#endif