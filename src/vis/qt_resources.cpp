#include "qt_resources.hpp"

inline void _InitResources() {
	Q_INIT_RESOURCE(shaders);
	Q_INIT_RESOURCE(textures);
}

namespace panoramix {
    namespace vis {

        QApplication* InitGui(int argc, char ** argv) {
            QApplication* app = new QApplication(argc, argv);
        	_InitResources();
            app->setQuitOnLastWindowClosed(true);
        	return app;
        }

    }
}