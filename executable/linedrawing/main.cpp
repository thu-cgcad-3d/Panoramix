#include "cache.hpp"
#include "singleton.hpp"

#include "tools.hpp"

using namespace pano;
using namespace pano::core;
using namespace pano::experimental;


void RoutineEditImageAnno();
void RoutineReconstruct0();
void RoutineReconstruct1();

int main(int argc, char **argv, char **env) {
  gui::Singleton::SetCmdArgs(argc, argv, env);
  gui::Singleton::InitGui(argc, argv);
  misc::SetCachePath(PANORAMIX_CACHE_DATA_DIR_STR "\\LineDrawing\\");

  //RoutineEditImageAnno();
  //RoutineReconstruct0();
	RoutineReconstruct1();


  return 0;
}