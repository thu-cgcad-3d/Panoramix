#include "panorama_reconstruction.hpp"

int main_label(int argc, char **argv,
               const std::map<std::string, std::string> &options);
int main_run(int argc, char **argv,
             const std::map<std::string, std::string> &options);

int main(int argc, char **argv) {
  std::map<std::string, std::string> options;
  for (int i = 1; i < argc; i += 2) {
    std::string name = argv[i];
    if (name.empty() || name[0] != '-') {
      std::cerr << "name " << name << " should starts with -" << std::endl;
      break;
    }
    std::string value = argv[i + 1];
    options[name.substr(1)] = value;
  }

  auto routine = &main_run;
  if (options.count("mode") && options.at("mode") == "label") {
    routine = &main_label;
  }
  return routine(argc, argv, options);
}