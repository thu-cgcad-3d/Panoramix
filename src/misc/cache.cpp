#include "../gui/qttools.hpp"
#include "cache.hpp"

namespace pano {
namespace misc {

std::string Tagify(const std::string &path) {
  auto tag = path;
  std::replace(tag.begin(), tag.end(), '.', '_');
  std::replace(tag.begin(), tag.end(), '\\', '.');
  std::replace(tag.begin(), tag.end(), '/', '.');
  std::replace(tag.begin(), tag.end(), ':', '.');
  return tag;
}

std::string CachePath() { return PROJECT_CACHE_DIR_STR "/"; }

std::string FolderOfFile(const std::string &filepath) {
  QFileInfo finfo(QString::fromStdString(filepath));
  if (!finfo.exists())
    return std::string();
  return finfo.absolutePath().toStdString();
}

std::string NameOfFile(const std::string &filepath) {
  QFileInfo finfo(QString::fromStdString(filepath));
  if (!finfo.exists())
    return std::string();
  return finfo.fileName().toStdString();
}
}
}