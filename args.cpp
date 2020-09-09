#include "args.hpp"
#include <sys/stat.h>
#include "config.hpp"

namespace fsys {
  static bool exists(const char* path) {
    struct stat st;
    return stat(path, &st);
  }

  static bool is_file(const char* path) {
    struct stat st;
    return stat(path, &st) && S_ISREG(st.st_mode);
  }
}

std::string idt::runjam::runjam_context::resodata() const {
  std::string file = this->get_config<std::string>("runjam_resodata", "ResonanceJam.dat");

  if (fsys::is_file(file.c_str())) return file;

  if (file[0] != '/') {
    std::string path = "data/" + file;
    if (fsys::is_file(path.c_str())) return path;

    if (std::strlen(PACKAGE_PREFIX)) {
      path = PACKAGE_PREFIX;
      path += "/data/";
      path += file;
      if (fsys::is_file(path.c_str())) return path;
    }

    path = PACKAGE_BUILD;
    path += "/";
    path += file;
    if (fsys::is_file(path.c_str())) return path;
  }

  return file;
}
