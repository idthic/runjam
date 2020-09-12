#include <sys/stat.h>
#include <iostream>
#include "args.hpp"
#include "config.hpp"

namespace fsys {
  static bool exists(const char* path) {
    struct stat st;
    return stat(path, &st) == 0;
  }

  static bool is_file(const char* path) {
    struct stat st;
    return stat(path, &st) == 0 && S_ISREG(st.st_mode);
  }
}

using namespace idt::runjam;

std::string runjam_context::resodata() const {
  std::string file = this->get_config<std::string>("runjam_resodata", "ResonanceJam.dat");

  if (fsys::is_file(file.c_str())) return file;

  if (file[0] != '/') {
    std::string path = "data/" + file;
    if (fsys::is_file(path.c_str())) return path;

    if (std::strlen(PACKAGE_PREFIX)) {
      path = PACKAGE_PREFIX;
      path += "/share/runjam/";
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

namespace {

  static void print_version() {
    std::cout << "runjam (idt) " << PACKAGE_VERSION << PACKAGE_HASH << std::endl;
  }

  static void print_help() {
    std::printf(
      "usage: runjam [SUBCOMMAND] [OPTIONS|VAR=VALUE]\n"
      "\n"
      "SUBCOMMAND\n"
      "  cascade (default)\n"
      "  decay\n"
      "  sample\n"
      "\n"
      "OPTIONS and VARIABLES\n"
      "      runjam_mode=SUBCOMMAND\n"
      "  -n, runjam_nevent=INT [1]              number of events to process\n"
      "  -s, runjam_seed=INT [18371]            seed for random numbers in runjam\n"
      "      runjam_jamseed=INT [runjam_seed]   seed for random numbers in JAM\n"
      "  -t, runjam_oversampling_factor=NUM [1] number of test particles\n"
      "  -w, runjam_switch_weak_decay=BOOL [false] enable weak decays\n"
      "      runjam_phi_decays=BOOL [true]\n"
      "  --resodata, runjam_resodata=FILE [ResonanceJam.dat] resodata\n"
      "\n"
      " Output options\n"
      "  -o,        runjam_output_directory=DIR [out]   directory of output files\n"
      "  --fphase,  runjam_fname_phdat=FILE []     output filename\n"
      "  --fphase0, runjam_fname_phdat0=FILE []    output filename\n"
      "  -d INT                                         (set the next two options)\n"
      "             runjam_output_phdat=BOOL [true]     output phasespace data\n"
      "             runjam_output_phdat0=BOOL [true]    output phasespace0 data\n"
      "             runjam_output_phbin=BOOL [false]    output binary phasespace\n"
      "             runjam_output_phbin0=BOOL [false]   output binary phasespace0\n"
      "             runjam_output_phdat_indexed=BOOL [false]\n"
      "             runjam_output_phdat0_indexed=BOOL [false]\n"
      "             runjam_output_index_start=INT [0]\n"
      "\n"
      " Initialization options\n"
      "  -i ICSPEC       specify initial condition\n"
      "    c0lrf:FILE    sample particles from the hypersurface data from\n"
      "                  rfh c0lrf format \"hypersurface_v1.txt\"\n"
      "    hydrojet:DIR  sample particles using the hypersurface data from\n"
      "                  hydrojet (DIR/freezeout.dat, DIR/position.dat)\n"
      "    phase:FILE    load particle lists from the text format \"phasespace.dat\".\n"
      "    phase1:FILE   load a particle list from FILE. This performs an additional\n"
      "                  check to require that FILE contains only a single event.\n"
      "    phbin:FILE    load particle lists from the binary format \"ph000k.bin\".\n"
      "    psample:FILE  read a particle list from the file in the format of\n"
      "                  runjam \"particlesample_pos.dat\"\n"
      "\n"
      " Options for viscous sampler\n"
      "  --switching-temperature, runjam_switching_temperature=TEMP [155]\n"
      "                    an advice to switching temperature in MeV\n"
      "  runjam_turnsOffViscousEffect=INT\n"
      "\n"
      " Options for hydrojet hypersurface\n"
      "  --hydrojet-ftemp, hydrojet_kintmp=INT [5]       freezeout temperature type\n"
      "  --hydrojet-pce,   hydrojet_eospce=INT [6]       eos_pce\n"
      "  --hydrojet-bfree, hydrojet_baryonfree=INT [1]   baryonfree\n"
      "  --hydrojet-dir,   hydrojet_directory=DIR [test] directory of freezeout.dat\n"
      "  --hydrojet-dt,    hydrojet_deltat=NUM [0.3]     delta tau\n"
      "  --hydrojet-dx,    hydrojet_deltax=NUM [0.3]     delta x\n"
      "  --hydrojet-dy,    hydrojet_deltay=NUM [0.3]     delta y\n"
      "  --hydrojet-dh,    hydrojet_deltah=NUM [0.3]     delta eta\n"
      "  hydrojet_reverse_particles=BOOL [false]         (debug) perform z-reflection\n"
      "  hydrojet_shuffle_particles=BOOL [false]         (debug) shuffle particles\n"
      "\n"
      " Other options\n"
      "  --help          show this help\n"
      "  --version       show version information\n"
      "\n"
      "SAMPLE\n"
      "\n"
      "$ ./runjam cascade -s 12345 -o jam --hydrojet-dir hydro\n"
      "$ ./runjam cascade -s 12345 -o jam -i phase:phasespace0.in\n"
      "$ ./runjam decay -s 12345 -o jam -i phase:phasespace0.in\n"
      "\n"
    );
  }

  struct runjam_commandline_reader {
  public:
    runjam_commandline_reader() {}

  private:
    int argc;
    char** argv;
    runjam_commandline_arguments* args;
    runjam_context* ctx;

    int i;
    const char* arg;
    int arg_index;
    const char* arg_optarg;
    bool flag_error;
    bool flag_version;
    bool flag_help;

    const char* get_optarg() {
      if (arg_optarg) {
        const char* ret = arg_optarg;
        arg_optarg = nullptr;
        return ret;
      } else if (i < argc) {
        return argv[i++];
      } else {
        std::cerr << "runjam:$" << arg_index << " (" << arg << "): missing optional argument." << std::endl;
        flag_error = true;
        return NULL;
      }
    }

    void assign_optarg(const char* key) {
      if (const char* optarg = get_optarg())
        ctx->set_value(key, optarg);
    }

    void assign_optarg_double(const char* key) {
      if (const char* optarg = get_optarg()) {
        if (*optarg) {
          ctx->set_value(key, std::atof(optarg));
        } else {
          std::cerr << "runjam:option(" << arg << "): the argument of the option is empty." << std::endl;
          flag_error = true;
        }
      }
    }

    void assign_optarg_int(const char* key) {
      if (const char* optarg = get_optarg())
        ctx->set_value(key, std::atoi(optarg));
    }

    void assign_optarg_input() {
      const char* optarg = get_optarg();
      if (!optarg) return;

      std::string spec = optarg;

      std::size_t const pos = spec.find(':');
      if (pos == std::string::npos) {
        std::cerr << "runjam: unrecognized input '-i " << optarg << "'" << std::endl;
        std::exit(EXIT_FAILURE);
      }

      args->initType = spec.substr(0, pos);
      args->initPath = spec.substr(pos + 1);
    }

  private:
    void read_longname_option() {
      std::string longname = arg + 2;
      std::size_t ia;
      if ((ia = longname.find('=', 0)) != std::string::npos) {
        arg_optarg = &arg[2 + ia + 1]; // = の右辺
        longname = longname.substr(0, ia);
      }

      if (longname == "help") {
        flag_help = true;
      } else if (longname == "version") {
        flag_version = true;
      } else if (longname == "resodata") {
        assign_optarg("runjam_resodata");
      } else if (longname == "fphase") {
        assign_optarg("runjam_fname_phdat");
      } else if (longname == "fphase0") {
        assign_optarg("runjam_fname_phdat0");

      } else if (longname == "hydrojet-dir") {
        assign_optarg("hydrojet_directory");
      } else if (longname == "hydrojet-ftemp")  {
        assign_optarg_int("hydrojet_kintmp");
      } else if (longname == "hydrojet-pce") {
        assign_optarg_int("hydrojet_eospce");
      } else if (longname == "hydrojet-bfree") {
        assign_optarg_int("hydrojet_baryonfree");
      } else if (longname == "hydrojet-dt") {
        assign_optarg_double("hydrojet_deltat");
      } else if (longname == "hydrojet-dx") {
        assign_optarg_double("hydrojet_deltax");
      } else if (longname == "hydrojet-dy") {
        assign_optarg_double("hydrojet_deltay");
      } else if (longname == "hydrojet-dh") {
        assign_optarg_double("hydrojet_deltah");

      } else if (longname == "switching-temperature") {
        assign_optarg_double("runjam_switching_temperature");
      } else {
        std::cerr << "unknown option '" << arg << "'" << std::endl;
        flag_error = true;
      }
    }

    bool process_option(int optchar) {
      switch (optchar) {
      case 's': assign_optarg_int("runjam_seed"); break;
      case 'n': assign_optarg_int("runjam_nevent"); break;
      case 't': assign_optarg_int("runjam_oversampling_factor"); break;
      case 'o': assign_optarg("runjam_output_directory"); break;
      case 'd':
        if (const char* optarg = get_optarg()) {
          int const value = std::atoi(optarg);
          ctx->set_value("runjam_output_phdat", value);
          ctx->set_value("runjam_output_phdat0", value);
        }
        break;
      case 'w': assign_optarg_int("runjam_switch_weak_decay"); break;
      case 'i': assign_optarg_input(); break;
      default: return false;
      }
      return true;
    }
    void read_option() {
      char c;
      for (const char* p = arg + 1; p && (c = *p); p++) {
        if (p[1]) arg_optarg = p + 1;
        if (!process_option(c)) {
          std::cerr << "runjam:$" << arg_index << ": unknown option '-" << c << "'" << std::endl;
          flag_error = true;
        }
      }
    }

    bool read_assign() {
      const char* p = arg;
      while (std::isalnum(*p) && *p == '_') p++;
      if (p == arg || *p != '=') return false;

      std::string const name(arg, p);
      const char* const value = p + 1;
      ctx->set_value(name.c_str(), value);
      return true;
    }

  public:
    int read(int argc, char** argv, runjam_commandline_arguments& args, runjam_context& ctx) {
      args.subcommand = ctx.get_config<std::string>("runjam_mode", "cascade");
      args.initType = "c0lrf";
      args.initPath = "hypersurface_v1.txt";

      this->argc = argc;
      this->argv = argv;
      this->args = &args;
      this->ctx = &ctx;
      this->flag_error = false;
      this->flag_help = false;
      this->flag_version = false;

      for (i = 1; i < argc; ) {
        arg = argv[arg_index = i++];
        arg_optarg = nullptr;
        if (arg[0] == '-') {
          if (arg[1] == '-') {
            this->read_longname_option();
          } else {
            this->read_option();
          }
        } else {
          // option 以外の文字列
          if (arg_index == 1) {
            this->args->subcommand = arg;
          } else if (!read_assign()) {
            std::cerr << "runjam:$" << arg_index << ": unrecognized argument '" << arg << "'" << std::endl;
            flag_error = true;
          }
        }
      }

      bool flag_exit = false;
      int exit_status = EXIT_SUCCESS;
      if (flag_error) exit_status = 2;

      if (flag_version) {
        print_version();
        flag_exit = 1;
      }
      if (flag_help) {
        print_help();
        flag_exit = 1;
      }
      if (flag_exit)
        std::exit(exit_status);

      return exit_status;
    }
  };

}

int runjam_commandline_arguments::read(int argc, char** argv, runjam_context& ctx) {
  runjam_commandline_reader reader;
  return reader.read(argc, argv, *this, ctx);
}
