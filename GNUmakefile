# -*- mode: makefile-gmake -*-

all:
.PHONY: all clean

PRECMD := $(shell ./mktool.sh update-commit-hash)

#------------------------------------------------------------------------------
# Compile config

configure := ./configure
config.mk: configure
	$(configure)

-include config.mk
user_CPPFLAGS := $(CPPFLAGS)
user_CXXFLAGS := $(CXXFLAGS)
user_LDFLAGS  := $(LDFLAGS)
user_LIBS     := $(LIBS)

config_CPPFLAGS :=
config_LDFLAGS  :=
ifneq ($(use_libjam1),)
  libjam1_LIBDIR  := $(libjam1_PREFIX)/lib
  config_LDFLAGS  += -L$(libjam1_LIBDIR) -Wl,-rpath,$(libjam1_LIBDIR)
endif
ifneq ($(use_libjam2),)
  libjam2_LIBDIR  := $(libjam2_PREFIX)/lib
  pythia8_LIBDIR  := $(pythia8_PREFIX)/lib
  config_CPPFLAGS +=  -isystem $(libjam2_PREFIX)/include -isystem $(pythia8_PREFIX)/include
  config_LDFLAGS  += -L$(libjam2_LIBDIR) -L$(pythia8_LIBDIR) -Wl,-rpath,$(libjam2_LIBDIR) -Wl,-rpath,$(pythia8_LIBDIR)
endif

INSDIR := $(DESTDIR)$(PREFIX)

CXXFLAGS := $(user_CXXFLAGS) -march=native -O3 -std=gnu++17
CPPFLAGS  =  $(user_CPPFLAGS) -I . $(config_CPPFLAGS) -MD -MP -MF $(@:.o=.dep)
LDFLAGS  := $(user_LDFLAGS) $(config_LDFLAGS)
LIBS     := $(user_LIBS)

#------------------------------------------------------------------------------

-include $(shell find $(OBJDIR) -name \*.dep 2>/dev/null)

OBJDIR := obj

runjam_OBJS := \
  $(OBJDIR)/config.o \
  $(OBJDIR)/util.o \
  $(OBJDIR)/args.o \
  $(OBJDIR)/jamimpl.o \
  $(OBJDIR)/ksh/integrator.o \
  $(OBJDIR)/ParticleSample.o \
  $(OBJDIR)/ResonanceList.o \
  $(OBJDIR)/cmd_resolist_feeddown_factor.o \
  $(OBJDIR)/main.o \
  $(patsubst %.cpp,$(OBJDIR)/%.o,$(wildcard ParticleSample*.cpp))

runjam_LIBS := $(LIBS)
ifneq ($(use_libjam1),)
  runjam_OBJS += $(OBJDIR)/libjam1.o
  runjam_LIBS := -ljam $(runjam_LIBS)
endif
ifneq ($(use_libjam2),)
  runjam_OBJS += $(OBJDIR)/libjam2.o
  runjam_LIBS := -ljam2 -lpythia8 $(runjam_LIBS)
endif

directories += $(OBJDIR) $(OBJDIR)/spectra $(OBJDIR)/ksh
$(OBJDIR)/%.o: %.cpp | $(OBJDIR) $(OBJDIR)/spectra $(OBJDIR)/ksh
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) -c -o $@ $<

all: runjam.exe
runjam.exe: $(runjam_OBJS)
	$(CXX) $(CXXFLAGS) $(LDFLAGS) -o $@ $^ $(runjam_LIBS)


#---------------------------------------
# Install

ifneq ($(INSDIR),)

$(INSDIR)/bin/runjam.exe: runjam.exe | $(INSDIR)/bin
	cp $< $@
$(INSDIR)/share/runjam/%: data/% | $(INSDIR)/share/runjam
	cp $< $@
directories += $(INSDIR)/bin $(INSDIR)/share/runjam
install-files += \
  $(INSDIR)/bin/runjam.exe \
  $(patsubst data/%,$(INSDIR)/share/runjam/%,$(wildcard data/Resonance*.dat))
install: $(install-files)
.PHONY: install

endif

#---------------------------------------
# Clean

clean:
	-rm -rvf obj

$(directories):
	mkdir -p $@
