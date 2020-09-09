# -*- mode: makefile-gmake -*-

all:
.PHONY: all clean

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

INSDIR := $(DESTDIR)$(PREFIX)
libjam_LIBDIR := $(libjam_PREFIX)/lib

PACKAGE_HASH := $(shell git show -s --format=%h 2>/dev/null)
ifneq ($(PACKAGE_HASH),)
  PACKAGE_HASH := +$(PACKAGE_HASH)
endif

CXXFLAGS := $(user_CXXFLAGS) -march=native -O3 -std=gnu++11
CPPFLAGS =  $(user_CPPFLAGS) -I . -MD -MP -MF $(@:.o=.dep) -D'PACKAGE_HASH="$(PACKAGE_HASH)"'
LDFLAGS  := $(user_LDFLAGS)  -L $(libjam_LIBDIR) -Wl,-rpath,$(libjam_LIBDIR)
LIBS     := $(user_LIBS)

#------------------------------------------------------------------------------

-include $(shell find $(OBJDIR) -name \*.dep 2>/dev/null)

OBJDIR := obj

runjam_OBJS := \
  $(OBJDIR)/main.o \
  $(OBJDIR)/args.o \
  $(OBJDIR)/jam/Jam1.o \
  $(OBJDIR)/ksh/integrator.o \
  $(OBJDIR)/spectra/ElementReso.o \
  $(OBJDIR)/spectra/HydroSpectrum.o \
  $(OBJDIR)/spectra/IParticleSample.o \
  $(OBJDIR)/spectra/IResonanceList.o \
  $(OBJDIR)/spectra/IntegratedCooperFrye.o \
  $(OBJDIR)/spectra/ParticleSampleHydrojet.o \
  $(OBJDIR)/spectra/ParticleSampleRead.o \
  $(OBJDIR)/spectra/ParticleSamplePhasespace.o \
  $(OBJDIR)/spectra/ParticleSampleViscous.o \
  $(OBJDIR)/RunJam.o \
  $(OBJDIR)/util/Math.o \
  $(OBJDIR)/util/PyRand.o \
  $(OBJDIR)/util/Random.o
runjam_LIBS := -ljam $(LIBS)

directories += $(OBJDIR) $(OBJDIR)/util $(OBJDIR)/spectra $(OBJDIR)/jam $(OBJDIR)/ksh
$(OBJDIR)/%.o: %.cpp | $(OBJDIR) $(OBJDIR)/util $(OBJDIR)/spectra $(OBJDIR)/jam $(OBJDIR)/ksh
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) -c -o $@ $<

all: runjam.exe
runjam.exe: $(runjam_OBJS)
	$(CXX) $(LDFLAGS) -o $@ $^ $(runjam_LIBS)


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
  $(INSDIR)/share/runjam/ResonanceEosqJam.dat \
  $(INSDIR)/share/runjam/ResonanceJam.dat
install: $(install-files)
.PHONY: install

endif

#---------------------------------------
# Clean

clean:
	-rm -rvf obj

$(directories):
	mkdir -p $@
