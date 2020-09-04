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

CXXFLAGS := $(user_CXXFLAGS) -march=native -O3 -std=gnu++11
CPPFLAGS =  $(user_CPPFLAGS) -I . -MD -MP -MF $(@:.o=.dep)
LDFLAGS  := $(user_LDFLAGS)  -L $(libjam_LIBDIR) -Wl,-rpath,$(libjam_LIBDIR)
LIBS     := $(user_LIBS)

#------------------------------------------------------------------------------

-include $(shell find $(OBJDIR) -name \*.dep 2>/dev/null)

OBJDIR := obj

hydro2jam_OBJS := \
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
  $(OBJDIR)/Hydro2Jam.o \
  $(OBJDIR)/util/Math.o \
  $(OBJDIR)/util/PyRand.o \
  $(OBJDIR)/util/Random.o
hydro2jam_LIBS := -ljam $(LIBS)

directories += $(OBJDIR) $(OBJDIR)/util $(OBJDIR)/spectra $(OBJDIR)/jam $(OBJDIR)/ksh
$(OBJDIR)/%.o: %.cpp | $(OBJDIR) $(OBJDIR)/util $(OBJDIR)/spectra $(OBJDIR)/jam $(OBJDIR)/ksh
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) -c -o $@ $<

all: hydro2jam.exe
hydro2jam.exe: $(hydro2jam_OBJS)
	$(CXX) $(LDFLAGS) -o $@ $^ $(hydro2jam_LIBS)


#---------------------------------------
# Install

ifneq ($(INSDIR),)

$(INSDIR)/bin/hydro2jam.exe: hydro2jam.exe | $(INSDIR)/bin
	cp $< $@
$(INSDIR)/share/hydro2jam/%: data/% | $(INSDIR)/share/hydro2jam
	cp $< $@
directories += $(INSDIR)/bin $(INSDIR)/share/hydro2jam
install-files += \
  $(INSDIR)/bin/hydro2jam.exe \
  $(INSDIR)/share/hydro2jam/ResonanceEosqJam.dat \
  $(INSDIR)/share/hydro2jam/ResonanceJam.dat
install: $(install-files)
.PHONY: install

endif

#---------------------------------------
# Clean

clean:
	-rm -rvf obj

$(directories):
	mkdir -p $@
