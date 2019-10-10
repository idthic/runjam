# -*- mode: makefile-gmake -*-

all:
.PHONY: all clean

#------------------------------------------------------------------------------
# Configuration

libjam_LIBDIR := $(HOME)/opt/jam-1.820/lib
libjam_MXV    := 200000

#CXX      := g++
CXXFLAGS := -march=native -O3
LDFLAGS  :=
LIBS     :=

#------------------------------------------------------------------------------

-include $(shell find $(OBJDIR) -name \*.dep 2>/dev/null)

OBJDIR := obj

CPPFLAGS = -DJAM_MXV=$(libjam_MXV) -I . -MD -MP -MF $(@:.o=.dep)
LDFLAGS += -L $(libjam_LIBDIR) -Wl,-rpath,$(libjam_LIBDIR)
hydro2jam_OBJS := \
  $(OBJDIR)/hydro2jam.o \
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

$(directories):
	mkdir -p $@

all: hydro2jam.exe
hydro2jam.exe: $(hydro2jam_OBJS)
	$(CXX) $(LDFLAGS) -o $@ $^ $(hydro2jam_LIBS)

clean:
	-rm -rvf obj
