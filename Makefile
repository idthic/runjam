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

-include $(wildcard */*.dep) hydro2jam.dep

CPPFLAGS = -DUSE_JAM -DJAM_MXV=$(libjam_MXV) -I . -MD -MP -MF $(@:.o=.dep)
LDFLAGS += -L $(libjam_LIBDIR) -Wl,-rpath,$(libjam_LIBDIR)
hydro2jam_OBJS := hydro2jam.o \
  jam/Jam1.o \
  ksh/integrator.o \
  spectra/ElementReso.o \
  spectra/HydroSpectrum.o \
  spectra/IParticleSample.o \
  spectra/IResonanceList.o \
  spectra/IntegratedCooperFrye.o \
  spectra/ParticleSampleHydrojet.o \
  spectra/ParticleSampleRead.o \
  spectra/ParticleSamplePhasespace.o \
  spectra/ParticleSampleViscous.o \
  Hydro2Jam.o \
  util/Math.o \
  util/PyRand.o \
  util/Random.o
hydro2jam_LIBS := -ljam $(LIBS)

jam/%.o: jam/%.cpp
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) -c -o $@ $<
ksh/%.o: ksh/%.cpp
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) -c -o $@ $<
spectra/%.o: spectra/%.cpp
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) -c -o $@ $<
user/%.o: user/%.cpp
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) -c -o $@ $<
uty/%.o: uty/%.cpp
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) -c -o $@ $<

all: hydro2jam.exe
hydro2jam.exe: $(hydro2jam_OBJS)
	$(CXX) $(LDFLAGS) -o $@ $^ $(hydro2jam_LIBS)

clean:
	-rm -f *.o */*.o *.dep */*.dep
