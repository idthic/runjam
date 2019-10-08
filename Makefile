#

all:
.PHONY: all clean

libjam_LIBDIR := $(HOME)/opt/jam-1.820/lib

CPPFLAGS = -DUSE_JAM -I . -MD -MP -MF $(@:.o=.dep)
CXXFLAGS := -march=native -O3
LDFLAGS := -L $(libjam_LIBDIR) -Wl,-rpath,$(libjam_LIBDIR)
hydro2jam_OBJS := hydro2jam.o \
  jam/Jam1.o \
  ksh/integrator.o \
  spectra/ElementReso.o \
  spectra/HydroSpectrum.o \
  spectra/IParticleSample.o \
  spectra/IResonanceList.o \
  spectra/IntegratedCooperFrye.o \
  spectra/ParticleSample.o \
  spectra/ParticleSampleDebugJamAsymmetry.o \
  spectra/ParticleSampleFromOversampledPhasespace.o \
  spectra/ParticleSampleFromPhasespaceDat.o \
  spectra/ParticleSampleViscous.o \
  user/Hydro2Jam.o \
  uty/Math.o \
  uty/PyRand.o \
  uty/Random.o \
  uty/Vector4.o
hydro2jam_LIBS := -ljam

jam/%.o: jam/%.cxx
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) -c -o $@ $<
ksh/%.o: ksh/%.cpp
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) -c -o $@ $<
spectra/%.o: spectra/%.cxx
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) -c -o $@ $<
spectra/%.o: spectra/%.cpp
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) -c -o $@ $<
user/%.o: user/%.cxx
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) -c -o $@ $<
uty/%.o: uty/%.cxx
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) -c -o $@ $<

all: hydro2jam.exe
hydro2jam.exe: $(hydro2jam_OBJS)
	$(CXX) $(LDFLAGS) -o $@ $^ $(hydro2jam_LIBS)

clean:
	-rm -f *.o */*.o
