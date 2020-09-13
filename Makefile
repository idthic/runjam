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

INSDIR := $(DESTDIR)$(PREFIX)
libjam_LIBDIR := $(libjam_PREFIX)/lib

CXXFLAGS := $(user_CXXFLAGS) -march=native -O3 -std=gnu++11
CPPFLAGS =  $(user_CPPFLAGS) -I . -MD -MP -MF $(@:.o=.dep)
LDFLAGS  := $(user_LDFLAGS)  -L $(libjam_LIBDIR) -Wl,-rpath,$(libjam_LIBDIR)
LIBS     := $(user_LIBS)

#------------------------------------------------------------------------------

-include $(shell find $(OBJDIR) -name \*.dep 2>/dev/null)

OBJDIR := obj

runjam_OBJS := \
  $(OBJDIR)/config.o \
  $(OBJDIR)/util.o \
  $(OBJDIR)/args.o \
  $(OBJDIR)/main.o \
  $(OBJDIR)/libjam.o \
  $(OBJDIR)/ksh/integrator.o \
  $(OBJDIR)/spectra/ParticleSample.o \
  $(OBJDIR)/spectra/ResonanceList.o \
  $(OBJDIR)/spectra/ParticleSampleHydrojet.o \
  $(OBJDIR)/spectra/ParticleSampleRead.o \
  $(OBJDIR)/spectra/ParticleSamplePhasespace.o \
  $(OBJDIR)/spectra/ParticleSampleViscous.o
runjam_LIBS := -ljam $(LIBS)

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
