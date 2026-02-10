/* This file is a part of runjam <https://github.com/idthic/runjam>.

   Copyright (C) 2022, Koichi Murase <myoga.murase at gmail.com>

   SPDX-License-Identifier: GPL-2.0-or-later

   This program is free software; you can redistribute it and/or
   modify it under the terms of the GNU General Public License as
   published by the Free Software Foundation; either version 2 of the
   License, or any later version.

   This program is distributed in the hope that it will be useful, but
   WITHOUT ANY WARRANTY; without even the implied warranty of
   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
   General Public License for more details.

   You should have received a copy of the GNU General Public License
   along with this program; if not, write to the Free Software
   Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA 02111-1307
   USA  */

#include <memory>
#include <vector>
#include <string>
#include <Pythia8/Settings.h>
#include "ParticleSample.hpp"

namespace jam2 { struct JAM; }

namespace libjam2 {

  class irunner {
  public:
    virtual ~irunner() {}

    virtual Pythia8::Settings* settings() = 0;
    virtual Pythia8::Settings const* settings() const = 0;
    virtual void initialize() = 0;
    virtual void run(std::vector<idt::runjam::Particle> const& initial_state, std::vector<idt::runjam::Particle>& final_state) = 0;

    virtual double get_particle_mass(int pdg) const = 0;
    virtual int get_particle_stable_code(int pdg) const = 0;

    //! The expected number of collisions in the last event, i.e., ncoll / oversample
    virtual double get_event_collision_number() const = 0;
  };

  std::unique_ptr<irunner> create_runner(idt::runjam::runjam_context const&, std::string const& input_filename = "/dev/null");
  std::string version_string();
}
