/*******************************************************************************
 * This file is part of Nonedf.
 * Copyright (c) 2021 Willem Elbers (whe@willemelbers.com)
 *
 * This program is free software: you can redistribute it and/or modify
 * it under the terms of the GNU Lesser General Public License as published
 * by the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU Lesser General Public License
 * along with this program.  If not, see <http://www.gnu.org/licenses/>.
 *
 ******************************************************************************/

#include <math.h>
#include "../include/phase_space.h"

double fermi_dirac_density(struct units *us, struct cosmology *cosmo, double a,
                           double *V, double m_eV) {

  /* Convert temperature to eV (to prevent overflows)*/
  const double k_b = us->kBoltzmann;
  const double eV = us->ElectronVolt;
  const double T_nu = cosmo->T_nu;
  const double T_eV = k_b * T_nu / eV;  // temperature in eV

  /* Calculate the momentum in eV */
  double p_eV = fermi_dirac_momentum(us, a, V, m_eV);

  printf("%e\n", T_eV);

  return 1.0 / (exp(p_eV / T_eV) + 1.0);
}

double fermi_dirac_momentum(struct units *us, double a, double *V, double m_eV) {
  double c = us->SpeedOfLight;
  double v = hypot(V[0], hypot(V[1], V[2]));
  double p0_eV = v * m_eV / c; // present-day momentum in eV

  return p0_eV;
}
