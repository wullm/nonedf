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

#ifndef PHASE_SPACE_H
#define PHASE_SPACE_H

#include "input.h"
#include "cosmology.h"

double fermi_dirac_density(struct units *us, struct cosmology *cosmo, double a,
                           double *V, double m_eV);
double fermi_dirac_momentum(struct units *us, double a, double *V, double m_eV);

#endif
