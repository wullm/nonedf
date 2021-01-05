/*******************************************************************************
 * This file is part of Mitos.
 * Copyright (c) 2020 Willem Elbers (whe@willemelbers.com)
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

#include <stdlib.h>
#include <string.h>
#include <hdf5.h>
#include <assert.h>
#include <math.h>
#include "../include/input.h"

int readParams(struct params *pars, const char *fname) {
     pars->Seed = ini_getl("Random", "Seed", 1, fname);

     pars->GridSize = ini_getl("Box", "GridSize", 64, fname);


     pars->MaxParticleTypes = ini_getl("Simulation", "MaxParticleTypes", 1, fname);
     pars->NumParticleTypes = 0; //should not be read, but inferred
     pars->Homogeneous = ini_getbool("Simulation", "Homogeneous", 0, fname);
     pars->MergeDarkMatterBaryons = ini_getbool("PerturbData", "MergeDarkMatterBaryons", 0, fname);

     pars->Snapshots = ini_getl("Read", "Snapshots", 0, fname);
     pars->SlabSize = ini_getl("Read", "SlabSize", 10000, fname);

     /* Read strings */
     int len = DEFAULT_STRING_LENGTH;
     pars->OutputDirectory = malloc(len);
     pars->Name = malloc(len);
     pars->InputFilename = malloc(len);
     pars->OutputFilename = malloc(len);
     ini_gets("Output", "Directory", "./output", pars->OutputDirectory, len, fname);
     ini_gets("Simulation", "Name", "No Name", pars->Name, len, fname);
     ini_gets("Output", "Filename", "particles.hdf5", pars->OutputFilename, len, fname);
     ini_gets("Read", "Filename", "", pars->InputFilename, len, fname);

     return 0;
}

int readUnits(struct units *us, const char *fname) {
    /* Internal units */
    us->UnitLengthMetres = ini_getd("Units", "UnitLengthMetres", 1.0, fname);
    us->UnitTimeSeconds = ini_getd("Units", "UnitTimeSeconds", 1.0, fname);
    us->UnitMassKilogram = ini_getd("Units", "UnitMassKilogram", 1.0, fname);
    us->UnitTemperatureKelvin = ini_getd("Units", "UnitTemperatureKelvin", 1.0, fname);
    us->UnitCurrentAmpere = ini_getd("Units", "UnitCurrentAmpere", 1.0, fname);

    /* Some physical constants */
    us->SpeedOfLight = SPEED_OF_LIGHT_METRES_SECONDS * us->UnitTimeSeconds
                        / us->UnitLengthMetres;
    us->GravityG = GRAVITY_G_SI_UNITS * us->UnitTimeSeconds * us->UnitTimeSeconds
                    / us->UnitLengthMetres / us->UnitLengthMetres / us->UnitLengthMetres
                    * us->UnitMassKilogram; // m^3 / kg / s^2 to internal
    us->hPlanck = PLANCK_CONST_SI_UNITS / us->UnitMassKilogram / us->UnitLengthMetres
                    / us->UnitLengthMetres * us->UnitTimeSeconds; //J*s = kg*m^2/s
    us->kBoltzmann = BOLTZMANN_CONST_SI_UNITS / us->UnitMassKilogram / us->UnitLengthMetres
                    / us->UnitLengthMetres * us->UnitTimeSeconds * us->UnitTimeSeconds
                    * us->UnitTemperatureKelvin; //J/K = kg*m^2/s^2/K
    us->ElectronVolt = ELECTRONVOLT_SI_UNITS / us->UnitMassKilogram / us->UnitLengthMetres
                    / us->UnitLengthMetres * us->UnitTimeSeconds
                    * us->UnitTimeSeconds; // J = kg*m^2/s^2

    return 0;
}


int cleanParams(struct params *pars) {
    free(pars->OutputDirectory);
    free(pars->Name);
    free(pars->InputFilename);
    free(pars->OutputFilename);

    return 0;
}


// hid_t openFile_MPI(MPI_Comm comm, const char *fname) {
//     /* Property list for MPI file access */
//     hid_t prop_faxs = H5Pcreate(H5P_FILE_ACCESS);
//     H5Pset_fapl_mpio(prop_faxs, comm, MPI_INFO_NULL);
//
//     /* Open the hdf5 file */
//     hid_t h_file = H5Fopen(fname, H5F_ACC_RDWR, prop_faxs);
//     H5Pclose(prop_faxs);
//
//     return h_file;
// }
