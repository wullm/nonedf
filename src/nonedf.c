/*******************************************************************************
 * This file is part of Nonedf.
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

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <hdf5.h>
#include <assert.h>
#include <math.h>

#include "../include/nonedf.h"

int main(int argc, char *argv[]) {
    if (argc == 1) {
        printf("No parameter file specified.\n");
        return 0;
    }

    /* Read options */
    const char *fname = argv[1];
    printf("The parameter file is %s\n", fname);

    struct params pars;
    struct units us;
    struct cosmology cosmo = {0};

    readParams(&pars, fname);
    readUnits(&us, fname);

    /* Seed the random number generator */
    rng_state seed = rand_uint64_init(pars.Seed + 1);

    printf("The parameter file is %s %d\n", pars.InputFilename, pars.Snapshots);

    /* Load the source particles */
    int snaps = pars.Snapshots;
    struct particle **sources = malloc(sizeof(struct particle*) * snaps);

    for (int i=0; i<snaps; i++) {
        char fname_snap[100];
        sprintf(fname_snap, "%s_%04d.hdf5", pars.InputFilename, i);
        printf("%s\n", fname_snap);

        /* Open the file */
        hid_t h_file = H5Fopen(fname_snap, H5F_ACC_RDONLY, H5P_DEFAULT);

        /* Read cosmology if this has not happened yet */
        if (cosmo.h == 0) {
            readCosmology(&cosmo, &us, h_file);
        }

        /* Close the file */
        H5Fclose(h_file);


        //
        // /* Open the Header group */
        // hid_t h_grp = H5Gopen(h_file, "Header", H5P_DEFAULT);
        //
        // /* Read the physical dimensions of the box */
        // double boxlen[3];
        // hid_t h_attr = H5Aopen(h_grp, "BoxSize", H5P_DEFAULT);
        // hid_t h_err = H5Aread(h_attr, H5T_NATIVE_DOUBLE, boxlen);
        // H5Aclose(h_attr);
        // assert(h_err >= 0);
        //
        // /* Read the numbers of particles of each type */
        // hsize_t numer_of_types;
        // h_attr = H5Aopen(h_grp, "NumPart_Total", H5P_DEFAULT);
        // hid_t h_atspace = H5Aget_space(h_attr);
        // H5Sget_simple_extent_dims(h_atspace, &numer_of_types, NULL);
        // H5Sclose(h_atspace);
        // H5Aclose(h_attr);
        //
        // /* Close the Header group again */
        // H5Gclose(h_grp);
        //
        // /* Check if the Cosmology group exists */
        // hid_t h_status = H5Eset_auto1(NULL, NULL);  //turn off error printing
        // h_status = H5Gget_objinfo(h_file, "/Cosmology", 0, NULL);
        //
        // /* If the group exists. */
        // if (h_status == 0) {
        //     /* Open the Cosmology group */
        //     h_grp = H5Gopen(h_file, "Cosmology", H5P_DEFAULT);
        //
        //     /* Read the redshift attribute */
        //     double redshift;
        //     h_attr = H5Aopen(h_grp, "Redshift", H5P_DEFAULT);
        //     h_err = H5Aread(h_attr, H5T_NATIVE_DOUBLE, &redshift);
        //     H5Aclose(h_attr);
        //     assert(h_err >= 0);
        //
        //     printf("The redshift was %f\n", redshift);
        //
        //     /* Close the Cosmology group */
        //     H5Gclose(h_grp);
        // }
        //
        // /* Try to open the desired import group */
        //
        // char ImportName[50] = "PartType1";
        //
        // /* Open the corresponding group */
        // h_grp = H5Gopen(h_file, ImportName, H5P_DEFAULT);
        //
        // /* Open the coordinates dataset */
        // hid_t h_dat = H5Dopen(h_grp, "Coordinates", H5P_DEFAULT);
        //
        // /* Find the dataspace (in the file) */
        // hid_t h_space = H5Dget_space (h_dat);
        //
        // /* Get the dimensions of this dataspace */
        // hsize_t dims[2];
        // H5Sget_simple_extent_dims(h_space, dims, NULL);
        //
        // /* How many particles do we want per slab? */
        // hid_t Npart = dims[0];
        // hid_t max_slab_size = pars.SlabSize;
        // int slabs = Npart/max_slab_size;
        // hid_t counter = 0;
        //
        // /* Close the data and memory spaces */
        // H5Sclose(h_space);
        //
        // /* Close the dataset */
        // H5Dclose(h_dat);
        //
        // double total_mass = 0; //for this particle type
        //
        // int slab_counter = 0;
        //
        // for (int k=0; k<slabs+1; k+=1) {
        //     /* All slabs have the same number of particles, except possibly the last */
        //     hid_t slab_size = fmin(Npart - k * max_slab_size, max_slab_size);
        //     counter += slab_size; //the number of particles read
        //
        //     /* Define the hyperslab */
        //     hsize_t slab_dims[2], start[2]; //for 3-vectors
        //     hsize_t slab_dims_one[1], start_one[1]; //for scalars
        //
        //     /* Slab dimensions for 3-vectors */
        //     slab_dims[0] = slab_size;
        //     slab_dims[1] = 3; //(x,y,z)
        //     start[0] = k * max_slab_size;
        //     start[1] = 0; //start with x
        //
        //     /* Slab dimensions for scalars */
        //     slab_dims_one[0] = slab_size;
        //     start_one[0] = k * max_slab_size;
        //
        //     /* Open the coordinates dataset */
        //     h_dat = H5Dopen(h_grp, "Coordinates", H5P_DEFAULT);
        //
        //     /* Find the dataspace (in the file) */
        //     h_space = H5Dget_space (h_dat);
        //
        //     /* Select the hyperslab */
        //     hid_t status = H5Sselect_hyperslab(h_space, H5S_SELECT_SET, start,
        //                                        NULL, slab_dims, NULL);
        //     assert(status >= 0);
        //
        //
        //     /* Create a memory space */
        //     hid_t h_mems = H5Screate_simple(2, slab_dims, NULL);
        //
        //     /* Create the data array */
        //     double data[slab_size][3];
        //
        //     status = H5Dread(h_dat, H5T_NATIVE_DOUBLE, h_mems, h_space, H5P_DEFAULT,
        //                      data);
        //
        //     /* Close the memory space */
        //     H5Sclose(h_mems);
        //
        //     /* Close the data and memory spaces */
        //     H5Sclose(h_space);
        //
        //     /* Close the dataset */
        //     H5Dclose(h_dat);
        //
        //
        //     /* Open the masses dataset */
        //     h_dat = H5Dopen(h_grp, "Masses", H5P_DEFAULT);
        //
        //     /* Find the dataspace (in the file) */
        //     h_space = H5Dget_space (h_dat);
        //
        //     /* Select the hyperslab */
        //     status = H5Sselect_hyperslab(h_space, H5S_SELECT_SET, start_one, NULL,
        //                                         slab_dims_one, NULL);
        //
        //     /* Create a memory space */
        //     h_mems = H5Screate_simple(1, slab_dims_one, NULL);
        //
        //     /* Create the data array */
        //     double mass_data[slab_size];
        //
        //     status = H5Dread(h_dat, H5T_NATIVE_DOUBLE, h_mems, h_space, H5P_DEFAULT,
        //                      mass_data);
        //
        //     /* Close the memory space */
        //     H5Sclose(h_mems);
        //
        //     /* Close the data and memory spaces */
        //     H5Sclose(h_space);
        //
        //     /* Close the dataset */
        //     H5Dclose(h_dat);
        // }
        //
        // /* Close the file */
        // H5Fclose(h_file);

    }

    /* Random sampler used for thermal species */
    struct sampler thermal_sampler;

    printf("The neutrino temperature is %f\n", cosmo.T_nu);

    /* Initialize the sampler */
    double T_eV = cosmo.T_nu * us.kBoltzmann / us.ElectronVolt;
    double thermal_params[2] = {T_eV, 0.0};
    /* Rescale the domain */
    double xl = THERMAL_MIN_MOMENTUM * T_eV; //units of kb*T
    double xr = THERMAL_MAX_MOMENTUM * T_eV; //units of kb*T

    int err = initSampler(&thermal_sampler, fd_pdf, xl, xr, thermal_params);
    if (err > 0) {
        printf("Error initializing the thermal motion sampler.\n");
        exit(1);
    }

    double p0_eV = samplerCustom(&thermal_sampler, &seed); //present-day momentum
    printf("%f\n", p0_eV);

    E_z(1.0, &cosmo);

    double fac = get_kick_factor(&cosmo, log(0.55), log(0.57));
    printf("The kick factor is %e\n", fac);

    /* Clean the random sampler */
    cleanSampler(&thermal_sampler);

    /* Free the source particle data */
    free(sources);

    /* Clean up */
    cleanParams(&pars);
    cleanCosmology(&cosmo);
}
