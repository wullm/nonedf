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


    /* Size of the box */
    double BoxLen = 0.;

    /* Seed the random number generator */
    rng_state seed = rand_uint64_init(pars.Seed + 1);

    printf("The parameter file is %s %d\n", pars.InputFilename, pars.Snapshots);

    /* Allocate memory for source particles */
    int snaps = pars.Snapshots;
    long long int NPartTot = 0;
    struct particle **sources = malloc(sizeof(struct particle*) * snaps);
    double *a_snapshots = calloc(snaps, sizeof(double));

    /* Load the source particles */
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

        /* Open the Header group */
        hid_t h_grp = H5Gopen(h_file, "Header", H5P_DEFAULT);

        /* Retrieve and store the physical dimensions of the box */
        if (BoxLen == 0) {
            double boxlen[3];
            hid_t h_attr = H5Aopen(h_grp, "BoxSize", H5P_DEFAULT);
            hid_t h_err = H5Aread(h_attr, H5T_NATIVE_DOUBLE, boxlen);
            H5Aclose(h_attr);
            assert(h_err >= 0);
            BoxLen = boxlen[0];
        }

        /* Close the Header group again */
        H5Gclose(h_grp);

        /* Check if the Cosmology group exists */
        hid_t h_status = H5Eset_auto1(NULL, NULL);  //turn off error printing
        h_status = H5Gget_objinfo(h_file, "/Cosmology", 0, NULL);

        /* If the group exists. */
        if (h_status == 0) {
            /* Open the Cosmology group */
            h_grp = H5Gopen(h_file, "Cosmology", H5P_DEFAULT);

            /* Read the redshift attribute */
            double redshift;
            hid_t h_attr = H5Aopen(h_grp, "Redshift", H5P_DEFAULT);
            hid_t h_err = H5Aread(h_attr, H5T_NATIVE_DOUBLE, &redshift);
            H5Aclose(h_attr);
            assert(h_err >= 0);

            printf("The redshift was %f\n", redshift);

            double a = 1./(1+redshift);
            a_snapshots[i] = a;

            /* Close the Cosmology group */
            H5Gclose(h_grp);
        }

        /* Try to open the desired import group */
        char ImportName[50] = "PartType1";

        /* Open the corresponding group */
        h_grp = H5Gopen(h_file, ImportName, H5P_DEFAULT);

        /* Open the coordinates dataset */
        hid_t h_dat = H5Dopen(h_grp, "Coordinates", H5P_DEFAULT);

        /* Find the dataspace (in the file) */
        hid_t h_space = H5Dget_space (h_dat);

        /* Get the dimensions of this dataspace */
        hsize_t dims[2];
        H5Sget_simple_extent_dims(h_space, dims, NULL);

        /* How many particles do we want per slab? */
        hid_t Npart = dims[0];
        Npart = (hid_t) pars.SlabSize*.5;
        NPartTot = Npart;
        hid_t max_slab_size = pars.SlabSize;
        int slabs = Npart/max_slab_size;
        hid_t counter = 0;

        printf("Met %d slabs\n", slabs);

        /* Allocate memory for the particles */
        sources[i] = malloc(Npart * sizeof(struct particle));
        printf("Allocated memory for %lu parts\n", Npart);

        /* Close the data and memory spaces */
        H5Sclose(h_space);

        /* Close the dataset */
        H5Dclose(h_dat);

        int slab_counter = 0;

        for (int k=0; k<slabs+1; k+=1) {
            /* All slabs have the same number of particles, except possibly the last */
            hid_t slab_size = fmin(Npart - k * max_slab_size, max_slab_size);
            counter += slab_size; //the number of particles read

            /* Define the hyperslab */
            hsize_t slab_dims[2], start[2]; //for 3-vectors
            hsize_t slab_dims_one[1], start_one[1]; //for scalars

            /* Slab dimensions for 3-vectors */
            slab_dims[0] = slab_size;
            slab_dims[1] = 3; //(x,y,z)
            start[0] = k * max_slab_size;
            start[1] = 0; //start with x

            /* Slab dimensions for scalars */
            slab_dims_one[0] = slab_size;
            start_one[0] = k * max_slab_size;

            /* Open the coordinates dataset */
            h_dat = H5Dopen(h_grp, "Coordinates", H5P_DEFAULT);

            /* Find the dataspace (in the file) */
            h_space = H5Dget_space (h_dat);

            /* Select the hyperslab */
            hid_t status = H5Sselect_hyperslab(h_space, H5S_SELECT_SET, start,
                                               NULL, slab_dims, NULL);
            assert(status >= 0);


            /* Create a memory space */
            hid_t h_mems = H5Screate_simple(2, slab_dims, NULL);

            /* Create the data array */
            double data[slab_size][3];

            status = H5Dread(h_dat, H5T_NATIVE_DOUBLE, h_mems, h_space, H5P_DEFAULT,
                             data);

            /* Close the memory space */
            H5Sclose(h_mems);

            /* Close the data and memory spaces */
            H5Sclose(h_space);

            /* Close the dataset */
            H5Dclose(h_dat);


            /* Open the masses dataset */
            h_dat = H5Dopen(h_grp, "Masses", H5P_DEFAULT);

            /* Find the dataspace (in the file) */
            h_space = H5Dget_space (h_dat);

            /* Select the hyperslab */
            status = H5Sselect_hyperslab(h_space, H5S_SELECT_SET, start_one, NULL,
                                                slab_dims_one, NULL);

            /* Create a memory space */
            h_mems = H5Screate_simple(1, slab_dims_one, NULL);

            /* Create the data array */
            double mass_data[slab_size];

            status = H5Dread(h_dat, H5T_NATIVE_DOUBLE, h_mems, h_space, H5P_DEFAULT,
                             mass_data);

            /* Close the memory space */
            H5Sclose(h_mems);

            /* Close the data and memory spaces */
            H5Sclose(h_space);

            /* Close the dataset */
            H5Dclose(h_dat);

            /* Store the particles */
            for (int l=0; l<slab_size; l++) {
                sources[i][counter - slab_size + l].x[0] = data[l][0];
                sources[i][counter - slab_size + l].x[1] = data[l][1];
                sources[i][counter - slab_size + l].x[2] = data[l][2];
                sources[i][counter - slab_size + l].mass = mass_data[l];
            }
        }

        /* Close the file */
        H5Fclose(h_file);

    }

    /* Make gravity meshes */
    const int N = pars.GridSize;
    double **grav_meshes = malloc(snaps * sizeof(double*));
    for (int i=0; i<snaps; i++) {
        /* Allocate the grid */
        grav_meshes[i] = calloc(N * N * N, sizeof(double));

        double total_mass = 0.0;
        double grid_cell_vol = (BoxLen * BoxLen * BoxLen) / (N * N * N);

        /* Assign particles with CIC */
        for (int l=0; l<NPartTot; l++) {
            double X = sources[i][l].x[0] / (BoxLen/N);
            double Y = sources[i][l].x[1] / (BoxLen/N);
            double Z = sources[i][l].x[2] / (BoxLen/N);
            double M = sources[i][l].mass;
            total_mass += M;

            int iX = (int) floor(X);
            int iY = (int) floor(Y);
            int iZ = (int) floor(Z);

            double shift = 0;

            //The search window with respect to the top-left-upper corner
    		int lookLftX = (int) floor((X-iX) - 1.5 + shift);
    		int lookRgtX = (int) floor((X-iX) + 1.5 + shift);
    		int lookLftY = (int) floor((Y-iY) - 1.5 + shift);
    		int lookRgtY = (int) floor((Y-iY) + 1.5 + shift);
    		int lookLftZ = (int) floor((Z-iZ) - 1.5 + shift);
    		int lookRgtZ = (int) floor((Z-iZ) + 1.5 + shift);

            //Do the mass assignment
    		for (int x=lookLftX; x<=lookRgtX; x++) {
    			for (int y=lookLftY; y<=lookRgtY; y++) {
    				for (int z=lookLftZ; z<=lookRgtZ; z++) {
                        double xx = fabs(X - (iX+x+shift));
                        double yy = fabs(Y - (iY+y+shift));
                        double zz = fabs(Z - (iZ+z+shift));

                        double part_x = xx < 0.5 ? (0.75-xx*xx)
                                                : (xx < 1.5 ? 0.5*(1.5-xx)*(1.5-xx) : 0);
        				double part_y = yy < 0.5 ? (0.75-yy*yy)
                                                : (yy < 1.5 ? 0.5*(1.5-yy)*(1.5-yy) : 0);
        				double part_z = zz < 0.5 ? (0.75-zz*zz)
                                                : (zz < 1.5 ? 0.5*(1.5-zz)*(1.5-zz) : 0);

                        grav_meshes[i][row_major(iX+x, iY+y, iZ+z, N)] += M/grid_cell_vol * (part_x*part_y*part_z);
    				}
    			}
    		}
        }

        /* The average density */
        double avg_density = total_mass / (BoxLen * BoxLen * BoxLen);

        /* Multiply by the density grid by Newton's constant */
        for (int j=0; j<N*N*N; j++) {
            grav_meshes[i][j] *=  (4.0 * M_PI) * us.GravityG;
        }

        /* Allocate memory for the Fourier transform */
        fftw_complex *fbox = malloc(N*N*(N/2+1) * sizeof(fftw_complex));

        /* Fourier transform */
        fftw_plan r2c = fftw_plan_dft_r2c_3d(N, N, N, grav_meshes[i], fbox, FFTW_ESTIMATE);
        fft_execute(r2c);
        fft_normalize_r2c(fbox, N, BoxLen);
        fftw_destroy_plan(r2c);

        /* Apply the inverse Poisson kernel -1/k^2 */
        fft_apply_kernel(fbox, fbox, N, BoxLen, kernel_inv_poisson, NULL);

        /* Fourier transform back */
        fftw_plan c2r = fftw_plan_dft_c2r_3d(N, N, N, fbox, grav_meshes[i], FFTW_ESTIMATE);
        fft_execute(c2r);
        fft_normalize_c2r(grav_meshes[i], N, BoxLen);
        fftw_destroy_plan(c2r);

        /* Free the complex grid */
        free(fbox);
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


    /* Generate random neutrino particles */
    for (int i=0; i<1; i++) {
        double a = cosmo.a_begin;
        double x[3];
        double v[3];

        /* Generate a random velocity from the Fermi-Dirac distribution */
        generateVelocity(&thermal_sampler, &us, &cosmo, &seed, a, v);

        /* Compute initial phase space density */
        double f_i = fermi_dirac_density(&us, &cosmo, a, v, cosmo.M_nu[0]);
        double f = f_i;

        /* The initial particle mass */
        double m_i = 1.0;
        double w = 0; //delta-f weight
        double m = 1e-20;


        /* Generate a random point in the central ball with radius R_nu */
        double R_nu = pars.CentralRadius;

        /* We first generate a random point on the sphere using Gaussians */
        x[0] = sampleNorm(&seed);
        x[1] = sampleNorm(&seed);
        x[2] = sampleNorm(&seed);

        /* And normalize */
        const double r_length = hypot(x[0], hypot(x[1], x[2]));
        if (r_length > 0) {
            x[0] /= r_length;
            x[1] /= r_length;
            x[2] /= r_length;
        }

        /* Next, use the last uniform random variate for the radial coordinate */
        double r = cbrt(sampleUniform(&seed)) * R_nu;

        /* Apply the radial coordinate and map to the centre of the rectangle */
        x[0] = x[0] * r + BoxLen * 0.5;
        x[1] = x[1] * r + BoxLen * 0.5;
        x[0] = x[2] * r + BoxLen * 0.5;


        /* Integrate the particle forward */
        double a_factor = 1.03;
        for (int step=0; step<40; step++) {
            double p = fermi_dirac_momentum(&us, a, v, cosmo.M_nu[0]);
            printf("%d %f %f %f %f %e %e\n", step, a, x[0], x[1], x[2], p, w);

            /* Find the bounding snapshots */
            int j;
            for (j=0; j<snaps; j++) {
                if (a_snapshots[j] > a) break;
            }
            int ind_prev = j-1;
            int ind_next = j;
            double a_prev = a_snapshots[ind_prev];
            double a_next = a_snapshots[ind_next];
            double delta = (a - a_prev) / (a_next - a_prev);

            /* Get the accelerations at the current location */
            double acc_prev[3];
            double acc_next[3];
            accelCIC(grav_meshes[ind_prev], N, BoxLen, x[0], x[1], x[2], acc_prev);
            accelCIC(grav_meshes[ind_next], N, BoxLen, x[0], x[1], x[2], acc_next);

            /* Interpolate the acceleration */
            double acc[3] = {acc_prev[0] + delta * (acc_next[0] - acc_prev[0]),
                             acc_prev[1] + delta * (acc_next[1] - acc_prev[1]),
                             acc_prev[2] + delta * (acc_next[2] - acc_prev[2])};

            /* Fetch the kick and drift factors */
            double kick_factor = get_kick_factor(&cosmo, log(a), log(a * a_factor));
            double drift_factor = get_drift_factor(&cosmo, log(a), log(a * a_factor));

            /* Fetch the relativistic correction factors */
            double relat_kick_correction = relativity_kick(v, a, &us);
            double relat_drift_correction = relativity_drift(v, a, &us);

            kick_factor *= relat_kick_correction;
            drift_factor *= relat_drift_correction;

            /* Execute kick */
            v[0] += acc[0] * kick_factor;
            v[1] += acc[1] * kick_factor;
            v[2] += acc[2] * kick_factor;

            /* Apply delta-f step */
            f = fermi_dirac_density(&us, &cosmo, a, v, cosmo.M_nu[0]);
            w = (f - f_i) / f_i;
            m = w * m_i;

            /* Execute drift */
            x[0] += v[0] * drift_factor;
            x[1] += v[1] * drift_factor;
            x[2] += v[2] * drift_factor;

            /* Ensure that particles wrap */
            x[0] = fwrap(x[0], BoxLen);
            x[1] = fwrap(x[1], BoxLen);
            x[2] = fwrap(x[2], BoxLen);

            /* If the particle leaves the central sphere, map it around */
            double r_x = x[0] - 0.5 * BoxLen;
            double r_y = x[1] - 0.5 * BoxLen;
            double r_z = x[2] - 0.5 * BoxLen;
            double r2 = r_x * r_x + r_y * r_y + r_z * r_z;
            if (r2 > R_nu * R_nu) {
                /* Map to its antipodal point */
                x[0] = 0.5 * BoxLen - r_x;
                x[1] = 0.5 * BoxLen - r_y;
                x[2] = 0.5 * BoxLen - r_z;

                /* Sample a new velocity */
                generateVelocity(&thermal_sampler, &us, &cosmo, &seed, a, v);
                /* Re-compute the initial phase space density */
                f_i = fermi_dirac_density(&us, &cosmo, a, v, cosmo.M_nu[0]);
                f = f_i;
                w = 0;
            }

            /* Step forward */
            a *= a_factor;

        }
    }

    /* Clean the random sampler */
    cleanSampler(&thermal_sampler);

    /* Free the source particle data */
    for (int i=0; i<snaps; i++) {
        free(sources[i]);
        free(grav_meshes[i]);
    }
    free(sources);
    free(grav_meshes);
    free(a_snapshots);

    /* Clean up */
    cleanParams(&pars);
    cleanCosmology(&cosmo);
}
