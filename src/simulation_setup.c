/*
 * Copyright (c) 2018 Maxime Schmitt <max.schmitt@unistra.fr>
 *
 * Redistribution and use in source and binary forms, with or without
 * modification, are permitted provided that the following conditions are met:
 *
 * 1. Redistributions of source code must retain the above copyright notice,
 * this list of conditions and the following disclaimer.
 *
 * 2. Redistributions in binary form must reproduce the above copyright notice,
 * this list of conditions and the following disclaimer in the documentation
 * and/or other materials provided with the distribution.
 *
 * 3. Neither the name of the copyright holder nor the names of its contributors
 * may be used to endorse or promote products derived from this software without
 * specific prior written permission.
 *
 * THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
 * AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
 * IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
 * ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE
 * LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
 * CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF
 * SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS
 * INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN
 * CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)
 * ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
 * POSSIBILITY OF SUCH DAMAGE.
 */

#include <stdio.h>
#include <stdlib.h>

#include "simulation_setup.h"
#include "solver_data_structure.h"

static const char err_read_from_file[] =
    "Syntax error in configuration file at argument %zu\n";

#define scan_or_fail(a, b, c)                                                  \
  {                                                                            \
    if (fscanf(config_file, a, b) != 1) {                                      \
      fprintf(stderr, err_read_from_file, c);                                  \
      exit(EXIT_FAILURE);                                                      \
    } else {                                                                   \
      c++;                                                                     \
    }                                                                          \
  }

static void initialize_from_file(
    const char *config_path, bool two_dims, size_t *M, size_t *N, size_t *O,
    size_t *sim_steps, FLOAT_TYPE *dt, FLOAT_TYPE *diff, FLOAT_TYPE *visc,
    size_t *num_dens_source, size_t *num_vel_source, size_t **dens_sourceX,
    size_t **dens_sourceY, size_t **dens_sourceZ, size_t **dens_frequency,
    FLOAT_TYPE **dens_add, size_t **vel_sourceX, size_t **vel_sourceY,
    size_t **vel_sourceZ, FLOAT_TYPE **vel_angleXY, FLOAT_TYPE **vel_angleXZ,
    FLOAT_TYPE **vel_force, size_t **vel_frequency) {

  FILE *config_file = fopen(config_path, "r");
  if (config_file == NULL) {
    fprintf(stderr, "Configuration file (%s) error:\n", config_path);
    perror("fopen");
    exit(1);
  }
  size_t curr_arg = 1;

  scan_or_fail("%zu", sim_steps, curr_arg);
  scan_or_fail("%zu", M, curr_arg);
  scan_or_fail("%zu", N, curr_arg);
  if (!two_dims)
    scan_or_fail("%zu", O, curr_arg);
  scan_or_fail("%f", dt, curr_arg);
  scan_or_fail("%f", diff, curr_arg);
  scan_or_fail("%f", visc, curr_arg);

  scan_or_fail("%zu", num_dens_source, curr_arg);
  *dens_sourceX = calloc(*num_dens_source, sizeof(**dens_sourceX));
  *dens_sourceY = calloc(*num_dens_source, sizeof(**dens_sourceY));
  if (!two_dims)
    *dens_sourceZ = calloc(*num_dens_source, sizeof(**dens_sourceZ));
  *dens_frequency = calloc(*num_dens_source, sizeof(**dens_frequency));
  *dens_add = calloc(*num_dens_source, sizeof(**dens_add));
  for (size_t i = 0; i < *num_dens_source; ++i) {
    scan_or_fail("%zu", &(*dens_sourceX)[i], curr_arg);
    scan_or_fail("%zu", &(*dens_sourceY)[i], curr_arg);
    if (!two_dims)
      scan_or_fail("%zu", &(*dens_sourceZ)[i], curr_arg);
    scan_or_fail("%zu", &(*dens_frequency)[i], curr_arg);
    scan_or_fail("%f", &(*dens_add)[i], curr_arg);
    fprintf(stderr, "Source at %zu %zu, frequency %zu value %f\n",
            (*dens_sourceX)[i], (*dens_sourceY)[i], (*dens_frequency)[i],
            (*dens_add)[i]);
  }

#define M_PI 3.14159265358979323846f /* pi */
  scan_or_fail("%zu", num_vel_source, curr_arg);
  *vel_sourceX = calloc(*num_vel_source, sizeof(**vel_sourceX));
  *vel_sourceY = calloc(*num_vel_source, sizeof(**vel_sourceY));
  if (!two_dims)
    *vel_sourceZ = calloc(*num_vel_source, sizeof(**vel_sourceZ));
  *vel_frequency = calloc(*num_vel_source, sizeof(**vel_frequency));
  *vel_angleXY = calloc(*num_vel_source, sizeof(**vel_angleXY));
  if (!two_dims)
    *vel_angleXZ = calloc(*num_vel_source, sizeof(**vel_angleXZ));
  *vel_force = calloc(*num_vel_source, sizeof(**vel_force));
  for (size_t i = 0; i < *num_vel_source; ++i) {
    scan_or_fail("%zu", &(*vel_sourceX)[i], curr_arg);
    scan_or_fail("%zu", &(*vel_sourceY)[i], curr_arg);
    if (!two_dims)
      scan_or_fail("%zu", &(*vel_sourceZ)[i], curr_arg);
    scan_or_fail("%zu", &(*vel_frequency)[i], curr_arg);
    scan_or_fail("%f", &(*vel_angleXY)[i], curr_arg);
    (*vel_angleXY)[i] *= M_PI / 180.f;
    if (!two_dims) {
      scan_or_fail("%f", &(*vel_angleXZ)[i], curr_arg);
      (*vel_angleXZ)[i] *= M_PI / 180.f;
    }
    scan_or_fail("%f", &(*vel_force)[i], curr_arg);
    fprintf(stderr, "Force at %zu %zu, frequency %zu angle %f value %f\n",
            (*vel_sourceX)[i], (*vel_sourceY)[i], (*vel_frequency)[i],
            (*vel_angleXY)[i], (*vel_force)[i]);
  }
  fclose(config_file);
#undef M_PI
}

static void alloc_cons_var_array(struct solverData *data) {
  if (data->dim == solver_dimension_twoD) {
    consVarLocals2D(data);
    data->cons_var.twoD.dens = calloc(1, sizeof(FLOAT_TYPE[M][N]));
    data->cons_var.twoD.velX = calloc(1, sizeof(FLOAT_TYPE[M][N]));
    data->cons_var.twoD.velY = calloc(1, sizeof(FLOAT_TYPE[M][N]));
    data->cons_var.twoD.dens_prev = calloc(1, sizeof(FLOAT_TYPE[M][N]));
    data->cons_var.twoD.velX_prev = calloc(1, sizeof(FLOAT_TYPE[M][N]));
    data->cons_var.twoD.velY_prev = calloc(1, sizeof(FLOAT_TYPE[M][N]));
  } else {
    consVarLocals3D(data);
    data->cons_var.threeD.dens = calloc(1, sizeof(FLOAT_TYPE[M][N][O]));
    data->cons_var.threeD.velX = calloc(1, sizeof(FLOAT_TYPE[M][N][O]));
    data->cons_var.threeD.velY = calloc(1, sizeof(FLOAT_TYPE[M][N][O]));
    data->cons_var.threeD.velZ = calloc(1, sizeof(FLOAT_TYPE[M][N][O]));
    data->cons_var.threeD.dens_prev = calloc(1, sizeof(FLOAT_TYPE[M][N][O]));
    data->cons_var.threeD.velX_prev = calloc(1, sizeof(FLOAT_TYPE[M][N][O]));
    data->cons_var.threeD.velY_prev = calloc(1, sizeof(FLOAT_TYPE[M][N][O]));
    data->cons_var.threeD.velZ_prev = calloc(1, sizeof(FLOAT_TYPE[M][N][O]));
  }
}

static void free_cons_var(struct solverData *data) {
  if (data->dim == solver_dimension_twoD) {
    consVarLocals2D(data);
    free(velX);
    free(velY);
    free(dens);
    free(velX_prev);
    free(velY_prev);
    free(dens_prev);
  } else {
    consVarLocals3D(data);
    free(velX);
    free(velY);
    free(velZ);
    free(dens);
    free(velX_prev);
    free(velY_prev);
    free(velZ_prev);
    free(dens_prev);
  }
}

static void free_simumation_specific(struct solverData *data) {
  if (data->dim == solver_dimension_twoD) {
    simDeclLocals2D(data);
    free(dens_add);
    free(dens_sourceX);
    free(dens_sourceY);
    free(vel_sourceX);
    free(vel_sourceY);
    free(vel_frequency);
    free(dens_frequency);
    free(vel_angleXY);
    free(vel_force);
  } else {
    simDeclLocals3D(data);
    free(dens_add);
    free(dens_sourceX);
    free(dens_sourceY);
    free(dens_sourceZ);
    free(vel_sourceX);
    free(vel_sourceY);
    free(vel_sourceZ);
    free(vel_frequency);
    free(dens_frequency);
    free(vel_angleXY);
    free(vel_angleXZ);
    free(vel_force);
  }
}

void initialize_solver(const char *config_path, size_t *num_steps,
                       struct solverData *data) {
  if (data->dim == solver_dimension_twoD) {
    initialize_from_file(
        config_path, true, &data->cons_var.twoD.M, &data->cons_var.twoD.N, 0,
        num_steps, &data->dt, &data->diff, &data->visc,
        &data->sim_decl.twoD.num_dens, &data->sim_decl.twoD.num_vel,
        &data->sim_decl.twoD.dens_sourceX, &data->sim_decl.twoD.dens_sourceY,
        NULL, &data->sim_decl.twoD.dens_frequency,
        &data->sim_decl.twoD.dens_add, &data->sim_decl.twoD.vel_sourceX,
        &data->sim_decl.twoD.vel_sourceY, NULL,
        &data->sim_decl.twoD.vel_angleXY, NULL, &data->sim_decl.twoD.vel_force,
        &data->sim_decl.twoD.vel_frequency);
  } else {
    initialize_from_file(
        config_path, false, &data->cons_var.threeD.M, &data->cons_var.threeD.N,
        &data->cons_var.threeD.O, num_steps, &data->dt, &data->diff,
        &data->visc, &data->sim_decl.threeD.num_dens,
        &data->sim_decl.threeD.num_vel, &data->sim_decl.threeD.dens_sourceX,
        &data->sim_decl.threeD.dens_sourceY,
        &data->sim_decl.threeD.dens_sourceZ,
        &data->sim_decl.threeD.dens_frequency, &data->sim_decl.threeD.dens_add,
        &data->sim_decl.threeD.vel_sourceX, &data->sim_decl.threeD.vel_sourceY,
        &data->sim_decl.threeD.vel_sourceZ, &data->sim_decl.threeD.vel_angleXY,
        &data->sim_decl.threeD.vel_angleXZ, &data->sim_decl.threeD.vel_force,
        &data->sim_decl.threeD.vel_frequency);
  }
  alloc_cons_var_array(data);
}

void free_solver_data(struct solverData *data) {
  free_cons_var(data);
  free_simumation_specific(data);
}

void print_data_to_file(const char *file_path, struct solverData *data) {
  FILE *out = fopen(file_path, "w");
  if (out == NULL) {
    perror("fopen:");
    return;
  }
  if (data->dim == solver_dimension_twoD) {
    consVarLocals2D(data);
    for (size_t i = 0; i < M; ++i) {
      for (size_t j = 0; j < N; ++j) {
        fprintf(out,
                "%zu %zu %." FLOAT_PRINT_PRECISION "e %." FLOAT_PRINT_PRECISION
                "e %." FLOAT_PRINT_PRECISION "e\n",
                i, j, dens[i][j], velX[i][j], velY[i][j]);
      }
    }
  } else {
    consVarLocals3D(data);
    for (size_t i = 0; i < M; ++i) {
      for (size_t j = 0; j < N; ++j) {
        for (size_t k = 0; k < O; ++k) {
          fprintf(out,
                  "%zu %zu %zu %." FLOAT_PRINT_PRECISION
                  "e %." FLOAT_PRINT_PRECISION "e %." FLOAT_PRINT_PRECISION
                  "e %" FLOAT_PRINT_PRECISION "e\n",
                  i, j, k, dens[i][j][k], velX[i][j][k], velY[i][j][k],
                  velZ[i][j][k]);
        }
      }
    }
  }
}
