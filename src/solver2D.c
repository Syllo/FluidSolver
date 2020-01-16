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

#include <stdbool.h>
#include <stdlib.h>
#include <string.h>
#include <tgmath.h>

#include "solver2D.h"

#define max(a, b) (((a) > (b)) ? (a) : (b))

enum border_condition {
  border_copy_value,
  border_negate_x,
  border_negate_y,
  border_negate_z,
};

static void add_source(size_t M, size_t N, float x[restrict M][N],
                       float s[restrict M][N], float dt) {
  for (size_t i = 0; i < M; ++i)
    for (size_t j = 0; j < N; ++j)
      x[i][j] += dt * s[i][j];
}

static void set_bnd(size_t M, size_t N, enum border_condition b,
                    float x[M][N]) {

  if (b == border_negate_x) {
    for (size_t i = 1; i < N - 1; ++i) {
      x[0][i] = -x[1][i];
      x[M - 1][i] = -x[M - 2][i];
    }
  } else {
    for (size_t i = 1; i < N - 1; ++i) {
      x[0][i] = x[1][i];
      x[M - 1][i] = x[M - 2][i];
    }
  }
  if (b == border_negate_y) {
    for (size_t i = 1; i < M - 1; ++i) {
      x[i][0] = -x[i][1];
      x[i][N - 1] = -x[i][N - 2];
    }
  } else {
    for (size_t i = 1; i < M - 1; ++i) {
      x[i][0] = x[i][1];
      x[i][N - 1] = x[i][N - 2];
    }
  }

  const float AHalf = 1.f / 2.f;

  x[0][0] = AHalf * (x[0][1] + x[1][0]);
  x[0][N - 1] = AHalf * (x[0][N - 2] + x[1][N - 1]);

  x[M - 1][0] = AHalf * (x[M - 1][1] + x[M - 2][0]);
  x[M - 1][N - 1] = AHalf * (x[M - 1][N - 2] + x[M - 2][N - 1]);
}

extern double rand_skip_percent;
extern bool sort_skip;

struct dataPosDeviation {
  double previous_val;
  double error;
  size_t posX, posY;
};

static int compareDeviation(const void *a, const void *b) {
  const struct dataPosDeviation *data1 = (const struct dataPosDeviation *)a;
  const struct dataPosDeviation *data2 = (const struct dataPosDeviation *)b;
  if (data1->error < data2->error) {
    return -1;
  } else {
    if (data1->error > data2->error) {
      return 1;
    } else {
      return 0;
    }
  }
}

static void lin_solve(size_t M, size_t N, enum border_condition b,
                      float x[restrict M][N], float x0[restrict M][N], float a,
                      float c) {
  struct dataPosDeviation(*dpd)[N - 2] =
      malloc(sizeof(struct dataPosDeviation[M - 2][N - 2]));

  for (size_t k = 0; k < 20; k++) {
    for (size_t i = 1; i < M - 1; ++i) {
      for (size_t j = 1; j < N - 1; ++j) {
        if (!sort_skip && drand48() < rand_skip_percent)
          continue;
        if (sort_skip) {
          dpd[i - 1][j - 1].previous_val = x[i][j];
        }
        x[i][j] = (x0[i][j] + a * (x[i - 1][j] + x[i + 1][j] + x[i][j - 1] +
                                   x[i][j + 1])) /
                  c;
        if (sort_skip) {
          dpd[i - 1][j - 1].error = dpd[i - 1][j - 1].previous_val - x[i][j];
          dpd[i - 1][j - 1].error *= dpd[i - 1][j - 1].error;
          dpd[i - 1][j - 1].posX = i;
          dpd[i - 1][j - 1].posY = j;
        }
      }
    }
    if (sort_skip) {
      qsort(dpd, (M - 2) * (N - 2), sizeof(struct dataPosDeviation),
            compareDeviation);
      double threshold_d = ceil(((M - 2) * (N - 2)) * rand_skip_percent);
      size_t threshold = (size_t)threshold_d;
      size_t whereAmI = 0;
      for (size_t i = 0; whereAmI < threshold && i < M - 2; ++i) {
        for (size_t j = 0; whereAmI < threshold && j < N - 2; ++j, whereAmI++) {
          // simulate skipping the computation for the values that have lowest
          // update derivative
          x[dpd[i][j].posX][dpd[i][j].posY] = dpd[i][j].previous_val;
        }
      }
    }
    set_bnd(M, N, b, x);
  }
  free(dpd);
}

static void diffuse(size_t M, size_t N, enum border_condition b,
                    float x[restrict M][N], float x0[restrict M][N], float diff,
                    float dt) {
  size_t max_size = max(M - 2, N - 2);
  float a = dt * diff * max_size * max_size;
  lin_solve(M, N, b, x, x0, a, 1 + 4 * a);
}

static void advect(size_t M, size_t N, enum border_condition b,
                   float d[restrict M][N], float d0[M][N], float u[M][N],
                   float v[M][N], float dt) {

  float dt0 = dt * (M - 2);
  float dt1 = dt * (N - 2);
  const float AHalf = 0.5f;

  for (size_t i = 1; i < M - 1; ++i) {
    for (size_t j = 1; j < N - 1; ++j) {

      float x = i - dt0 * u[i][j];
      if (x < AHalf)
        x = AHalf;
      if (x > (M - 2) + AHalf)
        x = (M - 2) + AHalf;
      size_t i0 = (size_t)x;
      size_t i1 = i0 + 1;

      float y = j - dt1 * v[i][j];
      if (y < AHalf)
        y = AHalf;
      if (y > (N - 2) + AHalf)
        y = (N - 2) + AHalf;
      size_t j0 = (size_t)y;
      size_t j1 = j0 + 1;
      float s1 = x - i0;
      float s0 = 1.f - s1;
      float t1 = y - j0;
      float t0 = 1.f - t1;
      d[i][j] = s0 * (t0 * d0[i0][j0] + t1 * d0[i0][j1]) +
                s1 * (t0 * d0[i1][j0] + t1 * d0[i1][j1]);
    }
  }
  set_bnd(M, N, b, d);
}

static void project(size_t M, size_t N, float u[restrict M][N],
                    float v[restrict M][N], float p[restrict M][N],
                    float div[restrict M][N]) {

  const float MinusAHalfDividedByM = -0.5f / (M - 2);
  const float MinusAHalfDividedByN = -0.5f / (N - 2);
  for (size_t i = 1; i < M - 1; ++i) {
    for (size_t j = 1; j < N - 1; ++j) {
      div[i][j] = MinusAHalfDividedByM * (u[i + 1][j] - u[i - 1][j]) +
                  MinusAHalfDividedByN * (v[i][j + 1] - v[i][j - 1]);
    }
  }
  set_bnd(M, N, border_copy_value, div);
  memset(p, 0, sizeof(float[M][N]));

  lin_solve(M, N, 0, p, div, 1, 4);

  const float MinusAHalfTimesM = -0.5f * (M - 2);
  const float MinusAHalfTimesN = -0.5f * (N - 2);
  for (size_t i = 1; i < M - 1; ++i) {
    for (size_t j = 1; j < N - 1; ++j) {
      u[i][j] += MinusAHalfTimesM * (p[i + 1][j] - p[i - 1][j]);
      v[i][j] += MinusAHalfTimesN * (p[i][j + 1] - p[i][j - 1]);
    }
  }
  set_bnd(M, N, border_negate_x, u);
  set_bnd(M, N, border_negate_y, v);
}

void dens_step2D(size_t M, size_t N, float x[M][N], float x0[M][N],
                 float u[M][N], float v[M][N], float diff, float dt) {

  add_source(M, N, x, x0, dt);

  diffuse(M, N, border_copy_value, x0, x, diff, dt);

  advect(M, N, border_copy_value, x, x0, u, v, dt);
}

void vel_step2D(size_t M, size_t N, float u[M][N], float v[M][N],
                float u0[M][N], float v0[M][N], float visc, float dt) {

  add_source(M, N, u, u0, dt);
  add_source(M, N, v, v0, dt);

  diffuse(M, N, border_negate_x, u0, u, visc, dt);
  diffuse(M, N, border_negate_y, v0, v, visc, dt);

  project(M, N, u0, v0, u, v);

  advect(M, N, border_negate_x, u, u0, u0, v0, dt);
  advect(M, N, border_negate_y, v, v0, u0, v0, dt);

  project(M, N, u, v, u0, v0);
}

void add_force_source2D(
    size_t sim_step, size_t M, size_t N, float d[restrict M][N],
    float u[restrict M][N], float v[restrict M][N], size_t num_dens,
    float dens_add[restrict num_dens], size_t dens_sourceX[restrict num_dens],
    size_t dens_sourceY[restrict num_dens],
    size_t dens_frequency[restrict num_dens], size_t num_vel,
    float vel_force[restrict num_vel], size_t vel_sourceX[restrict num_vel],
    size_t vel_sourceY[restrict num_vel], float vel_angleXY[restrict num_vel],
    size_t vel_frequency[restrict num_vel]) {

  for (size_t i = 0; i < M; ++i) {
    for (size_t j = 0; j < N; ++j) {
      u[i][j] = v[i][j] = d[i][j] = 0.0f;
    }
  }

  for (size_t i = 0; i < num_dens; ++i) {
    if (sim_step % dens_frequency[i] == 0) {
      d[dens_sourceX[i]][dens_sourceY[i]] = dens_add[i];
    }
  }
  for (size_t i = 0; i < num_vel; ++i) {
    if (sim_step % vel_frequency[i] == 0) {
      float cosVal = cos(vel_angleXY[i]);
      float sinVal = sin(vel_angleXY[i]);

      u[vel_sourceX[i]][vel_sourceY[i]] = cosVal * vel_force[i];
      v[vel_sourceX[i]][vel_sourceY[i]] = sinVal * vel_force[i];
    }
  }
}
