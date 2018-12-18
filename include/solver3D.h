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

#ifndef SOLVER_3D_H_
#define SOLVER_3D_H_
#include <stddef.h>

void vel_step3D(size_t M, size_t N, size_t O, float u[M][N][O],
                float v[M][N][O], float w[M][N][O], float u0[M][N][O],
                float v0[M][N][O], float w0[M][N][O], float visc, float dt);

void dens_step3D(size_t M, size_t N, size_t O, float x[M][N][O],
                 float x0[M][N][O], float u[M][N][O], float v[M][N][O],
                 float w[M][N][O], float diff, float dt);

void add_force_source3D(
    size_t sim_step, size_t M, size_t N, size_t O, float d[restrict M][N][O],
    float u[restrict M][N][O], float v[restrict M][N][O],
    float w[restrict M][N][O], size_t num_dens,
    float dens_add[restrict num_dens], size_t dens_sourceX[restrict num_dens],
    size_t dens_sourceY[restrict num_dens],
    size_t dens_sourceZ[restrict num_dens],
    size_t dens_frequency[restrict num_dens], size_t num_vel,
    float vel_force[restrict num_vel], size_t vel_sourceX[restrict num_vel],
    size_t vel_sourceY[restrict num_vel], size_t vel_sourceZ[restrict num_vel],
    float vel_angleXY[restrict num_vel], float vel_angleXZ[restrict num_vel],
    size_t vel_frequency[restrict num_vel]);

#endif // SOLVER_3D_H_
