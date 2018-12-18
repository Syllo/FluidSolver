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

#include "solver.h"
#include "solver2D.h"
#include "solver3D.h"

void solver_step(size_t step, struct solverData *data) {
  if (data->dim == solver_dimension_twoD) {
    simDeclLocals2D(data);
    consVarLocals2D(data);
    add_force_source2D(step, M, N, dens_prev, velX_prev, velY_prev, num_dens,
                       dens_add, dens_sourceX, dens_sourceY, dens_frequency,
                       num_vel, vel_force, vel_sourceX, vel_sourceY,
                       vel_angleXY, vel_frequency);
    vel_step2D(M, N, velX, velY, velX_prev, velY_prev, data->visc, data->dt);
    dens_step2D(M, N, dens, dens_prev, velX, velY, data->diff, data->dt);

  } else {
    simDeclLocals3D(data);
    consVarLocals3D(data);
    add_force_source3D(step, M, N, O, dens_prev, velX_prev, velY_prev,
                       velZ_prev, num_dens, dens_add, dens_sourceX,
                       dens_sourceY, dens_sourceZ, dens_frequency, num_vel,
                       vel_force, vel_sourceX, vel_sourceY, vel_sourceZ,
                       vel_angleXY, vel_angleXZ, vel_frequency);
    vel_step3D(M, N, O, velX, velY, velZ, velX_prev, velY_prev, velZ_prev,
               data->visc, data->dt);
    dens_step3D(M, N, O, dens, dens_prev, velX, velY, velZ, data->diff,
                data->dt);
  }
}
