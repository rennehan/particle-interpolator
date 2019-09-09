cdef extern from "interpolator.h":
    double interpolate_to_grid(double *pos_x,
                               double *pos_y,
                               double *pos_z,
                               double *radii,
                               double *quantities,
                               double *weights,
                               int N_cell,
                               int N_particles);

def interpolate(double[:] pos_x, 
                double[:] pos_y, 
                double[:] pos_z, 
                double[:] radii,
                double[:] quantities,
                double[:] weights,
                N_cell):
    return interpolate_to_grid(&pos_x[0], &pos_y[0], &pos_z[0], &radii[0], &quantities[0], &weights[0], N_cell, len(pos_x))
