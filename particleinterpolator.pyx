cdef extern from "interpolator.h":
    void interpolate_to_grid(char tmp_file_name[],
                             double *pos_x,
                             double *pos_y,
                             double *pos_z,
                             double *radii,
                             double *quantities,
                             double *weights,
                             int N_cell,
                             int N_particles);

def interpolate(tmp_file_name,
                double[:] pos_x, 
                double[:] pos_y, 
                double[:] pos_z, 
                double[:] radii,
                double[:] quantities,
                double[:] weights,
                N_cell):
    interpolate_to_grid(tmp_file_name,
                        &pos_x[0],
                        &pos_y[0],
                        &pos_z[0],
                        &radii[0],
                        &quantities[0],
                        &weights[0],
                        N_cell,
                        len(pos_x))
