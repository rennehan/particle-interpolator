double cubic_kernel(double u, double hinv3);
int real_idx(int idx, int N_cell);
double interpolate_to_grid(double *pos_x,
                           double *pos_y,
                           double *pos_z, 
                           double *radii,
                           double *quantities,
                           double *weights,
                           int N_cell,
                           int N_particles);
