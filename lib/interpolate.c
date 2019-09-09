#include "interpolate.h"
#include <math.h>
#include <stdlib.h>

double cubic_kernel(double u, double hinv3)
{
    double wk;

    if (u < 0.5)
    {
        wk = (1.0 + 6.0 * (u - 1.0) * u * u);
    } else {
        double t1 = 1.0 - u;
        double t2 = t1 * t1;
        wk = 2.0 * t2 * t1;
    }

    return (8.0 / M_PI) * wk * hinv3;
}

int real_idx(int idx, int N_cell)
{
    if (idx >= N_cell)
    {
        return idx - N_cell;
    }

    if (idx < 0)
    {
        return idx + N_cell;
    }

    return idx;
}

void interpolate_to_grid(char tmp_file_name[],
                         double *pos_x,
                         double *pos_y,
                         double *pos_z, 
                         double *radii,
                         double *quantities,
                         double *weights,
                         int N_cell,
                         int N_particles)
{
    int i, j, k, l;
    double U[N_cell];

    double size = 1;
    double min = 0;
    double half = 0.5;
    double delta = size / N_cell;

    double map[N_cell][N_cell][N_cell], map_weights[N_cell][N_cell][N_cell];
 
    for (i = 0; i < N_cell; i++)
    {
        U[i] = min + (i + 0.5) * delta;

        for (j = 0; j < N_cell; j++)
        {
            for (k = 0; k < N_cell; k++)
            {
                map[i][j][k] = 0;
                map_weights[i][j][k] = 0;
            }
        }
    }
 
    /* Important that we assume all coordinates are within [0,1] */
    for (i = 0; i < N_particles; i++)
    {
        if (radii[i] < delta)
        {
            radii[i] = 2.0 * delta;
        }

        /* These are the indices in the U vector */
        int max_x = (int)((pos_x[i] + radii[i]) / delta);
        int min_x = (int)((pos_x[i] - radii[i]) / delta);

        int max_y = (int)((pos_y[i] + radii[i]) / delta);
        int min_y = (int)((pos_y[i] - radii[i]) / delta);
    
        int max_z = (int)((pos_z[i] + radii[i]) / delta);
        int min_z = (int)((pos_z[i] - radii[i]) / delta);

        int total_cell_chunk = (max_x - min_x) + (max_y - min_y) + (max_z - min_z);

        int deposit_x[total_cell_chunk], deposit_y[total_cell_chunk], deposit_z[total_cell_chunk];
        double temp_weights[total_cell_chunk]
        double sum_weights = 0;

        double x_diff = y_diff = z_diff = 0;
        int x_idx = y_idx = z_idx = 0;

        int running_idx = 0;

        for (j = min_x; j < max_x; j++)
        {
            x_idx = real_idx(j, N_cell);
            x_diff = pos_x[i] - U[x_idx];

            if (x_diff > half)
            {
                x_diff = pos_x[i] - size - U[x_idx];
            }
            if (x_diff < -half)
            {
                x_diff = pos_x[i] + size - U[x_idx];
            }

            for (k = min_y; k < max_y; k++)
            {
                y_idx = real_idx(k, N_cell);
                y_diff = pos_y[i] - U[y_idx];

                if (y_diff > half)
                {
                    y_diff = pos_y[i] - size - U[y_idx];
                }
                if (y_diff < -half)
                {
                    y_diff = pos_y[i] + size - U[y_idx];
                }

                for (l = min_z; l < max_z; l++)
                {
                    z_idx = real_idx(l, N_cell);
                    z_diff = pos_z[i] - U[z_idx];

                    if (z_diff > half)
                    {
                        z_diff = pos_z[i] - size - U[z_idx];
                    }
                    if (z_diff < -half)
                    {
                        z_diff = pos_z[i] + size - U[z_idx];
                    }

                    distance = sqrt(x_diff * x_diff + y_diff * y_diff + z_diff * z_diff);

                    temp_weights[running_idx] = 0;
                    deposit_x[running_idx] = 0;
                    deposit_y[running_idx] = 0;
                    deposit_z[running_idx] = 0;

                    if (distance < radii[i])
                    {
                        temp_weights[running_idx] = cubic_kernel(distance / radii[i], 1.0 / (radii[i] * radii[i] * radii[i]));
                        sum_weights += temp_weights[running_idx];
                        deposit_x[running_idx] = U[x_idx] / delta;
                        deposit_y[running_idx] = U[y_idx] / delta;
                        deposit_z[running_idx] = U[z_idx] / delta;
                    }

                    running_idx++;
                }
            }
        }

        /* Add quantites to the chunk of cells where they belong */
        for (j = 0; j < running_idx; j++)
        {
            /* It's fine if it's adding to 0 by default because temp_weights == 0 in that case */
            x_idx = deposit_x[j];
            y_idx = deposit_y[j];
            z_idx = deposit_z[j];

            map[x_idx][y_idx][z_idx] += weights[i] * quantities[i] * temp_weights[j];
            map_weights[x_idx][y_idx][z_idx] += weights[i] * temp_weights[j]
        }    
    }

    /* Now we find all of the non-zero entries and write them to a temporary file */
    FILE * tmp_file;

    if ((tmp_file = fopen(tmp_file_name, "w")) == NULL)
    {
        printf("Cannot open temporary file %s\n", tmp_file_name);
        exit(-1);
    }

    for (i = 0; i < N_cell; i++)
    {
        for (j = 0; j < N_cell; j++)
        {
            for (k = 0; k < N_cell; k++)
            {
                if (map[i][j][k] == 0)
                {
                    continue;
                }

                fprintf(tmp_file, "%d\t%d\t%d\t%g\t%g\n", i, j, k, map[i][j][k], map_weights[i][j][k]);
            }
        }
    }

    fclose(tmp_file);

    return;
}

