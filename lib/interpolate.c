#include "interpolate.h"
#include <math.h>
#include <stdlib.h>
#include <stdio.h>

#define SIZE 1.0
#define HALF 0.5
#define MIN 0.0
#define DIMS 2

void out_of_memory(char var_name[], long long bytes)
{
    printf("Out of memory for %s and %lli bytes.\n", var_name, bytes);
    exit(-1);
}

double adjusted_difference(double difference)
{
    if (difference > HALF)
    {
        difference -= SIZE;
    }
    if (difference < -HALF)
    {
        difference += SIZE;
    }

    return difference;
}

/* D. Rennehan: hinv3 could be hinv2 if DIMS == 2 */
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

#if (DIMS == 3)
    return (8.0 / M_PI) * wk * hinv3;
#else
    return (40.0 / (7.0 * M_PI)) * wk * hinv3;
#endif
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
    long long map_idx;
    double U[N_cell];

    long long N_cell2 = N_cell * N_cell;
#if (DIMS == 3)
    long long N_cell3 = N_cell2 * N_cell;
#else
    long long N_cell3 = N_cell2;
#endif

    double delta = SIZE / N_cell;

    long long map_memory_bytes = N_cell3 * sizeof(double);

    double *map = (double *)malloc(map_memory_bytes);
    if (map == NULL)
    {
        out_of_memory("map", map_memory_bytes);
    }

    double *map_weights = (double *)malloc(map_memory_bytes);
    if (map_weights == NULL)
    {
        out_of_memory("map_weights", map_memory_bytes);
    }

    printf("Allocated 2*%lli bytes.\n", map_memory_bytes);

    for (i = 0; i < N_cell; i++)
    {
        U[i] = MIN + (i + 0.5) * delta;
    }

    for (map_idx = 0; map_idx < N_cell3; map_idx++)
    {
        map[map_idx] = 0;
        map_weights[map_idx] = 0;
    }

    printf("Loop over %d particles.\n", N_particles);

    /* Important that we assume all coordinates are within [0,1] */
    for (i = 0; i < N_particles; i++)
    {
        /* These are the indices in the U vector */
        int max_x = (int)((pos_x[i] + radii[i]) / delta);
        int min_x = (int)((pos_x[i] - radii[i]) / delta);

        int max_y = (int)((pos_y[i] + radii[i]) / delta);
        int min_y = (int)((pos_y[i] - radii[i]) / delta);
   
#if (DIMS == 3) 
        int max_z = (int)((pos_z[i] + radii[i]) / delta);
        int min_z = (int)((pos_z[i] - radii[i]) / delta);
#endif

        int num_x_cells = max_x - min_x + 1;
        int num_y_cells = max_y - min_y + 1;
#if (DIMS == 3)
        int num_z_cells = max_z - min_z + 1;
#else
        int num_z_cells = 1;
#endif

        long total_cell_chunk = (long)(num_x_cells * num_y_cells * num_z_cells);
        long long chunk_memory_bytes = total_cell_chunk * sizeof(int);

        int *deposit_x = (int *)malloc(chunk_memory_bytes);
        if (deposit_x == NULL)
        {
            out_of_memory("deposit_x", chunk_memory_bytes);
        }

        int *deposit_y = (int *)malloc(chunk_memory_bytes);
        if (deposit_y == NULL)
        {
            out_of_memory("deposit_y", chunk_memory_bytes);
        }

#if (DIMS == 3)
        int *deposit_z = (int *)malloc(chunk_memory_bytes);
        if (deposit_z == NULL)
        {
            out_of_memory("deposit_z", chunk_memory_bytes);
        }
#endif

        double *temp_weights = (double *)malloc(total_cell_chunk * sizeof(double));
        if (temp_weights == NULL)
        {
            out_of_memory("temp_weights", (long long)(total_cell_chunk * sizeof(double)));
        }

        long long x_int_bytes = (long long)(num_x_cells * sizeof(int));
        long long x_double_bytes = (long long)(num_x_cells * sizeof(double));
        long long y_int_bytes = (long long)(num_y_cells * sizeof(int));
        long long y_double_bytes = (long long)(num_y_cells * sizeof(double));
#if (DIMS == 3)
        long long z_int_bytes = (long long)(num_z_cells * sizeof(int));
        long long z_double_bytes = (long long)(num_z_cells * sizeof(double));
#endif

        int *x_indices = (int *)malloc(x_int_bytes);
        if (x_indices == NULL)
        {
            out_of_memory("x_indices", x_int_bytes);
        }

        int *y_indices = (int *)malloc(y_int_bytes);
        if (y_indices == NULL)
        {
            out_of_memory("y_indices", y_int_bytes);
        }

#if (DIMS == 3)
        int *z_indices = (int *)malloc(z_int_bytes);
        if (z_indices == NULL)
        {
            out_of_memory("z_indices", z_int_bytes);
        }
#endif

        double *x_diffs = (double *)malloc(x_double_bytes);
        if (x_diffs == NULL)
        {
            out_of_memory("x_diffs", x_double_bytes);
        }

        double *y_diffs = (double *)malloc(y_double_bytes);
        if (y_diffs == NULL)
        {
            out_of_memory("y_diffs", y_double_bytes);
        }

#if (DIMS == 3)
        double *z_diffs = (double *)malloc(z_double_bytes);
        if (z_diffs == NULL)
        {
            out_of_memory("z_diffs", z_double_bytes);
        }
#endif

        double distance = 0, sum_weights = 0;
        int running_idx = 0;

        for (j = min_x; j <= max_x; j++)
        {
            x_indices[running_idx] = real_idx(j, N_cell);
            x_diffs[running_idx] = adjusted_difference(pos_x[i] - U[x_indices[running_idx]]);
            running_idx++;
        }

        running_idx = 0;

        for (j = min_y; j <= max_y; j++)
        {
            y_indices[running_idx] = real_idx(j, N_cell);
            y_diffs[running_idx] = adjusted_difference(pos_y[i] - U[y_indices[running_idx]]);
            running_idx++;
        }

#if (DIMS == 3)
        running_idx = 0;

        for (j = min_z; j <= max_z; j++)
        {
            z_indices[running_idx] = real_idx(j, N_cell);
            z_diffs[running_idx] = adjusted_difference(pos_z[i] - U[z_indices[running_idx]]);
            running_idx++;
        }
#endif

        running_idx = 0;

        /**
         * Loop over the 3D cube that could have possible contributions from the particle.
         *
         * We are saving computation by calculating the indices and differences above so
         * here we only have to keep track of the running index in the one dimensional
         * arrays above.
         *
         * For each cell, we calculate the distance from the particle to that cell. If the
         * distance to the cell is less than the radius of the particle then we calculate
         * the smoothed quantity at that grid cell and set the weight based on the kernel.
         * Otherwise, the weight will be zero and nothing will be added.
         */
        if (num_x_cells == 1 && num_y_cells == 1 && num_z_cells == 1)
        {
            /* All of the contribution in the one cell */
            temp_weights[0] = 1.0;
            deposit_x[0] = x_indices[0];
            deposit_y[0] = y_indices[0];
#if (DIMS == 3)
            deposit_z[0] = z_indices[0];
#endif

            sum_weights += 1.0;
            running_idx++;

            /* Skip the loop, it only contributes to one cell completely */
            num_x_cells = 0;
            num_y_cells = 0;
            num_z_cells = 0;
        }

        for (j = 0; j < num_x_cells; j++)
        {
            for (k = 0; k < num_y_cells; k++)
            {
#if (DIMS == 3)
                for (l = 0; l < num_z_cells; l++)
                {
                    distance = sqrt(x_diffs[j] * x_diffs[j] + y_diffs[k] * y_diffs[k] + z_diffs[l] * z_diffs[l]);
#else
                    distance = sqrt(x_diffs[j] * x_diffs[j] + y_diffs[k] * y_diffs[k]);
#endif
                    temp_weights[running_idx] = 0;
                    deposit_x[running_idx] = x_indices[j];
                    deposit_y[running_idx] = y_indices[k];
#if (DIMS == 3)
                    deposit_z[running_idx] = z_indices[l];
#endif

                    if (distance < radii[i])
                    {
#if (DIMS == 3)
                        temp_weights[running_idx] = cubic_kernel(distance / radii[i], 1.0 / (radii[i] * radii[i] * radii[i]));
#else
                        temp_weights[running_idx] = cubic_kernel(distance / radii[i], 1.0 / (radii[i] * radii[i]));
#endif
                        sum_weights += temp_weights[running_idx];
                    }

                    running_idx++;
#if (DIMS == 3)
                }
#endif
            }
        }

        double wk;
        double true_weight;

        /* Add quantites to the chunk of cells where they belong */
        for (j = 0; j < running_idx; j++)
        {
            wk = temp_weights[j] / sum_weights;
            true_weight = weights[i] * wk;

#if (DIMS == 3)
            map_idx = N_cell2 * deposit_x[j] + N_cell * deposit_y[j] + deposit_z[j];
#else
            map_idx = N_cell * deposit_x[j] + deposit_y[j];
#endif
   
            map[map_idx] += quantities[i] * true_weight;
            map_weights[map_idx] += true_weight;
        }

#if (DIMS == 3)
        free(z_diffs);
#endif
        free(y_diffs);
        free(x_diffs);
#if (DIMS == 3)
        free(z_indices);
#endif
        free(y_indices);
        free(x_indices);
        free(temp_weights);
#if (DIMS == 3)
        free(deposit_z);
#endif
        free(deposit_y);
        free(deposit_x); 
    }

    printf("Write to file.\n");

    /* Now we find all of the non-zero entries and write them to a temporary file */
    FILE * tmp_file;

    if ((tmp_file = fopen(tmp_file_name, "wb")) == NULL)
    {
        printf("Cannot open temporary file %s\n", tmp_file_name);
        exit(-1);
    }

    /* Write a single number indicating the grid size */
    fwrite(&N_cell, sizeof(int), 1, tmp_file);
    fwrite(map, sizeof(double), N_cell3, tmp_file);
    fwrite(map_weights, sizeof(double), N_cell3, tmp_file);

    fclose(tmp_file);

    free(map_weights);
    free(map);

    return;
}

