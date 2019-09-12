import particleinterpolator
import numpy as np
import matplotlib.pyplot as plt
import argparse as ap


parser = ap.ArgumentParser()
parser.add_argument('--generate', help = 'Generate the projection.', action = 'store_true')
parser.add_argument('--plot', help = 'Plot the resulting slice as a test.', action = 'store_true')
args = parser.parse_args()

# N_cell in the grid
N = 128
N_part = 3

if args.generate:
    h = 0.4

    x = np.array([0.25, 0.75], dtype = float)
    y = np.array([0.5, 0.5], dtype = float)
    z = np.array([0.5, 0.5], dtype = float)

    # Radii, quantity to deposit, weights
    r = h * np.ones(N_part, dtype = float)
    q = np.ones(N_part, dtype = float)
    w = np.ones(N_part, dtype = float)

    particleinterpolator.interpolate(b'test.dat', x, y, z, r, q, w, N)

if args.plot:
    # Plot result
    slice_num = 64

    data = np.loadtxt('test.dat')

    values = data[:, 0]

    indices = np.arange(0, len(values))

    x_indices = np.array(indices / N**2, dtype = int)
    y_indices = np.array((indices % N**2) / N, dtype = int)
    z_indices = np.array(indices - y_indices * N - x_indices * N**2, dtype = int)

    matrix = np.zeros((N, N, N), dtype = float)

    for i, map_value in enumerate(values):
        matrix[x_indices[i]][y_indices[i]][z_indices[i]] = map_value

    plt.figure(figsize = (8, 8), dpi = 90)
    plt.imshow(matrix[:, slice_num, :], cmap = 'afmhot')
    plt.savefig('test.png')
    plt.close()


