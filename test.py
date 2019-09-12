import particleinterpolator
import numpy as np
import matplotlib.pyplot as plt

# N_cell in the grid
N = 128

# Number of fake particles in 1D
N_part_1d = 64

# Smoothing length of particles
h = 5.0 / N_part_1d

N_part = N_part_1d**3

# Coordinates
x = np.array(np.random.uniform(size = N_part), dtype = float)
y = np.array(np.random.uniform(size = N_part), dtype = float)
z = np.array(np.random.uniform(size = N_part), dtype = float)

# Radii, quantity to deposit, weights
r = h * np.ones(N_part, dtype = float)
q = np.ones(N_part, dtype = float)
w = np.ones(N_part, dtype = float)

x = x.flatten()
y = y.flatten()
z = z.flatten()

particleinterpolator.interpolate(b'test.dat', x, y, z, r, q, w, N)


# Plot result
slice_num = 0

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
plt.imshow(matrix[:, :, slice_num], cmap = 'afmhot')
plt.savefig('test.png')
plt.close()

