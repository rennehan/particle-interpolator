import array
import particleinterpolator


x = array.array('d', [0.5])
y = array.array('d', [0.5])
z = array.array('d', [0.5])

r = array.array('d', [0.1])
q = array.array('d', [1])
w = array.array('d', [1])
N = 100

particleinterpolator.interpolate(b'test.dat', x, y, z, r, q, w, N)
