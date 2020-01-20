import math
import gauss

def k(x,y):
    return 1

def phi(i, x, y):
    if i == 1:
        return (1 + x)*(1 + y)
    if i == 2:
        return -x * (y + 1)
    if i == 3:
        return x * y
    if i == 4:
        return -(x + 1) * y
    return 0

def phi_x_derivative(i, x, y):
    if phi(i, x, y) != 0:
        if i == 1 or i == 4: return 0.5
        if i == 2 or i == 3: return -0.5
    return 0

def phi_y_derivative(i, x, y):
    if phi(i, x, y) != 0:
        if i == 1 or i == 2: return 0.5
        if i == 3 or i == 4: return -0.5
    return 0

def e(i, x, y, func):
    if -1 <= x <= 0 and 0 <= y <= 1:
        if i == 0: return func(2, x, y - 1)
        if i == 1: return func(3, x, y - 1)
        return 0
    if -1 <= x <= 0 and -1 <= y <= 0:
        if i == 1: return func(2, x, y)
        if i == 2: return func(3, x, y)
        if i == 3: return func(4, x, y)
        return 0
    if 0 <= x <= 1 and -1 <= y <= 0:
        if i == 3: return func(3, x - 1, y)
        if i == 4: return func(4, x - 1, y)
        return 0
    return 0

def g(x, y):
    return math.pow(x*x, 1/3)

def B_matrix():
    B = [[0 for _ in range(5)] for _ in range(5)]
    centers = [(-0.5, 0.5), (-0.5, -0.5), (0.5, -0.5)]
    for i in range(5):
        for j in range(5):
            B[i][j] = sum([k(x, y) * e(i, x, y, phi_x_derivative) * e(j, x, y, phi_x_derivative) + k(x, y) * e(i, x, y, phi_y_derivative) * e(j, x, y, phi_y_derivative) for (x,y) in centers])
    return B

def L_matrix():
    L = [0 for _ in range(5)]
    centers = [(-0.5, 1), (-1, 0.5), (-1, -0.5), (-0.5, -1), (0.5, -1), (1, -0.5)]
    for i in range(5):
        L[i] = sum([k(x,y) * e(i, x, y, phi) * g(x, y) for (x, y) in centers])
    return L

def solve():
    B = B_matrix()
    L = L_matrix()
    return gauss.gauss_solve(B, L)

def u_h(w, x, y):
    return sum([w[i] * e(i, x, y, phi) for i in range(5)])

def ex_u_values(probes):
    matrix = [[(None, None) for _ in range(2*probes)] for _ in range(2*probes)]
    probe_length = 1/probes
    w = solve()
    for (x_start, y_start, i_start, j_start) in [(-1, 1, 0, 0), (-1, 0, probes, 0), (0, 0, probes, probes)]:
        y = y_start - probe_length/2
        for i in range(i_start, i_start + probes):
            x = x_start + probe_length/2
            for j in range(j_start, j_start + probes):
                matrix[i][j] = (x, y)
                x += probe_length
            y -= probe_length

    result = [[0 for _ in range(2*probes)] for _ in range(2*probes)]
    for (i, row) in enumerate(matrix):
        for (j, (x, y)) in enumerate(row):
            if x is not None:
                result[i][j] = u_h(w, x, y)
    return result
