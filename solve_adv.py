import math
import gauss

width = 0.5
square_center = [(-0.75, 0.75), (-0.25, 0.75),
          (-0.75, 0.25), (-0.25, 0.25),
          (-0.75, -0.25), (-0.25, -0.25), (0.25, -0.25), (0.75, -0.25),
          (-0.75, -0.75), (-0.25, -0.75), (0.25, -0.75), (0.75, -0.75)]

neumann_center = [(-0.25, 1), (-0.75, 1), (-1, 0.75), (-1, 0.25), (-1, -0.25), (-1, -0.75), (-0.75, -1), (-0.25, -1), (0.25, -1), (0.75, -1), (1, -0.75), (1, -0.25)]

max_e = [(-1, 1), (-0.5, 1),
         (-1, 0.5), (-0.5, 0.5),
         (-1, 0), (-0.5, 0),
         (-1, -0.5), (-0.5, -0.5), (0, -0.5), (0.5, -0.5), (1, -0.5),
         (-1, -1), (-0.5, -1), (0, -1), (0.5, -1), (1, -1),]

test_functions_for_corners = [(1, 0, 2, 3), (-1, 1, 3, -1),
                              (3, 2, 4, 5), (-1, 3, 5, -1),
                              (5, 4, 6, 7), (-1, 5, 7, 8), (-1, -1, 8, 9), (-1, -1, 9, 10),
                              (7, 6, 11, 12), (8, 7, 12, 13), (9, 8, 13, 14), (10, 9, 14, 15)]

def in_any_element(x,y):
    if -1 <= x <= 0 and 0 <= y <= 1: return True
    if -1 <= x <= 1 and -1 <= y <= 0: return True
    return False

def k(x,y):
    if not in_any_element(x, y):
        return 0
    if y > 0.5:
        return 2
    if y == 0.5:
        return 1.5
    return 1

def phi(i, x, y):
    if i == 1:
        return x*y
    if i == 2:
        return (1-x) * y
    if i == 3:
        return (x - 1) * (y - 1)
    if i == 4:
        return x * (1 - y)
    return 0

def phi_x_derivative(i, x, y):
    if phi(i, x, y) != 0:
        if i == 1 or i == 4: return 0.5/width
        if i == 2 or i == 3: return -0.5/width
    return 0

def phi_y_derivative(i, x, y):
    if phi(i, x, y) != 0:
        if i == 1 or i == 2: return 0.5/width
        if i == 3 or i == 4: return -0.5/width
    return 0

def e(i, x, y, func):
    x_max = max_e[i][0]
    y_max = max_e[i][1]
    if not in_any_element(x,y):
        return 0
    if abs(x_max - x) >= width or abs(y_max - y) >= width:
        return 0
    else:
        if x >= x_max and y >= y_max:
            return func(3, (x-x_max)/width, (y-y_max)/width)
        if x >= x_max and y <= y_max:
            return func(2, (x-x_max)/width, (y-(y_max - width))/width)
        if x <= x_max and y >= y_max:
            return func(4, (x-(x_max - width))/width, (y-y_max)/width)
        if x <= x_max and y <= y_max:
            return func(1, (x-(x_max - width))/width, (y-(y_max - width))/width)

def g(x, y):
    return math.pow(x*x, 1/3)

def B_matrix():
    B = [[0 for _ in range(16)] for _ in range(16)]
    for (sq, (x, y)) in enumerate(square_center):
        for i in test_functions_for_corners[sq]:
            for j in test_functions_for_corners[sq]:
                if i != -1 and j != -1:
                    B[i][j] += k(x, y) * e(i, x, y, phi_x_derivative) * e(j, x, y, phi_x_derivative) * width * width + k(x,y) * e(i, x, y, phi_y_derivative) * e(j, x, y, phi_y_derivative) * width * width
    return B

def L_matrix():
    L = [0 for _ in range(16)]
    for i in range(16):
        L[i] = sum([k(x,y) * e(i, x, y, phi) * g(x, y) * width for (x, y) in neumann_center])
    return L

def solve():
    B = B_matrix()
    L = L_matrix()
    return gauss.gauss_solve(B, L)

def u_h(w, x, y):
    return sum([w[i] * e(i, x, y, phi) for i in range(16)])

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