def print_matrix(A):
    for row in A:
        print(row)

def print_equation_system(system):
    (A, B) = system
    for i, row in enumerate(A):
        print(str(row) + " " + str([B[i]]))

def gaussian_matrix(A, B):
    l = len(A)
    for current_index in range(l-1):
        max_row = max(range(current_index, l-1), key=lambda i: abs(A[i][current_index]))
        (A[current_index], A[max_row]) = (A[max_row], A[current_index])
        (B[current_index], B[max_row]) = (B[max_row], B[current_index])
        for row in range(current_index + 1, l):
            multiplier = A[row][current_index] / A[current_index][current_index]
            A[row] = [0 for _ in range(current_index + 1)] + [A[row][col] - multiplier * A[current_index][col] for col in range(current_index + 1, l)]
            B[row] -= multiplier * B[current_index]
    return (A, B)

def gauss_solve(A, B):
    (A, B) = gaussian_matrix(A, B)
    l = len(A)
    for current_index in range(l-1, 0, -1):
        for row_index in range(current_index - 1, -1, -1):
            B[row_index] -= A[row_index][current_index] / A[current_index][current_index] * B[current_index]
            A[row_index][current_index] = 0
    return [B[i] / A[i][i] for i in range(l)]
