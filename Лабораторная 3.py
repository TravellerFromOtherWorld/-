def m_Jacobi(A, b, eps):
    cur_x1 = b[0] / A[0][0]
    cur_x2 = b[1] / A[1][1]
    cur_x3 = b[2] / A[2][2]

    next_x1 = (b[0] + (-1) * A[0][1] * cur_x2 + (-1) * A[0][2] * cur_x3) / A[0][0]
    next_x2 = (b[1] + (-1) * A[1][0] * cur_x1 + (-1) * A[1][2] * cur_x3) / A[1][1]
    next_x3 = (b[2] + (-1) * A[2][0] * cur_x1 + (-1) * A[2][1] * cur_x2) / A[2][2]
    n = 1

    while abs(next_x1 - cur_x1) > eps or abs(next_x2 - cur_x2) > eps or abs(next_x3 - cur_x3) > eps:
        n += 1
        cur_x1 = next_x1
        cur_x2 = next_x2
        cur_x3 = next_x3

        next_x1 = (b[0] + (-1) * A[0][1] * cur_x2 + (-1) * A[0][2] * cur_x3) / A[0][0]
        next_x2 = (b[1] + (-1) * A[1][0] * cur_x1 + (-1) * A[1][2] * cur_x3) / A[1][1]
        next_x3 = (b[2] + (-1) * A[2][0] * cur_x1 + (-1) * A[2][1] * cur_x2) / A[2][2]

    return next_x1, next_x2, next_x3, n


def m_Gaussa_Zeidelya(A, b, eps):
    cur_x1 = 0
    cur_x2 = b[1] / A[1][1]
    cur_x3 = b[2] / A[2][2]

    next_x1 = (b[0] + (-1) * A[0][1] * cur_x2 + (-1) * A[0][2] * cur_x3) / A[0][0]
    next_x2 = (b[1] + (-1) * A[1][0] * next_x1 + (-1) * A[1][2] * cur_x3) / A[1][1]
    next_x3 = (b[2] + (-1) * A[2][0] * next_x1 + (-1) * A[2][1] * next_x2) / A[2][2]
    n = 1

    while abs(next_x1 - cur_x1) > eps or abs(next_x2 - cur_x2) > eps or abs(next_x3 - cur_x3) > eps:
        n += 1
        cur_x1 = next_x1
        cur_x2 = next_x2
        cur_x3 = next_x3

        next_x1 = (b[0] + (-1) * A[0][1] * cur_x2 + (-1) * A[0][2] * cur_x3) / A[0][0]
        next_x2 = (b[1] + (-1) * A[1][0] * next_x1 + (-1) * A[1][2] * cur_x3) / A[1][1]
        next_x3 = (b[2] + (-1) * A[2][0] * next_x1 + (-1) * A[2][1] * next_x2) / A[2][2]

    return next_x1, next_x2, next_x3, n


if __name__ == '__main__':
    A = [
        [0.94, 0.12, -0.61],
        [-0.26, -0.6, -0.24],
        [-0.14, 1.21, 1.56]
    ]

    b = [1.63, -2.56, 4.91]

    print('Метод Якоби:\n', m_Jacobi(A, b, 0.5 * 10**(-4)))
    print('Метод Гаусса-Зейделя:\n', m_Gaussa_Zeidelya(A, b, 0.5 * 10**(-4)))
