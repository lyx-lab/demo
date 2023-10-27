import numpy as np
import sympy as sy

np.set_printoptions(suppress=True)

# 角度n转化为弧度
def rad(n):
    return ((n) * np.pi / 180)

# 弧度n转化为角度
def degree(n):
    return ((n) * 180 / np.pi)

# X轴的旋转，默认输入单位为弧度
# numpy版本，用于数值计算
def mat_rotate_x(theta):
    s_theta = np.sin(theta)
    c_theta = np.cos(theta)
    R = np.identity(3)
    R[1,1] = c_theta
    R[2,1] = s_theta
    R[1,2] = -s_theta
    R[2,2] = c_theta
    return R

# sympy版本，用于符号计算
def sy_mat_rotate_x(theta):
    s_theta = sy.sin(theta)
    c_theta = sy.cos(theta)
    R = sy.eye(3)
    R[1,1] = c_theta
    R[2,1] = s_theta
    R[1,2] = -s_theta
    R[2,2] = c_theta
    return R

# Y轴的旋转，默认输入单位为弧度
# numpy版本，用于数值计算
def mat_rotate_y(theta):
    s_theta = np.sin(theta)
    c_theta = np.cos(theta)
    R = np.identity(3)
    R[0,0] = c_theta
    R[2,0] = -s_theta
    R[0,2] = s_theta
    R[2,2] = c_theta
    return R

# sympy版本，用于符号计算
def sy_mat_rotate_y(theta):
    s_theta = sy.sin(theta)
    c_theta = sy.cos(theta)
    R = sy.eye(3)
    R[0,0] = c_theta
    R[2,0] = -s_theta
    R[0,2] = s_theta
    R[2,2] = c_theta
    return R

# Z轴的旋转，默认输入单位为弧度
# numpy版本，用于数值计算
def mat_rotate_z(theta):
    s_theta = np.sin(theta)
    c_theta = np.cos(theta)
    R = np.identity(3)
    R[0,0] = c_theta
    R[1,0] = s_theta
    R[0,1] = -s_theta
    R[1,1] = c_theta
    return R

# sympy版本，用于符号计算
def sy_mat_rotate_z(theta):
    s_theta = sy.sin(theta)
    c_theta = sy.cos(theta)
    R = sy.eye(3)
    R[0,0] = c_theta
    R[1,0] = s_theta
    R[0,1] = -s_theta
    R[1,1] = c_theta
    return R

# 生成一个n维向量
def vector(*ele):
    return np.array(ele, dtype="float64").reshape(-1, 1)
# 生成一个3维矩阵
def mat3(*ele):
    assert len(ele) == 9
    return np.array(ele, dtype="float64").reshape(3, 3)

# 变换算子（平移及旋转）
def mat_trans(R, P):
    T = np.identity(4)
    T[:-1, :-1] = R
    T[:-1, -1:] = P
    return T
# 逆变换
def mat_trans_inverse(mat):
    R = mat[:-1, :-1]
    P = mat[:-1, -1:]
    R_inv = R.T
    P_inv = - R_inv @ P
    return mat_trans(R_inv, P_inv)


# Z-Y-Z欧拉角
# def angles_euler_zyz(mat):
#     if (mat[0, 2]**2 + mat[1,2]**2) < 1e-5:
#         if mat[2, 2] > 0:
#             beta = 0
#         else:
#             beta = np.pi
#         alpha = 0
#         gamma = np.arctan(-mat[0, 1] / mat[0, 0])
#     else:
#         sin_beta = np.sqrt(mat[2, 0]**2 + mat[2, 1]**2)
#         beta = np.arctan(sin_beta / mat[2, 2])
#         alpha = np.arctan(mat[1, 2] / mat[0, 2])
#         gamma = np.arctan(mat[2, 1] / -mat[2, 0])
#
#     return np.array((alpha, beta, gamma))
def angles_euler_zyz(mat, interval=None):
    if interval is None:
        interval = ((-np.pi,np.pi),(-np.pi,np.pi),(-np.pi,np.pi))

    def rotate(alpha, beta, gamma):
        return mat_rotate_z(alpha) @ mat_rotate_y(beta) @ mat_rotate_z(gamma)

    def arctan_ls(angle):
        if angle > 0:
            return [angle, angle - np.pi]
        else:
            return [angle, angle + np.pi]

    if (mat[0, 2]**2 + mat[1,2]**2) < 1e-5:
        if mat[2, 2] > 0:
            beta = 0
        else:
            beta = np.pi
        alpha = 0
        gamma = np.arctan(-mat[0, 1] / mat[0, 0])
    else:
        sin_beta = np.sqrt(mat[2, 0]**2 + mat[2, 1]**2)
        beta = np.arctan(sin_beta / mat[2, 2])
        alpha = np.arctan(mat[1, 2] / mat[0, 2])
        gamma = np.arctan(mat[2, 1] / -mat[2, 0])

    beta_ls = arctan_ls(beta) + arctan_ls(-beta)
    beta_ls = [n for n in beta_ls if (interval[0][0] < n < interval[0][1])]

    alpha_ls = arctan_ls(alpha)
    alpha_ls = [n for n in alpha_ls if (interval[1][0] < n < interval[1][1])]

    gamma_ls = arctan_ls(gamma)
    gamma_ls = [n for n in gamma_ls if (interval[2][0] < n < interval[2][1])]

    cnt = 0
    for a in alpha_ls:
        for b in beta_ls:
            for g in gamma_ls:
                mat1 = rotate(a, b, g)
                if (abs(mat - mat1) < 1e-5).all():
                    if(cnt == 0):
                        a1 = a
                        b1 = b
                        g1 = g
                        cnt+=1
                    else:
                        a2 = a
                        b2 = b
                        g2 = g
                        return np.array([(a1, b1, g1), (a2, b2, g2)])

    return np.array([(np.inf, np.inf, np.inf), (np.inf, np.inf, np.inf)])

# Z-Y-X欧拉角
def angles_euler_zyx(mat, interval=None):
    if interval is None:
        interval = ((-np.pi,np.pi),(-np.pi,np.pi),(-np.pi,np.pi))

    def rotate(alpha, beta, gamma):
        return mat_rotate_z(alpha) @ mat_rotate_y(beta) @ mat_rotate_x(gamma)

    def arctan_ls(angle):
        if angle > 0:
            return [angle, angle - np.pi]
        else:
            return [angle, angle + np.pi]

    if (mat[2, 1]**2 + mat[2,2]**2) < 1e-5:
        if mat[2, 0] > 0:
            beta  = -np.pi/2
            gamma = -np.arctan(mat[0, 1] / mat[1, 1])
        else:
            beta  = np.pi/2
            gamma = np.arctan(mat[0, 1] / mat[1, 1])
        alpha = 0
    else:
        cos_beta = np.sqrt(mat[0, 0]**2 + mat[1, 0]**2)
        beta = np.arctan(mat[2, 0] / cos_beta)
        alpha = np.arctan(mat[1, 0] / mat[0, 0])
        gamma = np.arctan(mat[2, 1] / mat[2, 2])

    beta_ls = arctan_ls(beta) + arctan_ls(-beta)
    beta_ls = [n for n in beta_ls if (interval[0][0] < n < interval[0][1])]

    alpha_ls = arctan_ls(alpha)
    alpha_ls = [n for n in alpha_ls if (interval[1][0] < n < interval[1][1])]

    gamma_ls = arctan_ls(gamma)
    gamma_ls = [n for n in gamma_ls if (interval[2][0] < n < interval[2][1])]

    cnt = 0
    for a in alpha_ls:
        for b in beta_ls:
            for g in gamma_ls:
                mat1 = rotate(a, b, g)
                if (abs(mat - mat1) < 1e-5).all():
                    if(cnt == 0):
                        a1 = a
                        b1 = b
                        g1 = g
                        cnt+=1
                    else:
                        a2 = a
                        b2 = b
                        g2 = g
                        if ((a1 == a2) & (b1 == b2) & (g1 == g2)):
                            continue
                        else:
                            cnt+=1
                            break

    if(cnt == 1):
        return np.array([(g1, b1, a1), (np.inf, np.inf, np.inf)])
    elif(cnt == 2):
        return np.array([(g1, b1, a1), (g2, b2, a2)])
    else:
        return np.array([(np.inf, np.inf, np.inf), (np.inf, np.inf, np.inf)])

# Press the green button in the gutter to run the script.
# if __name__ == '__main__':
    # matrix_x = mat_rotate_x(rad(30))
    # print(matrix_x)
    # x = sy.symbols("x")
    # matrix_sy_x = sy_mat_rotate_x(x)
    # print(matrix_sy_x)
    #
    # matrix_y = mat_rotate_y(rad(30))
    # print(matrix_y)
    # y = sy.symbols("y")
    # matrix_sy_y = sy_mat_rotate_y(y)
    # print(matrix_sy_y)
    #
    # matrix_z = mat_rotate_z(rad(30))
    # print(matrix_z)
    # z = sy.symbols("z")
    # matrix_sy_z = sy_mat_rotate_z(z)
    # print(matrix_sy_z)


    # matrix_rot_z = mat_rotate_z(rad(30)) @ vector(-1,1.732,0)
    # print(matrix_rot_z)

    # zyz = (-100, 10, 120)
    # zyz = [rad(n) for n in zyz]
    # R = mat_rotate_z(zyz[0]) @ mat_rotate_y(zyz[1]) @ mat_rotate_z(zyz[2])
    # print(R)
    # print(degree(angles_euler_zyz(R)))

    # zyx = (10, 20, 30)
    # zyx = [rad(n) for n in zyx]
    # R = mat_rotate_z(zyx[0]) @ mat_rotate_y(zyx[1]) @ mat_rotate_x(zyx[2])
    # print("The rotation matrix of","x=",(int)(degree(zyx[0])), "y=",(int)(degree(zyx[1])), "z=",(int)(degree(zyx[2])))
    # print(R)
    # print(degree(angles_euler_zyx(R)))


    zyx = (30, 90, -55)
    zyx = [rad(n) for n in zyx]
    R = mat_rotate_z(zyx[0]) @ mat_rotate_y(zyx[1]) @ mat_rotate_x(zyx[2])
    print("The rotation matrix of","x=",(int)(degree(zyx[0])), "y=",(int)(degree(zyx[1])), "z=",(int)(degree(zyx[2])))
    print(R)
    print(degree(angles_euler_zyx(R)))

    vector_a = mat_rotate_y(rad(20)) @ vector(1, 0, 1)
    print(vector_a)
    vector_a = mat_rotate_y(rad(-20)) @ vector_a
    print(vector_a)
