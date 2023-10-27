
import roboticstoolbox as rtb
from spatialmath.base import *
from main import *
import numpy as np
# import swift
# import spatialmath as sm

if __name__ == '__main__':
    # robot = rtb.models.DH.Panda()
    # print(robot)


    deg_xyz = (10, 20, 30)
    deg_xyz = [rad(n) for n in deg_xyz]
    R = mat_rotate_z(deg_xyz[2]) @ mat_rotate_y(deg_xyz[1]) @ mat_rotate_x(deg_xyz[0])
    print("my rotation matrix of", "x=", (int)(degree(deg_xyz[0])), "y=", (int)(degree(deg_xyz[1])), "z=", (int)(degree(deg_xyz[2])))
    print(R)
    print("robotics tool rotation matrix of", "x=", (int)(degree(deg_xyz[0])), "y=", (int)(degree(deg_xyz[1])), "z=", (int)(degree(deg_xyz[2])))
    T = rpy2tr(deg_xyz[0], deg_xyz[1], deg_xyz[2], 'rad',order="zyx")
    print(T)

    print("my euler transform:")
    print(degree(angles_euler_zyx(R)))
    print("robotics tool euler transform:")
    print(tr2rpy(T, 'deg',"zyx"))

    deg_xyz = (0, 20, 0)
    B_ORG = np.array([3, 0, 1])

    print("my homogeneous transformation:")
    deg_xyz = [rad(n) for n in deg_xyz]
    R = mat_rotate_z(deg_xyz[2]) @ mat_rotate_y(deg_xyz[1]) @ mat_rotate_x(deg_xyz[0])
    T = mat_trans(R, B_ORG.reshape(3,1))

    vector_B = np.array([1, 0, 1])
    vector_A = R @ vector_B + B_ORG
    print(vector_A)



    print("robotics tool homogeneous transformation:")
    T = rpy2tr(deg_xyz[0], deg_xyz[1], deg_xyz[2], 'rad',order="zyx")

    vector_B = np.array([1, 0, 1, 1])
    vector_A = T @ vector_B + np.append(transl(transl(3, 0, 1)), [0])
    vector_A = np.delete(vector_A, 3)
    print(vector_A)




