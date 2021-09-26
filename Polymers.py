import math
import random
import numpy as np
import random
import matplotlib.pyplot as plt
import os
import time

"""
https://stackoverflow.com/questions/3838329/how-can-i-check-if-two-segments-intersect
https://math.stackexchange.com/questions/270767/find-intersection-of-two-3d-lines/271366
https://dwightreid.com/blog/2015/09/21/python-how-to-solve-simultaneous-equations/
https://docs.python.org/3/tutorial/errors.html
https://numpy.org/doc/stable/reference/generated/numpy.zeros.html
https://www.coursera.org/lecture/intro-to-numerical-analysis/convergence-criteria-for-simple-iteration-67Mme
https://digitalcommons.unl.edu/cgi/viewcontent.cgi?article=1114&context=usgsstaffpub
"""

# Global run parameters
PolymerLengths = [20, 60, 80]
MonomerLength = 10
EdgeDiameter = 2
Dim = 2
AngleInterval = 20
MeanErr = 0.000001
MaxTries = 10000

# initializing

file_name_radii_3d = 'radii_3d_N{}_l{}.txt'
file_name_radii_2d = 'radii_2d_N{}_l{}.txt'
Rs_Vs_N = []
fig_num = [0]
max_tries = 8
theta_list = np.arange(0, np.pi, np.deg2rad(AngleInterval))
phi_list = np.arange(0, 2*np.pi, np.deg2rad(AngleInterval))


def calc_dist(x1, y1, z1, x2=0, y2=0, z2=0):
    return math.sqrt((x2 - x1) ** 2 + (y2 - y1) ** 2 + (z2 - z1) ** 2)


def check_ball_dist():
    """
    check that distance btn 2 molecules is smaller than d
    """
    for i in range(3, len(coor_x) + 1):
        if Dim == 2:
            dist = calc_dist(coor_x[-1], coor_y[-1], 0, coor_x[-i], coor_y[-i], 0)
        elif Dim == 3:
            dist = calc_dist(coor_x[-1], coor_y[-1], coor_z[-1], coor_x[-i], coor_y[-i], coor_z[-i])
        if dist <= EdgeDiameter:
            return False
    return True


def check_crossing():
    """
    check that molecule doesn't cross another
    """
    return True

    if len(coor_x) < 4:
        return True

    # what if matrix is singular
    if Dim == 3:
        X1, Y1, Z1, X2, Y2, Z2 = coor_x[-1], coor_y[-1], coor_z[-1], coor_x[-2], coor_y[-2], coor_z[-2]
        for i in range(3, len(coor_x)):
            X3, Y3, Z3, X4, Y4, Z4 = coor_x[-i], coor_y[-i], coor_z[-i], coor_x[-i - 1], coor_y[-i - 1], coor_z[-i - 1]

            # solving equations to check if lines are crossing
            M1 = np.zeros((2, 2), dtype=np.int)
            M1[0, 0] = X1 - X2
            M1[1, 0] = Y1 - Y2
            M1[0, 1] = X4 - X3
            M1[1, 1] = Y4 - Y3
            M2 = np.zeros(2, dtype=np.int)
            M2[0] = X4 - X2
            M2[1] = Y4 - Y2

            try:  # check that matrix is not singular (parallel) in the x-y plain
                M3 = np.linalg.inv(M1)
            except np.linalg.LinAlgError:
                print('Parallel lines x-y')
                M1[1, 0] = Z1 - Z2
                M1[1, 1] = Z4 - Z3
                M2[1] = Z4 - Z2
                try:  # check that matrix is not singular (parallel) in the x-z plain
                    M3 = np.linalg.inv(M1)
                except np.linalg.LinAlgError:
                    print('Parallel lines x-z')
                    continue  # if parallel in x-y and x-z - no crossing
                else:
                    M4 = np.dot(M3, M2)
                    if np.all(M4 >= 0) and np.all(M4 <= 1):  # check s and t <1 and >0
                        if np.abs(((Y1 - Y2) * M4[0] + (Y4 - Y3) * M4[1]) - (Y4 - Y2)) < StickDiameter:
                            # check lines distance is less then StickDiameter in the remaining axis (Y)
                            return False

            M4 = np.dot(M3, M2)
            if np.all(M4 >= 0) and np.all(M4 <= 1):  # check s and t <1 and >0
                if np.abs(((Z1 - Z2) * M4[0] + (Z4 - Z3) * M4[1]) - (Z4 - Z2)) < StickDiameter:
                    # check lines distance is less then StickDiameter in the remaining axis (Z)
                    return False
            else:
                continue
        return True

    elif Dim == 2:
        X1, Y1, X2, Y2 = coor_x[-1], coor_y[-1], coor_x[-2], coor_y[-2]
        for i in range(3, len(coor_x)):
            X3, Y3, X4, Y4 = coor_x[-i], coor_y[-i], coor_x[-i - 1], coor_y[-i - 1]
            # solving equations to check if lines are crossing
            M1 = np.zeros((2, 2))
            M1[0, 0] = X1 - X2
            M1[1, 0] = Y1 - Y2
            M1[1, 0] = Y1 - Y2
            M1[0, 1] = X4 - X3
            M1[1, 1] = Y4 - Y3
            M2 = np.zeros(2)
            M2[0] = X4 - X2
            M2[1] = Y4 - Y2
            try:  # check that matrix is not singular (parallel) in the x-y
                M3 = np.linalg.inv(M1)
            except np.linalg.LinAlgError:
                print('Parallel lines x-y')
                continue  # if parallel in x-y - no crossing
            M4 = np.dot(M3, M2)
            if np.all(M4 >= 0) and np.all(M4 <= 1):  # check s and t <1 and >0
                return False
            else:
                continue
        return True


def polymer_sim(PolymerLength, dim, file_name='coordinates.txt'):

    global coor_x, coor_y, coor_z, fp, R
    file_object = open(file_name, 'w').close()
    coor_x = [0]
    coor_y = [0]
    coor_z = [0]

    if dim == 3:
        with open(file_name, mode='w') as fp:  # write first coordinate to file
            fp.write("X COOR\tY COOR\tZ COOR\n")
        with open(file_name, mode='a') as fp:  # write first coordinate to file
            fp.write("%+05.5f\t%+05.5f\t%+05.5f\n" % (coor_x[-1], coor_y[-1], coor_z[-1]))
        for step in range(PolymerLength):  # start progressing the chain
            try_num = 0  # if chain gets stuck then abort
            while True and try_num < max_tries:
                # theta = random.random() * math.pi  # picking a direction with angle theta between [0,pi]
                # phi = random.random() * 2 * math.pi  # and phi between [0,2pi]
                theta = np.random.choice(theta_list)
                phi = np.random.choice(phi_list)
                coor_x.append(coor_x[-1] + MonomerLength * math.sin(theta) * math.cos(phi))  # calc coordinates
                coor_y.append(coor_y[-1] + MonomerLength * math.sin(theta) * math.sin(phi))
                coor_z.append(coor_z[-1] + MonomerLength * math.cos(theta))
                if check_ball_dist() and check_crossing():  # check that distance and crossing is OK
                    with open(file_name, mode='a') as fp:
                        fp.write("%+05.5f\t%+05.5f\t%+05.5f\n" % (coor_x[-1], coor_y[-1], coor_z[-1]))
                    break
                else:  # if not OK, calc again
                    coor_x.pop()
                    coor_y.pop()
                    coor_z.pop()
                    try_num += 1

    elif dim == 2:
        with open(file_name, mode='w') as fp:  # write first coordinate to file
            fp.write("X COOR\tY COOR\tZ COOR\n")
        with open(file_name, mode='a') as fp:  # write first coordinate to file
            fp.write('{:+05.5f}\t{:+05.5f}\n'.format(coor_x[-1], coor_y[-1]))
        for step in range(PolymerLength):
            try_num = 0
            while True and try_num < max_tries:
                # phi = random.random() * 2 * math.pi  # and phi between [0,2pi]
                phi = np.random.choice(phi_list)
                coor_x.append(coor_x[-1] + MonomerLength * math.cos(phi))  # calc coordinates
                coor_y.append(coor_y[-1] + MonomerLength * math.sin(phi))
                if check_ball_dist() and check_crossing():  # check that distance and crossing is OK
                    with open(file_name, mode='a') as fp:  # write coordinate to file
                        fp.write('{:+05.5f}\t{:+05.5f}\n'.format(coor_x[-1], coor_y[-1]))
                    break
                else:  # if not OK, calc again
                    coor_x.pop()
                    coor_y.pop()
                    try_num += 1
    R = calc_dist(coor_x[-1], coor_y[-1], coor_z[-1])  # calc R
    return R


for PolymerLength in PolymerLengths:  # for all lengths of chain

    start=time.time()
    if Dim == 2:
        radii_file_name = file_name_radii_2d.format(PolymerLength, MonomerLength)
    elif Dim == 3:
        radii_file_name = file_name_radii_3d.format(PolymerLength, MonomerLength)

    filesize = os.path.exists(radii_file_name)

    if filesize == 0:
        with open(radii_file_name, mode='w') as fp:  # write R to file
            fp.write("Radii\n")

    Rs=[]
    R_avg_old = None
    find_mean = False
    while not find_mean:

        # initialize chain
        R=polymer_sim(PolymerLength, Dim)
        Rs.append(R)
        R_avg_new=np.mean(Rs)

        with open(radii_file_name, mode='a') as fp:  # write R to file
            fp.write("%.5f\n" % (R))

        with open(radii_file_name, mode='r') as fp:  # read file of R
            str = fp.read().split()[1:-1]
        # calculate mean

        # print(R_avg_new)
        if R_avg_old is None:
            R_avg_old = R_avg_new
            continue
        if np.abs(R_avg_new - R_avg_old)/R_avg_new < MeanErr:  # check if average R converged
            find_mean = True
            Rs_Vs_N.append(R_avg_new)
        else:
            R_avg_old = R_avg_new

    # plot histogram of Rs
    plt.figure(fig_num[0])
    fig_num[0] += 1
    plt.hist(Rs)
    plt.title('histogram for : dim-{} N-{} l-{}'.format(Dim, PolymerLength, MonomerLength))
    plt.xlabel('R')
    plt.ylabel('')
    plt.show()

    print('LENNNNNNNNN', len(Rs))
    print(PolymerLength, np.mean(Rs))
    print(time.time()-start)

# plot R Vs. Chain length
plt.figure(fig_num[0])
fig_num[0] += 1
plt.plot(PolymerLengths, Rs_Vs_N)
plt.title('Radii Vs. N : dim-{} l-{}'.format(Dim, PolymerLength, MonomerLength))
plt.xlabel('N')
plt.ylabel('R')
plt.show()
