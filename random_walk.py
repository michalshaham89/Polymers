import random
import timeit
import numpy as np
import operator
from matplotlib import pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
import pickle
import os
import re
import time

WALKS_DIR_PATH = "C:\Users\Michal\Documents\michuli\physics_2\cm_project\walks"
EE_RADII_DIR_PATH = "C:\Users\Michal\Documents\michuli\physics_2\cm_project\end_to_end_radii"
save = True
ERR = 0.01  # 1 percent


def plot_walk(walk, d):
    if type(walk) is dict:
        walk = zip(*sorted(walk.items(), key=operator.itemgetter(1)))[0]
    if d == 2:
        plt.plot(*zip(*walk))
    elif d == 3:
        fig = plt.figure()
        ax = Axes3D(fig)
        ax.plot(*zip(*walk))


def concatenating(w1, w2):
    last_vertex, n1 = max(w1.iteritems(), key=operator.itemgetter(1))
    for vertex, i in w2.items():
        if i != 0:
            new_vertex = tuple(np.add(vertex, last_vertex))
            if new_vertex in w1:
                return False
            else:
                w1[new_vertex] = n1 + i
    return True


def std_err(array):
    return array.std(ddof=1) / np.sqrt(len(array))


def load_data(n, d):
    print('LOADING DATA', n, d)
    dir_path = os.path.join(WALKS_DIR_PATH, 'saw_%dd_%dlen' % (d, n))
    ee_radii = []
    if os.path.exists(os.path.join(dir_path)):
        for walk_file in os.listdir(dir_path):
            try:
                walk = pickle.load(open(os.path.join(dir_path, walk_file), 'rb'))
            except EOFError as e:
                print(e.message)
                print(os.path.join(dir_path, walk_file))
                continue
            ee_radii.append(np.sqrt(sum([np.square(coordinate) for coordinate in walk[-1]])))
    return np.array(ee_radii)


def get_saw_simple(n, d, to_plot=False):
    axis_options = range(d)
    dir_options = (-1, 1)

    while True:
        walk = {}
        tmp_vertex = [0] * d
        walk[tuple(tmp_vertex)] = 0

        chosen_axis = 0
        chosen_dir = 0
        start_again = False

        for i in range(1, n + 1):
            while True:
                tmp_axis, tmp_dir = random.choice(axis_options), random.choice(dir_options)
                if tmp_axis == chosen_axis and tmp_dir == -chosen_dir:
                    continue
                else:
                    chosen_axis, chosen_dir = tmp_axis, tmp_dir
                    break
            tmp_vertex[chosen_axis] += chosen_dir
            if tuple(tmp_vertex) in walk:
                start_again = True
                break
            else:
                walk[tuple(tmp_vertex)] = i
        if start_again:
            continue
        else:
            break

    if to_plot:
        walk = zip(*sorted(walk.items(), key=operator.itemgetter(1)))[0]
        plot_walk(walk, d)

    return walk


def get_saw_dimerization(n, n_0, d):
    if n <= n_0:
        return get_saw_simple(n, d)
    else:
        while True:
            w1 = get_saw_dimerization(n / 2, n_0, d)
            w2 = get_saw_dimerization(n / 2, n_0, d)
            if concatenating(w1, w2):
                return w1


def get_end_to_end_radius(n, n_0, d):
    walk = get_saw_dimerization(n, n_0, d)
    walk = zip(*sorted(walk.items(), key=operator.itemgetter(1)))[0]
    if save:
        dir_path = os.path.join(WALKS_DIR_PATH, 'saw_%dd_%dlen' % (d, n))
        if not os.path.exists(os.path.join(dir_path)):
            os.makedirs(dir_path)
        nums = map(int, re.findall(r'\d+', ''.join(os.listdir(dir_path))))
        num = max(nums) + 1 if nums else 1
        pickle.dump(walk, open(os.path.join(dir_path, 'saw_' + str(num)), 'wb'))

    return np.sqrt(sum([np.square(coordinate) for coordinate in walk[-1]]))


def get_saw_samplings(n, n_0, d):
    tmp_err = 1
    ee_radii = load_data(n, d)
    ee_radii = ee_radii if ee_radii.size else np.array([get_end_to_end_radius(n, n_0, d)])
    start = time.time()
    i = 0
    print('STARTING', n, n_0, d, np.mean(ee_radii))
    while tmp_err > ERR:
        ee_radii = np.concatenate((ee_radii, np.array([get_end_to_end_radius(n, n_0, d)])))
        tmp_err = std_err(ee_radii) / ee_radii.mean()
        # plt.plot(e_to_e_radii)
        # plt.show()
        # print len(e_to_e_radii), e_to_e_radii.std(ddof=1), tmp_err, e_to_e_radii.mean()
        i += 1
        if i % 10 == 0:
            print(ee_radii.size, int(time.time() - start), n, n_0, d, ee_radii.std(ddof=1), \
                std_err(ee_radii) / np.mean(ee_radii), std_err(ee_radii), np.mean(ee_radii))

    print('RESULTS', ee_radii.size, int(time.time() - start), n, n_0, d, ee_radii.std(ddof=1), \
        std_err(ee_radii) / np.mean(ee_radii), std_err(ee_radii), np.mean(ee_radii))
    print()

    if save:
        dir_path = os.path.join(EE_RADII_DIR_PATH, 'ee_radii_%dd_%dlen' % (d, n))
        if not os.path.exists(os.path.join(dir_path)):
            os.makedirs(dir_path)
        nums = map(int, re.findall(r'\d+', ''.join(os.listdir(dir_path))))
        num = max(nums) + 1 if nums else 1
        pickle.dump(ee_radii, open(os.path.join(dir_path, 'ee_radii_' + str(num)), 'wb'))


if __name__ == '__main__':
    N_0 = 32
    D = 2
    for deg in range(11, 15):
        N = pow(2, deg)
        get_saw_samplings(N, N_0, D)
