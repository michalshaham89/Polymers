import os
import shutil
from matplotlib import pyplot as plt
import glob
import numpy as np

mean_err_filee=0.000001
err=0.01
Rs_list={2:{20:44, 60:79, 80:90}, 3:{100:1081,140:1502,180:1930}}
proj_path_name= '/Users/michal/PycharmProjects/Polymers'
proj_path_finish='/Users/michal/PycharmProjects/Polymers/tested_codes'
input_file_name= 'input.txt'
grading_file = 'grading_file.txt'
file_name_radii=['radii_{}d_N{}_l{}.txt', 'radii_{}d_N{}_I{}.txt', 'rdii_{}d_N{}.txt', 'radii_{}d_N{}',
                 'radii_{}d_N{}.txt', 'radii_{}d_N {}.txt', 'radii_{}d_N{}_l{}.0.txt']
wrong_file_name = False
dims=[2,3]
# dims=[3]
id_list=['209020114']


def test_coordiante():
    # f = open("coordinates.txt", "r")
    coor_files = ["coordinates.txt", "coordinate.txt", "coordinates", "cordi.txt", "3dpoints.txt", "2dpoints.txt"]
    for coor_file in coor_files:
        try:
            f = open(coor_file, "r")
            break
        except FileNotFoundError as err:
            print(err)

    X=[]
    Y=[]
    Z=[]
    ind=0
    for line in f:
        if ind==0:
            ind+=1
            continue
        line=line.split()
        X.append(float(line[0]))
        Y.append(float(line[1]))
        if len(line)==3:
            try:
                Z.append(float(line[2]))
            except ValueError as err:
                print(err)

    f.close()
    os.remove(os.path.join(proj_path_name, coor_file))
    print(len(X))
    if Z:
        fig = plt.figure()
        ax = fig.add_subplot(111, projection='3d')
        X=np.array(X)
        Y = np.array(Y)
        Z = np.array(Z)
        ax.plot3D(X, Y, Z)
        ax.scatter(X, Y, Z)
        max_range = np.array([X.max() - X.min(), Y.max() - Y.min(), Z.max() - Z.min()]).max()
        Xb = 0.5 * max_range * np.mgrid[-1:2:2, -1:2:2, -1:2:2][0].flatten() + 0.5 * (X.max() + X.min())
        Yb = 0.5 * max_range * np.mgrid[-1:2:2, -1:2:2, -1:2:2][1].flatten() + 0.5 * (Y.max() + Y.min())
        Zb = 0.5 * max_range * np.mgrid[-1:2:2, -1:2:2, -1:2:2][2].flatten() + 0.5 * (Z.max() + Z.min())
        # Comment or uncomment following both lines to test the fake bounding box:
        for xb, yb, zb in zip(Xb, Yb, Zb):
            ax.plot([xb], [yb], [zb], 'w')
    else:
        plt.figure()
        plt.plot(X, Y,'.')
        plt.plot(X, Y,'-')
        plt.gca().set_aspect('equal', adjustable='box')
    plt.grid()
    plt.show()


def test_radii(dim, grading_file):
    i=0
    files_rad=[]
    for path, dirs, files in os.walk(os.getcwd()):
        if path == '/Users/michal/PycharmProjects/Polymers':
            if dim==2:
                files_rad=glob.glob(file_name_radii[0].format(dim, '[0-9][0-9]', '[0-9][0-9]'))
            elif dim==3:
                files_rad=glob.glob(file_name_radii[0].format(dim, '[0-9][0-9][0-9]', '[0-9][0-9]'))
            if not files_rad:
                wrong_file_name = True
                for file_name_rad in files:
                    if 'radii' in file_name_rad.lower():
                        files_rad.append(file_name_rad)
    print(files_rad)

    i = 1
    for file in files_rad:
        n=int(file.split('N')[1].split('_')[0].split('.')[0].split(' ')[0])
        with open(file, "r") as f:
            text = f.read()
        text = text.split('\n')[1:]
        Rs = []
        for string in text:
            try:
                Rs.append(float(string))
            except ValueError as e:
                print(string, e)
        Rs_mean1=np.mean(Rs)
        Rs_mean2 = np.mean(Rs[:-1])
        if np.abs(Rs_mean1-Rs_mean2)/np.max([Rs_mean2,Rs_mean1])<=mean_err_filee:
            with open(grading_file, 'a') as f:
                f.write('\t{}:YES'.format(i))
        else:
            with open(grading_file, 'a') as f:
                print('Wrong condition: ', file, Rs_mean1, Rs_mean2)
                f.write('\t{}:NO {:.2f}/{:.2f}/{:.2e}'.format(i, round(Rs_mean1,2), round(Rs_mean2,2),
                        np.abs(Rs_mean1-Rs_mean2)/np.max([Rs_mean2,Rs_mean1])))
        if np.abs(Rs_mean2-Rs_list[dim][n])/np.min([Rs_mean2,Rs_list[dim][n]])<=err:
            with open(grading_file, 'a') as f:
                f.write(' {}:YES'.format(i))
        else:
            print('Wrong radius: ', file, Rs_mean2, Rs_list[dim][n])
            with open(grading_file, 'a') as f:
                f.write(' {}:NO {:.2f}/{:.2f}'.format(i, round(Rs_mean1, 2), Rs_list[dim][n]))
        i+=1

        os.remove(os.path.join(proj_path_name, file))

    if i!=3:
        wrong_file_name = True
    return wrong_file_name


def check_input(string, grading_file):
    user_input = input(string+' y/n/comment')
    if user_input == 'y':
        with open(grading_file, 'a') as f:
            f.write('\tYES')
    elif user_input == 'n':
        with open(grading_file, 'a') as f:
            f.write('\tNO')
    else:
        with open(grading_file, 'a') as f:
            f.write('\t'+user_input)


filesize = os.path.exists(grading_file)
if filesize == 0:
    with open(grading_file, 'w') as f:
        f.write('ID\tRuns_2d\tgraphs_2d\tcoordinate_2d\tmean_R_2d\tcondition_2d\tRuns_3d\tgraphs_3d\tcoordinate_3d\t'
                'mean_R_3d\tcondition_3d\t')
        #0 A2 B2 C2 D2 A3 B3 C3 D3 C

test_file=False
if dims==[2,3]:
    if not id_list:
        id_list=[1]
    for id in id_list:
        for path, dirs, files in os.walk(os.getcwd()):
            if path=='/Users/michal/PycharmProjects/Polymers/codes':
                if id==1:
                    test_file=files[0]
                else:
                    if str(id)+'.py' in files:
                        test_file = str(id)+'.py'
                    else:
                        break
                # test_file = '206061012.py'

                shutil.move(os.path.join(path, test_file), os.path.join(proj_path_name, test_file))
                with open(grading_file, 'a') as f:
                    f.write('\n' + test_file.split('.')[-2])
                break
        if test_file:
            break
    if not test_file:
        test_file = str(id) + '.py'
        shutil.move(os.path.join('/Users/michal/PycharmProjects/Polymers/tested_codes', test_file), os.path.join(proj_path_name, test_file))


elif dims==[3]:
    if id_list:
        test_file = id_list[0]+'.py'
    else:
        test_file = glob.glob('[0-9][0-9][0-9][0-9][0-9][0-9][0-9][0-9][0-9].py')[0]
print(test_file)


for dim in dims:
    print('TEST DIM', dim)
    filesize = os.path.exists(os.path.join(proj_path_name, input_file_name))
    if filesize != 0:
        os.remove(os.path.join(proj_path_name, input_file_name))
    shutil.copyfile(os.path.join(proj_path_name, input_file_name.replace('.', str(dim) + 'd.')), os.path.join(proj_path_name, input_file_name))

    # execfile(test_file)
    try:
        exec(open(test_file).read())
    except NameError as err:
        print(err)
        with open(grading_file, 'a') as f:
            f.write('\tError while running')
        check_input(string='what happened?', grading_file=grading_file)
        continue
    except SyntaxError as err:
        print(err)
        with open(grading_file, 'a') as f:
            f.write('\tError while running')
        check_input(string='what happened?', grading_file=grading_file)
        continue
    except ValueError as err:
        print(err)
        with open(grading_file, 'a') as f:
            f.write('\tError while running')
        check_input(string='what happened?', grading_file=grading_file)
        continue
    plt.show()

    with open(grading_file, 'a') as f:
        f.write('\tA{}:'.format(dim))
    check_input(string='does it run?', grading_file=grading_file)
    with open(grading_file, 'a') as f:
        f.write('\tB{}:'.format(dim))
    check_input(string='do the graphs make sense?', grading_file=grading_file)
    test_coordiante()
    with open(grading_file, 'a') as f:
        f.write('\tC{}:'.format(dim))
    check_input(string='are the coordinate OK?', grading_file=grading_file)
    plt.show()
    with open(grading_file, 'a') as f:
        f.write('\tD{}:'.format(dim))
    wrong_file_name = test_radii(dim, grading_file)

if wrong_file_name:
    with open(grading_file, 'a') as f:
        f.write('\tWrong_file_name')

with open(grading_file, 'a') as f:
    f.write('\tCOMMENTS:')
check_input(string='are there any comments?', grading_file=grading_file)


finish_test=True
if finish_test:
    os.rename(os.path.join(proj_path_name, test_file), os.path.join(proj_path_finish, test_file))
