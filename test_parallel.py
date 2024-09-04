import argparse

from multiprocessing import Pool, get_context


def return_x(x, y):

    suma = 0

    for i in range(200):
        for j in range(200):
            for k in range(200):
                suma += i+j+k

    return suma

if __name__ == '__main__':

    parser = argparse.ArgumentParser(description='Run paralle test')

    # required arguments
    parser.add_argument('-p','--parallel', type=int, help='Number of cpus', required=True)

    args = parser.parse_args()

    x_list = [(x) for x in range(100)]

    if args.parallel == 1:

        extended = [return_x(x, x) for x in x_list]

    else:

        with Pool(args.parallel) as p:

            extended = p.starmap(return_x, x_list)

    print(extended)