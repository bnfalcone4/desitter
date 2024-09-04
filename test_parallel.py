from multiprocessing import Pool, get_context


def return_x(x, y):

    suma = 0

    for i in range(200):
        for j in range(200):
            for k in range(200):
                suma += i+j+k

    return suma

if __name__ == '__main__':

    x_list = [(x, x) for x in range(100)]
    args = x_list
    print(x_list)

    with Pool(8) as p:

        extended = p.starmap(return_x, args)

    print(extended)
