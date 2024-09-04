from multiprocessing import Pool, get_context


def return_x(x):

    suma = 0

    for i in range(100):
        for j in range(100):
            for k in range(100):
                suma += i+j+k

    return suma

if __name__ == 'name':

    x_list = [(x) for x in range(100)]
    args = x_list
    print(x_list)

    with Pool(8) as p:

        extended = p.starmap(return_x, args)

    print(extended)
