from multiprocessing import Pool, get_context


def return_x(x,y):

	return x

x_list = [(x,x) for x in range(100)]
args = x_list

with Pool(8) as p:

    extended = p.starmap(return_x, args)

print(extended)


"""
#!/usr/bin/env python3
from functools import partial
from itertools import repeat
from multiprocessing import Pool, freeze_support

def func(a, b):
    return a + b

def main():
    a_args = [1,2,3]
    second_arg = 1
    with Pool() as pool:
        L = pool.starmap(func, [(1, 1), (2, 1), (3, 1)])
        M = pool.starmap(func, zip(a_args, repeat(second_arg)))
        N = pool.map(partial(func, b=second_arg), a_args)
        assert L == M == N

    print(L)

if __name__=="__main__":
    #freeze_support()
    main()
"""
