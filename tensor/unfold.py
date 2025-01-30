import numpy as np


def unfold(C, k): 
    """Mode k (>= 1) unfolding of tensor C. Following conventional mode, not NumPy mode"""
    ndims = len(C.shape)
    permutation = tuple( list(range(0, ndims-2)) + [ndims-1, ndims-2] )  # column --> row fiber
    C = C.transpose(permutation)
    k = ndims - k  # conventional --> NumPy fiber mode
    permutation = tuple( list(range(0, k)) + list(range(k+1, ndims)) + [k] )  # as row-wise fiber
    c = C.transpose(permutation)
    c = c.reshape( (-1, C.shape[k]) )  # unfold as row fiber
    c = c.transpose()  # row --> column fiber of unfolded matrix
    return c


if __name__ == '__main__':
    C = np.array([[[[1, 4, 7, 10], [2, 5, 8, 11], [3, 6, 9, 12]],
                   [[13, 16, 19, 22], [14, 17, 20, 23], [15, 18, 21, 24]]],

                  [[[25, 28, 31, 34], [26, 29, 32, 35], [27, 30, 33, 36]],
                   [[37, 40, 43, 46], [38, 41, 44, 47], [39, 42, 45, 48]]]])  # 4D tensor
    print(C)
    print('\n')
    
    for k in range(1, 5):
        print(unfold(C, k))
        print()
    
