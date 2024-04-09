# %%
import numpy as np
from scipy.sparse import coo_matrix

def creat_matrix(input):
    '''
    Read Hi-C sparse martix from a text file with the format is
    1 3 130
    2 2 500
    2 3 231
    The sparse will be switch to dense matrix in this function.
    '''
    data = np.loadtxt(input, delimiter="\t", dtype=int)
    sparse_matrix = coo_matrix((data[:, 2], (data[:, 0], data[:, 1])))
    dense_matrix = sparse_matrix.toarray()
    return dense_matrix



def findGap(matrix, threshold):
    # This function is finding Candidates for Centremeres and they will be return
    # with dict.
    degrees = []
    # Calculate degrees of nodes
    for i in range(4, len(matrix)-3):
        ddq = 0
        for j in range(i-3, i):
            ddq += matrix[j][i]
        for j in range(i,i+3):
            ddq += matrix[i][j]
        degrees.append(ddq)
    avg = np.mean(degrees)

    variance = []
    for i in degrees:
        variance.append((i-avg)/avg)
    gap = {}
    errToler = {}
    singal = None
    for i in range(0, len(variance)):
        if singal and (-variance[i] > threshold):
            gap [singal] += 1
        if (singal == None) and (-variance[i] > threshold):
            gap [i] = 1
            singal = i
            errToler[i] = 0
        if singal and (-variance[i] <= threshold):
            if errToler[singal] < 2:
                errToler[singal] += 1
                gap [singal] += 1
            else:
                k = i-1
                while -variance[k] <= threshold:
                    gap[singal] -= 1
                    errToler[singal] -= 1
                    k -= 1
                singal = None
    candidates = {key+4: value for key, value in gap.items()}
    return candidates


def findPosition(gap, finerMatrix):
    pass


if __name__ == "__main__":
    input = "chr1.matrix"
    matrix = creat_matrix(input)
    # for i in range(1, len(matrix)):
    #     print(matrix[i][i])
    gap = findGap(matrix, 0.7)
    print(gap)