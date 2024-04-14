# %%
import numpy as np
import pandas as pd
from scipy.sparse import coo_matrix
import argparse



MINGAP = 2
SIGNAL_THRESHOLD = 0.7

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

def read_finer_matrix(input, bed):
    
    # get diag weight and return it as  list format
    finer_diag = np.zeros(len(bed.index)+1)
    with open(input, 'r') as file:
        i = 1
        for line in file:
            arr = line.strip().split("\t")
            if arr[0] == arr[1]:
                finer_diag[int(arr[0])] = arr[2]
    return finer_diag

def findGap(matrix):
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
    avg = np.mean(sorted(degrees)[:-1000])

    variance = []
    for i in degrees:
        variance.append((i-avg)/avg)
    gap = {}
    errToler = {}
    singal = None
    for i in range(0, len(variance)):
        if singal and (-variance[i] > SIGNAL_THRESHOLD):
            gap [singal] += 1
        if (singal == None) and (-variance[i] > SIGNAL_THRESHOLD):
            gap [i] = 1
            singal = i
            errToler[i] = 0
        if singal and (-variance[i] <= SIGNAL_THRESHOLD):
            if errToler[singal] < 2:
                errToler[singal] += 1
                gap [singal] += 1
            else:
                k = i-1
                while -variance[k] <= SIGNAL_THRESHOLD:
                    gap[singal] -= 1
                    errToler[singal] -= 1
                    k -= 1
                singal = None
    candidates = {key+4: value for key, value in gap.items()}
    # Filter out the too short candidated gaps that are likely to be noise. 
    filter_candidates = {key: value for key, value in 
                         candidates.items() if value > MINGAP}
    return filter_candidates


def findPosition(gap, selfweights, bed1, bed2):
    centremere = {}
    mean_weight = np.mean(sorted(selfweights)[:-1000]) 
    # print(len(selfweights))
    # print(mean_weight)
    # print(selfweights)
    for key, item in gap.items():
        if bed1[bed1["id"] == key]["chr"].iloc[0] != bed1[bed1["id"] == (key+item)]["chr"].iloc[0]:
            continue
        chr = bed1[bed1["id"] == key]["chr"].iloc[0]
        start = bed1[bed1["id"] == key]["start"].iloc[0]
        end = bed1[bed1["id"] == (key+item)]["start"].iloc[0]

        Lstart = int(bed2[(bed2["chr"] == chr) & (bed2["start"]==start)]["id"].iloc[0])
        Rstart = int(bed2[(bed2["chr"] == chr) & (bed2["start"]==end)]["id"].iloc[0])
        # print(bed1[bed1["id"]==key])
        # print(start, end ,Lstart, Rstart)
        # Find right bound
        lbound, rbound = None, None
        
        for i in range(Lstart, Lstart + 10):
            # print(i)
            if selfweights[i]/mean_weight < 0.4:
                # print(selfweights[i])
                lbound = bed2[bed2["id"]==i]["start"].iloc[0]
                break
        # Find right bound
        for i in range(Rstart, Rstart+10):
            if selfweights[i]/mean_weight > 0.7:
                # print(selfweights[i])
                rbound = bed2[bed2["id"]==i]["end"].iloc[0]
                break

            
        if chr not in centremere:
            centremere[chr] = []
        centremere[chr].append([lbound, rbound])
    return centremere

# %%
def parse_arguments():
    parser = argparse.ArgumentParser(description='A program for find centremere position based on Hi-C.')
    parser.add_argument('matrix1', type=str, help='Path to matrix of 100000 file')
    parser.add_argument('matrix2', type=str, help='Path to matrix of 200000 file')
    parser.add_argument('bed1', type=str, help='Path to bed file of matrix1 ')
    parser.add_argument('bed2', type=str, help='Path to bed file of matrix2')
    parser.add_argument('--MINGAP', type=int, default=2, help='Minimum gap value n*100000 (default: 2)')
    parser.add_argument('--SIGNAL_THRESHOLD', type=float, default=0.7, help='Signal threshold value (default: 0.7)')
    parser.add_argument('-o', type=str, help='Output file')
    args = parser.parse_args()
    return args

if __name__ == "__main__":
    args = parse_arguments()
    # input = "../hicpro/matrix/1_100000.matrix"
    # finer_input = "../hicpro/matrix/1_20000.matrix"
    # bed_file1 = "../hicpro/matrix/1_100000_abs.bed"
    # bed_file2 = "../hicpro/matrix/1_20000_abs.bed"
    # print("reading bed")
    bed1 = pd.read_csv(args.bed1, sep = "\t", header= None, 
                       names=['chr', 'start', 'end', 'id'])
    bed2 = pd.read_csv(args.bed2, sep = "\t", header= None,
                       names=['chr', 'start', 'end', 'id'])
    # print("read bed finished")
    gap = findGap(creat_matrix(args.matrix1))
    # print("Candidated centremeres:")
    # for i in gap:
    #     tmp = bed1[bed1["id"]==i]
    #     print(tmp["chr"].iloc[0],tmp["start"].iloc[0], i,gap[i])
    finer_diag = read_finer_matrix(args.matrix2, bed2)
    centremere = findPosition(gap, finer_diag, bed1, bed2)
    print(centremere)