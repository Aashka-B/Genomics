from pprint import pprint
import random

distanceMatrix = [[0, 12, 12, 13, 15, 15],
                  [12, 0, 2, 6, 8, 8],
                  [12, 2, 0, 6, 9, 9],
                  [13, 6, 6, 0, 8, 8],
                  [15, 8, 9, 8, 0, 4],
                  [15, 8, 9, 8, 4, 0]]

speciesList = ["M_Spacii", "T_Pain", "G_Unit", "Q_Doba", "R_Mani", "A_Finch"]

def findSmallest(dM):
    min_val = float('inf')
    candidates = []
    for i in range(len(dM)):
        for j in range(i+1, len(dM[i])):
            if 0 < dM[i][j] < min_val:
                min_val = dM[i][j]
                candidates = [(i, j)]
            elif dM[i][j] == min_val:
                candidates.append((i, j))
    # Randomly choose among the smallest if there are multiple
    return random.choice(candidates) if candidates else (None, None)

def updateMatrix(dM, row, col):
    newMat = []
    n = len(dM)
    for i in range(n):
        if i == row or i == col:
            continue
        newRow = [dM[i][j] for j in range(n) if j != row and j != col]
        newRow.insert(row, (dM[row][i] + dM[col][i]) / 2)
        newMat.append(newRow)
    newRow = [(dM[row][i] + dM[col][i]) / 2 for i in range(n) if i != row and i != col]
    newRow.insert(row, 0)
    newMat.insert(row, newRow)
    return newMat

def updateSpecies(sp, r, c, branch_length):
    merged_species = f"({sp[r]},{sp[c]}):{branch_length}"
    sp[r] = merged_species
    del sp[c]
    return sp

def UPGMA(dM, sp):
    while len(dM) > 1:
        leastRow, leastCol = findSmallest(dM)
        branch_length = dM[leastRow][leastCol]
        dM = updateMatrix(dM, leastRow, leastCol)
        sp = updateSpecies(sp, leastRow, leastCol, branch_length)
        print(f"Updated Distance Matrix:")
        pprint(dM)
        print(f"Updated Species List: {sp}")

UPGMA(distanceMatrix, speciesList)
