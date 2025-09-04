import math
import numpy as np
from sklearn.metrics import roc_curve, auc


modelCodons = ['TTT', 'TTC', 'TTA', 'TTG', 'CTT', 'CTC', 'CTA',
               'CTG', 'ATT', 'ATC', 'ATA', 'ATG', 'GTT', 'GTC',
               'GTA', 'GTG', 'TCT', 'TCC', 'TCA', 'TCG', 'AGT',
               'AGC', 'CCT', 'CCC', 'CCA', 'CCG', 'ACT', 'ACC',
               'ACA', 'ACG', 'GCT', 'GCC', 'GCA', 'GCG', 'TAT',
               'TAC', 'CAT', 'CAC', 'CAA', 'CAG', 'AAT', 'AAC',
               'AAA', 'AAG', 'GAT', 'GAC', 'GAA', 'GAG', 'TGT',
               'TGC', 'CGT', 'CGC', 'CGA', 'CGG', 'AGA', 'AGG',
               'GGT', 'GGC', 'GGA', 'GGG', 'TGG', 'TAA', 'TAG',
               'TGA']


def getProbs(f1):
    f = open(f1)
    pMatrix = []
    for line in f:
        tmp = line.rstrip().split("\t")
        tmp = [float(i) for i in tmp]
        pMatrix.append(tmp)
    f.close()
    return pMatrix


def getSeq(filename):
    f = open(filename)
    id2seq = {}
    currkey = ""
    for line in f:
        if line.startswith(">"):
            currkey = line[1:].strip().split("|")[0]
            id2seq[currkey] = ""
        else:
            id2seq[currkey] += line.rstrip()
    f.close()
    return id2seq


def scoreModels():
    codingMatrix = getProbs("./codingModel.tab")
    noncodingMatrix = getProbs("./noncodingModel.tab")
    id2ancestorSeq = getSeq("./Ancestor.fa")
    id2spaciiSeq = getSeq("./Spacii_2100.fa")
    allID = list(id2ancestorSeq.keys())

    true_labels = []
    likelihood_ratios = []

    for ID in allID:
        cScore = 0
        nScore = 0

        ancestorSeq = id2ancestorSeq[ID]
        spaciiSeq = id2spaciiSeq[ID]

        seqLength = min(len(ancestorSeq), len(spaciiSeq))
        seqLength -= seqLength % 3

        for i in range(0, seqLength, 3):
            ancestorCodon = ancestorSeq[i:i+3]
            spaciiCodon = spaciiSeq[i:i+3]

            if ancestorCodon in modelCodons and spaciiCodon in modelCodons:
                aIndex = modelCodons.index(ancestorCodon)
                sIndex = modelCodons.index(spaciiCodon)
                cScore += math.log(codingMatrix[aIndex][sIndex])
                nScore += math.log(noncodingMatrix[aIndex][sIndex])

        likelihood_ratio = cScore - nScore
        likelihood_ratios.append(likelihood_ratio)
        true_labels.append(1 if '_c_' in ID else 0)

    return true_labels, likelihood_ratios


def printFPRandTPR(true_labels, likelihood_ratios):
    fpr, tpr, thresholds = roc_curve(true_labels, likelihood_ratios)
    for i in range(len(thresholds)):
        print(f"Threshold: {thresholds[i]}, FPR: {fpr[i]}, TPR: {tpr[i]}")


# Main execution
if __name__ == "__main__":
    true_labels, likelihood_ratios = scoreModels()
    printFPRandTPR(true_labels, likelihood_ratios)
