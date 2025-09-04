import math

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


def scoreModels():
    codingMatrix = getProbs("./codingModel.tab")  # Load coding model
    noncodingMatrix = getProbs("./noncodingModel.tab")  # Load non-coding model
    id2ancestorSeq = getSeq("./Ancestor.fa")  # Load ancestor sequences
    id2spaciiSeq = getSeq("./Spacii.fa")  # Load M. Spacii sequences
    allID = list(id2ancestorSeq.keys())  # Get all sequence IDs

    for ID in allID:
        cScore = 0  # Log likelihood score for the coding model
        nScore = 0  # Log likelihood score for the non-coding model

        ancestorSeq = id2ancestorSeq[ID]
        spaciiSeq = id2spaciiSeq[ID]

        # Ensure sequences are of equal length and divisible by 3
        seqLength = min(len(ancestorSeq), len(spaciiSeq))
        seqLength -= seqLength % 3  # Adjust to ensure divisibility by 3

        for i in range(0, seqLength, 3):
            ancestorCodon = ancestorSeq[i:i+3]
            spaciiCodon = spaciiSeq[i:i+3]

            if ancestorCodon in modelCodons and spaciiCodon in modelCodons:
                aIndex = modelCodons.index(ancestorCodon)
                sIndex = modelCodons.index(spaciiCodon)

                # Sum log probabilities for both models
                cScore += math.log(codingMatrix[aIndex][sIndex])
                nScore += math.log(noncodingMatrix[aIndex][sIndex])

        # Classification based on which score is higher
        if cScore > nScore:
            print(ID + " is likely coding. Scores: Coding = " + str(cScore) + ", Non-Coding = " + str(nScore))
        else:
            print(ID + " is likely NOT coding. Scores: Coding = " + str(cScore) + ", Non-Coding = " + str(nScore))


def getProbs(f1):
    f = open(f1)
    pMatrix = []
    for line in f:
        tmp = line.rstrip().split("\t")
        tmp = [float(i) for i in tmp]
        pMatrix.append(tmp)
    return pMatrix


def getSeq(filename):
    f = open(filename)
    id2seq = {}
    currkey = ""
    for line in f:
        if line.find(">") == 0:
            currkey = (line[1:].split("|")[0])
            id2seq[currkey] = ""
        else:
            id2seq[currkey] = id2seq[currkey] + line.rstrip()
    f.close()
    return id2seq


scoreModels()
