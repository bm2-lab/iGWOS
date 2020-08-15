hitScoreM = [0, 0, 0.014, 0, 0, 0.395, 0.317, 0, 0.389, 0.079, 0.445, 0.508, 0.613, 0.851, 0.732, 0.828, 0.615, 0.804,
             0.685, 0.583]

def calcMitScore(string1, string2, startPos=0):
    """
    The MIT off-target score
    see 'Scores of single hits' on http://crispr.mit.edu/about
    startPos can be used to feed sequences longer than 20bp into this function
    the most likely off-targets have a score of 100
    >>> int(calcMitScore("GGGGGGGGGGGGGGGGGGGG","GGGGGGGGGGGGGGGGGGGG"))
    100
    # mismatches in the first three positions have no effect
    >>> int(calcMitScore("AGGGGGGGGGGGGGGGGGGG","GGGGGGGGGGGGGGGGGGGG"))
    100
    # less likely off-targets have lower scores
    >>> int(calcMitScore("GGGGGGGGGGGGGGGGGGGG","GGGGGGGGGGGGGGGGGGGA"))
    41
    """
    # The Patrick Hsu weighting scheme
    # print string1, string2
    if len(string1) == len(string2) == 23:
        string1 = string1[:20]
        string2 = string2[:20]

    assert (len(string1) == len(string2) == 20)

    dists = []  # distances between mismatches, for part 2
    mmCount = 0  # number of mismatches, for part 3
    lastMmPos = None  # position of last mismatch, used to calculate distance

    score1 = 1.0
    for pos in range(0, len(string1)):
        if string1[pos] != string2[pos]:
            mmCount += 1
            if lastMmPos != None:
                dists.append(pos - lastMmPos)
            score1 *= 1 - hitScoreM[pos]
            lastMmPos = pos

    # 2nd part of the score - distribution of mismatches
    if mmCount < 2:  # special case, not shown in the paper
        score2 = 1.0
    else:
        avgDist = sum(dists) / len(dists)
        score2 = 1.0 / (((19 - avgDist) / 19.0) * 4 + 1)

    # 3rd part of the score - mismatch penalty
    if mmCount == 0:  # special case, not shown in the paper
        score3 = 1.0
    else:
        score3 = 1.0 / (mmCount ** 2)

    score = score1 * score2 * score3 * 100
    return score
