def findRuns(lst):
    """ yield (start, end) tuples for all runs of ident. numbers in lst
    >>> list(findRuns([1,1,1,0,0,1,0,1,1,1]))
    [(0, 3), (5, 6), (7, 10)]
    """
    start, end = False, False

    for i, x in enumerate(lst):
        if x and start is False:
            start = i
        if x == 0 and start is not False and end is False:
            end = i - 1
        if start is not False and end is not False:
            yield start, end + 1  # and len is (end-start)
            start, end = False, False

    if start is not False:
        yield start, i + 1  # and len is (end-start)
def calcCropitScore(guideSeq, otSeq):
    """
    see http://www.ncbi.nlm.nih.gov/pmc/articles/PMC4605288/ PMID 26032770
    >>> int(calcCropitScore("GGGGGGGGGGGGGGGGGGGG","GGGGGGGGGGGGGGGGGGGG"))
    650
    # mismatch in 3' part
    >>> int(calcCropitScore("GGGGGGGGGGGGGGGGGGGA","GGGGGGGGGGGGGGGGGGGG"))
    575
    # mismatch in 5' part
    >>> int(calcCropitScore("AGGGGGGGGGGGGGGGGGGG","GGGGGGGGGGGGGGGGGGGG"))
    642
    # only mismatches -> least likely offtarget
    >>> int(calcCropitScore("AAAAAAAAAAAAAAAAAAAA","GGGGGGGGGGGGGGGGGGGG"))
    -27
    """
    if len(guideSeq) == 23:
        guideSeq = guideSeq[:20]
        otSeq = otSeq[:20]

    assert (len(guideSeq) == len(otSeq) == 20)

    penalties = [5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 70, 70, 70, 70, 70, 50, 50, 50, 50, 50]
    score = 0.0

    # do the score only for the non-mism positions
    misList = []
    score = 0.0
    for i in range(0, 20):
        if guideSeq[i] != otSeq[i]:
            misList.append(1)
        else:
            misList.append(0)
            score += penalties[i]

    # get the runs of mismatches and update score for these positions
    consecPos = set()
    singlePos = set()
    for start, end in findRuns(misList):
        if end - start == 1:
            score += -penalties[start] / 2.0
        else:
            # mean if they happen to fall into different segments
            startScore = penalties[start]
            endScore = penalties[end - 1]
            score += -((startScore + endScore) / 2.0)
    return score