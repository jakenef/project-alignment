def align(
        seq1: str,
        seq2: str,
        match_award=-3,
        indel_penalty=5,
        sub_penalty=1,
        banded_width=-1,
        gap_open_penalty=0,
        gap='-',
) -> tuple[float, str | None, str | None]:
    """
        Align seq1 against seq2 using Needleman-Wunsch
        Put seq1 on left (j) and seq2 on top (i)
        => matrix[i][j]
        :param seq1: the first sequence to align; should be on the "left" of the matrix
        :param seq2: the second sequence to align; should be on the "top" of the matrix
        :param match_award: how many points to award a match
        :param indel_penalty: how many points to award a gap in either sequence
        :param sub_penalty: how many points to award a substitution
        :param banded_width: banded_width * 2 + 1 is the width of the banded alignment; -1 indicates full alignment
        :param gap_open_penalty: how much it costs to open a gap. If 0, there is no gap_open penalty
        :param gap: the character to use to represent gaps in the alignment strings
    """

    print(
        f"align params -> seq1={seq1!r}, seq2={seq2!r}, "
        f"match_award={match_award}, indel_penalty={indel_penalty}, "
        f"sub_penalty={sub_penalty}, banded_width={banded_width}, "
        f"gap_open_penalty={gap_open_penalty}, gap={gap!r}"
    )

    n = len(seq2) # rows
    m = len(seq1) # cols

    #1. Setup tables

    editCostTable = [[0]*(m+1) for _ in range(n+1)]

    backPointerTable = [[""]*(m+1) for _ in range(n+1)]
    
    #2. Initialize tables with starter values

    for i in range(1, n + 1):
        editCostTable[i][0] = editCostTable[i-1][0] + indel_penalty
        backPointerTable[i][0] = "U"

    for j in range(1, m + 1):
        editCostTable[0][j] = editCostTable[0][j - 1] + indel_penalty
        backPointerTable[0][j] = "L"

    # 3. Build matrix

    for i in range(1, n + 1):
        for j in range(1, m + 1):
            diagonal = editCostTable[i - 1][j - 1] + (match_award if seq2[i - 1] == seq1[j - 1] else sub_penalty)
            left = editCostTable[i][j - 1] + indel_penalty
            up = editCostTable[i - 1][j] + indel_penalty

            candidates = [
                (diagonal, 0, "D"),
                (left, 1, "L"),
                (up, 2, "U")
            ]

            min_cost, _, direction = min(candidates)

            editCostTable[i][j] = min_cost
            backPointerTable[i][j] = direction

    score = editCostTable[n][m]

    # 4. Trace back to build aligned strings

    i = n
    j = m 

    aligned1 = ""
    aligned2 = ""

    while (i > 0 or j > 0):
        if backPointerTable[i][j] == 'D':
            aligned1 = seq1[j - 1] + aligned1
            aligned2 = seq2[i - 1] + aligned2
            i -= 1
            j -= 1
        elif backPointerTable[i][j] == 'L':
            aligned1 = seq1[j - 1] + aligned1
            aligned2 = gap + aligned2
            j -= 1
        else: # == "U"
            aligned1 = gap + aligned1
            aligned2 = seq2[i - 1] + aligned2
            i -= 1

    printSequenceTable(editCostTable, seq1, seq2)
    print("-----------------------------------------------")
    printSequenceTable(backPointerTable, seq1, seq2)



    return score, aligned1, aligned2

def printSequenceTable(table: list[list], seq1: str, seq2: str) -> None:
    # table is assumed to be (len(seq1)+1) x (len(seq2)+1)
    if not table:
        print("Empty table")
        return

    rows = len(table)
    cols = len(table[0])

    # first row: "-" followed by seq1 characters
    header = [' '] + ['-'] + list(seq1)
    print('\t'.join(str(x) for x in header))

    # subsequent rows: first column is "-" then seq1 characters, then table values
    for i in range(rows):
        label = '-' if i == 0 else seq2[i - 1]
        values = [str(table[i][j]) for j in range(cols)]
        print('\t'.join([str(label)] + values))

if __name__ == "__main__":
    score, aligned1, aligned2 = align("A", "C", match_award=-1, sub_penalty=5, indel_penalty=1)
    print("Score:", score)
    print("Aligned 1:", aligned1)
    print("Aligned 2:", aligned2)