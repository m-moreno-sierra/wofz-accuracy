#!/bin/env python

"""
Enumerates maximal polyominoes in contact with circumscribing circles of increasing diameter.
Prints a LaTeX table.
Used to generate the table in the cgt paper.
"""

from math import sqrt
import sys
sys.path.insert(0, '../shared')
import functool as fut

def polyomino_rows(s,t):
    d2 = s**2 + t**2
    rows = []
    nj2 = int(sqrt(d2)) + 1
    for k in range(2-s%2, nj2, 2):
        dx = sqrt(d2 - k**2)/2 - (t%2)/2
        row = 2 * int(dx) + (t%2)
        if row > 0:
            if k > s%2:
                rows.insert(0, row)
            rows.append(row)
    return rows

def sorted_polyominoes(N):
    Out = []
    for s in range(1,N):
        for t in range(s,N):
            d2 = s**2 + t**2
            if d2 > N**2:
                continue
            rows = polyomino_rows(s,t)
            area = 0
            for n in rows:
                area += n

            Out.append((s, t, d2, area, rows))

    return sorted(Out, key=lambda e: e[2])

if __name__ == '__main__':
    fut.print_provenience()
    P = sorted_polyominoes(10)
    nout = 1
    oldrows = []
    for n in range(len(P)):
        p = P[n]
        i, j, d2, area, rows = p
        sym = 2
        if i%2 == j%2:
            sym = 4
        if d2 > 100:
            break
        if rows == oldrows:
            print("& %i & %i & & & & & & \\\\" % (i, j))
            continue

        str_rows = ', '.join([str(n) for n in rows])
        rows2 = polyomino_rows(j, i)
        name2 = ""
        str_rows2 = ""
        if rows2 != rows:
            name2 = "$\\Pi'_{%i}$" % (nout)
            str_rows2 = ', '.join([str(n) for n in rows2])

        #    print("%2i %2i %2i %3i %3i %s" % (n, i, j, d2, area, rows))
        print("$\\Pi_{%i}$ & %i & %i & $D_%i$ & %i & %i & %s & %s & %s\\\\" %
              (nout, j, i, sym, d2, area, str_rows, name2, str_rows2))
        nout += 1
        oldrows = rows
    print()
