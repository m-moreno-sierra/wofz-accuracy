#!/bin/env python

"""
Generates coloring for subdomain map.
"""
# Using code by Divyanshu Mehta from https://www.geeksforgeeks.org/m-coloring-problem/


class Graph():

    def __init__(self, vertices):
        self.V = vertices
        self.graph = [[0 for column in range(vertices)]
                      for row in range(vertices)]

    # A utility function to check
    # if the current color assignment
    # is safe for vertex v
    def isSafe(self, v, color, c):
        for i in range(self.V):
            if self.graph[v][i] == 1 and color[i] == c:
                return False
        return True

    # A recursive utility function to solve m
    # coloring  problem
    def graphColorUtil(self, m, color, v):
        if v == self.V:
            return True

        for c in range(1, m + 1):
            if self.isSafe(v, color, c) == True:
                color[v] = c
                if self.graphColorUtil(m, color, v + 1) == True:
                    return True
                color[v] = 0

    def graphColoring(self, m):
        color = [0] * self.V
        if self.graphColorUtil(m, color, 0) == None:
            raise Exception("No solution exists")

        return color

def process_neighbor(C, A, j, i, dj, di):
    if C[j][i] != C[j+dj][i+di] and C[j][i]!=-1 and C[j+dj][i+di]!=-1:
        A[C[j][i]][C[j+dj][i+di]] = 1
        A[C[j+dj][i+di]][C[j][i]] = 1

# Driver Code
if __name__ == '__main__':
    with open("../../lib/w_taylor_cover.c") as f:
        ta = f.readlines()
    C = []
    nmax = -1
    for t in ta[2:-2]:
        c = []
        for w in t.rstrip(",\n").split(","):
            n = int(w.strip())
            c.append(n)
            nmax = max(n, nmax)
        C.append(c)
    A = [[0 for i in range(nmax+1)] for j in range(nmax+1)]
    for j in range(len(C)):
        for i in range(len(C[j])):
            if j>0:
                process_neighbor(C, A, j, i, -1, 0)
            if j<len(C)-1:
                process_neighbor(C, A, j, i, +1, 0)
            if i>0:
                process_neighbor(C, A, j, i, 0, -1)
            if i<len(C[j])-1:
                process_neighbor(C, A, j, i, 0, +1)

    g = Graph(nmax+1)
    g.graph = A
    color = g.graphColoring(4)

    for c in C:
        for n in c:
            if n==-1:
                out = 0
            else:
                out = color[n]
            print(f"{out} ", end="")
        print()
