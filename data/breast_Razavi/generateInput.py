#!/usr/bin/python
import sys
import glob
import os

# determine germline mutation

patientIdx = 0
for filename in glob.glob("processed_1/*.trees"):
    with open(filename) as f:
        f.readline()
        line = f.readline().rstrip("\n")
        if line == "":
            continue
        numTrees = int(line.split()[0])
        if numTrees > 10000:
            continue
        if numTrees == 0:
            continue

        line2 = line
        ok = True
        for treeIdx in range(numTrees):
            line = f.readline().rstrip("\n")
            numEdges = int(line.split()[0])
            if numEdges == 0:
                ok = False
                continue
            if treeIdx == 0:
                print line2, "for patient", os.path.basename(filename).rstrip(".trees")

            print numEdges + 1, "#edges, tree", treeIdx
            tree = []
            vertices = set([])
            inner_vertices = set([])
            for edgeIdx in range(numEdges):
                line = f.readline().rstrip("\n")
                s = line.split()
                tree.append([s[0].split(":")[-1], s[1].split(":")[-1]])
                vertices.add(tree[-1][0])
                vertices.add(tree[-1][1])
                inner_vertices.add(tree[-1][1])
            assert len(set.difference(vertices, inner_vertices)) == 1
            root = list(set.difference(vertices, inner_vertices))[0]
            for edge in tree:
                print edge[0], edge[1]
            print "GL", root
        if ok:
            patientIdx += 1
