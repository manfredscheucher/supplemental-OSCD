import itertools
from scipy.special import binom
from sys import argv
from ast import literal_eval
from os import path

# We use the following variable names:
#  n for the dimension of the cube
#  F for a family of SCDs
#  D for an SCD
#  C for a chain

def is_symmetric_chain(C):
    if len(C) == 0:
        print("empty chains are not allowed")
        return False
    n = len(C[0])

    # test if all vertices are bitstrings
    for v in C:
        for b in v:
            if b not in [0, 1]:
                print("bitstring representation of vertex is corrupted")
                return False

    # test if all bitstrings have the same length
    for v in C:
        if len(v) != n:
            print("vertex of wrong length")
            return False

    # test if consecutive vertices along the chain differ in a single bit
    for i in range(1, len(C)):
        diff = 0
        for j in range(len(C[0])):
            if C[i][j] < C[i - 1][j]:
                print("chain is not level-increasing")
                return False
            else:
                diff += C[i][j] - C[i - 1][j]
        if diff != 1:
            print("chain is not saturated")
            return False

    # test if start and end vertex lie on symmetric levels
    k = (n + 1 - len(C)) / 2
    if sum(C[0]) != k or sum(C[-1]) != n - k:
        print("chain is not symmetric")
        return False

    return True

def is_scd(D, n):
    if len(D[0][0]) != n:
        print("SCD has wrong dimension")
        return False
    if len(D) == 0:
        print("SCD must not be empty")
        return False

    # test if all chains are symmetric
    for C in D:
        if not is_symmetric_chain(C):
            return False

    # test if the chains cover the entire cube
    vertices = set().union(*D)
    for v in itertools.product([0, 1], repeat=n):
        if v not in vertices:
            print("chains do not cover the entire cube")
            return False

    # test if the number of chains is correct. As the chains cover the entire
    # cube, this implies disjointness.
    if len(D) != binom(n, n // 2):
        print("chains are not disjoint")
        return False

    return True

def are_scds(F,n):
    # test if all SCDs in the family are indeed SCDs in the n-cube
    for D in F:
        if not is_scd(D, n):
            return False

    return True

def is_edge_disjoint_family(F, n):
  
    # test edge-disjointness between any pair of SCDs
    for D1, D2 in itertools.combinations(F, 2):
        for C1 in D1:
            for C2 in D2:
                for i in range(1, len(C1)):
                    if C1[i - 1] in C2 and C1[i] in C2:
                        print(C1, C2)
                        print("chains are not edge-disjoint")
                        return False

    return True

def is_almost_orthogonal_family(F, n):

    # test almost-orthogonality
    for D1, D2 in itertools.combinations(F, 2):
        for C1 in D1:
            for C2 in D2:
                intersection = (set(C1) - {tuple([0] * n)}) & set(C2)
                if len(intersection) > 1:
                    print("chains are not almost-orthogonal")
                    return False

    return True

def is_good_family(F, n):
    # if n is odd, test if union of all edges of 2-element chains is unicyclic
    if n == 1:
        return True
    if n % 2 == 1:
        adj = {}
        for v in itertools.permutations([0] * ((n + 1) // 2) + [1] * ((n - 1) // 2)):
            adj[v] = []
        for v in itertools.permutations([0] * ((n - 1) // 2) + [1] * ((n + 1) // 2)):
            adj[v] = []
        for D in F:
            for C in D:
                if len(C) == 2:
                    adj[C[0]].append(C[1])
                    adj[C[1]].append(C[0])
        return is_unicyclic(adj)
	# if n is even, check if all length 1 chains are distinct
	# (goodness is not defined and visited for even dimensions in our paper, but in Spink's paper)
    else:
        for Da, Db in itertools.combinations(F, 2):
            for Ca in Da:
                if len(Ca) == 1 and Ca in Db:
                    print("SCDs have colliding 1-element chains")
                    return False

        return True

def is_unicyclic(adj):
    visited = {v: False for v in adj}

    # test whether connected component is unicyclic via depth-first-search
    def is_unicyclic_component(v):
        cycle = False
        visited[v] = True
        stack = [(v, w) for w in adj[v]]
        while stack:
            v, w = stack.pop()
            if visited[w]:
                if not cycle:
                    stack.remove((w, v))
                    cycle = True
                else:
                    return False
            else:
                visited[w] = True
                stack.extend([(w, x) for x in adj[w] if x != v])
        return True

    # test all connected components
    for v in adj:
        if not visited[v]:
            if not is_unicyclic_component(v):
                print("graph of 2-element chains is not unicyclic")
                return False

    return True

def complementary_chains(C1, C2):
    if len(C1) != len(C2):
        return False

    for i in range(len(C1)):
        for j in range(len(C1[i])):
            if C1[i][j] + C2[-1 - i][j] != 1:
                return False

    return True

def complementary_scds(D1, D2):
    if len(D1) != len(D2):
        return False

    for C1 in D1:
        found_compl_chain = False
        for C2 in D2:
            if complementary_chains(C1, C2):
                found_compl_chain = True
                break
        if not found_compl_chain:
            return False

    return True

def check_scds(path_to_file, check_ortho, check_edge, check_compl):
    with open(path_to_file) as file:
        line = file.readline()
		# read entire family of SCDs from the file
        F = literal_eval(line)
        F = [[[tuple([int(b) for b in v]) for v in C] for C in D] for D in F]

        n = len(F[0][0][0])
        k = len(F)
        print ("dimension:", n)
        print ("number of SCDs:", k)

        print("All decompositions are SCDs: ", end="")
        print(are_scds(F,n))

        if check_ortho:
            print("SCDs are almost-orthogonal: ", end="")
            print(is_almost_orthogonal_family(F, n))
            print("SCDs are good: ", end="")
            print(is_good_family(F, n))

        if check_edge:
            print("SCDs are edge-disjoint: ", end="")
            print(is_edge_disjoint_family(F, n))

        if check_compl:
            print("complementary pairs:")
            compl = False
            for i in range(k):
                for j in range(i + 1, k):
                    if complementary_scds(F[i], F[j]):
                        print("SCDs", i + 1, "and", j + 1, "are complementary")
                        compl = True
            if not compl:
                print("None")

def help():
    print("verification tool for almost-orthogonal and edge-disjoint SCDs")
    print("usage: python3 verify.py file")

if __name__ == "__main__":
    if len(argv) < 2:
        help()
        exit()

    path_to_file = path.abspath(argv[1])

    check_scds(path_to_file, True, True, True)
