import RNA
import random as rd
import math
import numpy as np
import matplotlib.pyplot as plt

#print RNA.fold('ACUGUACUGGCGGUAGCUGCUAGUUUACUAUCGGACGAUCU')

def boltzman(D, beta):
    return sum(math.exp(-beta*d) for d in D)

def reproRate(d, beta, Z):
    return math.exp(-beta*d)/Z

def mutate(s, mu, A):
    return ''.join(x if rd.random() > mu else rd.choice(list(set(A)-set(x))) for x in s)

def RNAevo(T, N, G, mu):
    L = len(T)
    A = ['A','U','G','C']
    beta = 2/float(L)


    # Generate N random sequences of length L uniformly
    seqs = [''.join(rd.choice(A) for _ in range(L)) for _ in range(N)]
    avgD = []
    for g in range(G):
        # Fold sequences
        folds = [RNA.fold(s)[0] for s in seqs]
        # Calculate distances
        dists = [RNA.bp_distance(t, T) for t in folds]
        avgD.append(np.mean(dists))
        # Calculate reproduction rates
        Z = boltzman(dists, beta)
        R = [reproRate(d, beta, Z) for d in dists]
        # Replicate population according to rates
        seqs = np.random.choice(seqs, size=N, p=R)
        # Mutate sequences
        seqs = [mutate(s, mu, A) for s in seqs]

    return avgD

gens = 500
for t in ['((((((((....))))))))', '((((..(((....)))))))', '(((....)))(((....)))']:
    avg1 = RNAevo(t, 100, gens, 0.01)
    avg2 = RNAevo(t, 100, gens, 0.02)
    avg3 = RNAevo(t, 100, gens, 0.05)
    avg4 = RNAevo(t, 100, gens, 0.1)
    x = range(gens)
    plt.plot(x, avg1, 'b.', x, avg2, 'r.', x, avg3, 'g.', x, avg4, 'y.')
    plt.show()
