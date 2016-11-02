import RNA
import random as rd
import math
import numpy as np
import matplotlib.pyplot as plt

# ============================================ Student info methods================================================
def get_student_name():
	# @TO_STUDENT: Write your name here
	student_name = "Jeremy Georges-Filteau"
	if not student_name:
		raise ValueError("Error: you forgot to add your name in get_student_name method.")
	return student_name

def get_student_id():
	# @TO_STUDENT: Write your student id here
	student_id = "260713547"
	if not student_id:
		raise ValueError("Error: you forgot to add your student id in get_student_id method.")
	return student_id
# =================================================================================================================

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

print("This is a solution of %s, student_id is %s" % (get_student_name(), get_student_id()) )

gens = 500
for t in ['((((((((....))))))))', '((((..(((....)))))))', '(((....)))(((....)))']:
    avg1 = RNAevo(t, 100, gens, 0.01)
    avg2 = RNAevo(t, 100, gens, 0.02)
    avg3 = RNAevo(t, 100, gens, 0.05)
    avg4 = RNAevo(t, 100, gens, 0.1)
    x = range(gens)
    plt.plot(x, avg1, 'b.', label='0.01')
    plt.plot(x, avg2, 'r.', label='0.02')
    plt.plot(x, avg3, 'g.', label='0.05')
    plt.plot(x, avg4, 'y.', label='0.1')
    plt.ylabel('Average distance from the target')
    plt.xlabel('Generations')
    plt.title(t)
    plt.legend()
    plt.show()
