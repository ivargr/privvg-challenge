import numpy as np
import random
import math

def pick_random(weights):
    # implemented like https://github.com/pangenome/odgi/blob/630fa431bb747a7f645fddaddfe02979dc346597/src/algorithms/diffpriv.cpp#L66
    x = 0
    d = random.random() * sum(weights)
    for i, weight in enumerate(weights):
        if x + weight >= d:
            print(x + weight, d, sum(weights))
            return i
        x += weight


def exponential_mechanism_weight(utility, epsilon=0.01):
    u = math.log1p(utility)
    du = u - math.log1p(utility-1)
    #du = max(u - math.log1p(utility-1), math.log1p(utility+1)-u)
    #u = utility
    #du = u - (utility-1)
    weight = math.exp(epsilon*u / (2*du))
    return weight


epsilons = [0.00000000000001, 0.000001, 0.01, 0.1]
for u in range(1, 50):
    r = [exponential_mechanism_weight(u, epsilon) for epsilon in epsilons]
    print(u, r)

print()
for e in epsilons:
    e3 = exponential_mechanism_weight(3, e)
    e30 = exponential_mechanism_weight(30, e)

    print(e, e3, e30, e30 / e3)


"""
# simulate 5 random weights
weights = [2, 1, 0.99, 0.01]  # np.random.randint(0, 5, 5)
results = []
for i in range(10000):
    results.append(pick_random(weights))

for i in range(len(weights)):
    print(results.count(i) / len(results))
"""