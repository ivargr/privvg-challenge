import numpy as np
from scipy.stats import poisson
import matplotlib.pyplot as plt
COVERAGE = 30


def get_counts():
    sample_names = {i: line.strip() for i, line in enumerate(open(snakemake.input.sample_names).readlines())}
    counts = np.load(snakemake.input.counts)
    priv_counts = np.load(snakemake.input.priv_counts)
    return sample_names, counts, priv_counts


def predict_index(counts, priv_counts):
    lower_than_expected = priv_counts == 0

    scores = np.sum((counts >= 1) * lower_than_expected,axis=-1) / np.sum((counts > 0), axis=-1)
    scores[0] = 0  # never predict chm13
    print(scores)
    return np.argmax(scores)


def predict_new(counts, priv_counts):
    total_counts = counts.sum(axis=0)
    total_dim, priv_dim = ((int(total_counts.max()+1), (int(priv_counts.max()+1))))
    heatmap = np.zeros((total_dim, priv_dim), dtype=int)
    for count in range(3, total_dim):
        local_counts = priv_counts[total_counts == count]
        hist = np.bincount(local_counts, minlength=int(priv_dim))
        plt.plot(hist[1:])
        plt.plot(poisson.pmf(np.arange(1, priv_dim), count*COVERAGE/counts.shape[0])*(local_counts != 0).sum())
        plt.show()
        heatmap[count, :] = hist
    plt.imshow(heatmap[:, 1:])
    plt.show()
    return 0


def get_static_counts():
    counts = np.load("data/test5000/marker_kmer_counts.npy")
    priv_counts = np.load("data/test5000/priv_1_e0.01_marker_kmer_counts.npy")
    sample_names = {i: line.strip() for i, line in enumerate(open("data/test5000/sample_names.txt").readlines())}
    return sample_names, counts, priv_counts


try:
    snakemake
    func = get_counts
except Exception:
    func = get_static_counts


sample_names, counts, priv_counts = func()
# i = predict_index(counts, priv_counts)
i = predict_new(counts, priv_counts)

#i = predict_index(counts, priv_counts)
#print(i)
predicted_individual = sample_names[i]
# print(predicted_individual)

try:
    snakemake
    with open(snakemake.output.prediction, "w") as f:
        f.write("%s\n" % predicted_individual)
except:
    pass
