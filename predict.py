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
    scores = np.sum(counts * lower_than_expected, axis=-1) / np.sum(counts, axis=-1)
    scores[0] = 0  # never predict chm13
    print(scores)
    return np.argmax(scores), scores


def new_predict_index(counts, priv_counts):
    total_counts = counts.sum(axis=0)
    max_count = total_counts.max()+1
    h = np.array([(priv_counts[total_counts==i]==0).sum()/(total_counts==i).sum() for i in range(max_count)])
    zero_observed = priv_counts == 0
    print(total_counts)
    rate_vector = h[total_counts.astype(int)]
    expected_zeros = np.array([rate_vector[row>0].sum() for row in counts])
    scores = np.sum((counts>0) * zero_observed, axis=-1) / expected_zeros
    scores[[0, 1]] = 0
    print(scores)
    return np.argmax(scores), scores

def estimate_privacy_rate_per_count(counts, priv_counts):
    total_counts = np.sum(counts, axis=0)
    rates = np.zeros(np.max(total_counts)+1)
    for count in range(3, len(rates)):
        lambda_0 = count*COVERAGE/(counts.shape[0]-1)
        h = np.bincount(priv_counts[total_counts==count])
        predicted_zeros = poisson.pmf(0, lambda_0)*(total_counts==count).sum()
        rates[count] = ((priv_counts[total_counts==count]==0).sum()-predicted_zeros)/(total_counts==count).sum()
    plt.plot(rates); plt.show()
    
    return rates


def get_log_pmf(total_counts, priv_counts, privacy_rates, n_individs):
    log_pmf = 0
    total_dim, priv_dim = ((int(total_counts.max()+1), (int(priv_counts.max()+1))))
    for count in range(3, 6):
        local_counts = priv_counts[total_counts == count]
        theoretical = poisson.logpmf(np.arange(priv_dim), count*COVERAGE/(n_individs))+np.log(1-privacy_rates[count])
        theoretical[0] = np.logaddexp(theoretical[0], np.log(privacy_rates[count]))
        #plt.plot(np.exp(theoretical)*local_counts.size)
        #plt.plot(np.bincount(local_counts))
        # plt.show()
        log_pmf += (np.bincount(local_counts, minlength=int(priv_dim))*theoretical).sum()
    return log_pmf


def predict_new(counts, priv_counts):
    counts = counts.astype(int)
    priv_counts = priv_counts.astype(int)
    privacy_rates = estimate_privacy_rate_per_count(counts, priv_counts)
    total_counts = counts.sum(axis=0)
    total_dim, priv_dim = ((int(total_counts.max()+1), (int(priv_counts.max()+1))))
    heatmap = np.zeros((total_dim, priv_dim), dtype=int)
    log_pmfs = [get_log_pmf(total_counts-row, priv_counts, privacy_rates, counts.shape[0]-1)
                for row in counts]
    log_pmfs[0] = -np.inf
    # print(np.sort(log_pmfs))
    return np.argmax(log_pmfs)

    for count in range(3, total_dim):
        local_counts = priv_counts[total_counts == count]
        hist = np.bincount(local_counts, minlength=int(priv_dim))
        plt.plot(hist[1:])
        plt.plot(poisson.pmf(np.arange(1, priv_dim),
                             count*COVERAGE/counts.shape[0])*(local_counts != 0).sum())

        plt.show()
        heatmap[count, :] = hist
        # plt.imshow(heatmap[:, 1:])
        # plt.show()
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
counts = counts.astype(int)
priv_counts = priv_counts.astype(int)
# i = predict_index(counts, priv_counts)
i, scores = new_predict_index(counts, priv_counts)
print(i)
print(sample_names[i])
#i = predict_index(counts, priv_counts)
#print(i)
predicted_individual = sample_names[i]
print("Highest score:", np.max(scores), " at index ", np.argmax(scores))
if np.max(scores) < 1.15 and snakemake.wildcards.dataset == "real":
    print("Highest scores is quite low. Predicting no individual has been removed")
    predicted_individual = "No individual"
# print(predicted_individual)

try:
    snakemake
    with open(snakemake.output.prediction, "w") as f:
        f.write("%s\n" % predicted_individual)

    np.save(snakemake.output.prediction_scores, scores)
except:
    pass
