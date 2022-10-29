import gzip
import bionumpy as bnp
import numpy as np
import npstructures as nps
from shared_memory_wrapper import from_file, to_file


def get_sample_kmers(sample_id):
    f = bnp.open("samples/" + sample_id + ".fa.gz")
    chunk = f.read()
    kmers = bnp.kmers.fast_hash(bnp.as_encoded_array(chunk.sequence, bnp.DNAEncoding), 31).ravel()
    print(kmers)
    return kmers


def get_all_kmers(sample_names):

    for sample_name in sample_names:
        print(sample_name)
        kmers = get_sample_kmers(sample_name)
        np.save("samples/" + sample_name + ".kmers", kmers)


def merge_all_kmers(sample_names, out_file="all_kmers"):
    kmers = (np.load("samples/" + sample_name + ".kmers.npy") for sample_name in sample_names)
    all = []
    for k in kmers:
        print("Got %d kmers" % (len(k)))
        all.append(k)

    all = np.concatenate(all)
    np.save(out_file, all)


def find_unique_kmers(sample_names):
    # ignore kmers on linear ref
    linear_ref_kmers = np.load("samples/chm13.kmers.npy")
    print("Making counter")
    unique_linear_ref_kmers = nps.HashSet(np.unique(linear_ref_kmers), mod=200000033)

    for sample_name in sample_names:
        kmers = np.load("samples/" + sample_name + ".kmers.npy")
        unique, counts = np.unique(kmers, return_counts=True)
        unique_kmers = unique[counts == 1]
        print("%s has %d kmers with frequency 1" % (sample_name, len(unique_kmers)))
        # only keep those that are not on linear ref
        #if sample_name != "chm13":
        unique_kmers = unique_kmers[unique_linear_ref_kmers.contains(unique_kmers) == False]
        print("N kmers left after filtering with linear ref: %d" % len(unique_kmers))
        np.save("samples/" + sample_name + ".unique_kmers", unique_kmers)


def merge_all_unique_sample_kmers(sample_names):
    all = []
    for i, sample_name in enumerate(sample_names):
        print(i)
        all.append(np.load("samples/" + sample_name + ".unique_kmers.npy"))

    all = np.concatenate(all)
    print("%d kmers in total" % len(all))
    np.save("all_unique_kmers.npy", all)


def merge_all_unique_sample_kmers2(sample_names):
    current = None
    for i, sample_name in enumerate(sample_names):
        print(i)
        new_kmers = np.load("samples/" + sample_name + ".unique_kmers.npy")
        if current is None:
            current = new_kmers
        else:
            print("Now %d kmers in current" % len(current))
            # add kmers that are not already added
            current_set = nps.HashSet(current, mod=200000033)
            new_to_add = new_kmers[current_set.contains(new_kmers) == False]
            print("Found %d new kmers to add" % len(new_to_add))
            current = np.concatenate([current, new_to_add])

    np.save("all_unique_kmers_v2.npy", current)


def create_marker_kmers(max_frequency=6):
    print("Loading kmers")
    all = np.load("all_unique_kmers.npy")
    # find kmers with low frequency (that few individuals have)
    print("Finding unique kmers")
    unique, counts = np.unique(all, return_counts=True)
    print("N unique kmers: %d" % len(unique))
    keep = unique[counts <= max_frequency]
    keep_counts = counts[counts <= max_frequency]
    print("N kmers with frequency lower than %d: %d" % (max_frequency, len(keep)))
    kmers_with_counts = nps.HashTable(keep, keep_counts, mod=200000033)
    print("Made hash table")
    to_file(kmers_with_counts, "marker_kmers_with_counts")


def find_marker_kmers_in_individuals(sample_names):
    # for each individual, find only the kmers that are also marker kmers
    # these kmers will be used when checking whether the individual may have been removed
    marker_kmers = np.load("marker_kmers.npy")
    marker_kmers = nps.HashSet(marker_kmers, mod=200000033)


def find_unique_node_kmers(gfa_file_name, all_graph_kmers):
    # for every node, try to find a kmer that is unique in the rest of the graph
    f = bnp.open(gfa_file_name)
    n_skipped = 0
    n_nodes = 0
    for i, chunk in enumerate(f):
        print(i)
        print("N skipped: %d" % n_skipped)
        for node_sequence in chunk.sequence:
            if n_nodes % 1000 == 0:
                print("%d nodes processed" % n_nodes)

            if len(node_sequence) > 100:
                continue

            n_nodes += 1
            if np.sum(node_sequence == "N") > 0:
                n_skipped += 1
                continue

            kmers = bnp.kmers.fast_hash(
                bnp.as_encoded_array(node_sequence, bnp.DNAEncoding),
                window_size=31)


def count_all_unique_kmers(sample_names):
    all_unique = np.load("all_unique_kmers_v2.npy")
    print("Making counter")
    counter = nps.Counter(all_unique, mod=200000033)
    for i, sample_name in enumerate(sample_names):
        print(i, sample_name)
        sample_kmers = np.load("samples/" + sample_name + ".kmers.npy")
        counter.count(sample_kmers)

    to_file(counter, "kmer_counts")


if __name__ == "__main__":
    sample_names = [l.strip() for l in open("sample.names")]
    #get_sample_kmers("test50k")


    #get_all_kmers()
    #find_unique_kmers(sample_names)
    # all kmers that any individual can have:
    #merge_all_unique_sample_kmers2(sample_names)

    # count how many times each kmer occurs in any individual
    # we will not count kmers occuring on lineawr ref since these are not in counte
    count_all_unique_kmers(sample_names)