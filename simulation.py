import random
random.seed(1)
import numpy as np
np.random.seed(1)
from collections import defaultdict
import sys

nucleotides = ["A", "C", "T", "G"]

def random_sequence(n):
    return ''.join(random.choices(nucleotides, k=n))


def simulate_gfa(n_variants=100, n_individuals=44, out_file_name="graph.gfa"):
    # very simple:
    # simulate n SNPs with random sequences inbetween
    # spread individuals out with two haplotypes
    individual_ids = list(str(i) for i in range(n_individuals))
    path_ids = ["chm13#0", "grch38#0"] + \
        ["individual" + id + "#1" for id in individual_ids] + \
        ["individual" + id + "#2" for id in individual_ids]

    nodes = []
    links = []
    paths = defaultdict(list)
    marker_nodes = []

    allele_frequencies = []

    node_id = 1
    for variant in range(n_variants):
        # add a linear ref node
        # add two variant nodes
        nodes.append((node_id, random_sequence(50)))
        # all individuals will have this node
        for path_id in path_ids:
            paths[path_id].append(node_id)

        variant_nodes = [node_id+1, node_id+2]

        for variant_node in variant_nodes:
            links.append((node_id, variant_node))
            links.append((variant_node, node_id+3))
            nodes.append((variant_node, random_sequence(5)))

        marker_nodes.extend(variant_nodes)

        # some individuals should have each node, should be a bit skewed
        n_on_variant = np.random.randint(1, len(path_ids)//4)
        allele_frequencies.append(n_on_variant / len(path_ids))
        paths_on_variant = set(random.sample(path_ids, n_on_variant))
        for path_id in path_ids:
            if path_id in paths_on_variant:
                paths[path_id].append(variant_nodes[1])
            else:
                paths[path_id].append(variant_nodes[0])

        node_id += 3

    # add a final ref node
    nodes.append((node_id, random_sequence(30)))
    for path_id in path_ids:
        paths[path_id].append(node_id)

    with open(out_file_name, "w") as f:
        for node in nodes:
            f.write("S\t%d\t%s\n" % node)

        for link in links:
            f.write("L\t%d\t+\t%d\t+\t+0M\n" % link)

        for path_id, path in paths.items():
            path = ','.join([str(node) + "+" for node in path])
            f.write("P\t%s\t%s\n" % (path_id, path))

    print("Sample names:")
    print(",".join(["chm13", "grch38"] + ["individual" + str(i) for i in range(n_individuals)]))

    print(allele_frequencies)
    np.save("marker_nodes.npy", np.array(marker_nodes))

if __name__ == "__main__":
    simulate_gfa(5000, 30, sys.argv[1])

