import bionumpy as bnp
import numpy as np
from shared_memory_wrapper import from_file, to_file
import npstructures as nps

configfile: "config.yaml"


rule all:
    input:
        "data/test/node_counts.npy",
        "data/test/priv_1_node_counts.npy",
        #"data/test/priv_1_marker_kmer_counts.npy",
        #"data/test/priv_2_marker_kmer_counts.npy",
        #"data/test/marker_kmer_counts.npy",
        "data/test/sample_names.txt"


rule make_odgi_graph:
    input:
        "data/{dataset}/graph.gfa"
    output:
        "data/{dataset}/graph.og"
    conda:
        "envs/odgi.yml"
    shell:
        "odgi build -g {input} -o {output}"


rule get_paths:
    input:
        "data/{dataset}/graph.og"
    output:
        "data/{dataset}/paths.txt"
    conda:
        "envs/odgi.yml"
    shell:
        "odgi paths -i {input} -L > {output}"


rule get_sample_names:
    output:
        "data/{dataset}/sample_names.txt"
    input:
        "data/{dataset}/paths.txt"
    shell:
        "cut -d# -f 1 {input} | sort | uniq > {output}"


rule get_graph_paths_as_fasta:
    output:
        "data/{dataset}/graph.fa.gz"
    input:
        "data/{dataset}/graph.og"
    conda:
        "envs/odgi.yml"
    shell:
        """
        odgi paths -i {input} -f | bgzip -@ 24 > {output} && samtools faidx {output} 
        """
    

rule get_sample_sequence:
    output:
        "data/{dataset}/sample_{name}.fa.gz"
    input:
        graph_sequences="data/{dataset}/graph.fa.gz",
    shell:
        "samtools faidx {input.graph_sequences} $(grep {wildcards.name}'#' {input.graph_sequences}.fai | cut -f 1)  "
        "| bgzip -@ 24 > {output} && samtools faidx {output}"


rule get_sample_kmers:
    output:
        "data/{dataset}/{sample_name}.kmers.npy"
    input:
        "data/{dataset}/sample_{sample_name}.fa.gz"
    shell:
        "python3 scripts.py get-sample-kmers {input} {output}"


rule get_unique_sample_kmers:
    output:
        "data/{dataset}/{sample_name}.unique_kmers.npy"
    input:
        sample_kmers="data/{dataset}/{sample_name}.kmers.npy",
        linear_ref_kmers="data/{dataset}/chm13.kmers.npy"
    shell:
        "python3 scripts.py get-unique-kmers-not-on-linear-ref {input.sample_kmers} {input.linear_ref_kmers} {output}"


def merge_all_unique_kmers_input(wildcards):
    return ["data/" + wildcards.dataset + "/" + sample_name + ".unique_kmers.npy" for sample_name in config["datasets"][wildcards.dataset]["individuals"].split(",")]

def sample_kmers_input(wildcards):
    return ["data/" + wildcards.dataset + "/" + sample_name + ".kmers.npy" for sample_name in config["datasets"][wildcards.dataset]["individuals"].split(",")]


rule merge_all_unique_kmers:
    output:
        "data/{dataset}/all_unique_sample_kmers.npy"
    input:
        kmers=merge_all_unique_kmers_input
    shell:
        "python3 scripts.py merge-all-unique-sample-kmers {output} {input}"


rule count_all_kmers:
    output:
        "data/{dataset}/kmer_counts.npz"
    input:
        "data/{dataset}/all_unique_sample_kmers.npy",
        sample_kmers_input
    shell:
        "python3 scripts.py count-all-unique-kmers {output} {input}"


rule get_marker_kmers:
    output:
        marker_kmers="data/{dataset}/marker_kmers.npy"
    input:
        counter="data/{dataset}/kmer_counts.npz",
        sample_kmers="data/{dataset}/all_unique_sample_kmers.npy"
    run:
        counter = from_file(input.counter)
        kmers = np.load(input.sample_kmers)
        markers = kmers[(counter[kmers] < 15) & (counter[kmers] > 3)]
        print("Found ", len(markers), "marker kmers")
        np.save(output.marker_kmers, markers)


rule get_sample_counts_for_marker_kmers:
    input:
        "data/{dataset}/marker_kmers.npy",
        sample_kmers_input
    output:
        counts="data/{dataset}/marker_kmer_counts.npy"
    run:
        marker_kmers = np.load(input[0])
        counter = nps.Counter(marker_kmers)

        individuals = input[1:]
        output_matrix = np.zeros((len(individuals), len(marker_kmers)))
        for i, individual in enumerate(individuals):
            print("Counting kmers in %s" % individual)
            individual_kmers = np.load(individual)
            counter.count(individual_kmers)
            output_matrix[i,:] = counter[marker_kmers]
            counter.fill(0)  # reset before counting next individual

        np.save(output.counts, output_matrix)




rule get_random_individuals_to_be_removed:
    input:
        "data/{dataset}/paths.txt"
    output:
        "data/{dataset}/random_individuals.txt"
    shell:
        """
        cut -d "#" -f 1 {input} | sort | uniq | grep -v "grch38" | grep -v "chm13" | shuf | head -n 6 > {output}
        """


rule make_priv_graph:
    output:
        graph="data/{dataset}/{i,\d+}.og",
        sequences="data/{dataset}/{i,\d+}.fa.gz",
        gfa="data/{dataset}/{i,\d+}.gfa"
    input:
        random_individuals="data/{dataset}/random_individuals.txt",
        paths="data/{dataset}/paths.txt",
        graph="data/{dataset}/graph.og",
    conda:
        "envs/odgi.yml"
    shell:
        "INDIVIDUAL=$(head -n {wildcards.i} {input.random_individuals} | tail -n 1) && "
        "target_haplotype_length=500 && "
        """
        cat {input.paths} | grep -v "^$INDIVIDUAL#" > data/{wildcards.dataset}/keep.{wildcards.i}
        odgi paths -i {input.graph} -K data/{wildcards.dataset}/keep.{wildcards.i} -o - \
          | odgi priv -i - -d 30 -e 0.01 -c 3 -b $target_haplotype_length -t 16 -P -o {output.graph}
        odgi paths -i {output.graph} -f | bgzip -@ 48 > {output.sequences}
        odgi view -i {output.graph} -g  > {output.gfa}
        """



rule count_marker_kmers_in_priv_graph:
    output:
        counts="data/{dataset}/priv_{i}_marker_kmer_counts.npy"
    input:
        sequences="data/{dataset}/{i,\d+}.fa.gz",
        marker_kmers="data/{dataset}/marker_kmers.npy"
    run:
        marker_kmers = np.load(input.marker_kmers)
        counter = nps.Counter(marker_kmers)

        priv_sequences = bnp.open(input.sequences).read().sequence
        print("Read priv sequences")
        priv_kmers = bnp.kmers.fast_hash(bnp.as_encoded_array(priv_sequences, bnp.DNAEncoding), window_size=31).ravel()
        print("%d kmers in priv graph" % len(priv_kmers))

        counter.count(priv_kmers)
        counts = counter[marker_kmers]

        np.save(output.counts, counts)


rule get_node_counts:
    input:
        gfa="data/{dataset}/graph.gfa",
        sample_names="data/{dataset}/sample_names.txt"
    output:
        counts="data/{dataset}/node_counts.npy"
    run:
        max_node_id = 20000
        sample_names = {line.strip(): i for i, line in enumerate(open(input.sample_names).readlines())}
        out_counts = np.zeros((len(sample_names), max_node_id))

        for line in open(input.gfa):
            if not line.startswith("P"):
                continue

            sample_name = line.split()[1].split("#")[0]
            name_id = sample_names[sample_name]
            nodes = np.array([int(node[:-1]) for node in line.split()[2].split(",")])
            out_counts[name_id,nodes] += 1

        np.save(output.counts, out_counts)


rule get_priv_node_counts:
    input:
        graph="data/{dataset}/{i,\d+}.gfa"
    output:
        counts="data/{dataset}/priv_{i}_node_counts.npy"
    run:
        nodes = []
        f = open(input.graph)
        for line in f:
            if not line.startswith("P"):
                continue

            nodes.extend([int(node[:-1]) for node in line.split()[2].split(",")])

        counts = np.bincount(nodes)
        np.save(output.counts, counts)


rule predict:
    input:
        counts="data/{dataset}/node_counts.npy",
        priv_counts="data/{dataset}/priv_{i}_node_counts_from_alignments.npy",
        sample_names="data/{dataset}/sample_names.txt"
    output:
        prediction="data/{dataset}/prediction_{i}.npy"

    run:
        sample_names = {i: line.strip() for i, line in enumerate(open(input.sample_names).readlines())}

        node_counts = np.load(input.counts)
        priv_node_counts = np.load(input.priv_counts)
        sum_node_counts = np.sum(node_counts,axis=0)
        marker_nodes = np.where((sum_node_counts < 5) & (sum_node_counts >= 2))[0]

        print("%d marker nodes" % (len(marker_nodes)))
        print(marker_nodes)

        counts_on_marker_nodes = node_counts[:, marker_nodes]
        lower_than_expected = priv_node_counts[marker_nodes] <= 0
        prediction = np.sum((counts_on_marker_nodes > 0) * lower_than_expected,axis=-1) / np.sum((counts_on_marker_nodes > 0),axis=-1)
        prediction[0] = 0  # never predict chm13
        print(prediction)
        print(np.argmax(prediction))
        print(sample_names[np.argmax(prediction)])

        np.save(output.prediction, np.argmax(prediction))


rule align_priv_paths_to_original_graph:
    output:
        alignments="data/{dataset}/{i,\d+}_alignments.gaf"
    input:
        graph="data/{dataset}/graph.gfa",
        haplotype_sequences="data/{dataset}/{i}.fa.gz"
    conda:
        "envs/graphaligner.yml"
    shell:
        "GraphAligner -g {input.graph} -f {input.haplotype_sequences} -a {output.alignments} -x vg"


rule get_priv_node_counts_from_alignments:
    input:
        alignments="data/{dataset}/{i}_alignments.gaf"
    output:
        counts="data/{dataset}/priv_{i}_node_counts_from_alignments.npy"
    run:
        import re
        all_nodes = []
        max_node_id = config["datasets"][wildcards.dataset]["max_node_id"]

        for line in open(input.alignments):
            nodes = line.split()[5]
            regex = r"(<|>)(\d+)"
            matches = re.finditer(regex, nodes)
            nodes = [m.groups() for m in matches]
            nodes = [int(node) if dir == ">" else int(node)+max_node_id for dir, node in nodes]
            print()
            print(nodes)
            print(line.split()[5])
            all_nodes.extend(nodes)

        np.save(output.counts, np.bincount(np.array(all_nodes).astype(np.int64)))