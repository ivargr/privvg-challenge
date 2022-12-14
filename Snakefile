import bionumpy as bnp
import numpy as np
import random
from shared_memory_wrapper import from_file, to_file
import npstructures as nps
np.random.seed(1)
random.seed(1)
import matplotlib.pyplot as plt

configfile: "config.yaml"

def get_sample_names(dataset):
    if "test" in dataset:
        return config["test_individual_names"].split(",")
    else:
        assert dataset == "real"
        return config["real_data_individual_names"].split(",")


    with open("data/" + dataset + "/sample_names.csv") as f:
        names = f.read().strip().split(",")
        print(names)
        return names


def get_path_length_for_odgi(wildcards):
    if wildcards.dataset == "real":
        return 10000
    else:
        return 500 # too high can make odgi priv run forever on the small test graphs


rule all:
    input:
        "data/test5000/sample_names.txt"
    
    
rule make_test_graph:
    output:
        "data/test{n_variants,\d+}/graph.gfa",
        "data/test{n_variants,\d+}/sample_names.txt",
        "data/test{n_variants,\d+}/sample_names.csv"
    run:
        from simulation import simulate_gfa
        simulate_gfa(int(wildcards.n_variants), 44, "data/test" + wildcards.n_variants + "/")


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


rule get_graph_paths_as_fasta:
    output:
        "data/{dataset}/graph.fa.gz"
    input:
        "data/{dataset}/graph.og"
    conda:
        "envs/odgi.yml"
    shell:
        """
        odgi paths -i {input} -f | bgzip -@ 12 > {output} && samtools faidx {output} 
        """
    

rule get_sample_sequence:
    output:
        "data/{dataset}/sample_{name}.fa.gz"
    input:
        graph_sequences="data/{dataset}/graph.fa.gz",
    conda:
        "envs/samtools.yml"
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
    #return ["data/" + wildcards.dataset + "/" + sample_name + ".unique_kmers.npy" for sample_name in config["datasets"][wildcards.dataset]["individuals"].split(",")]
    return ["data/" + wildcards.dataset + "/" + sample_name + ".unique_kmers.npy" for sample_name in get_sample_names(wildcards.dataset)]

def sample_kmers_input(wildcards):
    #return ["data/" + wildcards.dataset + "/" + sample_name + ".kmers.npy" for sample_name in config["datasets"][wildcards.dataset]["individuals"].split(",")]
    return ["data/" + wildcards.dataset + "/" + sample_name + ".kmers.npy" for sample_name in get_sample_names(wildcards.dataset)]


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
        markers = kmers[(counter[kmers] < 8) & (counter[kmers] > 2)]
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
    output:
        random_individuals="data/{dataset}/random_individuals.txt"
    run:
        sample_names = get_sample_names(wildcards.dataset)  # config["datasets"][wildcards.dataset]["individuals"].split(",")
        sample_names.remove("grch38")
        sample_names.remove("chm13")
        random.shuffle(sample_names)
        print(sample_names)
        with open(output.random_individuals, "w") as f:
            f.write("\n".join(sample_names))


def get_haplotype_length(wildcards):
    return get_path_length_for_odgi(wildcards)


rule make_priv_graph:
    output:
        graph="data/{dataset}/{i,\d+}_e{epsilon}.og",
        sequences="data/{dataset}/{i,\d+}_e{epsilon}.fa.gz",
        gfa="data/{dataset}/{i,\d+}_e{epsilon}.gfa"
    input:
        random_individuals="data/{dataset}/random_individuals.txt",
        paths="data/{dataset}/paths.txt",
        graph="data/{dataset}/graph.og",
    conda:
        "envs/odgi.yml"
    resources:
        mem_gb=15
    threads:
        2
    params:
        haplotype_length=get_haplotype_length
    shell:
        "INDIVIDUAL=$(head -n {wildcards.i} {input.random_individuals} | tail -n 1) && "
        "target_haplotype_length={params.haplotype_length} && "
        """
        cat {input.paths} | grep -v "^$INDIVIDUAL#" > data/{wildcards.dataset}/keep.{wildcards.i}
        odgi paths -i {input.graph} -K data/{wildcards.dataset}/keep.{wildcards.i} -o - \
          | odgi priv -i - -d 30 -e {wildcards.epsilon} -c 3 -b $target_haplotype_length -t 2 -P -o {output.graph}
        odgi paths -i {output.graph} -f | bgzip -@ 48 > {output.sequences}
        odgi view -i {output.graph} -g  > {output.gfa}
        """



rule count_marker_kmers_in_priv_graph:
    output:
        counts="data/{dataset}/priv_{i,\w+}_e{epsilon}_marker_kmer_counts.npy"
    input:
        sequences="data/{dataset}/{i,\w+}_e{epsilon}.fa.gz",
        marker_kmers="data/{dataset}/marker_kmers.npy"
    run:
        marker_kmers = np.load(input.marker_kmers)
        counter = nps.Counter(marker_kmers, mod=20000033)

        for chunk in bnp.open(input.sequences).read_chunks(100000000):
            priv_sequences = chunk.sequence
            print("Read priv sequences")
            priv_kmers = bnp.get_kmers(bnp.as_encoded_array(priv_sequences, bnp.DNAEncoding), 31).raw().ravel()
            print("%d kmers processed" % len(priv_kmers))
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


rule predict_with_kmers:
    input:
        counts="data/{dataset}/marker_kmer_counts.npy",
        priv_counts="data/{dataset}/priv_{i}_e{epsilon}_marker_kmer_counts.npy",
        sample_names="data/{dataset}/sample_names.txt"
    output:
        prediction="data/{dataset}/prediction_{i,\w+}_e{epsilon}_with_kmers.txt",
        prediction_scores="data/{dataset}/prediction_{i,\w+}_e{epsilon}_with_kmers_scores.npy",
    script:
        "predict.py"



def all_predictions(wildcards):
    n_individuals = len(get_sample_names(wildcards.dataset))-2
    if wildcards.dataset == "real":
        n_individuals = 30
    #n_individuals = 20
    return ["data/" + wildcards.dataset + "/prediction_" + str(i) + "_e" + wildcards.epsilon + "_with_kmers.txt"
            for i in range(1, n_individuals+1)]


rule predict_all:
    input:
        "data/{dataset}/random_individuals.txt",
        all_predictions
    output:
        "data/{dataset}/all_predictions_e{epsilon}.txt",
        "data/{dataset}/accuracy_e{epsilon}.txt",
    run:
        truth = [line.strip() for line in open(input[0]).readlines()]
        predicted = [open(name).read().strip() for name in input[1:]]

        print(truth, len(truth))

        with open(output[0], "w") as f:
            n_correct = 0
            n_tot = 0
            for truth_individual, predicted in zip(truth, predicted):
                n_tot += 1
                if truth_individual == predicted:
                    n_correct += 1
                f.write("%s,%s\n" % (truth_individual, predicted))

        accuracy = n_correct / n_tot
        print("N correct: %d/%d" % (n_correct, n_tot))

        with open(output[1], "w") as f:
            f.write("%.5f\n"  % accuracy)


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
            #print()
            #print(nodes)
            #print(line.split()[5])
            all_nodes.extend(nodes)

        np.save(output.counts, np.bincount(np.array(all_nodes).astype(np.int64)))


rule get_challenge_dataset:
    output:
        "data/real/challenge{i,\d+}_e0.01.fa.gz"
    shell:
        "wget -O {output} --show-progress https://f004.backblazeb2.com/file/pangenome/privvg/{wildcards.i}.fa.gz"


rule challenge_solution:
    input:
        "data/real/prediction_challenge0_e0.01_with_kmers.txt",
        "data/real/prediction_challenge1_e0.01_with_kmers.txt",
        "data/real/prediction_challenge2_e0.01_with_kmers.txt",
        "data/real/prediction_challenge3_e0.01_with_kmers.txt",
        "data/real/prediction_challenge4_e0.01_with_kmers.txt",
        "data/real/prediction_challenge5_e0.01_with_kmers.txt",
    output:
        solution="challenge_solution.txt"
    shell:
        "cat {input} > {output}"


rule plot_accuracy_across_n_variants:
    input:
        expand("data/test{n_variants}/accuracy_e{{epsilon}}.txt", n_variants=[500, 1000, 2500, 5000, 7500, 10000, 12500, 15000, 20000, 30000])
    output:
        "plot_across_n_variants_e{epsilon}.png"
    run:
        n_variants = []
        accuracies = []
        for file_name in input:
            accuracies.append(float(open(file_name).read().strip()))
            n_variants.append(int(file_name.split("/")[1].replace("test", "")))

        plt.plot(n_variants, accuracies)
        plt.ylabel("Accuracy")
        plt.xlabel("Number of variants")
        plt.title("Accuracy of predicting which sample has been removed, epsilon = " + wildcards.epsilon)
        plt.savefig(output[0])


rule plot_accuracy_across_epsilons:
    input:
        expand("data/test{{n_variants}}/accuracy_e{epsilon}.txt", epsilon=[0.000001, 0.00001, 0.0001, 0.001, 0.01, 0.05, 0.1, 0.2, 0.5, 1, 5, 10])
        #expand("data/test{{n_variants}}/accuracy_e{epsilon}.txt", epsilon=[0.0001, 0.001, 0.01, 0.05, 0.1, 0.25, 0.5, 1])
    output:
        "plot_across_epsilon_{n_variants}variants.png"
    run:
        epsilons = []
        accuracies = []
        for file_name in input:
            accuracies.append(float(open(file_name).read().strip()))
            epsilons.append(float(file_name.split("accuracy_e")[1].replace(".txt", "")))

        plt.plot(epsilons, accuracies)
        print(epsilons, accuracies)
        plt.title("Accuracy on varying epsilons, number of variants = " + wildcards.n_variants)
        plt.xlabel("Epsilon (log-scale)")
        plt.ylabel("Accuracy")
        plt.xscale("log")
        plt.savefig(output[0])

