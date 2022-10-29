

rule main:
    output:
        "data/test/1.og",


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
        sequences="data/{dataset}/{i}.fa.gz"
    input:
        random_individuals="data/{dataset}/random_individuals.txt",
        paths="data/{dataset}/paths.txt",
        graph="data/{dataset}/graph.og"
    conda:
        "envs/odgi.yml"
    shell:
        "INDIVIDUAL=$(head -n {wildcards.i} {input.random_individuals} | tail -n 1) && "
        "target_haplotype_length=200 && "
        """
        cat {input.paths} | grep -v "^$INDIVIDUAL#" > data/{wildcards.dataset}/keep.{wildcards.i}
        odgi paths -i {input.graph} -K data/{wildcards.dataset}/keep.{wildcards.i} -o - \
          | odgi priv -i - -d 30 -e 0.01 -c 3 -b $target_haplotype_length -t 16 -P -o {output.graph}
        odgi paths -i {wildcards.i}.og -f | bgzip -@ 48 > {output.sequences}
        """
