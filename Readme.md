# Human Chromostome 6 privacy challenge, attempt at solution

This is an attempt at solving the [Human chromosome 6 privacy callenge](https://privvg.github.io/2022/09/30/human-chromosome-6-privacy-challenge.html).

We believe the following is the solution to the challenge:

```
0.fa.gz: 
1.fa.gz:
...
```

## Overview of solution
* We believe the competion creators might be wrong in their assumption about the "composition theorem", i.e. assuming that if there is epsilon-differential privacy at a single variant, then privacy is also preserved when querying multiple variants (more about this in the last section).
* If our assumption is true, then, given enough variants in the graph, there should be enough information leakage to infer wich sample is missing from each provided sample.
* An example of information leakage would be that a node from the original chr 6 graph is supported by fewer than expected haplotypes sampled from a differential private graph. This indicates that haplotypes supporting the node in the chr 6 graph might have been removed. 
* We don't have access to the nodes and hapotype paths in the differential private graph, but since the haplotypes are sampled without errors, we can quite easily map these sequences back to the original chr 6 graph to see which nodes are supported and not supported by haplotypes.

### Testing whether our idea could work
To test whether the above hypothesis is true, we simulated a simple graph in GFA format and generated differential-private graphs and sampled haplotypes according to the recipy given in the competition instructions.

We implemented a very simple method for determining haplotypes from the original graph that seems to be missing from the differential private graph. We do this by mapping the sampled haplotype sequence back to the original graph, and look at nodes with zero support (no reads mapped to them). We predict the missing haplotype by choosing the haplotype with most cases of having a node in the original graph that has zero support when mapping reads from the differential private sequences.

This simple method gives a 100% prediction accuracy with `epsilon=0.01` on a graph with `5000 variants` with random sequences. The accuracy decreased with fewer nodes, and also with lower accuracy, but appearantly `epsilon=0.01` is far too high even for a small graph like this:




### Solving the real case
The real case is slightly more tricky since mapping the reads back to the original graph takes time. To save some time, and to be able to experiment effectively with different settings (such as changing epsilon), we instead looked at kmers in the original graph and in the differential private samples.

The idea is the same as what presented above:
* If a sample is removed, we expect to see fewer kmers from this sample in the differential private data than if it is not removed.
* It is easiest to see such differences if looking at kmers that few individuals have.

We use [BioNumpy](https://github.com/bionumpy/bionumpy) to first scan all kmers in all original samples, and pick out "*marker kmers*", kmers that occur few times and that few indiviudals have.

We then count how many times these kmers occur in each individual in the original graph and in the provided private sample. It would probably be a lot more accurate to map the sequences to the original graph, but just looking at kmers seems to give a 100% accuracy when we simulate differential private graphs from the providedd chromosome 6 graph:



## Reproducing the findings
We have implemented the solution is a small Snakemake pipeline. Dependencies (BioNumpy, odgi, and a few Python packages) should be installed automatically through Conda. 

To run the prediction of the real data, run:

```python
snakemake --use-conda --cores 2 data/real/prediction.txt
```

The results can then be found in the file `data/real/prediction.txt`.

It takes a few hours or so to scan and find marker kmers in the original chromosome 6 graph. The hashtable we use for counting kmers is a bit memory demanding, so one should have around 30GB of memory available. Try running with fewer cores if you run out of memory, and more cores to make it faster if you have enough memory available.


Running on a small simulated graph is very quick and should only take a few minutes:
```python

```



``` 
snakemake --use-conda --cores 8 --resources mem_gb=30 -R count_marker_kmers_in_priv_graph data/real/all_predictions_e0.01.txt
```





* We believe there might be a problem in the implementation/interpretation of the "Composition theorem". According to [the description here](https://privvg.github.io/2022/06/13/Differential-Privacy.html), the privvg-developers say that:
```
It suffice to say that when one makes differentially private queries in succession, the result also preserves privacy.
```
We are no experts on differential privacy, but our intuition tells us that if there is epsilon-differential privacy at each single variant, this means that there should be less privacy and more information when looking at multiple variants. From our understanding, the priacy might even be as bad as `k * epsilon` for k variants ()

* We believe there might be an error/mistake in the implementation in that there is only epsilon-differential privacy at each single variant, but when variants are preserved  

## Running on simulated data
``` 

```


## How to reproduce
###
```python
# get data sets
wget ...

# run the prediction
snakemake --use-conda data/test/all_predictions.txt
```


### A few other notes
* In the challenge description, the developers say that `[the sequences sampled from the graph are] less likely to expose private information than the graph`. We believe this assumption might be a bit naive, and that an implementation of differential privacy should first and foremost make sure that the graph satifies differential privacy and not rely on increased privacy by sampling sequences. As long as we have the original graph that the differential private graph has been created from (which we have in our case), then in theory the differential private graph can be reconstructed by mapping the reads to this graph. Mapping these reads should be easy since they are sampled without errors. This is true as long as nodes with zero haplotypes are removed from the differential private graph (which we think they should).

## Work by
* Ivar Grytten (ivargry@ifi.uio.no)
* Knut Rand (knutdr@ifi.uio.no)


