# Human Chromostome 6 privacy challenge, attempt at solution

This is an attempt at solving the [Human chromosome 6 privacy callenge](https://privvg.github.io/2022/09/30/human-chromosome-6-privacy-challenge.html).

We believe the following is the solution to the challenge:

```
0.fa.gz: 
1.fa.gz:
...
```
## Overview of our solution
* Assuming we have an original non-private graph `A` and an epsilon-differentially private graph `B` which is created by removing an unknown individual `i` from `A`. Both graphs contain paths (in the challenge we only observe the sequences of paths in B). 
* We hypothesise that if something is wrong with the privvg-implementation, then for nodes covered by the individual `i` in `A`, we would observe lower than "expected" coverage of paths over these nodes in `B`. What is expected is a bit unclear due to the exponential mechanism, but for simplicity we could look at nodes with zero coverage (no paths) in `B`. If an individual is removed, then for nodes in `A` that the indivdual has, and that few other individual have, we would more often (than if the individual was not removed) see zero coverage on that node in `B`.
* We don't have access to the graph `B`, but we have access to the error-free sampled path sequences. This means that we should be able to quite easily reconstruct `B` by mapping these path sequences to `A`.


### Testing whether our idea could work
We did some simulations to test whether our approach could work. We simulated a simple graph in GFA format and generated differential-private graphs and sampled haplotypes according to the recipy given in the competition instructions.

A problem with mapping reads back to `A` to look at node coverage is that mapping 30x reads for multiple individuals at chromosome 6 will take some time (at least some hours, and maybe even a few days). To be able to test things quickly and play around with parameters, we think the above described approach should also work when looking at kmer frequencies instead of node coverage (a bit similar to the approach of [KAGE](https://github.com/ivargr/kage). I.e., instead of looking for nodes with lower than expected coverage, we look at kmers from `A` with lower than expected coverage in the sequences sampled from `B`. 

We use [BioNumpy](https://github.com/bionumpy/bionumpy) to first scan all kmers in the original graph `A`, and pick out "*marker kmers*", which are kmers that occur few times and that few indiviudals have in `A`. We then count how many times these kmers occur in `B`, and ...

On our simulated data, this approach gave a 100% prediction accuracy with `epsilon=0.01` on a graph with `5000 variants` with random sequences. The accuracy decreased with fewer nodes, but 5000 nodes (which is not much) seemed to be enough to break privacy, so something seems to be wrong with the implementation/concept of privvg. 


We also experimented with different epsilon-values, but were confused on the role of epsilon here (see notes further down). In general, lower epsilon did not make it more difficult to get correct predictions.



### Solving the real case
Before solving the real case, we made a bunch of samples to test on by following the procedure as given in the challenge description, to check the accuracy. Also here, we got a 100% accuracy.


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

## Some speculation and notes
* We have a feeling that the competion creators might be wrong in their assumption about the "composition theorem", i.e. assuming that if there is epsilon-differential privacy at a single variant, then privacy is also preserved when querying multiple variants. We don't know the theory here, but according to e.g. [this paper](https://arxiv.org/abs/1311.0776) it seems that the privacy when performing `k` samples in worst case can be `epsilon * k`, which intuitively also makes somewhat sense. Maybe this is why our approach works, but we're not in any way certain here. 
* We experience that it is more difficult to get correct predictions with larger epsilon. It seems that with large epsilon only major alleles are chosen, which makes sense since the probability-weight of choosing these alleles grows exponentially according to the exponential mechanism. We thus have a feeling that the exponential mechanism does not make so much sense for this problem, since we would rather would want the allele frequencies to be closer to the true allele frequencies with larger epsilon. This could be interesting to discuss further. 
* In the challenge description, the developers say that `[the sequences sampled from the graph are] less likely to expose private information than the graph`. We believe this assumption might be a bit naive, and that an implementation of differential privacy should first and foremost make sure that the graph satifies differential privacy and not rely on increased privacy by sampling sequences. As long as we have the original graph that the differential private graph has been created from (which we have in our case), then in theory the differential private graph can be reconstructed by mapping the reads to this graph. Mapping these reads should be easy since they are sampled without errors. This is true as long as nodes with zero haplotypes are removed from the differential private graph (which we think they should).

## Work by
* Ivar Grytten (ivargry@ifi.uio.no)
* Knut Rand (knutdr@ifi.uio.no)


