import gzip
import bionumpy as bnp
import numpy as np
import npstructures as nps
from shared_memory_wrapper import from_file, to_file
from typing import List
import typer
app = typer.Typer()

HASHTABLE_MODULO = 160000003


def get_sample_names():
    with open("sample.names") as f:
        return [l.strip() for l in f]


def get_kmers(sequence_chunk):
    sequences = sequence_chunk.sequence
    # hack, set all Ns to 0 to be able to hash
    sequences._data[sequences._data == "N"] = "A"
    kmers = bnp.kmers.fast_hash(bnp.as_encoded_array(sequences, bnp.DNAEncoding), 31).ravel()
    return kmers


@app.command()
def get_sample_kmers(input_fasta, output):
    f = bnp.open(input_fasta)
    np.save(output, get_kmers(f.read()))


@app.command()
def get_unique_kmers_not_on_linear_ref(input_file, linear_ref_kmer_file, output_file):
    linear_ref_kmers = np.load(linear_ref_kmer_file)
    print("Making counter")
    unique_linear_ref_kmers = nps.HashSet(np.unique(linear_ref_kmers), mod=HASHTABLE_MODULO)
    kmers = np.load(input_file)
    unique, counts = np.unique(kmers, return_counts=True)
    unique_kmers = unique[counts == 1]
    print("%s has %d kmers with frequency 1" % (input_file, len(unique_kmers)))
    unique_kmers = unique_kmers[unique_linear_ref_kmers.contains(unique_kmers) == False]
    print("N kmers left after filtering with linear ref: %d" % len(unique_kmers))
    np.save(output_file, unique_kmers)


def get_all_kmers():
    for sample_name in get_sample_names():
        print(sample_name)
        kmers = get_sample_kmers(sample_name)
        np.save(sample_name + ".kmers", kmers)


@app.command()
def merge_all_unique_sample_kmers(output_file: str, input_files: List[str]):
    current = None
    for i, sample_file_name in enumerate(input_files):
        print(i, sample_file_name)
        new_kmers = np.load(sample_file_name)
        if current is None:
            current = new_kmers
        else:
            print("Now %d kmers in current" % len(current))
            # add kmers that are not already added
            current_set = nps.HashSet(current, mod=HASHTABLE_MODULO)
            new_to_add = new_kmers[current_set.contains(new_kmers) == False]
            print("Found %d new kmers to add" % len(new_to_add))
            current = np.concatenate([current, new_to_add])

    np.save(output_file, current)
    print("Saved to %s" % output_file)


@app.command()
def count_all_unique_kmers(out_file_name, unique_kmers, sample_kmer_files: List[str]):
    all_unique = np.load(unique_kmers)
    print("Making counter")
    counter = nps.Counter(all_unique, mod=HASHTABLE_MODULO)
    for i, sample_kmers in enumerate(sample_kmer_files):
        print(i, sample_kmers)
        sample_kmers = np.load(sample_kmers)
        counter.count(sample_kmers)

    to_file(counter, out_file_name)


if __name__ == "__main__":
    app()