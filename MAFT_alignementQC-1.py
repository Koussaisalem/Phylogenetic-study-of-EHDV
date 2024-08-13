# quality_control.py
from __future__ import print_function
import sys
print("Python executable: {}".format(sys.executable))

from Bio import SeqIO, AlignIO
from Bio.Align.Applications import MafftCommandline
from Bio.Seq import Seq
import subprocess
import os

def initial_quality_check(input_file, min_length=100):
    """
    Perform initial quality checks on sequences.
    """
    good_sequences = []
    total_sequences = 0
    for record in SeqIO.parse(input_file, "fasta"):
        total_sequences += 1
        seq = str(record.seq).upper()
        if len(seq) >= min_length and 'N' not in seq:
            good_sequences.append(record)
    
    print("Total sequences: {}".format(total_sequences))
    print("Sequences passing initial QC: {}".format(len(good_sequences)))
    return good_sequences

def run_mafft(input_file, output_file):
    """
    Run MAFFT alignment.
    """
    mafft_path = "/home/koussai/miniconda3/envs/stage/bin/mafft"
    print("MAFFT path: {}".format(mafft_path))
    print("Input file exists: {}".format(os.path.exists(input_file)))
    print("MAFFT executable exists: {}".format(os.path.exists(mafft_path)))
    
    try:
        result = subprocess.run([mafft_path, input_file], stdout=subprocess.PIPE, stderr=subprocess.PIPE, universal_newlines=True)
        with open(output_file, "w") as handle:
            handle.write(result.stdout)
        print("MAFFT stderr: {}".format(result.stderr))
        
        # Count sequences in MAFFT output
        mafft_sequences = list(SeqIO.parse(output_file, "fasta"))
        print("Sequences after MAFFT: {}".format(len(mafft_sequences)))
    except Exception as e:
        print("Error running MAFFT: {}".format(str(e)))

def post_alignment_qc(aligned_file, gap_threshold=0.7, coverage_threshold=0.3):
    """
    Perform quality control on aligned sequences.
    """
    alignment = AlignIO.read(aligned_file, "fasta")
    good_sequences = []
    
    print("Total sequences in alignment: {}".format(len(alignment)))
    print("Alignment length: {}".format(alignment.get_alignment_length()))
    
    for record in alignment:
        seq = str(record.seq)
        gap_ratio = seq.count('-') / len(seq)
        non_gap_length = len(seq) - seq.count('-')
        coverage = non_gap_length / len(alignment[0])
        
        print("Sequence {}: gap_ratio={:.2f}, coverage={:.2f}".format(record.id, gap_ratio, coverage))
        
        if gap_ratio <= gap_threshold and coverage >= coverage_threshold:
            good_sequences.append(record)
        else:
            print("  Failed QC")
    
    print("Sequences passing post-alignment QC: {}".format(len(good_sequences)))
    return good_sequences


def main(input_file, output_file, gap_threshold=0.7, coverage_threshold=0.3):

    # Initial QC
    print("Performing initial quality check on ".format(input_file))
    good_seqs = initial_quality_check(input_file)
    SeqIO.write(good_seqs, "initial_qc.fasta", "fasta")

    # Check if initial_qc.fasta is created and not empty
    if os.path.exists("initial_qc.fasta") and os.path.getsize("initial_qc.fasta") > 0:
        print("initial_qc.fasta created successfully")
    else:
        print("Error: initial_qc.fasta is empty or not created")
        return
    
    # MAFFT alignment
    print("Running MAFFT alignment")
    run_mafft("initial_qc.fasta", "aligned.fasta")
    # Check mafft output
    mafft_sequences = list(SeqIO.parse("aligned.fasta", "fasta"))
    print("Sequences in MAFFT output: {}".format(len(mafft_sequences)))

    # Post-alignment QC
    print("Performing post-alignment quality control")
    final_seqs = post_alignment_qc("aligned.fasta", gap_threshold, coverage_threshold)
    
    # Write final sequences
    SeqIO.write(final_seqs, output_file, "fasta")
    print("Final curated sequences written to".format (output_file))

if __name__ == "__main__":
    main("VP5.fasta", "vp5_curated.fasta", gap_threshold=0.9, coverage_threshold=0.1)
    main("VP2.fasta", "vp2_curated.fasta", gap_threshold=0.9, coverage_threshold=0.1)

