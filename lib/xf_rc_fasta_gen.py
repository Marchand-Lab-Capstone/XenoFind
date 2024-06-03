import os
import sys
from Bio import SeqIO
from Bio.Seq import Seq

def reverse_complement_fasta(input_fasta, output_fasta):
    # Read the input fasta file and create the reverse complement
    records = []
    for record in SeqIO.parse(input_fasta, "fasta"):
        record.seq = record.seq.reverse_complement()
        record.id += "+RC"
        record.description = record.id
        records.append(record)
    
    # Write the reverse complemented sequences to the output file with sequences in one line
    with open(output_fasta, 'w') as output_handle:
        for record in records:
            output_handle.write(f">{record.id}\n{str(record.seq)}\n")

def main():
    input_fasta = sys.argv[1]
    output_fasta = sys.argv[2]

    # Check if the input fasta file exists
    if not os.path.isfile(input_fasta):
        print(f"Error: Input fasta file '{input_fasta}' not found.")
        sys.exit(1)

    # Copy the input fasta file to the output path
    output_copy = os.path.join(os.path.dirname(output_fasta), os.path.basename(input_fasta))
    os.system(f"cp {input_fasta} {output_copy}")
    
    # Create the reverse complement fasta file
    rev_comp_output_fasta = output_fasta.replace('.fa', '_rc.fa')
    reverse_complement_fasta(input_fasta, rev_comp_output_fasta)
    
    print(f"Reverse complemented fasta file saved as '{rev_comp_output_fasta}'")

if __name__ == "__main__":
    main()

