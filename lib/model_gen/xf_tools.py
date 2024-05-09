import pysam

def filter_primary_alignments(bam_path, out_name): #can edit this function to take in an input file name 
    """
    Filters a BAM file to retain only primary alignments and outputs the result
    to a new BAM file named 'primary_alignments.bam' in the same directory as
    the input file.
    
    Parameters: 
    bam_path: path to the BAM file as a string.
    
    Returns: 
    Path to the filtered BAM file containing only primary alignments.
    """
    
    # Generating the output BAM file path
    directory = os.path.dirname(bam_path)
    output_bam = os.path.join(directory, '{}.bam'.format(out_name))
    
    # Using pysam to filter for primary alignments
    with pysam.AlignmentFile(bam_path, "rb") as infile, \
         pysam.AlignmentFile(output_bam, "wb", header=infile.header) as outfile:

        for read in infile:
            if not read.is_secondary and not read.is_supplementary and not read.is_unmapped:
                outfile.write(read)
                
    print('XenoFind [STATUS] - Primary Only BAM file generated.')

    return output_bam