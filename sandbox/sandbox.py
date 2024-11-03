import pysam


def print_read_details(read):
    """Prints the details of a read in a nice format."""
    print("Read Name:", read.query_name)
    print("Sequence:", read.query_sequence)
    print("Flag:", read.flag)
    print("Reference Name:", read.reference_name)
    print("Reference Start Position:", read.reference_start)
    print("Mapping Quality:", read.mapping_quality)
    print("CIGAR String:", read.cigarstring)
    print("Alignment Length:", read.query_alignment_length)
    print("Is Read Paired:", read.is_paired)
    print("Is Read Mapped:", read.is_mapped)
    print("Is Read Reverse Complemented:", read.is_reverse)
    print("Is Read Secondary Alignment:", read.is_secondary)
    print("-" * 50)  # Separator for readability

def main():
    # Open the BAM file
    BAM_FP = ''
    with pysam.AlignmentFile(BAM_FP, "rb") as bamfile:
        # Iterate over the first 2 reads in the BAM file
        for i, read in enumerate(bamfile.fetch()):
            if i == 2:
                break
            print_read_details(read)


