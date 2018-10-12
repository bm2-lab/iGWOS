from __future__ import division
try:
    import pysam
except:
    print('*'*20+'\ncan not import pysam\n'+'*'*20)
import re
from collections import Counter


def find_initial_read_pileups(bam_file,
                              banned_chroms=["alt", "random", "Un"],
                              depth_thresh=5):
    """
    Identifies any pileups in aligned sequencing reads that exceed
    a user specified threshold. Returns the genomic locus (coordinates)
    of contiguous sequence above the pileup threshold.

    Args:
        bam_file (str):
            The full path to a bam file that has been sorted and
            indexed.
        depth_tresh (int):
            The minimum number of reads required to call a pileup.
            Default: 5
        banned_chroms (list):
            A list of strings that cannot exist within the chromosome
            name for a given pileup location.
            Default: ["alt", "random", "Un", "chrM"]
    Returns:
        peak_locations (list):
            A list of genomic coordinates in the form chr1:1234-5678
            that exceed the read pileup threshold.
    """
    peak_locations = []
    bam = pysam.AlignmentFile(bam_file, 'rb')
    pileup = bam.pileup()
    in_peak = False
    # Various parameters for the peak identification.
    peak_start = -1
    peak_end = -1
    last_chrom = -1
    peak_start = -1
    peak_end = -1
    last_pos = -1
    for col in pileup:
        # Ignore undesired chromosomes.
        bad_chroms = [b for b in banned_chroms if b in col.reference_name]
        if any(bad_chroms):
            continue
        if in_peak:
            # If any of these are true, terminate the peak.
            if(int(col.nsegments) < int(depth_thresh)) or \
                    (int(col.reference_pos) != (last_pos + 1)) or \
                    (col.reference_name != last_chrom):
                in_peak = False
                peak_end = last_pos
                # Ignore 1 length peaks.
                if peak_end == peak_start:
                    continue
                coord = "{c}:{s}-{e}".format(c=last_chrom,
                                             s=peak_start,
                                             e=peak_end)
                peak_locations.append(coord)
        else:
            if(int(col.nsegments) >= int(depth_thresh)):
                in_peak = True
                peak_start = col.reference_pos
        last_pos = col.reference_pos
        last_chrom = col.reference_name
    return peak_locations


def call_site_seq_features(initial_peaks, bam_file,
                           pyfaidx_fasta_obj, out_name,
                           min_dup=5, min_read_length=60,
                           flank_from_feature=20,
                           required_cliff_ratio=0.9,
                           banned_chroms=["alt", "random", "Un", "chrM"]):
    """
    This function takes peaks and determines whether or
    not they are site-seq features. The steps for feature
    detection are as follows:

    1) Remove peaks that lie within undesired chromosomes.
    2) For each peak, identify the count for locations where sequencing
       reads terminate (cliffs, either 5' or 3') at the same position.
       Reads that are not long enough (determined by min_read_length)
       are removed from the counts.
    3) Ensure that the number of cliff reads for that peak is high enough
       (determined by min_read_length). If the read count for the cliff is
       too low the potential feature is ignored.
    4) Ensure that the ratio of non-cliff reads to cliff reads is not
       too high. If the ratio is too high (determined by required_cliff_ratio)
       then the potential feature is ignored.
    5) All potential features that pass the above filters are reported.

    Args:
        initial_peaks (list):
            A list of peaks already found. These are strings formatted like
            'chr#:#-#'
        bam_file (str):
            The path and name of a sorted and indexed bam file.
        pyfaidx_fasta_obj (pyFaidx Fasta Object):
            A pyFaidx object created from a reference genome.
        out_name (str):
            The name/path of the output fasta file. The fasta file contains
            information about the peak in the header as well as the sequence
            around
        min_dup (int):
            The number of reads that terminate at the exact same location
            (cliff_reads) needed for a potential feature to be considered.
            Default: 5 reads
        min_read_length (int):
            The minimum read length required. Reads shorter than this will
            be removed.
            Default: 60.
        flank_from_feature (int):
            Sets the sequence around the SITE-Seq feature to return in the
            fasta file output by the function.
            Default: 20.
        required_cliff_ratio (float):
            The minimum ratio of cliff reads to cliff spanning reads
            required in order for a feature to be considered.
            Default: 0.9.
        banned_chroms (list):
            A list of chromosome names that should not be considered when
            attempting to identify features. This may need to be
            changed depending on the reference used.

    Returns:
        out_name (str):
            The path to the fasta file that contains the features.
    """
    bam = pysam.AlignmentFile(bam_file, 'rb')
    features = []

    for i, peak in enumerate(initial_peaks):
        r1_starts = []
        chrom, start, end = parse_coordinate(peak)
        if any(p for p in banned_chroms if p in chrom):
            continue
        try:
            reads = bam.fetch(chrom, start - 5, end + 5)
        except:
            continue
        for read in reads:
            read_start = read.reference_start
            read_end = read.reference_end
            read_on_minus = read.is_reverse
            if read_start is None or read_end is None:
                continue
            aligned_length = abs(int(read_start) - int(read_end))
            if aligned_length < min_read_length:
                continue
            read_append = read_end if read_on_minus else read_start
            r1_starts.append(read_append)
        terminations = Counter(r1_starts)
        if len(terminations) == 0:
            continue
        r1_start = terminations.most_common(1)[0][0]
        r1_count = terminations.most_common(1)[0][1]
        if r1_count < min_dup:
            continue
        non_cliff_reads_at_cliff = 0
        if not r1_start:
            continue
        for read in bam.fetch(chrom, r1_start - 1, r1_start + 1):
            if read.reference_start != r1_start and \
                    read.reference_end != r1_start:
                non_cliff_reads_at_cliff += 1
        try:
            cliff_ratio = non_cliff_reads_at_cliff / r1_count
        except ZeroDivisionError:
            cliff_ratio = 0
        if cliff_ratio >= required_cliff_ratio:
            continue
        fasta_dict = {}
        fasta_dict["feature_region"] = "{c}:{s}-{e}".format(
            c=chrom, s=r1_start - flank_from_feature,
            e=r1_start + flank_from_feature)
        fasta_dict["feature_cut"] = "{c}:{s}".format(c=chrom,
                                                     s=r1_start)
        fasta_dict["r1_start_count"] = r1_count
        fasta_dict["reads_near_r1"] = non_cliff_reads_at_cliff
        features.append(fasta_dict)
    fasta_from_features(out_name, features, pyfaidx_fasta_obj)
    return out_name


### HELPER FUNCTIONS ###


def fasta_from_features(out_name, features, pyfaidx_fasta_obj):
    """
    This function requires a reference genome in a fasta file.
    See https://github.com/mdshw5/pyfaidx for more information.
    """
    with open(out_name, "w+") as f:
        for feature in features:
            sequence = retrieve_sequence(feature["feature_region"],
                                         pyfaidx_fasta_obj)
            header = ">{region}|{cut}|{r1_start}|{start_count}".format(
                region=feature["feature_region"], cut=feature["feature_cut"],
                r1_start=feature["reads_near_r1"],
                start_count=feature["r1_start_count"])
            f.write(header + "\n")
            f.write(sequence + "\n")
    return out_name


def retrieve_sequence(coordinate, pyfaidx_fasta_obj):
    """
    Obtain sequence from a pyfaidx fasta object given a
    coordinate.
    """
    chrom, start, stop = parse_coordinate(coordinate)
    return pyfaidx_fasta_obj[chrom][start - 1:stop].seq


def parse_coordinate(coordinate):
    """
    Parse a coordinate formatted like :
    chrZ:#-# OR chrZ:# where Z and # are placeholders.
    """
    coord_range_pattern = r'(.*):\d+-\d+$'
    coord_site_pattern = r'(.*):\d+$'
    if re.match(coord_range_pattern, coordinate):
        chrom = coordinate.split(':')[0]
        start = int(coordinate.split(':')[1].split('-')[0])
        stop = int(coordinate.split(':')[1].split('-')[1])
        return chrom, start, stop
    elif re.match(coord_site_pattern, coordinate):
        chrom = coordinate.split(':')[0]
        start = int(coordinate.split(':')[1])
        return chrom, start
