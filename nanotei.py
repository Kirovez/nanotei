import os
import pandas as pd
import pysam
from Bio import SeqIO
from bin.SplitPosInfo import SplitPosInfo
from bin.CoverageEstimator import  CovEstimation
import random


def _getInsertion(align, min_indel=1000):
    """
    identify insertion positions on reference and return False if no insertions found.
    If insertion presents then list [ins position on reference, ins length] is returned
    """
    bases_before_insertion = 0
    for cigar in align.cigartuples:
        if cigar[0] == 1 and cigar[1] > min_indel:
            return (align.reference_start + bases_before_insertion, cigar[1])
        elif cigar[0] == 0 or cigar[0] == 2:
            bases_before_insertion += cigar[1]
    return False


def getSplittedPositions(bam_file, mapping_quality=40, min_len_clipped=1000):
    """
    read bam alignments and collect split coordinates into SplitPosInfo object
    """

    split_positions = {}  # chromosome:split SplitPosInfo

    cnt = 0
    cnt_all_reads = 0
    cnt_i = 0
    non_cnt = 0
    cov = 0
    for algns in pysam.AlignmentFile(bam_file, 'rb'):
        if algns.cigartuples and algns.reference_name != None:
            if algns.mapping_quality >= mapping_quality and not algns.is_supplementary:
                cnt_all_reads += 1
                cov += abs(int(algns.reference_end) - int(algns.reference_start))
                # add chromosome to the dictionary
                if algns.reference_name not in split_positions:
                    split_positions[algns.reference_name] = []
                ##check end of the reads
                if (algns.cigartuples[-1][0] == 5 or algns.cigartuples[-1][0] == 4) and algns.cigartuples[-1][
                    1] > min_len_clipped:
                    if algns.is_reverse:
                        split_positions[algns.reference_name].append(
                            SplitPosInfo(algns.reference_name, algns.reference_end, algns.query_name,
                                         algns.cigartuples[-1][1], 0))
                    else:
                        split_positions[algns.reference_name].append(
                            SplitPosInfo(algns.reference_name, algns.reference_end, algns.query_name, 0,
                                         algns.cigartuples[-1][1]))
                    cnt_i += 1

                ##check start of the reads
                if (algns.cigartuples[0][0] == 5 or algns.cigartuples[0][0] == 4) and algns.cigartuples[0][
                    1] > min_len_clipped:
                    if algns.is_reverse:
                        split_positions[algns.reference_name].append(
                            SplitPosInfo(algns.reference_name, algns.reference_start, algns.query_name, 0,
                                         algns.cigartuples[0][1]))
                    else:
                        split_positions[algns.reference_name].append(
                            SplitPosInfo(algns.reference_name, algns.reference_start, algns.query_name,
                                         algns.cigartuples[0][1], 0))
                    cnt_i += 1

                isInsPresent = _getInsertion(algns)
                if isInsPresent:
                    if algns.is_reverse:
                        inRead_coord = int(isInsPresent[0]) - algns.reference_start
                        split_positions[algns.reference_name].append(
                            SplitPosInfo(algns.reference_name, isInsPresent[0], algns.query_name,
                                         algns.infer_read_length() - inRead_coord - isInsPresent[1],
                                         algns.infer_read_length() - inRead_coord))
                    else:
                        inRead_coord = int(isInsPresent[0]) - algns.reference_start
                        split_positions[algns.reference_name].append(
                            SplitPosInfo(algns.reference_name, isInsPresent[0], algns.query_name, inRead_coord,
                                         inRead_coord + isInsPresent[1]))
                    # print("INSERTION:", algns.reference_name, isInsPresent[0], 'Size:', isInsPresent[1], algns.query_name)

        else:
            non_cnt += 1
    print('Number of aligns', cnt_all_reads)
    print('Number of splitted ends in the reads', cnt_i)
    print('Number of aligns without reference', non_cnt)
    print('Total number of bp from aligned fragments:', cov)
    return (split_positions)


def mergeSplitCount(split_positions, min_distance=100):
    """
    merge SplitPosInfo objects and count a number of split reads for each TEI
    """
    merged_split_positions = {}
    cnt = 0
    for chromosome in split_positions:
        merged_split_positions[chromosome] = []
        new_splits = []
        split_positions[chromosome] = sorted(split_positions[chromosome], key=lambda pos: pos.pos_start)
        for spl_pos in split_positions[chromosome]:
            if not new_splits:
                new_splits.append(spl_pos)
            else:
                if abs(new_splits[-1].pos_end - spl_pos.pos_start) <= min_distance:
                    new_splits[-1].mergePositions(spl_pos)
                else:
                    cnt += 1
                    # print("{4}\t{0}:{1}..{2}\t{3}".format(chromosome,new_splits[-1][0],new_splits[-1][1], new_splits[-1][-1], cnt))
                    new_splits.append(spl_pos)
        merged_split_positions[chromosome] = new_splits
        print("Sites for chromosome", chromosome, 'is', len(new_splits))
    return merged_split_positions

def collectReadClips(merged_splitCounts, fastq):
    index_seq = SeqIO.index(fastq, 'fastq')
    cnt_seqs = 0
    no_reads = 0
    clipFq = fastq + '_collectedSplits.fastq'
    with open(clipFq, 'w') as outFile:
        for chr in merged_splitCounts:
            for teis in merged_splitCounts[chr]:
                tei_id = "{0}:{1}..{2}".format(teis.chrom, teis.pos_start, teis.pos_end)
                for ind, seqs in enumerate(teis.reads):
                    read_id, start_in_read, end_in_read = seqs
                    #                     if read_id == 'fc9ba83b-72e9-4860-8696-40fba7f47537':
                    #                         print( read_id, start_in_read, end_in_read)
                    if read_id in index_seq:
                        seq_new = str(index_seq[read_id].seq)
                        if start_in_read != 0 and end_in_read != 0:  # insertion
                            seq_split = seq_new[start_in_read:end_in_read]
                        elif start_in_read == 0:  # clip in the end
                            seq_split = seq_new[-end_in_read:]
                        else:  # clip start read
                            seq_split = seq_new[:start_in_read]
                        seq_id = tei_id + '#' + str(read_id)
                        outFile.write(">{0}\n{1}\n".format(seq_id, seq_split))
                        cnt_seqs += 1
                    else:
                        no_reads += 1
    print("Total reads were collected:", cnt_seqs)
    print("No reads were found:", no_reads)
    return clipFq


def mapppingBamSort(reads, genome_fasta, mm2='minimap2', samtools_path="samtools"):
    os.system('{3} -ax map-ont -t 150 {0} {1} > {2}.sam'.format(genome_fasta, reads, reads, mm2))
    sam_file = '{0}.sam'.format(reads)
    bam_file = sam_file.rsplit('.', 1)[0] + ".bam"
    if r'/' in bam_file:
        sort_bam_file = bam_file.rsplit(r'/', 1)[0] + r"/sorted_" + bam_file.rsplit(r'/', 1)[1]
    else:
        sort_bam_file = r"./sorted_" + bam_file

    ##sam to bam
    samview = '{2} view -Sb {0} > {1}'.format(sam_file, bam_file, samtools_path)
    print(samview)
    os.system(samview)

    ##sort bam
    sortbam = '{2} sort -o {1} -@ 100 {0} '.format(bam_file, sort_bam_file, samtools_path)
    print(sortbam)
    os.system(sortbam)

    ##index
    os.system('{1} index -@ 100 {0}'.format(sort_bam_file, samtools_path))

    return (sort_bam_file)


def getBedFromClippedBam(bamfile_clips):
    originTEs_per_TEIpositions = {}  # tei_position:[]
    cnt_tei_with_origin = 0
    bedOut = '{}_clipped_reads_pos.bed'.format(bamfile_clips)
    clipped_reads_pos = open(bedOut, 'w')
    for algns in pysam.AlignmentFile(bamfile_clips, 'rb'):
        if algns.cigartuples and algns.reference_name != None:
            if not algns.is_supplementary and not algns.is_secondary:
                TEIposition = algns.query_name.split('#')[0]
                if TEIposition not in originTEs_per_TEIpositions:
                    originTEs_per_TEIpositions[TEIposition] = []
                clipped_reads_pos.write(
                    "{0}\t{1}\t{2}\t{3}\n".format(algns.reference_name, algns.reference_start, algns.reference_end,
                                                  algns.query_name))
                originTEs_per_TEIpositions[TEIposition].append(
                    "{0}:{1}..{2}".format(algns.reference_name, algns.reference_start, algns.reference_end))
    clipped_reads_pos.close()
    return bedOut


def getFinalTable(inter_in, GenCov, overlap_with_origTE=0.5, min_read_support = 2, minpvalue = 0.05, filter_by_pvalue = True):
    # select the best overlapping
    inter_in = pd.read_csv(inter_in, sep="\t", header=None)
    inter = inter_in.loc[(inter_in[2] - inter_in[1]) * overlap_with_origTE < inter_in[8]]  # .groupby([7]).greater()

    # add tei coord column and read id
    inter[['TEIcoord', 'read']] = inter[7].str.split('#', 1, expand=True)
    # find well supported TEIs
    well_supported_TEIs = inter.groupby('TEIcoord').agg({7: 'nunique'}).reset_index()  # .read.nunique()
    well_supported_TEIs.reset_index()

    ## select TEIs with 2 or more reads supported
    before_fil = len(well_supported_TEIs)
    well_supported_TEIs = well_supported_TEIs[well_supported_TEIs[7] >= min_read_support]
    print(f"Filtered by number of min supported reads ({min_read_support}) region coverage: ", before_fil - len(well_supported_TEIs), )

    ## add pvalue column where H0 is that TEI is not somatic
    well_supported_TEIs['pvalue'] = well_supported_TEIs[7].map(GenCov.getPvalueSmaller)

    ## evaluate the read coverage in TEI region
    well_supported_TEIs['isTEIregioIsOK'] = well_supported_TEIs['TEIcoord'].map(GenCov.isTEI_cov_OK)

    well_supported_TEIs_filtered = well_supported_TEIs #
    print("TEI with too high region coverage (not fit Poisson distribution): ", len(well_supported_TEIs) - len(well_supported_TEIs.loc[well_supported_TEIs.isTEIregioIsOK]))

    ## aggregate data by TEI coordintate and TE id
    oneTEI_manyTE = inter.groupby(['TEIcoord', 3]).agg({7: 'nunique'}).reset_index()

    ## leave only TEIs that are in the selected list
    oneTEI_manyTE = oneTEI_manyTE[oneTEI_manyTE.TEIcoord.isin(well_supported_TEIs_filtered['TEIcoord'])]

    # concatenate reads per each TE and TE ids
    oneTEI_manyTE["TE_reads"] = oneTEI_manyTE[3].astype(str) + ":" + oneTEI_manyTE[7].astype(str)

    # join all TE:reads per each TEI
    final_oneTEI_manyTE = oneTEI_manyTE.groupby('TEIcoord')['TE_reads'].apply(lambda x: ','.join(x)).reset_index()

    # merge table with #supported reads and pvalue and TE per each read
    merge_nreads_final_oneTEI_manyTE = pd.merge(well_supported_TEIs_filtered, final_oneTEI_manyTE, on="TEIcoord")

    ## choose one TE per each TEI
    tmp_tab = inter[inter.TEIcoord.isin(well_supported_TEIs_filtered['TEIcoord'])].groupby(['TEIcoord', 3]).agg(
        {7: 'nunique'}).reset_index()
    tmp_tab.groupby(['TEIcoord'])[7].max()
    tmp_tab = tmp_tab.loc[tmp_tab.groupby(['TEIcoord'])[7].idxmax()]
    merge_nreads_final_oneTEI_manyTE2 = pd.merge(merge_nreads_final_oneTEI_manyTE, tmp_tab, on="TEIcoord")

    # add tab with selected single TE
    merge_nreads_final_oneTEI_manyTE2 = merge_nreads_final_oneTEI_manyTE2.drop('7_y', axis=1)
    merge_nreads_final_oneTEI_manyTE2.set_axis(
        ['TEIcoord', '#supp.reads', 'pvalue', 'isTEIregioIsOK', 'TE_reads', 'chosenTE'], axis=1, inplace=True)

    # add TE coordinates
    unique_inter = inter[[0, 1, 2, 3]].drop_duplicates()

    addTE = pd.merge(merge_nreads_final_oneTEI_manyTE2, unique_inter, left_on='chosenTE', right_on=3)
    addTE.set_axis(
        ['TEIcoord', '#supp.reads', 'pvalue', 'isTEIregioIsOK', 'TE_reads', 'chosenTE', "TE.chr", "TE.start", 'TE.end',
         'TEid'], axis=1, inplace=True)
    addTE = addTE.drop('TEid', axis=1)

    if filter_by_pvalue:
        return addTE[addTE['pvalue'] > minpvalue]
    return addTE

def getbed(df):
    tmpbed = df['TEIcoord'].str.split(r':', 1, expand=True)
    tmpbed[['start', 'end']] = tmpbed[1].str.split('\..', 1, expand=True)
    tmpbed['TEIcoord'] = df['TEIcoord']

    return tmpbed[[0, 'start', 'end', 'TEIcoord']]

def run_nanotei(bam_file, fastq, genome_fasta, outfolder, bed_TE, outtab, mapping_quality=40, min_len_clipped=1000,
                 min_distance=20, minimap2_path = 'minimap2',
         samtools_path = "samtools", bedtools = 'bedtools', overlap_with_origTE = 0.5, minpvalue = 0.05,
                filter_by_pvalue = True, write_bed = False, min_read_support = 2):
    outFolder = outfolder
    if outFolder[-1] != '/':
        outFolder += "/"
    os.system(f'mkdir {outFolder} ')

    inter_file = f'{outFolder}refclipped_vs_TEbed.inter'
    print("1. Getting clipped reads from the bam file")
    split_positions = getSplittedPositions(bam_file, mapping_quality=mapping_quality, min_len_clipped=min_len_clipped)
    print("2. Merging neighbor TEIs into one TEIs")
    merged_splitCounts = mergeSplitCount(split_positions, min_distance=min_distance)
    print("3. Collecting clipped ends of the reads")
    cliped_fq = collectReadClips(merged_splitCounts, fastq)
    print("4. Mapping clipped ends of the reads to the reference genome")
    bamfile_clips = mapppingBamSort(cliped_fq, genome_fasta, samtools_path = samtools_path, mm2=minimap2_path)
    print("5. Getting mapping positions of the clipped reads and write bed file")
    bedfile_clips = getBedFromClippedBam(bamfile_clips)
    os.system(f"mv {bedfile_clips} {outFolder}")
    print(f"mv {bedfile_clips} {outFolder}")
    bedfile_clips = outFolder + bedfile_clips.rsplit("/")[-1]
    print(bedfile_clips)
    print("6. Intersecting the bed file with clipped reads and annotated TE positions")
    os.system(f"{bedtools} sort -i {bedfile_clips} > {bedfile_clips}_sorted.bed")
    ## intersect clipped reads mapping positions with TE from Slotkin
    os.system(f"{bedtools} intersect -a {bed_TE} -b {bedfile_clips}_sorted.bed -wo > {inter_file}")
    print("7. Parsing the table and filtering TEIs by number of supported reads and read coverage in the TEI region")
    ## estimate coverage of the genome by bam file
    GenCov = CovEstimation(bam_file, sample_per_chr=20)
    final_tab = getFinalTable(inter_file, GenCov, overlap_with_origTE=overlap_with_origTE, minpvalue=minpvalue, filter_by_pvalue = filter_by_pvalue, min_read_support = min_read_support)
    print('Total number of TEIs in the output table', len(final_tab))
    final_tab.to_csv(outtab, sep='\t', index=False)

    if write_bed:
        getbed(final_tab).to_csv(outtab + ".bed", sep='\t', header = None, index=False)




if __name__ == "__main__":
    import argparse
    parser = argparse.ArgumentParser(description='Script to find regions of insertions using Nanopore reads')
    parser.add_argument('bam_file', help='path to the sorted bam file with mapped reads to the genome')
    parser.add_argument('fastq', help='path to fastq file of reads')
    parser.add_argument('genome_fasta', help='path to target sequence in fasta format')
    parser.add_argument('outfolder', help='path to the output folder')
    parser.add_argument('bed_TE', help='path to the bed file with annotated TE')
    parser.add_argument('outtab', help='path to the output table')
    parser.add_argument('-q', '--map_q', help='minimum mapping quality', default=20,type=int)
    parser.add_argument('-mlc', '--min_len_clipped', help='minimum length of the clipped part', default=1000,type=int)
    parser.add_argument('-md', '--min_merge_distance', help='minimum distance between two reference positions of clipping to be merged into one', default=20,type=int)
    parser.add_argument('-samtp', '--samtools_path', help='path to samtools program', default='samtools')
    parser.add_argument('-mm2', '--minimap2_path', help='path to minimap2 program', default="minimap2")
    parser.add_argument('-bt', '--bedtools', help='path to minimap2 program', default="bedtools")
    parser.add_argument('-ovt', '--overlap_with_origTE', help='minimum portion of origin TE overlapped with clipped repa alignemnt', default=0.5,  type = float)
    parser.add_argument('-mpv', '--minpvalue', help='minimum pvalue to filter (it is used with --fpv)', default=0.05, type = float)
    parser.add_argument('--fpv', help='filter by pvalue', action='store_true')
    parser.add_argument('-mrs', '--min_read_support', help='minimum read support for filtering', default=2, type = int)
    parser.add_argument('--bed', help='filter by pvalue', action='store_true')



    args = parser.parse_args()

    print(args.fastq)
    run_nanotei(args.bam_file, args.fastq, args.genome_fasta, args.outfolder, args.bed_TE, args.outtab, mapping_quality=args.map_q,
                min_len_clipped=args.min_len_clipped,
                min_distance= args.min_merge_distance, minimap2_path = args.minimap2_path,
         samtools_path = args.samtools_path,  bedtools = args.bedtools, minpvalue = args.minpvalue, min_read_support = args.min_read_support,
                filter_by_pvalue = args.fpv, write_bed = args.bed)
