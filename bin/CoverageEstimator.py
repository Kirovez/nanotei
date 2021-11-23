import random
import subprocess
from statistics import median, mean
import numpy as np
import pysam


class CovEstimation:
    def __init__(self, bam_file, sample_per_chr=20):
        print('Estimation of genome coverage by random sampling')
        self.bam_file = bam_file
        self.sample_per_chr = sample_per_chr
        self.chr_indexes = self.__getIndexes(self.bam_file)
        self.gcov = self.__getCoverageGenomeBam(self.bam_file, sample_per_chr=self.sample_per_chr)
        print('Poisson distribution modelling')
        self.pois_distribution = np.random.poisson(round(self.gcov / 2), 1000)

    def getPvalueGreater(self, val):
        """
        return number of elements of 1000 from Poisson distribution greater than val
        """
        return (np.sum(self.pois_distribution >= val) / 1000)

    def getPvalueSmaller(self, val):
        """
        return number of elements of 1000 from Poisson distribution smaller than val
        """
        return (np.sum(self.pois_distribution <= val) / 1000)

    def __getIndexes(self, bam_file):
        """
        returns dictionary with indexed chromosomes
        """
        d = {}
        samfile = pysam.AlignmentFile(bam_file, "rb")
        for chr in samfile.references:
            chr_i = samfile.references.index(chr)
            print(chr, samfile.lengths[chr_i])
            d[chr] = chr_i
        return d

    def __getCoverageGenomeBam(self, bam_file, sample_per_chr=20):
        """
        estimate the read coverage over sample_per_chr random dots per chromosome
        returns:: mean genome coverage
        """
        samfile = pysam.AlignmentFile(bam_file, "rb")
        # print(samfile.pileup("NC_003070.9",8515945,8515946)[0].n)
        cov = []
        for chr in samfile.references:
            chr_i = self.chr_indexes[chr]
            random_pos = random.sample(range(0, samfile.lengths[chr_i], 1), sample_per_chr)
            for positions in random_pos:
                pos_cv = samfile.count_coverage(chr, positions, positions + 1, quality_threshold=0)
                cov += [i[0] for i in pos_cv if i[0] != 0]
        samfile.close()

        gcov = median(sorted(cov)) + 1
        print("Genome coverage median", gcov)
        return gcov

    ### estimate coverage in the TEI
    def isTEI_cov_OK(self, jbrpos):
        """
        it will check median coverage in the interval and return true if the coverage fits the Poisoon dis with lamda = mean genome coverage / 2
        returns false if the coverage fits the Poisoon dis
        """
        jbrpos = jbrpos.replace("..", '-')
        out = subprocess.check_output(["samtools", "depth", "-r", f'{jbrpos}', self.bam_file], encoding="utf-8")
        cov = []
        for lines in out.split("\n"):
            if lines:
                chr, pos, depth = lines.split('\t')
                cov.append(int(depth))
        if len(cov) < 2:
            return False

        m = mean(sorted(cov))
        if self.getPvalueGreater(round(m / 2)) > 0.05:
            return True
        else:
            print(f"Reads depth ({m}) in TEI {jbrpos} is outlier")
            return False