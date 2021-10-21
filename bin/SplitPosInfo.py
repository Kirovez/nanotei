
class SplitPosInfo:
    """

    """

    def __init__(self, chrom, pos, read_id, start, end):
        self.chrom = chrom
        self.pos_start = pos
        self.pos_end = pos
        self.read_id = read_id
        self.start_in_read = start
        self.end_in_read = end
        self.count_support = 0
        self.reads = [[self.read_id, self.start_in_read, self.end_in_read]]

    def mergePositions(self, SplitPosInfo_object):
        self.pos_end = SplitPosInfo_object.pos_end
        self.count_support += 1
        self.reads.append([SplitPosInfo_object.read_id,
                           SplitPosInfo_object.start_in_read,
                           SplitPosInfo_object.end_in_read])

    def getJbrowseCoords(self):
        return "{0}:{1}..{2}".format(self.chrom, self.pos_start, self.pos_end)

