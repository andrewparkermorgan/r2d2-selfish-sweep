from __future__ import print_function

import sys
from time import gmtime, strftime
import numpy as np

class MsRun(object):

    def __init__(self, segsites, positions, samples, command = None, seeds = None):
        self.command = command
        self.seeds = list(seeds)
        #assert(segsites == len(positions) == len(samples[0]))
        self.segsites = segsites
        self.positions = list(positions)
        self.samples = samples
        self.nchroms = len(samples)
        self.nind = int(self.nchroms/2)
        self._pos_format = 1.4 # force some trailing zeros
        self._parse_samples()

    def _parse_samples(self):
        """Parse the 'samples' (haplotypes) spit out by ms into a matrix of site x haps."""
        self.genotypes = np.ndarray((self.segsites, self.nchroms), dtype = int)
        for i in range(0, self.nchroms):
            geno = [ int(x) for x in self.samples[i] ]
            self.genotypes[:,i] = geno

    def write_vcf(self, stream = sys.stdout, length = 1e6, samples = None):
        """Write this simulation run as a VCF file."""
        if not samples:
            samples = [ "N1." + str(x+1) for x in range(0, self.nchroms/2) ]
        print("##fileformat=VCFv4.2", file = stream)
        print("##date={}".format(strftime("%Y-%m-%d %H:%M:%S", gmtime())), file = stream)
        print("##seeds={}".format(self.seeds), file = stream)
        print("##command={}".format(self.command), file = stream)
        print("#CHROM","POS","ID","REF","ALT","QUAL","FILTER","INFO","FORMAT",*samples, file = stream, sep = "\t")
        last_pos = 0
        for i in range(0, self.segsites):
            pos = int(length*self.positions[i])
            if pos == 0 or pos <= last_pos:
                pos = last_pos + 1
            last_pos = pos
            geno = [ "{}|{}".format(self.genotypes[i,2*j], self.genotypes[i,2*j+1]) for j in range(0, self.nind) ]
            print("chr0", pos, ".", "C", "T", ".", "PASS", ".", "GT", *geno, file = stream, sep = "\t")

    def __str__(self):
        out = (self.segsites, ' '.join(['{0:.4f}'.format(pos) for pos in self.positions]), '\n'.join(self.samples))
        # note: the trailing space after positions is in the original MS
        return "\n//\nsegsites: %d\npositions: %s \n%s\n" % out

class MsReader(object):
    def __init__(self, file_handle):
        """
        Initialize MS reader.
        """
        self._file_handle = file_handle
        self.command = None
        self.seeds = None
        self._get_header()
        self.replicates = list()

    def _get_header(self):
        """
        Get the MS command and seeds.
        """
        self.command = next(self._file_handle).strip()
        #try:
        #    self.seeds = map(int, next(self._file_handle).strip().split())
        #except:
        #    self.seeds = map(str, next(self._file_handle).strip().split())
        self.seeds = next(self._file_handle).strip().split()
        while True:
            line = next(self._file_handle)
            if line.startswith("//"):
                break

    def __iter__(self):
        return self

    def __next__(self):
        return self.next()

    def next(self):
        """
        Iterator over simulations.
        """
        segsites = int(next(self._file_handle).strip().split(": ")[1])
        positions = map(float, next(self._file_handle).strip().split(": ")[1].split(" "))
        samples = list()
        line = next(self._file_handle).strip()
        while not line.startswith("//"):
            if len(line) > 0:
                samples.append(line)
            try:
                line = next(self._file_handle).strip()
            except StopIteration:
                break
        return MsRun(segsites, positions, samples, self.command, self.seeds)

    def read_replicates(self):
        """
        Read all replicates into a list.
        """
        for rep in self:
            self.replicates.append(rep)
        return self.replicates

    @property
    def header(self):
        """
        Return a string representing the header of this MS object (e.g.
        the command and seeds).
        """
        # trailing space after command follows MS.
        return "%s \n%s\n" % (self.command, ' '.join(map(str, self  .seeds)))

if __name__ == "__main__":

    reader = MsReader(open("OUT", "r"))
    for sim in reader:
        print(sim.command)
        print(sim.segsites)
        sim.write_vcf()
