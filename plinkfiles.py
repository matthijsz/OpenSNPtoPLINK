import math
import gc
import zlib


# BedBlock class defines the binary .bed-block for a single variant
# Can be compressed whenever possible in order to save memory, but this will significantly impact performance
class BedBlock:
    def __init__(self, n_participants, compress):
        self.compress = compress
        self.n_subjects = n_participants
        self.full_participant_bytes = math.floor(len(n_participants)/4)
        block_length = math.ceil(len(n_participants)/4)
        # Some trickery is required because the content of the last byte depends on number of participants:
        if self.full_participant_bytes == block_length:
            self.bin = bytearray(b'\x55' * block_length)  # \x55 is 4 missing genotypes:  01010101
        else:
            self.bin = bytearray(b'\x55' * self.full_participant_bytes)
            n_in_last_block = n_participants % 4  # Number of missing genotypes in last block
            last_block = 0  # Last block is intially \x00
            for i in n_in_last_block:
                last_block += 1 << (i*2)  # And 01 is added for each participant that will still end up in this block
            self.bin += last_block.to_bytes(1, byteorder='little')
        if self.compress:
            self.bin = zlib.compress(self.bin)
        self.counter = 0

    def modify(self, n_participant, genotype):
        if self.compress:
            self.bin = bytearray(zlib.decompress(self.bin))
        bitshift = (n_participant % 4) * 2
        genotype = genotype << bitshift
        # First subtract the standard missing value, then add new genotype to existing byte
        new_byte = (self.bin[n_participant // 4] - (1 << bitshift)) + genotype
        self.bin[n_participant // 4] = new_byte
        if self.compress:
            self.bin = zlib.compress(self.bin)

    def __bytes__(self):
        if self.compress:
            return zlib.decompress(self.bin)
        else:
            return bytes(self.bin)

    def close(self):
        self.bin = b'\x00'
        gc.collect()


# Bed class defines a full bed file, with as much BedBlock classes as required
class Bed:
    def __init__(self, n_participants, compress=False):
        self.n_participants = n_participants
        self.compress = compress
        self.variants = {}

    def add(self, n_participant, variant, genotype):
        if variant not in self.variants.keys():
            self.variants[variant] = BedBlock(self.n_participants, compress=self.compress)
        self.variants[variant].modify(n_participant, genotype)

    def to_file(self, base_name):
        with open('{}.bed'.format(base_name), 'wb') as f:
            f.write(b'\x6c\x1b\x01')  # PLINK1.9 magic number
            for bedblock in self.variants.values():
                f.write(bytes(bedblock))
                # bedblock.close()  # This can be used to save memory while writing
            # self.variants = {}

    def remove(self, snp):
        del self.variants[snp]


# Bim class defines the .bim file
# This class is also responsible for converting string genotypes (like AG, CT) to numeric (0, 1, 2 or 3)
class Bim:
    def __init__(self):
        self.data = {}

    def update(self, variant, chromosome, bp, str_genotype):
        if len(str_genotype) == 2:
            a1, a2 = str_genotype[0], str_genotype[1]
        elif len(str_genotype) == 1:
            a1, a2 = str_genotype, None
        else:
            return False
        if variant in self.data.keys():
            if len(self.data[variant]) == 5:
                if a1 != self.data[variant][-1]:
                    self.data[variant].append(a1)
                elif (a2 is not None) and (a2 != self.data[variant][-1]):
                    self.data[variant].append(a2)
        else:
            self.data[variant] = [chromosome, variant, 0, bp]
            if a1 == a2:
                self.data[variant].append(a1)
            elif a2 is not None:
                self.data[variant].append(a2)
        return True

    def genotype_to_numeric(self, variant, str_genotype):
        a1 = self.data[variant][4]
        if str_genotype is None:
            return 1
        if len(str_genotype) != 2:
            return 1
        if str_genotype[0] == a1:
            return 0 if str_genotype[1] == a1 else 2
        else:
            return 2 if str_genotype[1] == a1 else 3

    def remove(self, snp):
        del self.data[snp]

    def __len__(self):
        return len(self.data)

    def to_file(self, base_name):
        with open('{}.bim'.format(base_name), 'w') as f:
            for variant, data in self.data.items():
                if len(data) < 6:
                    data.append('-')  # This should only happen with missing a2
                f.write('\t'.join([str(i) for i in data])+'\n')


# Wrapper for working with both Bed and Bim classes simultaneously
class PlinkFiles:
    def __init__(self, n_participants, compress=False):
        self.bed = Bed(n_participants, compress=compress)
        self.bim = Bim()

    def add(self, n_participant, variant, chromosome, bp, str_genotype):
        bim_check = self.bim.update(variant, chromosome, bp, str_genotype)
        if bim_check:
            genotype = self.bim.genotype_to_numeric(variant, str_genotype)
            self.bed.add(n_participant, variant, genotype)

    def save(self, base_name):
        self.bed.to_file(base_name)
        self.bim.to_file(base_name)

    def remove(self, remove_snps):
        for i in remove_snps:
            self.bed.remove(i)
            self.bim.remove(i)