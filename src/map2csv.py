"""
Convert a map file to a csv of the format:

base position, base, read count

"""
import sys
from utils import get_ecoli_genome,verbose_gen

def interpret_map_file(filename):
    chip_seq_data = [0 for _ in get_ecoli_genome()]
    lines_so_far = 0
    reads = 0
    read_length = 51
    n = len(get_ecoli_genome())
    with open(filename) as f:
        for line in f:
            lines_so_far += 1
            if line.startswith("UB-NGS"):
                reads += 1
                fields = line.split("\t")
                strand = fields[1]
                assert strand == "+" or strand == "-"
                strand_direction = 1 if strand == "+" else -1
                index = int(fields[3])
                if reads % 10000 == 0:
                    print index,lines_so_far/13460000.0
                for i in xrange(read_length):
                    chip_seq_data[(index + i*strand_direction)%n] += 1
    return chip_seq_data

def write_chip_seq_data(chip_seq_data,filename):
    genome = get_ecoli_genome()
    with open(filename,'w') as f:
        for (i,(base,val)) in verbose_gen(enumerate(zip(genome,chip_seq_data)),10**5):
            f.write("%s,%s,%s\n"%(i,base,val))

if __name__ == "__main__":
    map_filename,out_filename = sys.argv[1:3]
    chip_seq_data = interpret_map_file(map_filename)
    write_chip_seq_data(chip_seq_data,out_filename)
