"""
Patrick O'Neill
pon2@umbc.edu
Wed May 22 13:43:26 EDT 2013

The purpose of this program is to provide a simple generative model
for chIP-seq datasets.

"""

from utils import *
genome_length = 1000000
pre_genome = random_site(genome_length) # random string of genome length
binding_sites_num = 20
binding_site_width = 10
binding_site_ic = 10
binding_site_epsilon = 0.05 # accept a motif w/ ic within epsilon of
                            # binding_site_ic.
chip_seq_reads = 10**6

binding_sites = generate_greedy_motif_with_ic(binding_site_ic,
                                              binding_site_epsilon,
                                              binding_sites_num,
                                              binding_site_width,
                                              verbose=True)

binding_site_locations = [random.randrange(genome_length - binding_site_width + 1)
                          for _ in binding_sites]

# stich the binding sites into the genome...
genome = foldl(lambda g,(s,i):subst(g,s,i),pre_genome,zip(binding_sites,
                                                         binding_site_locations))

def make_fragments(G,b):
    """Return fragments of genome of length G, for b breaks"""
    fragments = [(0,G)]
    for i in xrange(b):
        breakage_location = random.randrange(G)
        fragment = head(fragments,
                        lambda(start,stop):start < breakage_location < stop)
        fragments.remove(fragment)
        start,stop = fragment
        fragments.extend([(start,breakage_location),(breakage_location,stop)])
    return fragments

def analyze_lengths(fragments,G,b):
    lengths = [stop - start for (start,stop) in fragments]
    _lambda = b/float(G)
    plt.hist(lengths)
    plt.plot(*pl(lambda x:len(lengths)*_lambda * exp(-x*_lambda),range(max(lengths))))
    plt.show()
