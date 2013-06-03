from utils import *
from sufficache import PSSM
from matplotlib import pyplot as plt
from motifs import *
from sklearn import svm
from sklearn.ensemble import RandomForestClassifier
from itertools import product,combinations

AraCa = "Sample_AMS-ara-gal-Ec_AraCa/AraCa.map"
AraCb = "Sample_AMS-ara-gal-Ec_AraCb/AraCb.map"        
    
# http://arep.med.harvard.edu/ecoli_matrices/dat/araC.dat
arac_motif = ["TATGCCATAGCATTTTTATCCATAAGATTAGCGGATCCTACCTGACGC",
              "AAACGCGTAACAAAAGTGTCTATAATCACGGCAGAAAAGTCCACATTG",
              "CATGTCGCAGCAATTTAATCCATATTTATGCTGTTTCCGACCTGACAC",
              "AGAATTACAGTGAGAACGTGCATAAATTTAGCGGGAAAAGACATAAGG",
              "GAATGCACAGCAGATTAATCCATAAGATTAGCCTGGAAATCCTTGTTG",
              "AACTATTCAGCAGGATAATGAATACAGAGGGGCGAATTATCTCTTGGC"]


print "loaded"

def call_peaks(chip_seq_data,n):
    """Find the n largest peaks in a chipseq dataset"""
    correlation_length = 1000
    avg = mean(chip_seq_data)
    sorted_is = list(reversed(sorted_indices(chip_seq_data)))
    peaks = []
    def i_in_peaks(i,peaks):
        return any([start <= i <= stop for (start,stop) in peaks])
    def i_near_peaks(i,peaks):
        near_peaks = [(start,stop) for (start,stop) in peaks
                      if start - correlation_length < i < stop + correlation_length]
        return near_peaks[0] if near_peaks else None
    for i in sorted_is:
        if i_in_peaks(i,peaks):
            print "found ",i,"in peaks"
            continue
        near_peak = i_near_peaks(i,peaks)
        if near_peak:
            start,stop = near_peak
            print "found ",i,"near peak",start,stop
            peaks.remove((start,stop))
            if i < start:
                new_start,new_stop = i,stop
            else:
                new_start,new_stop = start,i
            print "revising peak",new_start,new_stop
            peaks.append((new_start,new_stop))
        else:
            print "adding new peak:",(i,i)
            peaks.append((i,i))
        if len(peaks) == n:
            break
    return peaks
        
    
true_signal = [0] * 100
true_signal[25] = 1
true_signal[75] = 1
true_signal[50] = 1
true_signal[55] = 0

num_reads = 10000
read_length = 10
fragments = [(r,r+read_length)
         for r in [random.randrange(len(true_signal)) for __ in range(num_reads)]]
pulldowns = [(start,stop) for (start,stop) in fragments
             if any(true_signal[start:stop])]
signal = [len([(start,stop) for (start,stop) in pulldowns if start <= i < stop])
          for i in range(len(true_signal))]

triangle_kernel = normalize(range(9,0-1,-1) + [0] * 80 + range(0,10))

def deconvolve(hs,gs):
    """Assume h = f***g, where f is the true signal and g is a
    convolution filter.  (*** is convolution.)  Given h and g, recover
    f through the identity f = INV_FFT(H/G) where / is point-wise division """
    Hs = dft(hs)
    Gs = dft(gs)
    return inv_dft([H/G for H,G in zip(Hs,Gs)])

def write_arff(chip_seq_data,w,filename):
    n = len(chip_seq_data)
    genome = get_ecoli_genome()
    genome = genome[-w:] + genome + genome[:w]
    with open(filename,'w') as f:
        for i, score in enumerate(chip_seq_data):
            window = genome[i-w:i+w]
            window_string = str(list(window))[1:-1]
            f.write(window_string + "," + str(score) + "\n")
            if i % 100000 == 0:
                print i,window_string

def dummy_encode_base(b):
    return map(lambda c:int(b==c),"ACGT")

def dummy_encode(seq):
    return concat(map(dummy_encode_base,seq))

def cross_validate_svm(motif=Escherichia_coli.Crp,kernel="linear",
                       negative_set_factor=10**2,verbose=False):
    print "starting"
    n = len(motif)
    w = len(motif[0])
    negative_instances = negative_set_factor*n
    print "constructing data set"
    X = map(dummy_encode,motif + [random_site(w) for i in range(negative_instances)])
    y = [1 for i in range(n)] + [0 for i in range(negative_instances)]
    data = transpose([X,y])
    print "constructing cv dataset"
    cv_data = cv(data)
    precisions = []
    recalls = []
    accuracies = []
    f_scores = []
    print "entering for loop"
    for i,(train_set,test_set) in enumerate(cv_data):
        train_instances,train_labels = transpose(train_set)
        test_instances,test_labels = transpose(test_set)
        #classifier = svm.SVC(kernel=kernel)
        classifier = RandomForestClassifier(n_estimators=1000)
        print "fitting"
        classifier.fit(train_instances,train_labels)
        print "done fitting"
        predictions = classifier.predict(test_instances)
        p = precision(predictions,test_labels)
        r = recall(predictions,test_labels)
        a = accuracy(predictions,test_labels)
        f = f_score(predictions,test_labels)
        precisions.append(p)
        recalls.append(r)
        accuracies.append(a)
        f_scores.append(f)
        if verbose:
            print "Round %s of cross-validation" % i
            print "Precision:",p
            print "Recall:",r
            print "Accuracy:",a
            print "F-score:",f
    print "10-fold CV statistics: mean w/ 95% CI"
    print "Average precision:",mean(precisions), "+/-", 1.96*sd(precisions)
    print "Average recall:",mean(recalls), "+/-", 1.96*sd(recalls)
    print "Average accuracy:",mean(accuracies), "+/-", 1.96*sd(accuracies)
    print "Average f_score:",mean(f_scores), "+/-", 1.96*sd(f_scores)

def make_training_set(chip_data,n,window_length=750):
    """Construct a training set of the form [(window,count value)] of
    size n"""
    genome = get_ecoli_genome()
    genome_length = len(genome)
    #pre-empt problems with wrap-around...
    genome = genome[-window_length:] + genome + genome[:window_length] 
    data_set = [(i,genome[i-window_length/2:i+window_length/2],chip_data[i])
                for i in [random.randrange(genome_length)
                          for __ in xrange(n)]]
    for (i,seq,val) in data_set:
        assert len(seq) == 750, (i,seq,val)
    return [(seq,val) for (i,seq,val) in data_set]

lamb = 0.5

def Kp(s,t,i):
    if i == 0:
        return 1
    elif min(len(s),len(t)) < i:
        return 0
    else:
        s,x = s[:-1],s[-1]
        return lamb * Kp(s,t,i) + sum([Kp(s,t[:j],i-1) * lamb**(len(t) - (j + 1) + 2)
                                       for (j,t_j) in enumerate(t) if t_j == x])

def Kp_test(s,t,i):
    if i == 0:
        return 1
    elif min(len(s),len(t)) < i:
        return 0
    else:
        s,x = s[:-1],s[-1]
        return lamb * Kp_test(s,t,i) + Kpp(s+x,t,i)

def Kpp(s,t,i):
    if min(len(s),len(t)) < i:
        return 0
    #print "Kpp: ",s,t,i
    s,x = s[:-1],s[-1]
    u_index = t.rfind(x) #what is the last occurrence of x in t?
    if u_index == len(t) - 1: # if t ends with x...
        #print "iffing"
        t,x_prime = t[:-1],t[-1]
        assert(x_prime == x)
        return lamb * (Kpp(s+x,t,i) + lamb * Kp(s,t,i-1))
    elif u_index >= 0:
        #print "eliffing"
        t,u = t[:u_index+1],t[u_index+1:]
    else: # if x does not appear in t...
        #print "elsing"
        t,u = "",t
    return lamb **(len(u)) * Kpp(s+x,t,i)

def K(s,t,i):
    if min(len(s),len(t)) < i:
        return 0
    else:
        s,x = s[:-1],s[-1]
        return K(s,t,i) + sum([Kp(s,t[:j],i-1) * lamb**2
                               for (j,t_j) in enumerate(t) if t_j == x])

def K_test(s,t,i):
    if min(len(s),len(t)) < i:
        return 0
    else:
        s,x = s[:-1],s[-1]
        return K_test(s,t,i) + sum([Kp_test(s,t[:j],i-1) * lamb**2
                               for (j,t_j) in enumerate(t) if t_j == x])

def normK(s,t,i):
    return K(s,t,i)/(sqrt(K(s,s,i)*K(t,t,i)))
    
def sigma(alphabet,n):
    """Return all words of length n over alphabet"""
    return product(alphabet,repeat=n)

def l(seq):
    """implement slice length function (Lodhi et al. 2002 p. 423)"""
    return seq[-1] - seq[0] + 1

def slices_of(seq,n):
    return combinations(range(len(seq)),n)

def K_explicit(s,t,n):
    alphabet = set(s + t)
    words = sigma(alphabet,n)
    return sum([lamb**(l(I)+l(J))
                for u in words
                for I in slices_of(s,n)
                for J in slices_of(t,n)
                if "".join(rslice(s,I)) == "".join(rslice(t,J)) == "".join(u)])
    
def rolling_window_correlation_exp():
    """What is the radius of the rolling average window that maximizes
    the correlation between TRAP probabilities, chIP-seq read intensities?"""
    pssm = PSSM(arac_motif)
    genome = get_ecoli_genome()
    traps = pssm.slide_trap(genome)
    beta = 1.61 #approximately = 1/(kBT) in kcal/mol @ room temp
    z = sum(exp(-beta*ep) for ep in traps)
    probs = [exp(-beta*ep)/z for ep in traps]
    corrs = ([(k,spearmanr(probs[:10**5],
                           circular_rolling_average(chip_seq_data[:10**5],k))[0])
              for k in verbose_gen(range(250,750,25))])
    plt.plot(*transpose(corrs))
    plt.xlabel("Rolling Window Radius (bp)")
    plt.ylabel(r"Spearman $\rho$")
    def smartwrap(text):
        from textwrap import wrap
        return "\n".join(wrap(text))
    plt.subplots_adjust(top=0.8)
    plt.title(smartwrap("Spearman Correlation vs. Radius of Rolling Average for TRAP probabilities, ChIP-seq reads (AraC a) on bp 1-10**5 of E. coli genome"))
    plt.savefig("correlation_vs_rolling_average_radius.png",dpi=400)
    """Conclusion: correlation is maximized for window radius 400-500
    bp, which is consistent with a fragment size on the same length scale."""
