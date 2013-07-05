__version__ = "0.4"

import _khmer
from _khmer import new_ktable
from _khmer import new_hashtable
from _khmer import _new_counting_hash
from _khmer import _new_hashbits
from _khmer import new_readmask
from _khmer import new_minmax
from _khmer import consume_genome
from _khmer import forward_hash, forward_hash_no_rc, reverse_hash
from _khmer import set_reporting_callback

from filter_utils import filter_fasta_file_any, filter_fasta_file_all, filter_fasta_file_limit_n


from hashlib import sha1
import bisect
import math


###

def new_hashbits(k, starting_size, n_tables=2):
    primes = get_n_primes_above_x(n_tables, starting_size)
    
    return _new_hashbits(k, primes)

def new_counting_hash(k, starting_size, n_tables=2):
    primes = get_n_primes_above_x(n_tables, starting_size)
    
    return _new_counting_hash(k, primes)

def load_hashbits(filename):
    ht = _new_hashbits(1, [1])
    ht.load(filename)

    return ht

def load_counting_hash(filename):
    ht = _new_counting_hash(1, [1])
    ht.load(filename)
    
    return ht

def _default_reporting_callback(info, n_reads, other):
    print '...', info, n_reads, other

def reset_reporting_callback():
    set_reporting_callback(_default_reporting_callback)

reset_reporting_callback()

def calc_expected_collisions(ht):
    """
    A quick & dirty expected collision rate calculation
    """
    sizes = ht.hashsizes()
    n_ht = float(len(sizes))
    occupancy = float(ht.n_occupied())
    min_size = min(sizes)

    fp_one = occupancy / min_size
    fp_all = fp_one ** n_ht

    return fp_all

###

class KmerCount(object):
    def __init__(self, size, report_zero=False):
        self._kt = new_ktable(size)
        self.report_zero = report_zero

    def consume(self, seq):
        self._kt.consume(seq)

    def _get_pairs(self):
        kt = self._kt
        size = kt.n_entries()
        
        for i in range(0, size):
            count = kt.get(i)
            if count or self.report_zero:
                kmer = kt.reverse_hash(i)
                yield kmer, kt.get(i)

    pairs = property(_get_pairs)

    def __getitem__(self, k):
        return self._kt.get(k)

def is_prime(n):
   '''
   checks if a number is prime
   '''
   if n < 2:
      return False
   if n == 2:
      return True
   if n % 2 == 0:
      return False
   for x in range(3, int(n**0.5)+1, 2):
      if n % x == 0:
         return False
   return True

def get_n_primes_near_x(n, x):
   '''
   steps backward until n primes (other than 2) have been
   found that are smaller than x.
   '''
   primes = []
   i = x-1
   if i % 2 == 0:
      i -= 1
   while len(primes) != n and i > 0:
      if is_prime(i):
         primes.append(i)
      i -= 2
   return primes

def get_n_primes_above_x(n, x):
   '''
   steps forward until n primes (other than 2) have been
   found that are smaller than x.
   '''
   primes = []
   i = x+1
   if i % 2 == 0:
      i += 1
   while len(primes) != n and i > 0:
      if is_prime(i):
         primes.append(i)
      i += 2
   return primes


class hll(object):
  def __init__(self,error_rate):

    self.error_rate=error_rate
    self.b=int(math.ceil(math.log((1.04 / self.error_rate) ** 2, 2)))
    self.bits=self.b
    self.alpha=self.get_alpha(self.b)
    self.num_bins= 1 << self.bits
    self.bit_bins=[ 1L << i for i in range(160 - self.bits + 1) ]
    self.estimators = [0 for i in range(self.num_bins)]




  def get_alpha(self, b):


    if not (4 <= b <= 16):
      raise ValueError("b=%d should be in range [4 : 16]" % b)

    if b == 4:
        return 0.673

    if b == 5:
        return 0.697

    if b == 6:
        return 0.709

    return 0.7213 / (1.0 + 1.079 / (1 << b))

  def rho(self,w):
    return len(self.bit_bins) - bisect.bisect_right(self.bit_bins, w)

  def add(self,word):
    word = str(word)
    hash = long(sha1(word).hexdigest(), 16)


    # here, 'bin' is determined by the first 'bits' bits of hash


    bin = hash & ((1 << self.bits) - 1)


    # now count the number of 0s in the remaining bits
    remaining_bits = hash >> self.bits


    count = self.rho(remaining_bits)



    # take max of currently stored estimation & this one
    self.estimators[bin] = max(self.estimators[bin], count)

  def cardinality(self):
    E = self.alpha * float(len(self.estimators) ** 2) / sum(math.pow(2.0, -x) for x in self.estimators)


    if E <= 2.5 * self.bits:                      
      V = self.estimators.count(0)           
      return self.estimators * math.log(self.estimators/ float(V)) if V > 0 else E
    elif E <= float(1L << 160) / 30.0:
      return E
    else:
      return -(1L << 160) * math.log(1.0 - E / (1L << 160))




    
    
