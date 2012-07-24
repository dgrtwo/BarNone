"""
BarcodeCache

Performs caching of barcodes, identifying mismatched ones using flamingo
"""

import copy
import itertools
import collections
import warnings

import flamingo


### FUNCTIONS ###

def iterate_fasta(infile):
    """Iterate over the fastq sequences in a file"""
    inf = open(infile)
    for i, l in enumerate(inf):
        if i % 2 == 1:
            yield l
    inf.close()


def iterate_fastq(infile):
    """Iterate over the fastq sequences in a file"""
    inf = open(infile)
    for i, l in enumerate(inf):
        if i % 4 == 1:
            yield l
    inf.close()


def iterate_qseq(infile):
    """Iterate over the fastq sequences in a file"""
    inf = open(infile)

    for l in inf:
        yield l.split("\t")[8]

    inf.close()


def iterate_txt(infile):
    inf = open(infile)

    for l in inf:
        yield l[:-1]

    inf.close()


def closest_match(original, matches, unique=False):
    """
    Return the closest Levenshtein match if there is one. If unique, return
    None if there is a tie for closest
    """
    if len(matches) == 1:
        return matches[0]

    distances = [flamingo.distance(m, original) for m in matches]

    if unique and distances.count(min(distances)) > 1:
        return None

    return min(zip(distances, matches))[1]


### CLASSES ###

class BarcodeCache(object):
    def __init__(self, barcode_dict):
        """
        Given a dictionary mapping barcodes to values (none of which can
        be None)
        """
        self.barcode_dict = dict([(k, v)
                                    for k, v in barcode_dict.items()])

        # need to keep two caches- one for searches that are limited to unique
        # items, and one for searches that aren't
        self.cache_dicts = dict([(u, dict(self.barcode_dict.items()))
                                    for u in (True, False)])

        self.original = dict([(v, k) for k, v in barcode_dict.items()])
        self.strains = list(set(barcode_dict.values()))
        self.index = flamingo.WrapperSimpleEd(self.barcode_dict.keys())
        self.total_looked_up = 0
        self.total_cached = 0
        self.total = 0

    def search(self, barcode, distance, verbose=False, details=False,
                unique=False):
        """Search for the object mapping from a barcode"""
        self.total += 1

        cache_match = self.cache_dicts[unique].get(barcode)
        if cache_match != None:
            # this is the closest match, so if it doesn't fit with this
            # distance, it won't work
            original = self.original[cache_match]
            if flamingo.distance(barcode, original) > distance:
                return None

            self.total_cached += 1

            return (cache_match, original) if details else cache_match

        if distance == 0:
            # won't be able to find any inexact matches anyway
            return None

        # inexact match
        matches = self.index.search(barcode, distance)

        if len(matches) == 0:
            self.cache_dicts[unique][barcode] = None
            return None

        best = closest_match(barcode, matches, unique=unique)

        if best == None:
            return None

        ret = self.barcode_dict[best]
        self.cache_dicts[unique][barcode] = ret

        self.total_looked_up += 1

        return (ret, best) if details else ret


class BarcodeCacheMultipleLen(object):
    """Contains multiple BarcodeCaches, one for each possible length"""
    def __init__(self, barcode_dict):
        """Given a dictionary mapping barcodes to any kind of object"""
        # divide up the keys by length
        self.original = dict([(v, k) for k, v in barcode_dict.items()])

        lengths = [l for l in map(len, barcode_dict) if l != 0]
        self.common_lengths = list(set(lengths))
        self.common_lengths.sort(key=lambda l: -lengths.count(l))

        self.descending_lengths = sorted(self.common_lengths, reverse=True)

        if len(self.common_lengths) > 1:
            warnings.warn("multiple lengths (" +
                            ",".join(map(str, self.common_lengths)) +
                            ") in barcode dictionary")

        self.barcode_caches = dict([(l, BarcodeCache(dict([(k, v)
                                        for k, v in barcode_dict.items()
                                            if len(k) == l])))
                                        for l in self.common_lengths])

    def get_strains(self):
        return list(set(itertools.chain(*[c.strains
                        for c in self.barcode_caches.values()])))

    strains = property(get_strains)

    def search_all(self, barcode, distance, unique=False, details=False,
                   verbose=False):
        """
        Search all sublengths of this barcode. Return either the matching
        name, or, if details=True, a 3-tuple:
        (name, matching barcode, length)
        """
        #print "search_all distance", distance
        matches = []

        # the order depends on the distance- if we need an exact match,
        # do it in descending order of length, if we are looking for an
        # approximate one, do it in descending order of frequency
        order = (self.common_lengths if distance != 0
                    else self.descending_lengths)

        for l in order:
            m = self.barcode_caches[l].search(barcode[:l],
                                                distance, unique=unique,
                                                details=True)
            if m != None:
                if m[1] == barcode[:l]:
                    # perfect match- unnecessary to check others
                    return m + (l, ) if details else m[0]
                matches.append(m + (l, ))

        if len(matches) == 0:
            return None
        elif len(matches) == 1:
            return matches[0] if details else matches[0][0]
        else:
            # find the best one
            closest = closest_match(barcode, [o for n, o, l in matches],
                                    unique=unique)
            if closest == None:
                return None
            best_match = [(n, o, l) for n, o, l in matches if o == closest]

            return best_match[0] if details else best_match[0][0]

    def search(self, barcode, distance, unique=False, details=False,
               verbose=False):
        """
        Find match for a barcode. First, check all possible lengths for an
        exact match, then all for an approximate match (since searching for
        an exact match is much, much faster than an approximate)
        """
        exact_m = self.search_all(barcode, 0, details=details)
        if exact_m != None:
            return exact_m
        ret = self.search_all(barcode, distance, unique=unique,
                                    details=details, verbose=verbose)
        return ret

    def mismatch_table(self):
        return "".join([c.mismatch_table()
                            for c in self.barcode_caches.values()])


class SampleCounter(object):
    """Count barcodes for one up/down tag within one sample"""
    def __init__(self, cache, name=None):
        """
        Given a BarcodeCache and a name.
        """
        self.name = name
        self.cache = cache
        self.data = collections.defaultdict(int)

    def add(self, barcode, dist, verbose=False):
        """
        Add a barcode to the dictionary, and return the object and the barcode
        it was matched to
        """
        matched = self.cache.search(barcode, dist, unique=True, details=True,
                                    verbose=verbose)
        if matched != None:
            self.data[matched[0]] += 1

        return matched


class BarcodeCounter(object):
    """Can count barcodes based on a barcode file"""
    def __init__(self, infile, upcode, downcode, multiplex_file=None,
                 track_mismatches=False):
        """
        Given a tab-delimited barcode file, of the format
        Strain  Uptag   Downtag
        Keep track of tags and multiplex barcodes using a BarcodeCache, and
        the barcodes themselves using a BarcodeCacheMultipleLen.
        """
        self.tagcache = BarcodeCache({upcode: 0, downcode: 1})

        self.total = 0
        self.total_found = 0

        inf = open(infile)
        uptags = {}
        downtags = {}
        self.original = [l[:-1].split("\t") for l in inf]
        for strain, uptag, downtag in self.original:
            uptags[uptag] = downtags[downtag] = strain

        inf.close()

        self.upcache = BarcodeCacheMultipleLen(uptags)
        self.downcache = BarcodeCacheMultipleLen(downtags)

        if multiplex_file != None:
            # create SampleCounters as a dictionary of 2-tuples
            inf = open(multiplex_file)
            counters = {}
            for l in inf:
                strain, barcode = l[:-1].split("\t")
                counters[barcode] = (SampleCounter(self.upcache,
                                                    strain + "_UP"),
                                     SampleCounter(self.downcache,
                                                    strain + "_DOWN"))
            self.counters = BarcodeCache(counters)
            inf.close()
            self.multiplexed = True
            self.ordered_counters = list(itertools.chain(*
                                            self.counters.strains))
        else:
            # single pair of SampleCounter
            self.counter = (SampleCounter(self.upcache),
                             SampleCounter(self.downcache))
            self.multiplexed = False
            self.ordered_counters = self.counter[:]

        if track_mismatches:
            # mismatches are indexed by tag, then by name, then by mismatched
            # barcode
            self.mismatches = [collections.defaultdict(lambda:
                                        collections.defaultdict(int))
                                    for i in range(2)]
        else:
            self.mismatches = None

    def add(self, barcode, tagcode, dist, multiplex_code=None, verbose=False):
        """Add a barcode to the appropriate tag"""
        if self.multiplexed == False and multiplex_code != None:
            raise ValueError("Cannot use multiplex_code, BarcodeCounter " +
                             "was not given a multiplexing file")
        elif self.multiplexed == True and multiplex_code == None:
            raise ValueError("Need multiplexed code")

        self.total += 1

        whichtag = self.tagcache.search(tagcode, 1)
        if whichtag == None:
            return False

        if self.multiplexed:
            counter = self.counters.search(multiplex_code, 1)
            if counter == None:
                return False
        else:
            counter = self.counter

        found = counter[whichtag].add(barcode, dist, verbose=verbose)
        if found != None:
            self.total_found += 1

            name, original, length = found

            # add to mismatch dictionary
            if self.mismatches:
                self.mismatches[whichtag][name][barcode[:length]] += 1

        return found

    def mismatch_table(self, outfile=None):
        """return a string describing all mismatches that have occured"""
        if self.mismatches == None:
            raise Exception("Cannot create mismatch file; BarcodeCounter " +
                            "was initialized without track_mismatches")

        ret = ""

        originals = [c.original for c in [self.upcache, self.downcache]]

        for mm_dict, original in zip(self.mismatches, originals):
            for n, bc_dict in mm_dict.items():
                o = original[n]
                ret += "\t".join(map(str, [n, o, bc_dict[o],
                            "/".join(["%s (%d)" % (k, v)
                                        for k, v in bc_dict.items()
                                            if k != o])])) + "\n"

        if outfile != None:
            with open(outfile, "w") as outf:
                outf.write(ret)

        return ret

    def revised_catalog(self, outfile=None):
        """Return a string or write a file with a revised barcode catalog"""
        if self.mismatches == None:
            raise Exception("Cannot create revised catalog file; " +
                            "BarcodeCounter was initialized without " +
                            "track_mismatches")

        originals = [c.original for c in [self.upcache, self.downcache]]

        bests = [dict([(n, max(bc_dict, key=bc_dict.get))
                            for n, bc_dict in mm_dict.items()])
                                for mm_dict, original_dict
                                    in zip(self.mismatches, originals)]

        ret = "\n".join(["\t".join((strain, bests[0].get(strain, up),
                            bests[1].get(strain, dn)))
                            for strain, up, dn in self.original])

        if outfile != None:
            with open(outfile, "w") as outf:
                outf.write(ret)

        return ret

    def report(self):
        """return a one-line description"""
        if self.total == 0:
            return "-"
        return "%d\t%.5f" % (self.total, float(self.total_found) / self.total)

    def write_file(self, outfile):
        """Write to a tab-delimited table"""
        outf = open(outfile, "w")
        outf.write(str(self))
        outf.close()

    def __repr__(self):
        if self.multiplexed:
            ret = "Strain\t" + "\t".join([c.name
                                    for c in self.ordered_counters]) + "\n"

            for s in self.upcache.strains:
                ret += s + "\t" + "\t".join([str(c.data[s])
                        for c in self.ordered_counters]) + "\n"
        else:
            ret = "Strain\tUP\tDOWN\n"
            for s in self.upcache.strains:
                ret += s + "\t%s\t%s\n" % (self.counter[0].data[s],
                                           self.counter[1].data[s])
        return ret
