"""
Unit test suite for BarNone.
"""

import os
import itertools
import unittest
import random
import shutil

from BarNone import matching
import flamingo

NUCLEOTIDES = "ACGT"

random.seed(1)


def random_barcode(l):
    return "".join([random.choice(NUCLEOTIDES) for i in range(l)])


def brute_force_levenshtein(w, lst):
    """find the closest levenshtein matches by brute force"""
    matches = [(e, flamingo.distance(w, e)) for e in lst]
    best_distance = min((d for w, d in matches))
    return [w for w, d in matches if d == best_distance]


def all_length(n):
    if n == 0:
        return [""]
    return list(itertools.chain(*[[s + n for n in NUCLEOTIDES]
                    for s in all_length(n - 1)]))


class TestFlamingo(unittest.TestCase):
    """Tests the flamingo WrapperSimpleEd"""
    def check_consistent(self, query, results, dist, msg=None):
        """check that a query is within the distance of all results"""
        for e in results:
            self.assertTrue(dist >= flamingo.distance(query, e),
                                    msg=msg)

    def check_search(self, indexer, query, expected, dist):
        ret = indexer.search(query, dist)
        self.assertEqual(ret, expected)
        self.check_consistent(query, ret, dist,
                        "Got expected array, but expected had incorrect value")

    def test_basic(self):
        """Basic, hardcoded tests"""
        words = ["apple", "banana", "orange"]
        indexer = flamingo.WrapperSimpleEd(words)
        for w in words:
            self.check_search(indexer, w, [w], 0)
            self.check_search(indexer, w, [w], 1)

        words = ["ABCDEFGHIJ"[:i] for i in range(1, 10)]
        indexer = flamingo.WrapperSimpleEd(words)
        for d in range(5):
            self.assertEqual(indexer.search(words[5], d), words[5 - d:6 + d])

    def test_barcodes(self):
        """
        Run on all barcodes of a certain length and confirm that it finds all
        the right ones
        """
        for l in range(1, 6):
            barcodes = all_length(l)
            indexer = flamingo.WrapperSimpleEd(barcodes)
            for b in barcodes:
                self.check_search(indexer, b, [b], 0)
                results = indexer.search(b, 1)
                self.assertEqual(len(results), 1 + 3 * l)

    def test_distance(self):
        """
        Test the method for finding the Levenshtein distance between
        strings
        """
        self.assertEqual(flamingo.distance("apple", "aaple"), 1)
        self.assertEqual(flamingo.distance("A", "AAA"), 2)
        self.assertEqual(flamingo.distance("apple", "applesauce"), 5)

#     def test_time_distance(self):
#         """Time the flamingo.distance function"""
#         import timeit
#         import Levenshtein
#
#         t = timeit.Timer("flamingo.distance('ACAGACTGATAGACAGATA', " +
#                                             "'ACAGCCTGATAGACATATA')",
#                          setup="import flamingo")
#         print t.timeit()
#
#         t = timeit.Timer("Levenshtein.distance('ACAGACTGATAGACAGATA', " +
#                                             "'ACAGCCTGATAGACATATA')",
#                          setup="import Levenshtein")
#         print t.timeit()


class TestBarcodeCache(unittest.TestCase):
    """Test the BarcodeCache object"""
    def test_barcode_cache(self):
        """basic tests of BarcodeCache"""
        # check the basic functionality: return best match
        words = ["ABCDEFGHIJ"[:i] for i in range(1, 11)]
        cache = matching.BarcodeCache(dict([(w, w) for w in words]))
        self.assertEqual(cache.search("ABCDE", 0), "ABCDE")
        self.assertEqual(cache.search("ABCDE", 1), "ABCDE")
        for j in range(1, 10):
            self.assertEqual(cache.search("ABGDE", j), "ABCDE")
        self.assertEqual(cache.search("ABCDEFGHIJK", 0), None)
        self.assertEqual(cache.search("ABCDEFGHIJK", 1), "ABCDEFGHIJ")

        words = all_length(5)
        cache = matching.BarcodeCache(dict([(w, w) for w in words]))

        for query, t_cached, t_looked_up in [("ACGAT", 1, 0),
                ("ACGA", 1, 1), ("ACGA", 2, 1), ("ACGATA", 2, 2),
                ("ACGATAA", 2, 2), ("ACGA", 3, 2)]:
            cache.search(query, 1)
            self.assertEqual(cache.total_cached, t_cached)
            self.assertEqual(cache.total_looked_up, t_looked_up)

    def test_barcode_cache_multiple_len(self):
        """basic tests of BarcodeCacheMultipleLen"""
        words = list(itertools.chain(*[[l * 5, l * 6]
                        for l in "ABCDEFG"]))

        cache = matching.BarcodeCacheMultipleLen(dict([(w, w) for w in words]))

        for d in range(0, 4):
            for l in "ABCDEFG":
                self.assertEqual(cache.search(l * 5, d), l * 5)
                self.assertEqual(cache.search(l * 6, d), l * 6)
        for d in range(1, 4):
            for l in "ABCDEFG":
                self.assertEqual(cache.search(l * 7, d), l * 6)
                self.assertEqual(cache.search(l * 8, d), l * 6)

        # check negative results as well
        for i in range(10):
            self.assertEqual(cache.search("T" * i, 1), None)

        # create and test a random, unique catalog
        LENGTH = 9
        catalog = list(set([random_barcode(LENGTH) for i in range(2000)]))

        cache = matching.BarcodeCacheMultipleLen(dict([(w, w)
                                                    for w in catalog]))
        for w in catalog:
            for d in range(4):
                self.assertEqual(cache.search(w, d), w)

        # pick random barcodes and see that they match correctly
        for i in range(250):
            b = random_barcode(LENGTH)
            closest = brute_force_levenshtein(b, catalog)
            dists = [flamingo.distance(b, m) for m in closest]
            # sanity check on brute_force_levenshtein:
            self.assertEqual(len(set(dists)), 1)

            if len(closest) > 1:
                self.assertEqual(cache.search(b, dists[0], unique=True), None)
            else:
                # make sure the cache finds it
                self.assertEqual(cache.search(b, dists[0], unique=True),
                                    closest[0])
                self.assertEqual(cache.search(b, dists[0] + 1, unique=True),
                                    closest[0])

            self.assertTrue(cache.search(b, dists[0], unique=False) in closest)

            # make sure it can't find anything closer
            for d in range(dists[0]):
                self.assertEqual(cache.search(b, d), None)

        # test an unusual case- close to ones of multiple length
        words = ["BAAAA", "BAAAAA"]
        cache = matching.BarcodeCacheMultipleLen(dict([(w, w) for w in words]))
        self.assertEqual(cache.search("AAAAA", 1), "BAAAA")
        self.assertEqual(cache.search("AAAAA", 1, details=True),
                            ("BAAAA", "BAAAA", 5))


class TestBarcodeCounter(unittest.TestCase):
    """
    Test the BarcodeCounter
    """
    def setUp(self):
        """Create a test directory to work in"""
        self.test_directory = "test_directory"
        if not os.path.exists(self.test_directory):
            os.mkdir(self.test_directory)

        self.original_directory = os.getcwd()
        self.test_file = "test_barcodes.txt"

        os.chdir(self.test_directory)

    def tearDown(self):
        """Go back to original directory and remove test directory"""
        os.chdir(self.original_directory)
        shutil.rmtree(self.test_directory)

    def test_barcode_counter(self):
        """Simple tests of barcode counter"""
        names = list("0123456789")
        uptags = list(itertools.chain(*[(l * 5, l * 6) for l in "ABCDE"]))
        dntags = list(itertools.chain(*[(l * 5, l * 6) for l in "VWXYZ"]))

        with open(self.test_file, "w") as outf:
            for v in zip(names, uptags, dntags):
                outf.write("\t".join(v) + "\n")

        # make sure it's not counting mismatches unless we specify
        counter = matching.BarcodeCounter(self.test_file, "UPT", "DNT")
        self.assertEqual(counter.mismatches, None)

        counter = matching.BarcodeCounter(self.test_file, "UPT", "DNT",
                                         track_mismatches=True)

        for n, u, d in zip(names, uptags, dntags):
            self.assertEqual(counter.mismatches[0][n][u], 0)
            # add one uptag
            self.assertEqual(counter.add(u, "UPT", 0), (n, u, len(u)))
            self.assertEqual(counter.mismatches[0][n][u], 1)

            self.assertEqual(counter.add(u + "Q", "UPT", 1), (n, u, len(u)))

            # add close ones
            if len(u) == 5:
                self.assertEqual(counter.add(u[:4], "UPT", 1), (n, u, 5))
                self.assertEqual(counter.mismatches[0][n][u[:4]], 1)

        # basic checks of mismatch table
        mm_table = [l[:-1].split("\t")
                     for l in counter.mismatch_table().split("\n") if l != ""]
        self.assertEqual(len(mm_table), 10)

        for spl in mm_table:
            self.assertTrue(spl[0] in names)
            self.assertTrue(spl[1] in uptags)
            self.assertEqual(spl[2], "2")

        # another case, since there are some that can't be tested with the
        # above cache
        words = ["apple", "banana", "orange", "bacon", "tomato", "lettuce"]
        with open("test_barcodes.txt", "w") as outf:
            for w in words:
                outf.write("\t".join(["_" + w, w, w]) + "\n")

        counter = matching.BarcodeCounter("test_barcodes.txt", "UPT", "DNT",
                                         track_mismatches=True)

        for w in words:
            self.assertEqual(counter.add(w, "UPT", 1), ("_" + w, w, len(w)))
            self.assertEqual(counter.add(w + "Q", "UPT", 1),
                                ("_" + w, w, len(w)))
            self.assertEqual(counter.add("Q" + w, "UPT", 2),
                                ("_" + w, w, len(w)))
            self.assertEqual(counter.add("Q" + w[1:], "UPT", 1),
                                ("_" + w, w, len(w)))

        mm_table = [l[:-1].split("\t")
                     for l in counter.mismatch_table().split("\n") if l != ""]

        self.assertEqual(len(mm_table), 6)

        for spl in mm_table:
            self.assertEqual(spl[2], "2")
            self.assertEqual(len(spl[3].split("/")), 2)

        self.assertTrue("_apple\tapple\tapple" in
                            counter.revised_catalog().split("\n"))

        for i in range(10):
            self.assertEqual(counter.add("aaple", "UPT", 1),
                                ("_apple", "apple", 5))

        self.assertTrue("_apple\taaple\tapple" in
                            counter.revised_catalog().split("\n"))


def run_tests():
    unittest.main()


if __name__ == "__main__":
    run_tests()
