# -*- coding: utf-8 -*-
"""
A cython library to filter pairs by genome distance fast.
"""


def distance_filter(list1, list2, distance):
    """Filter lists by distance.

    All items in list2 are compared to list1 by the distance specified in
    distance.

    Note: All SNPs must be on the same chromosome.

    Parameters
    ----------
    list1/list2 : list
        The input lists must a list of tuples in the format (name, location)
    distance : int
        Integer number of base pairs to filter distance with

    Returns
    -------
    list_of_tuples
        A list of all SNPs in list1 with overlapping SNPs in list2:
            [(snp1, location, {(snp2, loc2), (snp3, loc3), ...}), ...]
    """
    list1 = list(list1)
    list2 = list(list2)
    if not isinstance(list1[0], tuple):
        raise TypeError("SNP lists must be lists of tuples, list1 is '{}'"
                        .format(type(list1[0])))
    if not isinstance(list2[0], tuple):
        raise TypeError("SNP lists must be lists of tuples, list2 is '{}'"
                        .format(type(list1[0])))
    return _filter_pairs(list(list1), list(list2), int(distance))


cdef list _filter_pairs(list list1, list list2, int distance):
    """Private cythonized function to do calculation, see filter_pairs()."""
    cdef list results = []
    cdef list matches = []
    for snp1, loc1 in list1:
        matches = []
        count = len(list2)
        for snp2, loc2 in list2:
            if abs(loc1-loc2) <= distance:
                matches.append((snp2, loc2))
        results.append((snp1, loc1, set(matches)))
    return results
