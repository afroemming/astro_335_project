 # -*- coding: utf-8 -*-
from sys import exit, argv
import requests
import numpy as np
import h5py
import illustris_python as il
import matplotlib.pyplot as plt

MILKY_WAY_MASS = 1.30e12 / 1e10 * 0.704 # Solar Masses. 1.30e12 +/- 0.30e12
# Source: McMillan, Paul J. (February 11, 2017). "The mass distribution and
# gravitational potential of the Milky Way". Monthly Notices of the Royal
#Astronomical Society. 465 (1): 76â€“94. arXiv:1608.00971â€¯

# We should determine what 'near milky way mass' means. For now, I've used one
# error away from the Milky Ways mass.
NEAR_MWM_DEVIATION = 0.30e12  / 1e10 * 0.704
API_KEY='874387f0fbf3b68439b727ae5cdde978'
BASE_URL = 'http://www.illustris-project.org/api/'
SIM_NAME = 'Illustris-3' #We'll probably want to use Illustris-1, but that
# dataset is superbig, so I'll use Illustris-3 for now
BASE_PATH = "."

if __name__ == "__main__":
    try: 
        SIM_NAME = argv[1]
    except IndexError:
        pass

def get(path, params=None):
    """Helper function to access API endpoints"""
    headers = {'api-key':API_KEY}
    # Make HTTP GET request to path
    r = requests.get(path, params=params, headers=headers)

    # Raise exception of response code is not HTTP SUCCESS (200)
    r.raise_for_status()

    if r.headers['content-type'] == 'application/json':
        return r.json() # Parse json responses automatically

    if 'content-disposition' in r.headers:
        filename = r.headers['content-disposition'].split("filename=")[1]
        with open(filename, 'wb') as f:
            f.write(r.content)
        return filename

    return r

def get_subhalos_in_mass_range(target_mass, deviation,
                               z = 0, limit=105, SIM_NAME=SIM_NAME):
    """
    Given a mass and deviation in group catalog units, return a search results for
    subhalos within that for a optionially specified z.
    """
    search_query = "?mass__gt=" + str(target_mass - deviation) + "&mass__lt=" +  \
                  str(target_mass + deviation)
    url = "http://www.illustris-project.org/api/" + '/' + SIM_NAME + '/' + \
    "snapshots/z=" + str(z) + "/subhalos/" + search_query + "&limit=" + str(limit)

    subhalos = get(url)

    return subhalos

def get_parent_subhalos_in_range(target_mass, deviations, z = 0, limit=105, SIM_NAME=SIM_NAME):
    parent_subhalos = []
    for subhalo in get_subhalos_in_mass_range['results']:
        if get(subhalo['url'])['parent'] ==  0:
            parent_subhalos.append(subhalo)

def get_subhalo_merger_tree(id, fields=[]):
    """
    Given a sublink id and a list of fields, return a dictionary with keys
    snapshots where a merger occured and with values tuples with the ID and
    mass of the subhalo along the mpb, followed bythe subhalo ID and mass
    of each subhalo that merged with the mpb of the given subhalo at that
    snapshot.
    """
    # TODO: Determine tree a subhalo is in given its sublink ID

    f = il.sublink.loadTree(BASE_PATH, 135, id)
    id0 = f['SubhaloID'][0] # ID of the first subhalo in the file, subtract
    # this from an ID to index in tree
    n = 0
    mpb_n = 0
    mergers = {}
    #Walk through tree until FirstProgenitorID == -1, i.e., we've reached the
    #earliest step where the subhalos mpb exists
    while mpb_n != -1:
        print(mpb_n)
        # Store mass and id of current subhalo along mpb
        if not f['SnapNum'][mpb_n] in mergers.keys():
            mergers[f['SnapNum'][mpb_n]] = [(f['SubhaloID'][mpb_n], f['Mass'][mpb_n])]
        else:
            mergers[f['SnapNum'][mpb_n]].insert(0, (f['SubhaloID'][mpb_n], f['Mass'][mpb_n]))

        #Find all next progenitors
        n = f['NextProgenitorID'][mpb_n]
        if n != -1:
            n -= id0
        while n != -1:
            if f['SnapNum'][n] in mergers.keys():
                mergers[f['SnapNum'][n]].append((f['SubhaloID'][n], f['Mass'][n]))
            else:
                mergers[f['SnapNum'][n]] = [(f['SubhaloID'][n], f['Mass'][n])]
            n = f['NextProgenitorID'][n]
            if n != -1:
                n -= id0
        # Move to the next subhalo along the mpb
        mpb_n = f['FirstProgenitorID'][mpb_n]
        if mpb_n != -1:
            mpb_n -= id0
    return mergers

#print(get_subhalo_merger_tree(900))

def get_merger_fractions(id):
    """ 
    Given a subhalo id, return a dict with keys the SnapNum where the 
    subhalo existed, and as values, the merger fractions of each merger
    """
    mtree = get_subhalo_merger_tree(id)
    mass_mtree = {}
    for snapNum in mtree:
        mass_mtree[snapNum] = []
        for merger in mtree[snapNum][1:]:
            mass_mtree[snapNum].append(merger[1]/mtree[snapNum][0][1])

    return mass_mtree

def test_get_merger_fraction():
    for i in range(600,1200):
        a = get_merger_fractions(i)
        for j in a:
            if len(a[j]) != 0:
                print(i, a[j])