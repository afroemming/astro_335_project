 # -*- coding: utf-8 -*-
from sys import exit
import requests
import numpy as np
import h5py
#import illustris_python as il
import matplotlib.pyplot as plt
from os.path import isfile
import os
import sys
import glob
import six

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
BASE_PATH = "E:/Astro 335 Project/Illustris-3"
TREE_NAME = "SubLink"

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
                               z = 0, limit=105):
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

def get_subhalo_merger_tree(id, fields=[]):
    """
    Given a sublink id and a list of fields, return a dictionary with keys
    snapshots where a merger occured and with values tuples with the ID and
    mass of the subhalo along the mpb, followed bythe subhalo ID and mass
    of each subhalo that merged with the mpb of the given subhalo at that
    snapshot.
    """
    # TODO: Determine tree a subhalo is in given its sublink ID
    chunkNum = return_tree_chunk_number(id)
    #print(chunkNum)
    #treeFile = "tree_extended.%d.hdf5" % chunkNum
    #path = os.path.join("E:", "Astro 335 Project", "Illustris-3", "trees", "SubLink", treeFile)
    f = h5py.File("E:/Astro 335 Project/Illustris-3/trees/SubLink/tree_extended.0.hdf5", 'r')
    #f = h5py.File(path, 'r')

    # Each index among the dataset goes with single subhalo. SubhaloID's are
    # contiguous, and the indexes are releated to SubhaloID's by the SubhaloID
    # of the first halo in the file + the index for a given subhalo
    first_sh_id = f['SubhaloID'][0]
    n = id
    # This dictionary will have as keys a Next Progenitor associated with a
    # a First Progenitor and
    mergers = {}
    first_progenitor_ID = f['FirstProgenitorID'][n]
    # Move through data file, following main branch until there earlist progenitor
    while f['FirstProgenitorID'][n] != -1:
        mergers[f['SnapNum'][n]] = [(n, f['Mass'][n])]
        # Move through next progenitor until we reach the last one
        while f['NextProgenitorID'][n] != -1:
            try:
                mergers[f['SnapNum'][n]].append((f['NextProgenitorID'][n], f['Mass'][f['NextProgenitorID'][n]]))
                n = f['NextProgenitorID'][n] - first_sh_id
            except:
                n = f['NextProgenitorID'][n] - first_sh_id
        # move along main branch
        n = f['FirstProgenitorID'][first_progenitor_ID] - first_sh_id
        first_progenitor_ID = n
    return mergers

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

def get_mass_hist_for_subhalo(id):

    ratio_mass_list = []
    mass_dict = get_merger_fractions(id)
    for snapNum in mass_dict:
        for merger in mass_dict[snapNum]:
            ratio_mass_list.append(merger)
    
    n, bins, patches = plt.hist(ratio_mass_list, bins = 10000, range = (0,.001), histtype = 'step')
    plt.show()

#print(get_subhalos_in_mass_range(MILKY_WAY_MASS, NEAR_MWM_DEVIATION))


#print(i_like_good_pussy_and_i_like_good_trees(1))
#print(get_subhalo_merger_tree(1))
#print(get_merger_fractions(1))



def treePath(basePath, treeName, chunkNum=0):
    """ Return absolute path to a SubLink HDF5 file (modify as needed). """
    # tree_path = '/trees/' + treeName + '/' + 'tree_extended.' + str(chunkNum) + '.hdf5'
    tree_path = os.path.join('trees', treeName, 'tree_extended.' + str(chunkNum) + '.hdf5')

    _path = os.path.join(basePath, tree_path)
    if len(glob.glob(_path)):
        return _path

    # new path scheme
    _path = os.path.join(basePath, os.path.pardir, 'postprocessing', tree_path)
    if len(glob.glob(_path)):
        return _path

    # try one or more alternative path schemes before failing
    _path = os.path.join(basePath, 'postprocessing', tree_path)
    if len(glob.glob(_path)):
        return _path

    raise ValueError("Could not construct treePath from basePath = '{}'".format(basePath))



offsetCache = dict()

def subLinkOffsets(basePath, treeName, cache=True):
    # create quick offset table for rows in the SubLink files
    if cache is True:
        cache = offsetCache

    if type(cache) is dict:
        path = os.path.join(basePath, treeName)
        try:
            return cache[path]
        except KeyError:
            pass

    search_path = treePath(basePath, treeName, '*')
    numTreeFiles = len(glob.glob(search_path))
    if numTreeFiles == 0:
        raise ValueError("No tree files found! for path '{}'".format(search_path))
    offsets = np.zeros(numTreeFiles, dtype='int64')

    for i in range(numTreeFiles-1):
        with h5py.File(treePath(basePath, treeName, i), 'r') as f:
            offsets[i+1] = offsets[i] + f['SubhaloID'].shape[0]

    if type(cache) is dict:
        cache[path] = offsets

    return offsets

offsets = subLinkOffsets(BASE_PATH, TREE_NAME)

#print(offsets)

def return_tree_chunk_number(id):
    for i in range(len(offsets)):
        if id >= offsets[i] and id < offsets[i+1]:
            return i
get_mass_hist_for_subhalo(2)
        

