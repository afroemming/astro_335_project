 # -*- coding: utf-8 -*-
from sys import exit
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
    snapshots where a merger occured and with values the subhalo ID and requested
    fields of each subhalo that merged with the mpb of the given subhalo at that
    snapshot.
    """
    # TODO: Determine tree a subhalo is in given its sublink ID
    # TODO: Given a list of fields,return those fields associated with
    # each

    f = h5py.File("trees/SubLink/tree_extended.0.hdf5", 'r')
    # Each index among the dataset goes with single subhalo. SubhaloID's are
    # contiguous, and the indexes are releated to SubhaloID's by the SubhaloID
    # of the first halo in the file + the index for a given subhalo
    first_sh_id = f['SubhaloID'][0]
    # n is would be the input for a function that takes this code
    n = 0
    # This dictionary will have as keys a Next Progenitor associated with a
    # a First Progenitor and
    mergers = {}
    first_progenitor_ID = f['FirstProgenitorID'][n]
    # Move through data file, following main branch until there earlist progenitor
    while f['FirstProgenitorID'][n] != -1:
        # Move through next progenitor until we reach the last one
        while f['NextProgenitorID'][n] != -1:
            # Try to add ID to dictionary, create list associated with ID if
            # key does not exist
            try:
                mergers[f['SnapNum']][n].append(f['NextProgenitorID'][n])
            except:
                print n
                mergers[f['SnapNum'][n]] = [f['NextProgenitorID'][n]]
            # move to next next progenitor
            n = f['NextProgenitorID'][n] - first_sh_id
        # move along main branch
        n = f['FirstProgenitorID'][first_progenitor_ID] - first_sh_id
        first_progenitor_ID = n
    return mergers
    