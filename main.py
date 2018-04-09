 # -*- coding: utf-8 -*-
from sys import exit
import requests
import numpy as np
import h5py
import illustris_python as il

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

def test_merger_explore():
    test_tree = il.sublink.loadTree(BASE_PATH, 135, 15, fields='NextProgenitorID')
    print test_tree.tolist()
    # The NextProgenitorID field for the selected z=0 subhalo gives the most
    # massive galaxy that merged with the selected subhalo for each snapshot.
    # We fetch the sublink tree for that subhalo and look and the same field
    # to find the next most massive subhalo and repeat until we get the id -1
    # indicatiing that there are no more subhalos that merged with the selected
    # subhalo in the current snapshot.  Each of these IDs are recorded in a
    # list, so that list gives the IDs of all subhalos that merged with the
    # primary subhalo in a snapshot.  This process is repeated for each snapshot
    # that the primary existed in, and all of these lists are arranged into a
    # dictionary that has as keys the snapshot number each merger occurred in


if __name__=="__main__":
    test_merger_explore()
