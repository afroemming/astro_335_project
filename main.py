from sys import exit
import requests
import numpy as np

MILKY_WAY_MASS = 1.30e12 / 1e10 * 0.704 # Solar Masses. 1.30e12 +/- 0.30e12
# Source: McMillan, Paul J. (February 11, 2017). "The mass distribution and 
# gravitational potential of the Milky Way". Monthly Notices of the Royal 
#Astronomical Society. 465 (1): 76â€“94. arXiv:1608.00971â€¯

# We should determine what 'near milky way mass' means. For now, I've used one 
# error away from the Milky Ways mass.
NEAR_MWM_DEVIATION = 0.30e12  / 1e10 * 0.704

API_KEY='874387f0fbf3b68439b727ae5cdde978'

BASE_URL = 'http://www.illustris-project.org/api/'
headers = {'api-key':API_KEY}

def get(path, params=None):
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

def get_subhalos_in_mass_range(target_mass, deviation, run_name = 'Illustris-1',
                               z = 0, limit=105):
  """ 
  Given a mass and deviation in group catalog units, return a search results for
  subhalos within that for a optionially specified run and z.
  """
  search_query = "?mass__gt=" + str(target_mass - deviation) + "&mass__lt=" +  \
                  str(target_mass + deviation)
  print(search_query)
  url = "http://www.illustris-project.org/api/" + '/' + run_name + '/' + \
  "snapshots/z=" + str(z) + "/subhalos/" + search_query + "&limit=" + str(limit)
  
  subhalos = get(url)
  
  return subhalos

test_halo = get(get_subhalos_in_mass_range(MILKY_WAY_MASS, NEAR_MWM_DEVIATION)['results'][0]['url'])

import h5py

mpb1 = get(test_halo['trees']['sublink_mpb'] )
f = h5py.File(mpb1,'r')
print(f.keys())