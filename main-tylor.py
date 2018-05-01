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
SIM_NAME = 'Illustris-1' #We'll probably want to use Illustris-1, but that
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
                               z, limit, SIM_NAME=SIM_NAME):
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

def get_parent_subhalos_in_range(target_mass, deviations, z, limit, SIM_NAME=SIM_NAME):
    parent_subhalos = []
    for subhalo in get_subhalos_in_mass_range(target_mass, deviations, z, limit, SIM_NAME=SIM_NAME)['results']:
        if get(subhalo['url'])['parent'] ==  0:
            parent_subhalos.append(subhalo)
    return parent_subhalos

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
        #print(mpb_n)
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

def get_merger_fractions_of_subhalos_in_mass_range(target_mass, deviation, z=0, limit=105):
    subhalo_id_list=[]
    subhalo_dict = get_parent_subhalos_in_range(target_mass, deviation, z, limit)
    merger_fraction_list = []
    for subhalo in subhalo_dict:
        subhalo_id_list.append(subhalo[u'id'])
    for subhalo_id in subhalo_id_list:
        for snapshot in get_merger_fractions(subhalo_id):
            for mass in get_merger_fractions(subhalo_id)[snapshot]:
                merger_fraction_list.append(mass)
    return merger_fraction_list

def merger_fractions_at_redshift_galaxy_search(target_mass, deviation, z, limit):

    merger_dictionary = {x+1: [] for x in range(135)}
    subhalo_id_list = []
    subhalo_dict = get_parent_subhalos_in_range(target_mass, deviation, z, limit)

    for subhalo in subhalo_dict:
        subhalo_id_list.append(subhalo[u'id'])
    
    
    for subhalo_id in subhalo_id_list:
        merger_fractions = get_merger_fractions(subhalo_id)
        for snapshot in merger_fractions:
            for mass in merger_fractions[snapshot]:
                merger_dictionary[snapshot].append(mass)
    return merger_dictionary

def major_merger_dict(target_mass, deviation, z, limit):
    merger_dictionary = merger_fractions_at_redshift_galaxy_search(target_mass, deviation, z, limit)
    major_merger_dictionary = {x+1: [] for x in range(135)}
    for snapshot in merger_dictionary:
        for mass in merger_dictionary[snapshot]:
            if mass > 0.25:
                major_merger_dictionary[snapshot].append(mass)
    return major_merger_dictionary

def minor_merger_dict(target_mass, deviation, z, limit):
    merger_dictionary = merger_fractions_at_redshift_galaxy_search(target_mass, deviation, z, limit)
    minor_merger_dictionary = {x+1: [] for x in range(135)}
    for snapshot in merger_dictionary:
        for mass in merger_dictionary[snapshot]:
            if mass <= 0.25 and mass > 0.1 :
                minor_merger_dictionary[snapshot].append(mass)
    return minor_merger_dictionary

def very_minor_merger_dict(target_mass, deviation, z, limit):
    merger_dictionary = merger_fractions_at_redshift_galaxy_search(target_mass, deviation, z, limit)
    very_minor_merger_dictionary = {x+1: [] for x in range(135)}
    for snapshot in merger_dictionary:
        for mass in merger_dictionary[snapshot]:
            if mass <= 0.1:
                very_minor_merger_dictionary[snapshot].append(mass)
    return very_minor_merger_dictionary

def fractional_major_mergers_dict(target_mass, deviation, z, limit):
    merger_dictionary = merger_fractions_at_redshift_galaxy_search(target_mass, deviation, z, limit)
    major_merger_dictionary = major_merger_dict(target_mass, deviation, z, limit)
    major_merger_fractional_dictionary = {x+1: 0.0 for x in range(135)}
    total_merger_dictionary = {x+1: [] for x in range(135)}
    for snapshot in merger_dictionary:
        for mass in merger_dictionary[snapshot]:
            total_merger_dictionary[snapshot].append(mass)
    for snapshot in major_merger_dictionary:
        if len(major_merger_dictionary[snapshot]) != 0 and len(total_merger_dictionary[snapshot]) != 0:
            fraction = float(len(major_merger_dictionary[snapshot])) / float(len(total_merger_dictionary[snapshot]))
            major_merger_fractional_dictionary[snapshot] = fraction
    return major_merger_fractional_dictionary

def fractional_minor_mergers_dict(target_mass, deviation, z, limit):
    merger_dictionary = merger_fractions_at_redshift_galaxy_search(target_mass, deviation, z, limit)
    minor_merger_dictionary = minor_merger_dict(target_mass, deviation, z, limit)
    minor_merger_fractional_dictionary = {x+1: 0.0 for x in range(135)}
    total_merger_dictionary = {x+1: [] for x in range(135)}
    for snapshot in merger_dictionary:
        for mass in merger_dictionary[snapshot]:
            total_merger_dictionary[snapshot].append(mass)
    for snapshot in minor_merger_dictionary:
        if len(minor_merger_dictionary[snapshot]) != 0 and len(total_merger_dictionary[snapshot]) != 0:
            fraction = float(len(minor_merger_dictionary[snapshot])) / float(len(total_merger_dictionary[snapshot]))
            minor_merger_fractional_dictionary[snapshot] = fraction
    return minor_merger_fractional_dictionary

def fractional_very_minor_mergers_dict(target_mass, deviation, z, limit):
    merger_dictionary = merger_fractions_at_redshift_galaxy_search(target_mass, deviation, z, limit)
    very_minor_merger_dictionary = very_minor_merger_dict(target_mass, deviation, z, limit)
    very_minor_merger_fractional_dictionary = {x+1: 0.0 for x in range(135)}
    total_merger_dictionary = {x+1: [] for x in range(135)}
    for snapshot in merger_dictionary:
        for mass in merger_dictionary[snapshot]:
            total_merger_dictionary[snapshot].append(mass)
    for snapshot in very_minor_merger_dictionary:
        if len(very_minor_merger_dictionary[snapshot]) != 0 and len(total_merger_dictionary[snapshot]) != 0:
            fraction = float(len(very_minor_merger_dictionary[snapshot])) / float(len(total_merger_dictionary[snapshot]))
            very_minor_merger_fractional_dictionary[snapshot] = fraction
    
    return very_minor_merger_fractional_dictionary
    

#print(fractional_very_minor_mergers_dict(MILKY_WAY_MASS, NEAR_MWM_DEVIATION, 0, 50))

def plot_of_merger_fraction_categories_vs_snapshot(target_mass, deviation, z, limit):
    very_minor_merger_list = []
    minor_merger_list = []    
    major_merger_list = []
    major_merger_fractional_dict = fractional_major_mergers_dict(target_mass, deviation, z, limit)
    minor_merger_fractional_dict = fractional_minor_mergers_dict(target_mass, deviation, z, limit)
    very_minor_merger_fractional_dict = fractional_very_minor_mergers_dict(target_mass, deviation, z, limit)
    # snapshot_list = range(1, 136)
    # snapshot_list.remove(53)
    # snapshot_list.remove(55)
    sup = get("http://www.illustris-project.org/api/Illustris-1/snapshots/")
    redshifts = []
    for i in range(1,134):
            redshifts.append(sup[i]["redshift"])
    redshifts = redshifts[11:]
    for snapshot in major_merger_fractional_dict:
        if snapshot != 53 and snapshot != 55:
            major_merger_list.append(major_merger_fractional_dict[snapshot])
    major_merger_list = major_merger_list[11:]
    for snapshot in minor_merger_fractional_dict:
        if snapshot != 53 and snapshot != 55:
            minor_merger_list.append(minor_merger_fractional_dict[snapshot])
    minor_merger_list = minor_merger_list[11:]
    for snapshot in very_minor_merger_fractional_dict:
        if snapshot != 53 and snapshot != 55:
            very_minor_merger_list.append(very_minor_merger_fractional_dict[snapshot])
    very_minor_merger_list = very_minor_merger_list[11:]
    fig, ax = plt.subplots(1)
    plt.rcParams['font.family'] = 'Georgia'
    ax.plot(redshifts, major_merger_list, color = '#eaadea', label='Fraction of Major Mergers')
    ax.plot(redshifts, minor_merger_list, color='#add8e6', label='Fraction of Minor Mergers')
    ax.set_title('Fraction of Mergers vs Redshift', fontsize=12, fontname = 'Georgia')
    ax.plot(redshifts, very_minor_merger_list, color="burlywood", label='Fraction of Very Minor Mergers')
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    ax.spines['left'].set_linewidth(.5)
    ax.tick_params(labelsize = 8, width = .5)
    ax.set_xlabel('Redshift', fontsize=12, fontname = 'Georgia')
    ax.spines['bottom'].set_linewidth(.5)
    box = ax.get_position()
    ax.set_position([box.x0, box.y0, box.width * 0.8, box.height])
    ax.legend(loc='center left', bbox_to_anchor=(.75,.5), fancybox=True, shadow=True)
    fig.savefig('merger fraction vs redshift.png', bbox_inches = 'tight')
    #plt.show()

def histogram_of_mergers_vs_mass(target_mass, deviation, z, limit):
    merger_dictionary = merger_fractions_at_redshift_galaxy_search(target_mass, deviation, z, limit)
    merger_mass_list = []
    for snapshot in merger_dictionary:
        for mass in merger_dictionary[snapshot]:
            merger_mass_list.append(mass)
    n, bins, patches = plt.hist(merger_mass_list, bins = 10000, range = (0,.1), histtype = 'step', color = '#03c03c', label = 'Number of Mergers Per Mass Ratio')
    plt.title('Number of Mergers Per Mass Ratio', fontsize = 12, fontname = "Georgia")
    plt.rcParams['font.family'] = 'Georgia'
    plt.xlabel('Merger Mass Ratio', fontsize = 12, fontname = 'Georgia')
    plt.ylabel('Number of Mergers', fontsize = 12, fontname = 'Georgia')
    ax = plt.gca()
    ax.set_facecolor('#e5e5e5')
    plt.yscale('log', nonposy='clip')
    plt.savefig('merger mass histogram.png', bbox_inches = 'tight')
    #plt.show()


def plot_total_mergers_vs_snapshot(target_mass, deviation, z, limit):
    merger_fractions = merger_fractions_at_redshift_galaxy_search(target_mass, deviation, z, limit)
    number_of_mergers = []
    normalized_mergers =[]
    for snapshot in merger_fractions:
        if snapshot != 53 and snapshot != 55:
            number_of_mergers.append(len(merger_fractions[snapshot]))
    normalize = float(max(number_of_mergers))
    for number in number_of_mergers:
        i = float(number)/normalize
        normalized_mergers.append(i)
    sup = get("http://www.illustris-project.org/api/Illustris-1/snapshots/")
    redshifts = []
    for i in range(1,134):
            redshifts.append(sup[i]["redshift"])
    fig, ax = plt.subplots(1,1)
    ax.plot(redshifts, normalized_mergers, color = '#fbb149')
    ax.set_xlabel('Redshift', fontsize = 12, fontname = 'Georgia')
    ax.set_title('Number of Mergers vs Redshift', fontsize = 12, fontname = 'Georgia')
    ax.set_ylabel('Number of Mergers', fontsize = 12, fontname = 'Georgia')
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    ax.spines['left'].set_linewidth(.5)
    ax.spines['bottom'].set_linewidth(.5)
    ax.tick_params(labelsize = 8, width = .5)
    #plt.show()
    fig.savefig('mergers vs redshift.png', bbox_inches='tight')
plot_total_mergers_vs_snapshot(MILKY_WAY_MASS, NEAR_MWM_DEVIATION, 0, 1000)
histogram_of_mergers_vs_mass(MILKY_WAY_MASS, NEAR_MWM_DEVIATION, 0, 1000)
plot_of_merger_fraction_categories_vs_snapshot(MILKY_WAY_MASS, NEAR_MWM_DEVIATION, 0, 1000)





