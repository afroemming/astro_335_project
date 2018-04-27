from os import mkdir,system, chdir 
from os.path import join, getcwd
from sys import exit, argv

try:
    SIM_NAME = argv[1]
except:
    print "Please supply a simulation name as a commandline argument (ex.: setup.py Illustris-3)"
    exit(1)

groupcat_URL = "http://www.illustris-project.org/api/" + SIM_NAME + \
               "/files/groupcat-135/?format=api" 
sublink_URL = "http://www.illustris-project.org/api/" + SIM_NAME + \
              "Illustris-3/files/sublink/?format=api" 

mkdir("trees")
mkdir("trees/SubLink")
mkdir("groups_135")

wget_path = 'wget'
if system(wget_path) != 256:
    wget_path = raw_input("Please type the full path to the wget excutable: \n")

chdir(join(getcwd(), "groups_135"))
system(wget_path + ' -nd -nc -nv -e robots=off -l 1 -r -A hdf5 --content-disposition --header="API-Key: 874387f0fbf3b68439b727ae5cdde978" ' + \
        groupcat_URL)
chdir(join('..', 'trees', 'SubLink'))
system(wget_path + '-nd -nc -nv -e robots=off -l 1 -r -A hdf5 --content-disposition --header="API-Key: 874387f0fbf3b68439b727ae5cdde978" ' + \
        sublink_URL)