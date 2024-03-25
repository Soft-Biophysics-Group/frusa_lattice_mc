import json
import os

def make_json_file(data,address):
    """
    Writes the data, provided as a python dictionary, to a .json file
    """
    with open(address,"w") as write_file:
        json.dump(data, write_file, indent=2)

def make_dir(name):
    """
    Creates a new directory if it doesn't already exist at the 
    specified address/name.
    
    Arguments:
    
    name - (str) name and location of the new directory.
    """
    
    if not os.path.exists(name):
        os.makedirs(name)

