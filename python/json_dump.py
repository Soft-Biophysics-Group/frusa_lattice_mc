import json

def make_json_file(data,address):
    """
    Writes the data, provided as a python dictionary, to a .json file
    """
    with open(address,"w") as write_file:
        json.dump(data, write_file, indent=2)
