#!/usr/bin/env python3

# Helpful functions for the io modules

def notfound(err):
    print(err)
    raise FileNotFoundError

def any_to_filelist(input_data=None):
    """ 
    Take input_data and convert it to a list of file paths
    
    Parameters
    ----------
    input_data (str, list-like, path): raw user input. 
        Can be str, list, str with *'s, etc. This is what we will work on
    
    Returns
    -------
    list of strings: Hopefully the files that the user was trying to input...

    Notes
    -----
    - Will error on unexpected data
    - Will error if no files can be found
    - Can handle:
        - list to one or multiple files
        - str with '*'
        - str to single file
        - str without '*', like "path/to/3DALL"
        - str to directory, we take all '*.bin' files there

    """
    
    from glob import glob
    import os

    if isinstance(input_data, (str, os.PathLike)):
        # Single file that exists:
        if os.path.isfile(input_data):
            return [input_data]
        # directory that exists, return all .bin files inside
        elif os.path.isdir(input_data):
            return sorted(glob(os.path.join(input_data, "*.bin")))
        # User gave us a string to glob:
        elif '*' in input_data:
            outfiles = sorted(glob(input_data))
            # make sure the files actually exist
            if len(outfiles) > 0:
                return outfiles
            else:
                notfound("glob pattern in input_data, but no files found\n"
                         f"Received:\n\t{input_data}")
        # maybe the user wants us to glob for them?
        else:
            # directory is already handled, it must be a path or something
            outfiles = sorted(glob(input_data + "*.bin"))
            if len(outfiles) > 0:
                return outfiles
            else:
                notfound(f"No *.bin files found within:\n\t{input_data}")

    # we were probably given a list of files already
    else:
        return input_data