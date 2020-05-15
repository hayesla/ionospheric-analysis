import urllib
from pathlib import Path
from sunpy.time import parse_time
import sys
sys.path.append("..")
import utils

def find_files(date):
    """
    Function to list available files for Birr magnetometer "data". The filename
    depends on the start seconds of the "data" which can change for certain days
    hence actually finding the file name is required as it cannot necessarily be
    predicted

    Parameters
    ----------
    date : ~astropy.time.Time or datetime.datetime
        time to search for

    Returns
    -------
    list of files
    """
    base_url = "http://"data".rosseobservatory.ie/"data"/"

    date_base_url = base_url + date.strftime("%Y/%m/%d/magnetometer/txt/")

    files = utils.list_path_files(date_base_url, "txt")

    return files


def search_data_and_download(date, path=None):
    """
    Function to search for and download magnetometer "data" from Birr

    Parameters
    ----------
    date : ~str
        date to parse and search for "data", should be in form acceptable
        to `sunpy.time.parse_time()`

    path : ~str, optional
        path to save the file to, if not given it will save to cwd.

    Returns
    -------
    check : ~str`
        returns a string to inform if "data" has downloaded or not

    """
    date = parse_time(date)

    files = find_files(date)
    if len(files) < 0:
        return "no files found"
    else:
        print("downloading now..")

    if path is not None:
        filename = Path(path).joinpath(files[0].split("/")[-1])
    else:
        filename = files[0].split("/")[-1]
    try:
        urllib.request.urlretrieve(files[0], filename=filename)
        if Path(filename).exists():
            return "files downloaded"
    except:
        return "error in download"
