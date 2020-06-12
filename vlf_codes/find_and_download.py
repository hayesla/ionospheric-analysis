import urllib
import sys
from pathlib import Path
sys.path.append("..")
import utils
from sunpy.time import parse_time

def find_sid_files(date):
    """
    Function to list available files for Birr magnetometer data. The filename
    depends on the start seconds of the data which can change for certain days
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
    base_url = "http://data.rosseobservatory.ie/data/"

    date_base_url = base_url + date.strftime("%Y/%m/%d/sid/txt/")

    files = utils.list_path_files(date_base_url, "txt")

    return files

def search_and_download_vlf(date, path=None):

    date = parse_time(date)

    files = find_sid_files(date)
    if len(files) < 0:
        return "no files found"
    else:
        print("downloading now..")

    if path is not None:
        filename = Path(path).joinpath(files[0].split("/")[-1])
    else:
        filename = Path(files[0].split("/")[-1])
    
    if filename.exists():
        return "file already exists!"
    try:
        urllib.request.urlretrieve(files[0], filename=filename)
        if Path(filename).exists():
            return "files downloaded"
    except:
        return "error in download"