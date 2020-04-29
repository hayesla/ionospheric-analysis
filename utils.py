import urllib
from bs4 import BeautifulSoup 

def list_path_files(url_path, file_ext):
    """
    Function to list all urls available given a certain path
    and search and return files with a given file extension.

    Parameters
    ----------
    url_path : ~str
        url path to search for urls 
    file_ext : ~str
        file extension, e.g. fits, txt, csv etc

    Returns
    -------
    links : ~list
        list of the files found

    Example
    -------
    >>> base_url = 'http://data.rosseobservatory.ie/data/2014/01/01/SID/txt/'
    >>> files = list_path_files(base_url, 'txt')
    >>> print(files)
    ... ['http://data.rosseobservatory.ie/data/2014/01/01/SID/txt/birr_SID_20140101_000000.txt']

    """
    test = urllib.request.urlopen(url_path)
    soup = BeautifulSoup(test, features="lxml")

    links = []
    for link in soup.findAll('a'):
        if link.get('href') is not None and link.get('href').find(file_ext) != -1:
            links.append(url_path + link.get('href').split('/')[-1]) 

    return links


