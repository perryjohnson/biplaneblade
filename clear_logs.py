"""Delete the contents of all the log files"""


def clear_log(path):
    """Deletes the contents of a log file located at 'path'."""
    f = open(path, 'w')
    f.close()

clear_log('blade.log')
print " Cleared 'blade.log'"
clear_log('station.log')
print " Cleared 'station.log'"
