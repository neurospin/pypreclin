##########################################################################
# NSAp - Copyright (C) CEA, 2013 - 2018
# Distributed under the terms of the CeCILL-B license, as published by
# the CEA-CNRS-INRIA. Refer to the LICENSE file or to
# http://www.cecill.info/licences/Licence_CeCILL-B_V1-en.html
# for details.
##########################################################################


"""
deprecated
"""

# System import
import subprocess
import fcntl
import os
import sys
import time

# Third party import
from filelock import FileLock


def non_block_read(output):
    """ Function to read a stream without blocking.
    """
    fd = output.fileno()
    fl = fcntl.fcntl(fd, fcntl.F_GETFL)
    fcntl.fcntl(fd, fcntl.F_SETFL, fl | os.O_NONBLOCK)
    try:
        return output.read()
    except:
        return ""


def search_align_windows():
    """ List the windows called align.
    
    Returns
    -------
    align_windows: list
        the align windows ids.
    """
    cmd = ["xdotool", "search", "--name", "align"]
    process = subprocess.Popen(
        cmd, env=os.environ, stdout=subprocess.PIPE, universal_newlines=True)
    align_windows = process.communicate()[0].split("\n")
    return align_windows


def send_key(key, win_id):
    """ Send a keyboard key to a windows.

    Parameters
    ----------
    key: str
        the keyboard event to be sent.
    win_id: int
        the window id to control.
    """
    cmd = ["xdotool", "windowactivate", "--sync", str(win_id), "key", key]
    subprocess.check_call(cmd)


def create_win(cmd, lock_file, default_windows):
    """ Create a JIP windows and return his id.

    Parameters
    ----------
    cmd: list of str
        the JIP command that will popup a window.
    lock_file: str
        the path to a file used to acquired a lock in order to identify
        properly the window id.
    default_windows: list
        the default windows ID.

    Returns
    -------
    process: object
        the JIP process.
    win_id: int
        the associated window id.
    """
    with FileLock(lock_file):
        process = subprocess.Popen(
                cmd, env=os.environ, stdout=subprocess.PIPE,
                stderr=subprocess.PIPE, universal_newlines=True)
        win_id = None
        while True:
            time.sleep(2)
            if win_id is None:
                new_windows = search_align_windows()
                active_windows = list(set(new_windows) - set(default_windows))
                assert len(active_windows) <= 1, "Parallel jobs not supported yet."
                if len(active_windows) == 1:
                    win_id = active_windows[0]
                    break
    return process, win_id


def call_auto(cmd, lock_file, refresh_winid=True):
    """ Execute a JIP command and control the graphical interface.

    Parameters
    ----------
    cmd: list of str
        the JIP command that will popup a window.
    lock_file: str
        the path to a file used to acquired a lock in order to identify
        properly the window id.
    refresh_winid: bool (optional, default True)
        if set check that only one JIP alignment is performed at a time. The
        windows ID is refreshed at each iteration.

    Returns
    -------
    mi: float
        the final optimizer cost function value (mutual information).
    """
    default_windows = search_align_windows()
    process, win_id = create_win(cmd, lock_file, default_windows)
    send_run = False
    while True:
        time.sleep(2)          
        if not send_run:
            send_key("a", win_id)
            send_run = True
        line = non_block_read(process.stdout)
        if line == "":
            continue
        if refresh_winid:
            new_windows = search_align_windows()
            active_windows = list(set(new_windows) - set(default_windows))
            assert len(active_windows) == 1, "Parallel jobs not supported yet."
            if active_windows[0] != win_id:
                print("Refreshing JIP windows ID {0} -> {1}.".format(
                    win_id, active_windows[0]))
                win_id = active_windows[0]
        if "MI =" in line:
            mi = float(line.split("=")[1])
            break
    send_key("q", win_id)
    send_key("y", win_id)
    stdout, stderr = process.communicate()

    return mi
