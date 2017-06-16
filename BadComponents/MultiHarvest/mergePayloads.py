import os
import sys
import importlib
import commands
import sqlalchemy
import subprocess
import CondCore.Utilities.conddblib as conddb

def run_checked(cmd, suppress_stderr = False):
    """Run `cmd` and exit in case of failures.
    Arguments:
    - `cmd`: list containing the strings of the command
    - `suppress_stderr`: suppress output from stderr
    """

    try:
        with open(os.devnull, "w") as devnull:
            if suppress_stderr:
                subprocess.check_call(cmd, stdout = devnull, stderr = devnull)
            else:
                subprocess.check_call(cmd, stdout = devnull)
    except subprocess.CalledProcessError as e:
        print "Problem in running the following command:"
        print " ".join(e.cmd)
        sys.exit(1)


def get_iovs(db, tag):
    """Retrieve the list of IOVs from `db` for `tag`.
    Arguments:
    - `db`: database connection string
    - `tag`: tag of database record
    """

    db = db.replace("sqlite_file:", "").replace("sqlite:", "")

    con = conddb.connect(url = conddb.make_url(db))
    session = con.session()
    IOV = session.get_dbtype(conddb.IOV)

    iovs = set(session.query(IOV.since).filter(IOV.tag_name == tag).all())
    if len(iovs) == 0:
        pass
        #print "No IOVs found for tag '"+tag+"' in database '"+db+"'."
        #sys.exit(1)

    session.close()

    return sorted([int(item[0]) for item in iovs])


merged_file="sqlite_file:merged.db"

files = commands.getstatusoutput("ls . | grep .db")[1].split("\n")
for file in files:
    print file
    packed_file="sqlite_file:"+file
    theIOVs = get_iovs(packed_file,"SiStripBadStrip_pcl")
    if(len(theIOVs)!=0):
        for iov in theIOVs:
            cmd = ("conddb_import",
                   "-f", packed_file,
                   "-c", merged_file,
                   "-i", "SiStripBadStrip_pcl",
                   "-t", "SiStripBadStrip_pcl",
                   "-b", str(iov),
                   "-e", str(iov))
            print cmd
            run_checked(cmd)
