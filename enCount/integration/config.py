from enCount.config import ENDATA_TEST_REPO, volume_root
import subprocess as sp
import os

def sync_test_data():
    """
    Update data from remote Git repository.
    :return:
    """

    if not os.path.exists(os.path.join(volume_root, ".git")):
        # Directory is not yet a git repo
        args = ["git", "lfs", "clone", ENDATA_TEST_REPO, volume_root]
        print(" ".join(args))
        r = sp.call(args)
        assert r == 0
    else:
        back = os.getcwd()
        os.chdir(volume_root)
        print("Working directory: %s" % os.getcwd())
        args = ["git", "lfs", "pull"]
        print(" ".join(args))
        r = sp.call(args)
        assert r == 0
        args = ["git", "pull", "--rebase"]
        print(" ".join(args))
        r = sp.call(args)
        assert r == 0
        os.chdir(back)
        print("Working directory: %s" % os.getcwd())