#!/usr/bin/env python

import os
import sys
import glob
import shutil

try:
    from setuptools import setup, Extension, Command
    from setuptools.command.build_ext import build_ext
except ImportError:
    from distutils.core import setup, Extension, Command
    from distutils.command.build_ext import build_ext

def invoke_f2py(files, flags=[], wd=None):
    """
    From python-fsps
    """
    from numpy.f2py import main

    olddir = os.path.abspath(os.curdir)
    oldargv = list(sys.argv)
    try:
        if wd is not None:
            os.chdir(wd)
        sys.argv = ['f2py']
        sys.argv.extend(files)
        sys.argv.extend(flags)

        main()
    finally:
        sys.argv = oldargv
        os.chdir(olddir)

# Hack to break
if '--break' in sys.argv:
    index = sys.argv.index('--break')
    sys.argv.pop(index)  # Removes the '--foo'
    BREAK_INSTALL = True
else:
    BREAK_INSTALL = False
        
class build_alf(build_ext):

    def run(self):
        # Generate the Fortran signature/interface.
        files = ['alf.f90']
        flags = " -m _alf -h alf.pyf --overwrite-signature".split()
        print("Running f2py on {0} with flags {1}".format(files, flags))
        invoke_f2py(['alf.f90'], flags, wd='alf')

        # Find the ALF source files.
        alf_dir = os.path.join(os.environ["ALF_HOME"], "src")
        fns = [f for f in glob.glob(os.path.join(alf_dir, "*.o"))
               if os.path.basename(f) not in ["spec_from_sum.o", "write_a_model.o", "alf.o"]]
               
        # Check to make sure that all of the required modules exist.
        flag = len(fns)
        # flag *= os.path.exists(os.path.join(alf_dir, "sps_utils.mod"))
        # flag *= os.path.exists(os.path.join(alf_dir, "sps_vars.mod"))
        if not flag:
            raise RuntimeError("You need to run make in $ALF_HOME/src first")

        # Add the interface source files to the file list.
        fns += ["alf.f90", "alf.pyf"]

        # Compile the library.
        #flags = '-c -I{0} --f90flags=-cpp --f90flags=-fPIC'.format(alf_dir)
        if BREAK_INSTALL:
            flags = '-c -I{0} -I/sw/include/openmpi -L/sw/lib/openmpi  --f90flags=-cpp --f90flags=-fno-strict-overflow --f90flags=-mcmodel=medium --f90flags=-fwrapv -lmpix --f90flags=-O3'.format(alf_dir)
        else:
            flags = '-c -I{0} -I/sw/include/openmpi -L/sw/lib/openmpi  --f90flags=-cpp --f90flags=-fno-strict-overflow --f90flags=-mcmodel=medium --f90flags=-fwrapv -lmpi --f90flags=-O3'.format(alf_dir)
            
        flags = flags.split()
        print("Running f2py on {0} with flags {1}".format(fns, flags))
        invoke_f2py(fns, flags, wd='alf')

        # Move the compiled library to the correct directory.
        infn = os.path.abspath(self.get_ext_filename("alf._alf"))
        outfn = os.path.abspath(self.get_ext_fullpath("alf._alf"))
        if infn != outfn:
            try:
                os.makedirs(os.path.dirname(outfn))
            except os.error:
                pass
            print("Copying {0} to {1}".format(infn, outfn))
            shutil.copyfile(infn, outfn)


if "publish" in sys.argv[-1]:
    os.system("python setup.py sdist upload")
    sys.exit()


# Hackishly inject a constant into builtins to enable importing of the
# package before the library is built.
if sys.version_info[0] < 3:
    import __builtin__ as builtins
else:
    import builtins
    
builtins.__ALF_SETUP__ = True
from alf import __version__  # NOQA

# This is a fake extension that is used to trick distutils into building our
# real library using the `build_fsps` function above even when `install` is
# called.
ext = Extension("alf._alf", sources=["alf/alf.f90"])

# The final setup command. Note: we override the `build_ext` command with our
# custom version from above.
setup(
    name="alf",
    url="https://github.com/gbrammer/python-alf",
    version=__version__,
    author="Gabriel Brammer",
    author_email="brammer@stsci.edu",
    description="Python bindings for Charlie Conroy's Alf.",
    long_description=open("README.md").read(),
    packages=["alf"],
    package_data={
        "": ["README.md"],
        "alf": ["_alf.so"],
    },
    include_package_data=True,
    ext_modules=[ext],
    scripts=glob.glob("scripts/*.py"),
    cmdclass={
        "build_ext": build_alf,
    },
    classifiers=[
        # "Development Status :: 5 - Production/Stable",
        "Intended Audience :: Developers",
        "Intended Audience :: Science/Research",
        "License :: OSI Approved :: MIT License",
        "Operating System :: OS Independent",
        "Programming Language :: Python",
    ],
)
