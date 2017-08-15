echo "Run mpifort"
python setup.py install --break 2>&1 >/dev/null | grep mpix |sed "s/mpix/mpi/" |sed "s/error: Command \"\/sw\/bin\/gfortran/mpifort/" | sed "s/\" failed with exit status 1//" > /tmp/alf_mpi.sh
cat /tmp/alf_mpi.sh
sh /tmp/alf_mpi.sh

echo ""
echo "cp _alf.cpython-35m-darwin.so /Users/brammer/anaconda3/lib/python3.5/site-packages/alf-0.0.0-py3.5-macosx-10.7-x86_64.egg/alf/"
cp _alf.cpython-35m-darwin.so /Users/brammer/anaconda3/lib/python3.5/site-packages/alf-0.0.0-py3.5-macosx-10.7-x86_64.egg/alf/

echo "Done!"

