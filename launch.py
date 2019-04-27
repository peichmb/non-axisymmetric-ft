echo "0100.0500.0025" 
cp params.py.0100.0500.0025 params.py
mkdir /scratch/belda/toroidal/decay/0100.0500.0025
mpiexec -n 96 python main.py > /scratch/belda/toroidal/decay/out_0100.0500.0025 

echo "0100.0500.0050" 
cp params.py.0100.0500.0050 params.py
mkdir /scratch/belda/toroidal/decay/0100.0500.0050
mpiexec -n 96 python main.py > /scratch/belda/toroidal/decay/out_0100.0500.0050 

echo "0100.0500.0100" 
cp params.py.0100.0500.0100 params.py
mkdir /scratch/belda/toroidal/decay/0100.0500.0100
mpiexec -n 96 python main.py > /scratch/belda/toroidal/decay/out_0100.0500.0100 

echo "0250.0150.0025" 
cp params.py.0250.0150.0025 params.py
mkdir /scratch/belda/toroidal/decay/0250.0150.0025
mpiexec -n 96 python main.py > /scratch/belda/toroidal/decay/out_0250.0150.0025 

echo "0250.0150.0050" 
cp params.py.0250.0150.0050 params.py
mkdir /scratch/belda/toroidal/decay/0250.0150.0050
mpiexec -n 96 python main.py > /scratch/belda/toroidal/decay/out_0250.0150.0050 

echo "0250.0150.0100" 
cp params.py.0250.0150.0100 params.py
mkdir /scratch/belda/toroidal/decay/0250.0150.0100
mpiexec -n 96 python main.py > /scratch/belda/toroidal/decay/out_0250.0150.0100 

echo "0100.0150.0025" 
cp params.py.0100.0150.0025 params.py
mkdir /scratch/belda/toroidal/decay/0100.0150.0025
mpiexec -n 96 python main.py > /scratch/belda/toroidal/decay/out_0100.0150.0025 

echo "0100.0150.0050" 
cp params.py.0100.0150.0050 params.py
mkdir /scratch/belda/toroidal/decay/0100.0150.0050
mpiexec -n 96 python main.py > /scratch/belda/toroidal/decay/out_0100.0150.0050 

echo "0100.0150.0100" 
cp params.py.0100.0150.0100 params.py
mkdir /scratch/belda/toroidal/decay/0100.0150.0100
mpiexec -n 96 python main.py > /scratch/belda/toroidal/decay/out_0100.0150.0100 
