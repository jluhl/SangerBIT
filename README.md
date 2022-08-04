SangerClone
====

SangerClone is designed to merge, trim and align Sanger sequences to reference genes.  


## SangerClone pipeline can process multiple Sanger clones with multiple reference genes by running only one single script.  


## System requirments  
### CAP3

```
# Download CAP3 from https://faculty.sites.iastate.edu/xqhuang/cap3-assembly-program
export PATH=PATH_TO_CAP3_DIR/CAP3:$PATH
```

### Python packages  
```
pip install biopython  
pip install pyfaidx  
pip install -U textwrap3  
conda install blast  
```


## Usage  
## Parameters  
```
python scripts/SangerClone.py  
--o: output directory path [default='sangerOut']  
--csv: a comma-delimited text description file with 4 columns: Clone Name, Forward Sanger Sequence File Name, Reverse Sanger Sequence File Name, Reference Gene Name [required]    
--ref: path to the reference gene(s), a single FASTA format file consist of input gene(s) sequences with ID(s) corresponding to the CSV input. [required]    
--seq: directory path, all Sanger Sequences shold be put in this directory with names corresponding to the CSV input. 
```
