# metaPhOrs API

## Installation

```bash
# install Python modules (MySQLdb and numpy), git and IPython
sudo apt-get install python-mysqldb python-numpy python-setuptools git ipython
sudo easy_install wordcloud

# clone metaphors_api repository
git clone https://github.com/lpryszcz/metaphors_api 
cd metaphors_api

# run IPython
ipython
```

## API functionality

```python
import dbClient
m = dbClient.metaphors()

# if you wish info about any function type its name finished by question mark ie.
m.get_orthologs_and_paralogs?

# show species
m.species

# see information about human
## TaxonomyDB id for human is 9606 http://www.ncbi.nlm.nih.gov/Taxonomy/Browser/wwwtax.cgi?id=9606
m.species[9606]

# get internal ID for TP63
protid = m.get_metaid('TP63')

# get entry description
m.get_description(protid)

# get external IDs for internal ID
extids = m.get_external_ids(protid)

# get FastA-formatted sequence 
fasta = m.get_fasta(protid)
## returns fasta of random sequence if no protid is provided

# get gene symbol
m.get_gene_name(protid)

# get GO 
go = m.get_GO(protid)
go['molecular_function']

# get orthologs at CS cut-off of 0.5
orthologs, paralogs = m.get_orthologs_and_paralogs(protid)
# orthologs = { taxid: [ [ protid, extid, consistency_score, evidence_level, no_of_tress, db_info, co-orthologs* ] ] }
## db_info: consistency_score and no_of_trees from PhylomeDB, Ensembl, EggNOG, OrthoMCL, TreeFAM and Hogenome
## co-orthologs: protid and consistency_score of each co-ortholog
* or co-paralogs for paralogs

# get protids of all orthologs
protids = [protid for taxid, odata in orthologs.iteritems() for protid, extid, CS, EL, noTrees, dbInfo, coOrthologs in odata]
# get fasta sequences
fastas = [m.get_fasta(protid) for protid in protids]

```
