# metaPhOrs API

## Installation

```bash
# install Python modules (MySQLdb and numpy), git and IPython
sudo apt-get install python-mysqldb python-numpy git ipython

# clone metaphors_api repository
git clone git@github.com:lpryszcz/metaphors_api.git
cd metaphors_api

# run IPython
ipython
```

## Connecting
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

# get external IDs for internal ID
extids = m.get_external_ids(protid)

# get orthologs at CS cut-off of 0.5
orthologs = m.get_orthologs_and_paralogs(protid, cs=0.5)

```
