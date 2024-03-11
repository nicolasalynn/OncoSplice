To install Oncosplice:

1. Download the .whl file

2. pip install oncosplice.whl

3.
poetry run datasetup -c /path/to/conservation_data.pkl -b /path/to/database/location -s /path/to/missplicing/location

4.
python
from oncosplice import *
mut_id = 'KRAS:12:25227343:G:T'
results = oncosplice(mut_id, sai_threshold=0.5, annotate=True)


