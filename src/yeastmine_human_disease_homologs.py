#!/usr/bin/env python

# This is an automatically generated script to run your query
# to use it you will require the intermine python client.
# To install the client, run the following command from a terminal:
#
#     sudo easy_install intermine
#
# For further documentation you can visit:
#     http://intermine.readthedocs.org/en/latest/web-services/

# The line below will be needed if you are running this script with python 2.
# Python 3 will ignore it.
from __future__ import print_function

# The following two lines will be needed in every python script:
from intermine.webservice import Service
service = Service("https://yeastmine.yeastgenome.org/yeastmine/service")

# Get a new query on the class (table) you will be querying:
query = service.new_query("Gene")

# The view specifies the output columns
query.add_view(
    "symbol", "secondaryIdentifier", "name", "homologues.homologue.symbol",
    "homologues.homologue.name", "homologues.homologue.organism.shortName",
    "homologues.homologue.primaryIdentifier",
    "homologues.homologue.crossReferences.identifier",
    "homologues.homologue.diseases.identifier",
    "homologues.homologue.diseases.name", "homologues.dataSets.dataSource.name"
)

# Uncomment and edit the line below (the default) to select a custom sort order:
# query.add_sort_order("Gene.symbol", "ASC")

# You can edit the constraint values below
query.add_constraint("homologues.homologue.crossReferences.source.name", "=", "MIM", code="D")
query.add_constraint("homologues.homologue.organism.shortName", "=", "H. sapiens", code="C")
query.add_constraint("organism.shortName", "=", "S. cerevisiae", code="B")
query.add_constraint("Gene", "IN", "ALL_Verified_Uncharacterized_Dubious_ORFs", code="A")

# Uncomment and edit the code below to specify your own custom logic:
# query.set_logic("A and B and C and D")

# Outer Joins
# (display properties of these relations if they exist,
# but also show objects without these relationships)
query.outerjoin("homologues.homologue.diseases")

for row in query.rows():
    print(row["symbol"], row["secondaryIdentifier"], row["name"], \
        row["homologues.homologue.symbol"], row["homologues.homologue.name"], \
        row["homologues.homologue.organism.shortName"], \
        row["homologues.homologue.primaryIdentifier"], \
        row["homologues.homologue.crossReferences.identifier"], \
        row["homologues.homologue.diseases.identifier"], row["homologues.homologue.diseases.name"], \
        row["homologues.dataSets.dataSource.name"])

