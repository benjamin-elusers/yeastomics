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
    "symbol", "featureType", "name", "primaryIdentifier",
    "phenotypes.publications.pubMedId", "phenotypes.experimentType",
    "phenotypes.mutantType", "phenotypes.allele.name",
    "phenotypes.alleleComment", "phenotypes.strainBackground",
    "phenotypes.qualifier", "phenotypes.observable", "phenotypes.chemical",
    "phenotypes.condition", "phenotypes.details", "phenotypes.reporter"
)

# Uncomment and edit the line below (the default) to select a custom sort order:
# query.add_sort_order("Gene.symbol", "ASC")

# Outer Joins
# (display properties of these relations if they exist,
# but also show objects without these relationships)
query.outerjoin("phenotypes.allele")

for row in query.rows():
    print(row["symbol"], row["featureType"], row["name"], row["primaryIdentifier"], \
        row["phenotypes.publications.pubMedId"], row["phenotypes.experimentType"], \
        row["phenotypes.mutantType"], row["phenotypes.allele.name"], \
        row["phenotypes.alleleComment"], row["phenotypes.strainBackground"], \
        row["phenotypes.qualifier"], row["phenotypes.observable"], row["phenotypes.chemical"], \
        row["phenotypes.condition"], row["phenotypes.details"], row["phenotypes.reporter"])

