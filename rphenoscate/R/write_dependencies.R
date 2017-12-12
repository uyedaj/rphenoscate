#!/usr/bin/r
library('ontologyIndex')
print("Loading ontology 'uberon.obo'")
uberon = get_ontology('uberon.obo',
                      propagate_relationships = c("is_a","part_of"),
                      extract_tags = 'minimal')

# 1. read list of uberon terms, one per line
filename = "input.txt"
input=file(filename,open="r")
lines=readLines(input)

# 2. for each pair of terms, check dependencies
for(i in 1:length(lines))
{
    for(j in 1:length(lines))
    {
        if (i == j) next;
        print(sprintf("i = %i j = %i",i,j))
        # OK, so does i depend on j?
    }
}
