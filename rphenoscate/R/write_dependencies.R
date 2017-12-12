#!/usr/bin/r
# ./write_dependencies.R ../../../uberon.obo ../../data/jackson_chars_nopolys.txt

# 0. Get arguments
oberon_filename = argv[1]
input_filename = argv[2]

# open obo file with DEX editor to see
part_of = "BFO:0000050"

library('ontologyIndex')
print("Loading ontology 'uberon.obo'")
uberon = get_ontology(oberon_filename,
                      propagate_relationships = c("is_a",part_of),
                      extract_tags = 'minimal')

# 1. read list of uberon terms, one per line
input=file(input_filename,open="r")
lines=readLines(input)

uterms = c()
for(i in 1:length(lines))
{
    uterms[i] = sprintf("UBERON:%s",lines[i])
}

# 2. for each pair of terms, check dependencies
D = get_term_descendancy_matrix(uberon, uterms);
for(i in 1:length(terms))
{
    for(j in 1:length(lines))
    {
        if (i == j) next;
        print(sprintf("i = %s j = %s",uterms[i],uterms[j]))
        # OK, so does i depend on j?
    }
}
