#!/usr/bin/r
# ./write_dependencies.R ../../../uberon.obo ../../data/jackson_chars_nopolys.txt

# 0. Get arguments
uberon_filename = argv[1]
input_filename = argv[2]

# open obo file with DEX editor to see
part_of = "BFO:0000050"

library('ontologyIndex')
print("Loading ontology 'uberon.obo'")
uberon = get_ontology(uberon_filename,
                      propagate_relationships = c("is_a", part_of),
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
for(i in 1:length(uterms))
{
    name1 = get_term_property(uberon,"name",uterms[i])
    for(j in 1:length(uterms))
    {
        if (i == j) next;

        name2 = get_term_property(uberon,"name",uterms[j])
        if (D[i,j])
        {
            cat(sprintf("%s %s\n",uterms[i],uterms[j]))
            cat(sprintf("'%s' -> '%s'\n",name1,name2))
        }
        # OK, so does i depend on j?
    }
}

# 2201586 = pectoral fin radial cartilege
# 2202027 = pectoral fin proximal radial cartilege 2
# UBERON:2201586 UBERON:2202027

