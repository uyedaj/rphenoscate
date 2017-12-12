#!/usr/bin/r
# ./write_dependencies.R ../../../uberon.obo ../../data/jackson_chars_nopolys.txt

# 0. Get arguments
uberon_filename = argv[1]
input_filename = argv[2]

# this does not work.
part_of = "BFO:0000050"

library('ontologyIndex')
print("Loading ontology 'uberon.obo'")
uberon_is_a= get_ontology(uberon_filename,
                      propagate_relationships = c("is_a"),
                      extract_tags = 'minimal')

uberon_part_of = get_ontology(uberon_filename,
                      propagate_relationships = c("part_of"),
                      extract_tags = 'minimal')

uberon = get_ontology(uberon_filename,
                      propagate_relationships = c("is_a", "part_of"),
                      extract_tags = 'minimal')

# 1. read list of uberon terms, one per line
input=file(input_filename,open="r")
lines=readLines(input)

uterms = c()
unames = c()
for(i in 1:length(lines))
{
    uterms[i] = sprintf("UBERON:%s",lines[i])
    unames[i] = get_term_property(uberon,"name",uterms[i])
}


# 2. for each pair of terms, check dependencies
D = get_term_descendancy_matrix(uberon, uterms);
Disa = get_term_descendancy_matrix(uberon_is_a, uterms);
Dpartof = get_term_descendancy_matrix(uberon_part_of, uterms);
outfile = file("dependencies.txt",open="w")
for(i in 1:length(uterms))
{
    for(j in 1:length(uterms))
    {
        if (i == j) next;

        if (D[i,j])
        {
            cat(sprintf("%s %s\n",uterms[i],uterms[j]))
            x = sprintf("'%s' -> '%s'",unames[i],unames[j])
            if (Disa[i,j]) { x = cat(x," [is_a]")}
            if (Dpartof[i,j]) { x = cat(x," [part_of]")}
            cat(x,"\n")
                                        # chr.id | state | chr.ancestor | state |
            cat(sprintf("%i,0,%i,1\n",j,i), file=outfile)
            cat(sprintf("%i,1,%i,1\n",j,i), file=outfile)
        }
        # OK, so does i depend on j?
    }
}

close(outfile)

# 2201586 = pectoral fin radial cartilege
# 2202027 = pectoral fin proximal radial cartilege 2
# UBERON:2201586 UBERON:2202027

