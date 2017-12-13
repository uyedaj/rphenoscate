#!/usr/bin/r
# ./write_dependencies.R ../../../uberon.obo ../../data/jackson_chars_nopolys.txt


# This function does logical matrix multiplication (|| and && instead of + and *)
# It is used to find out if i is connected to j through any k
logical_mult = function(x,y)
{
    I = dim(x)[1]
    K = dim(x)[2]
    J = dim(y)[2]

    z = rep(FALSE,I*J)
    dim(z) = c(I,J)
    for(i in 1:I)
    {
        for(j in 1:J)
        {
            for(k in 1:K)
            {
                z[i,j] = z[i,j] || (x[i,k] && y[k,j])
            }
        }
    }
    z
}

# connected(i,j) if i is connected to j but not through any other k
remove_indirect = function(x)
{
    x & (!logical_mult(x,x))
}


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

uberon_develops_from = get_ontology(uberon_filename,
                                    propagate_relationships = c("develops_from"),
                                    extract_tags = 'minimal')

uberon = get_ontology(uberon_filename,
                      propagate_relationships = c("is_a", "part_of","develops_from"),
                      extract_tags = 'minimal')

# FIXME: optionally handle develops_from ?

dotfile = file("dependencies.dot",open="w")
cat("digraph \"token0\" {\ngraph [ranksep=0.25, fontname=Arial,  nodesep=0.25, ranksep=0.5];\nnode [fontname=Arial, style=filled, height=0, width=0, shape=box];\nedge [style=\"setlinewidth(2)\"];\n",file=dotfile);

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
D = remove_indirect(get_term_descendancy_matrix(uberon, uterms));
Disa = remove_indirect(get_term_descendancy_matrix(uberon_is_a, uterms));
Dpartof = remove_indirect(get_term_descendancy_matrix(uberon_part_of, uterms));
Ddevelopsfrom = remove_indirect(get_term_descendancy_matrix(uberon_develops_from, uterms));
outfile = file("dependencies.txt",open="w")
for(i in 1:length(uterms))
{
    cat(sprintf("\"%s\" [label=\"%s\n(%s)\"]\n",uterms[i], unames[i], uterms[i]),file=dotfile)

    for(j in 1:length(uterms))
    {
        if (i == j) next;

        if (D[i,j])
        {
            cat(sprintf("%s %s\n",uterms[i],uterms[j]))

            labels = c()
            rgb=c("00","00","00")

            if (Disa[i,j])
            {
                labels = c(labels,"is_a")
                rgb[1] = "ff"
            }
            if (Dpartof[i,j])
            {
                labels = c(labels, "part_of")
                rgb[3] = "ff"
            }
            if (Ddevelopsfrom[i,j])
            {
                labels = c(labels, "develops_from")
                rgb[2] = "ff"
            }

            label = ""
            attributes = c(paste0("color=\"#",paste0(rgb,collapse=''),"\""))
            if (length(labels) > 0)
            {
                label = paste(labels,sep=",")
                attributes=c(attributes,paste0("label=\"",label,"\""))
            }


            x = sprintf("'%s' -> '%s'",unames[i],unames[j])
            cat(x,"\n")

            # chr.id | state | chr.ancestor | state |
            cat(sprintf("%s,0,%s,1\n",unames[j],unames[i]), file=outfile)
            cat(sprintf("%s,1,%s,1\n",unames[j],unames[i]), file=outfile)

            edge = sprintf("\"%s\" -> \"%s\"",uterms[j],uterms[i])
            attributes = paste0("[",paste0(attributes,collapse=','),"]")
            cat(edge, " ", attributes, "\n", file=dotfile)
        }
        # OK, so does i depend on j?
    }
}

close(outfile)

cat("}\n", file=dotfile)
close(dotfile)

# 2201586 = pectoral fin radial cartilege
# 2202027 = pectoral fin proximal radial cartilege 2
# UBERON:2201586 UBERON:2202027

