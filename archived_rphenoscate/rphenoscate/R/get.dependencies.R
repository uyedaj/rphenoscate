library('hashmap')
library('sets')
library('ontologyIndex')

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
remove.indirect = function(x)
{
    x & (!logical_mult(x,x))
}


empty.df = function(colnames)
{
    setNames(data.frame(matrix(ncol=length(colnames), nrow=0)),colnames)
}

                                        # this doesn't really work with lists...
addrow = function (df,row)
{
    colnames = names(df)
                                        #    print("addrow:")
                                        #    print(row)
                                        #    print(colnames)
    dfrow = data.frame(matrix(ncol=length(row),nrow=1,row), stringsAsFactors=FALSE)
    dfrow = setNames(dfrow,colnames)
    rbind(df,dfrow)
}

get.dependencies = function(uberon.filename, uterms, indirect=FALSE)
{
                                        # 1. Load uberon propagating different relationships
    print("Loading ontology 'uberon.obo'")
    uberon_is_a= get_ontology(uberon.filename,
                              propagate_relationships = c("is_a"),
                              extract_tags = 'minimal')

    uberon_part_of = get_ontology(uberon.filename,
                                  propagate_relationships = c("part_of"),
                                  extract_tags = 'minimal')

    uberon_develops_from = get_ontology(uberon.filename,
                                        propagate_relationships = c("develops_from"),
                                        extract_tags = 'minimal')

    uberon = get_ontology(uberon.filename,
                          propagate_relationships = c("is_a", "part_of","develops_from"),
                          extract_tags = 'minimal')

                                        # 2. Get names for the uberon terms
    unames = c()
    for(i in 1:length(uterms))
    {
        unames[i] = get_term_property(uberon,"name",uterms[i])
    }

                                        # 3. for each pair of terms, check dependencies
    colnames = c("chr.id","dependent.chr.id","subrels")
    df = empty.df(colnames)

    if (indirect)
    {
        maybe.remove.indirect = function(x) {x}
    }
    else
    {
        maybe.remove.indirect = remove.indirect
    }

    D = maybe.remove.indirect(get_term_descendancy_matrix(uberon, uterms));
    Disa = maybe.remove.indirect(get_term_descendancy_matrix(uberon_is_a, uterms));
    Dpartof = maybe.remove.indirect(get_term_descendancy_matrix(uberon_part_of, uterms));
    Ddevelopsfrom = maybe.remove.indirect(get_term_descendancy_matrix(uberon_develops_from, uterms));

    for(i in 1:length(uterms))
    {
        for(j in 1:length(uterms))
        {
            if (i == j) next;

            if (D[i,j])
            {
                cat(sprintf("%s %s\n",uterms[i],uterms[j]))

                labels = c()

                if (Disa[i,j])
                {
                    labels = c(labels,"is_a")
                }
                if (Dpartof[i,j])
                {
                    labels = c(labels, "part_of")
                }
                if (Ddevelopsfrom[i,j])
                {
                    labels = c(labels, "develops_from")
                }

                label = ""
                if (length(labels) > 0)
                {
                    label = paste(labels,sep=":")
                }

                cat(sprintf("'%s' -> '%s'",unames[i],unames[j]),"\n")

                # chr.id | chr.ancestor | label
                df = addrow(df,c(uterms[i], uterms[j], label))

                                        #                cat(sprintf("%s,0,%s,1\n",unames[j],unames[i]), file=outfile)
                                        #                cat(sprintf("%s,1,%s,1\n",unames[j],unames[i]), file=outfile)


                                        #                cat(edge, " ", attributes, "\n", file=dotfile)
            }
                                        # OK, so does i depend on j?
        }
    }

    tuple(df,uterms,hashmap(uterms,unames),uberon)
}

