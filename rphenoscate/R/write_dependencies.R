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

# Read list of uberon terms, one integer per line (separate into different function)
read.uberon.terms = function(filename)
{
    input=file(filename, open="r")
    lines=readLines(input)
    uterms = c()
    for(i in 1:length(lines))
    {
        uterms[i] = sprintf("UBERON:%s",lines[i])
    }
    close(input)
    uterms
}

                                        # 0. Get arguments

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

    triple(df,uterms,hashmap(uterms,unames))
}

add.names = function(deps,outfile)
{
    df = deps[[1]]
    id.to.name = deps[[3]]
    df$chr.name = id.to.name[[ df$chr.id ]]
    df$dependent.chr.name = id.to.name[[ df$dependent.chr.id ]]
    df2 = data.frame(df$chr.id,df$chr.name,df$dependent.chr.id,df$dependent.chr.name,df$subrels)
    triple(df2,df[[2]],df[[3]])
}

write.dependencies = function(deps,outfile)
{
    df = deps[[1]]
    write.csv(df, outfile, row.names=F, quote=F)
}

write.dot = function(deps.tuple,filename)
{
    deps = deps.tuple[[1]]
    uterms = deps.tuple[[2]]
    id.to.name = deps.tuple[[3]]

    dotfile = file(filename,open="w")

    cat("digraph \"token0\" {\ngraph [ranksep=0.25, fontname=Arial,  nodesep=0.25, ranksep=0.5];\nnode [fontname=Arial, style=filled, height=0, width=0, shape=box];\nedge [style=\"setlinewidth(2)\"];\n",file=dotfile)

    uterms1 = deps$chr.id
    unames1 = id.to.name[[uterms1]]

    uterms2 = deps$dependent.chr.id
    unames2 = id.to.name[[uterms2]]

    subrels = deps$subrels
    for(i in 1:length(uterms))
    {
        uterm = uterms[i]
        uname = id.to.name[[uterm]]
        cat(sprintf("\"%s\" [label=\"%s\\n(%s)\"]\n",uterm, uname, uterm),file=dotfile)
    }

    for(i in 1:nrow(deps))
    {
        rgb=c("00","00","00")

        cat("subrels[i] = '",subrels[i],"'\n")
        labels = strsplit( subrels[i], ':')
        
        if ("is_a" %in% labels)
        {
            rgb[1] = "ff"
        }

        if ("part_of" %in% labels)
        {
            rgb[3] = "ff"
        }

        if ("develops_from" %in% labels)
        {
            rgb[2] = "ff"
        }

        attributes = c(paste0("color=\"#",paste0(rgb,collapse=''),"\""))
        if (length(labels) > 0)
        {
            label = paste(labels,sep=",")
            attributes=c(attributes,paste0("label=\"",label,"\""))
        }

        attributes = paste0("[",paste0(attributes,collapse=','),"]")
        edge = sprintf("\"%s\" -> \"%s\"",uterms2[i],uterms1[i])
        cat("edge: ",edge," ",uterms2[i]," ",uterms2[i],"\n")
        cat(edge, " ", attributes, "\n", file=dotfile)
    }

    cat("}\n", file=dotfile)
    close(dotfile)
}

write.terms=function(deps,filename)
{
    uterms = deps[[2]]

    outfile = file(filename,open="w")
    for(i in 1:length(uterms))
    {
        cat(uterms[i],"\n",file=outfile)
    }
    close(outfile)
}

library('hashmap')
library('sets')
library('ontologyIndex')
uberon.terms = read.uberon.terms('../../data/jackson_chars_nopolys.txt')
deps = get.dependencies('../../../uberon.obo', uberon.terms)

uberon.filename = argv[1]
uberon.terms = read.uberon.terms(argv[2])

deps = get.dependencies(uberon.filename, uberon.terms)

write.dependencies(deps,"dependencies.txt")
write.dependencies(add.names(deps),"dependencies_with_names.txt")

write.dot(deps,"dependencies.dot")
write.terms(deps,"terms.txt")

indirect.deps = get.dependencies(uberon.filename, uberon.terms, indirect=TRUE)
write.dependencies(indirect.deps,"indirect_dependencies.txt")

# 2201586 = pectoral fin radial cartilege
# 2202027 = pectoral fin proximal radial cartilege 2
# UBERON:2201586 UBERON:2202027

# uberon.terms = read.uberon.terms('../../data/jackson_chars_nopolys.txt')
# deps = get.dependencies('../../../oberon.obo', uberon.terms)
