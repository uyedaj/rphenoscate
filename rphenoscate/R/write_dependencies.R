#!/usr/bin/r
# ./write_dependencies.R ../../../uberon.obo ../../data/jackson_chars_nopolys.txt

source("get.dependencies.R")

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

add.names = function(deps,outfile)
{
    df = deps[[1]]
    id.to.name = deps[[3]]
    df$chr.name = id.to.name[[ df$chr.id ]]
    df$dependent.chr.name = id.to.name[[ df$dependent.chr.id ]]
    df2 = data.frame(df$chr.id,df$chr.name,df$dependent.chr.id,df$dependent.chr.name,df$subrels)

    deps2 = deps
    deps2[[1]] = df2
    deps2
}

write.dependencies = function(deps,outfile)
{
    df = deps[[1]]
    write.csv(df, outfile, row.names=F, quote=F)
}

source("write.dot.R")

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
