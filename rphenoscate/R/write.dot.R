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

