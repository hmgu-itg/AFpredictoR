#!/usr/bin/env Rscript
library(data.table)
library(VarianceGamma)


read.frqx = function(frqx.file){
    frqx = fread(frqx.file)
    setnames(frqx, c(5:10), c("homa1", "het", "homa2", "hap1", "hap2", "missing"))

    frqx[,c("frqx_ac", "frqx_an", "frqx_a1f"):=list((2*homa1+het), (2*homa1+2*het+2*homa2+2*missing), (2*homa1+het)/(2*homa1+2*het+2*homa2+2*missing))]
    # setnames(frqx, c('A1', 'A2', 'SNP'), c('frqx_A1', 'frqx_A2', 'frqx_SNP'))
    frqx[,c("A1", "A2", "homa1", "het", "homa2", "hap1", "hap2", "missing", 'CHR'):=NULL]
    return(frqx)
}


read.bim = function(bim.file){
    bim = fread(bim.file)
    bim[,V3:=NULL]
    setnames(bim, c('frqx_chrom', 'SNP', 'frqx_pos', 'frqx_A1', 'frqx_A2'))
    return(bim)
}

# Assumes columns: SNP, chrom, pos, A1, A2, AF, N
read.sumstats = function(sumstats.file){
    sumstats = fread(sumstats.file)
    sumstats[,an:=N*2]
    sumstats[,ac:=as.integer(an*AF)] # This will round down
    sumstats[,id_a1_a2:=paste0('chr', chrom, ':', pos, ':', A1, ':', A2)]
    sumstats[,id_a2_a1:=paste0('chr', chrom, ':', pos, ':', A2, ':', A1)]
    return(sumstats)
}


main = function(sumstats.file,
                bim.file,
                frqx.file,
                output,
                threshold.pdiff=0.05,
                threshold.pmax=1e-10,
                debug=FALSE) {
    output.file = paste0(output, '.csv')
    output.debug.file = paste0(output, '.debug.csv')

    sumstats = read.sumstats(sumstats.file)
    bim = read.bim(bim.file)
    frqx = read.frqx(frqx.file)
    


    frqx2 = merge(bim, frqx, by = 'SNP')
    frqx2[,var_id:=paste0('chr', frqx_chrom, ':', frqx_pos, ':', frqx_A1, ':', frqx_A2)]
    setnames(frqx2, 'SNP', 'frqx_SNP')


    m1=merge(frqx2, sumstats, by.x="var_id", by.y="id_a1_a2")
    m1[,id_a2_a1:=NULL]

    m2=merge(frqx2, sumstats, by.x="var_id", by.y="id_a2_a1")
    m2[,id_a1_a2:=NULL]

    m = rbind(m1, m2)
    m = m[!duplicated(m)] # There were some sumstats with duplicate rows.
    m[,var_id:=NULL]

    # Run proportions test and add test statistics to the table
    binomp=apply(m[,.(ac, frqx_ac, an, frqx_an)], 1, function(x){
        te=prop.test(x=x[1:2], n=x[3:4])
        return(list(p=te$p.value, chisq=te$statistic))
    })
    binomp2=apply(m[,.(ac, frqx_ac, an, frqx_an)], 1, function(x){
        x2=x[4]-x[2]
        te=prop.test(x=c(x[1],x2), n=x[3:4])
        return(list(p=te$p.value, chisq=te$statistic))
    })



    addendum=cbind(rbindlist(binomp), rbindlist(binomp2))
    setnames(addendum, c("p1", "chi1", "p2", "chi2"))
    addendum[,c("pmax", "chidiff"):=list(pmax(p1,p2), -abs(chi1-chi2))]
    addendum[,AF_is_for_frqx_A1:=p1>p2]
    addendum[,pdiff:=pvg(chidiff, vgC=0, theta=0, nu=2, sigma=2)]

    m_addendum=cbind(m, addendum)

    if (debug==TRUE) fwrite(m_addendum, output.debug.file)

    m_addendum2 = m_addendum[pdiff<threshold.pdiff 
                            & pmax>threshold.pmax,
                            .(SNP, chrom, pos, A1, A2, AF, N, AF_is_for_frqx_A1, frqx_A1, frqx_A2, pdiff, pmax)]


    # Case 1 and 7
    m_addendum2[AF_is_for_frqx_A1==TRUE
                & A1==frqx_A2
                & A2==frqx_A1, AF:=1-AF]
    # Case 2 and 8
    m_addendum2[AF_is_for_frqx_A1==FALSE
                & A1==frqx_A1
                & A2==frqx_A2, AF:=1-AF]

    out = m_addendum2[,.(SNP, chrom, pos, A1, A2, AF, N)]

    fwrite(out, output.file)
}


args = commandArgs(trailingOnly=TRUE)

args.length = length(args)

if (args.length==4){
    main(
        sumstats.file = args[1],
        bim.file = args[2],
        frqx.file = args[3],
        output = args[4],
        debug=TRUE
    )
} else if (args.length==6){
    main(
        sumstats.file = args[1],
        bim.file = args[2],
        frqx.file = args[3],
        output = args[4],
        pdiff = args[5],
        pmax = args[6],
        debug = TRUE
    )
} else {
    print('incorrect run option')
}


