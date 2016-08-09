library(optparse)

if (!exists('opt'))
{
    option_list = list(
        make_option(c("-c", "--cov"), type = "character", help = "Path to .cov file"),
        make_option(c("-n", "--bin.red"), type = "integer", help = "Fold increase in size of coverage bins"),
        make_option(c("-s", "--sub.sample.fold"), type = "integer", help = "Fold by which to sub-sample .cov file"),
        make_option(c("-b", "--bin.num"), type = "integer", default = 50, help = "Number of segments to evaluate per array"),
        make_option(c("-p", "--prob.cutoff"), type = "integer", default = 0.01, help = "Cutoff purity ploidy probability value for inclusion in CN probability matrix"),
        make_option(c("-o", "--outdir"), type = "character", default = './', help = "output directory"),
        make_option(c("-l", "--libdir"), type = "character", help = "Directory containing this R file")
)

parseobj = OptionParser(option_list=option_list)
opt = parse_args(parseobj)

if (is.null(opt$libdir) | (is.null(opt$cov)))
        stop(print_help(parseobj))

print(opt)

print(.libPaths())
options(error=function()traceback(2))

## keep record of run
writeLines(paste(paste('--', names(opt), ' ', sapply(opt, function(x) paste(x, collapse = ',')), sep = '', collapse = ' '), sep = ''), paste(opt$outdir, 'cmd.args', sep = '/'))
    saveRDS(opt, paste(opt$outdir, 'cmd.args.rds', sep = '/'))
}

library(skitools)
library(DNAcopy)
library(data.table)
library(pscl)
library(gplots)
library(d3heatmap)

system(paste('mkdir -p', opt$outdir))

##############################################################
############## Generate input data ###########################
message('Generating input data')

.segment = function(tcov){
    ix = which(!is.na(tcov$ratio))
    cat('sending ', length(ix), ' segments\n')
    cna = CNA(log(tcov$ratio[ix]), as.character(seqnames(tcov))[ix], start(tcov)[ix], data.type = 'logratio')
    gc()
    cat('finished making cna\n')
    Seg = segment(smooth.CNA(cna), alpha = 1e-5, verbose = T)
    cat('finished segmenting\n')
    ## out = seg2gr(print(seg), new.sl) ## remove seqlengths that have not been segmented
    ## out = gr.fix(out, new.sl, drop = T)
    return(Seg)
}

set.bin.size <- function(n = opt$bin.red, data_path = opt$cov){
    CBS_cov<-readRDS(data_path)
    raw <- cbind(CBS_cov$ratio,as.data.table(ranges(CBS_cov))$start,as.data.table(seqnames(CBS_cov)))
    raw[V1 ==Inf,V1:=NA]
    new <- raw[,.(ratio = mean(V1,na.rm=T),start = min(V2)),by=.(value,(seq(nrow(raw)) - 1) %/% n)]
    new <- new[!is.na(ratio)]
    new <- new[ratio > 0]
    seq.names <- as.data.table(table(new$value))
                                        #Create new grange table which is required for .segment function
    grange_new <- GRanges(seqnames = Rle(seq.names$V1, seq.names$N), ranges = IRanges(new$start,end=new$start+n*200-1,width=NULL,names=NULL), ratio = new$ratio)
    new.seg<-.segment(grange_new)
    final <- cbind(new, rep(seq(1, length(new.seg$output$num.mark), 1), new.seg$output$num.mark))
    names(final) <- c("chr","seq","data","start","segment")
    return(final)
}

set.seed(1)
y = set.bin.size()

sub.sample <- function(num = opt$sub.sample.fold, z = y){
    x <- z[sample(nrow(z),nrow(z)/num)]
    return(x)
}

log.sum.exp<- function(x){
    offset <- max(x)
    log(sum(exp(x - offset))) + offset
}

x <- sub.sample()

###############################################################
############ Run this version do determine the SD #############
message('Calculating SD')

sloidy = mean(x$data)

combo_whole = data.table(expand.grid(ploidy = seq(0.5, 5, 0.5), purity = seq(0.1, 1.0, 0.1), SD = seq(sd(x$data/100),sd(x$data),sd(x$data/100)), CN = 0:10))
setkeyv(combo_whole,c("ploidy", "purity", "SD"))
num.ppv = nrow(unique(combo_whole[,list(ploidy, purity,SD)]))
ppv <- unique(combo_whole[,list(ploidy, purity, SD)])
num.ppv = ceiling(num.ppv/opt$bin.num)

my.list<-vector("list",num.ppv)
system.time(
    for(i in 1:num.ppv){
        combo = merge(combo_whole, ppv[((i * opt$bin.num) - opt$bin.num + 1) : (i * opt$bin.num),])
        input <- combo[,.(data = rep(x$data, each = nrow(combo)), purity = rep(purity, nrow(x)), ploidy = rep(ploidy, nrow(x)), SD = rep (SD, nrow(x)), CN = rep(CN, nrow(x)), segment = rep(x$segment, each = nrow(combo)))]
        input[, CN_prior:= dnorm(CN, ploidy, 1)]
        input[, CN_prior := CN_prior / sum(unique(CN_prior)), by=ploidy]
                                        #Added log=T to dnorm. Also now adding low of CN_prior
        input[, prob:= dnorm(data, (sloidy * (purity * CN + (2 * (1 - purity)))) / ((2 * (1 - purity)) + (purity * ploidy)), SD, log=T) + log( CN_prior)]
                                        #Removed the log from inside the sume here
        input<-input[,.(logsum = sum(prob)),by=.(purity, ploidy, SD, CN, segment)]
        input<- input[,.(lse = log.sum.exp(logsum)), by = .(purity, ploidy, SD, segment)]
        input<- input[,.(final = sum(lse)), by = .(purity, ploidy, SD)]
        my.list[[i]]<-input
        print(i)
    }
)

input <- rbindlist(my.list)
setkey(input,SD)
int = input[,.(lse = log.sum.exp(final)), keyby = list(SD)]
out = merge(input, int)[, prob := exp(final - lse)]
MAP = out[,. (MAP = max(lse) + log(densigamma(SD,1,1))),by=SD]
SD.MAP = MAP[MAP == max(MAP)]$SD

saveRDS(out, file = paste0(opt$outdir, "/SD_prob.Rds"))


##############################################################
########## Run this w/ max SD to get purity and ploidy########
message('Calculating purity and ploidy')

combo_whole = data.table(expand.grid(ploidy = seq(0.55, 5.5, 0.05), purity = seq(0.01, 1.0, 0.01), SD = SD.MAP, CN = 0:10))
setkeyv(combo_whole,c("ploidy", "purity", "SD"))
num.ppv = nrow(unique(combo_whole[,list(ploidy, purity,SD)]))
ppv <- unique(combo_whole[,list(ploidy, purity, SD)])
num.ppv = ceiling(num.ppv/opt$bin.num)

CN.list<-vector("list",num.ppv)
my.list<-vector("list",num.ppv)
system.time(
    for(i in 1:num.ppv){
        combo = merge(combo_whole, ppv[((i * opt$bin.num) - opt$bin.num + 1) : (i * opt$bin.num),])
        input <- combo[,.(data = rep(x$data, each = nrow(combo)), purity = rep(purity, nrow(x)), ploidy = rep(ploidy, nrow(x)), SD = rep (SD, nrow(x)), CN = rep(CN, nrow(x)), segment = rep(x$segment, each = nrow(combo)))]
        input[, CN_prior:= dnorm(CN, ploidy, 1)]
        input[, CN_prior := CN_prior / sum(unique(CN_prior)), by=ploidy]
        input[, prob:= dnorm(data, (sloidy * (purity * CN + (2 * (1 - purity)))) / ((2 * (1 - purity)) + (purity * ploidy)), SD, log=T) + log( CN_prior)]
        input<-input[,.(logsum = sum(prob)),by=.(purity, ploidy, SD, CN, segment)]
        CN.list[[i]] <- dcast(input,purity + ploidy + segment ~ CN, value.var = 'logsum')
        input<- input[,.(lse = log.sum.exp(logsum)), by = .(purity, ploidy, SD, segment)]
        input<- input[,.(final = sum(lse)), by = .(purity, ploidy, SD)]
        my.list[[i]]<-input
        print(i)
    }
)

input <- rbindlist(my.list)
CN.matrix <- rbindlist(CN.list)

setkey(input,SD)
int = input[,.(lse = log.sum.exp(final)), keyby = list(SD)]
out = merge(input, int)[, prob := exp(final - lse)]

saveRDS(out, file = paste0(opt$outdir, "/PP_prob.Rds"))


#Make heat map
hm <- acast(out,ploidy~purity,value.var="prob")

pdf(file = paste0(opt$outdir, '/PurityPloidy.pdf'), width = 12, height = 12)
heatmap.2(hm,Colv=NA,Rowv=NA,scale="none",dendrogram="none",revC = TRUE, trace='none',key=FALSE,ylab = "Ploidy", xlab = "Purity",colsep = c(1:100),rowsep=c(1:100),sepwidth = c(0.001,0.001),sepcolor="black", main = paste0("SD = ", SD.MAP))
dev.off()
saveRDS(hm, file = paste0(opt$outdir, "/heatmap.Rds"))
#wij(heatmap(hm,dendrogram=NULL,Rowv=FALSE,Colv=FALSE),filename=paste0(opt$outdir, '/PurityPloidy.html'))

#Make CN matrix
out <- out[prob > opt$prob.cutoff]
setkeyv(out,c("purity", "ploidy"))
setkeyv(CN.matrix,c("purity", "ploidy"))
high.prob.CN.matrix <- merge(out, CN.matrix)
CN<-high.prob.CN.matrix[,(exp(.SD)/rowSums(exp(.SD))), .SDcols = 8:18]
CN.output <- cbind(high.prob.CN.matrix[,.(purity, ploidy, segment)], CN)

saveRDS(CN.output, file = paste0(opt$outdir, "/CN_matrix.Rds"))

message('done')
