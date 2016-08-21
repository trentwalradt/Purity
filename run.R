### Saturday, Aug 20, 2016 06:55:02 PM twalradt ####
### Heavily annotated version
####################################################


### All of these optparse commands are included so that the R script can be run with Flow.  Each make option corresponds to an option in the .task file located at ~/task/PurityPloidy.task.  The first option with a single letter is in the .task file, and the second option in in this R script.  Reference those with opt$[name].  e.g. opt$cov calls the coverage file path.  NOTE, make sure to change the type for difference inputs.  I didn't change prob.cutoff's type to double and it was a huge pain.  Double for decimal points!

library(optparse)

if (!exists('opt'))
{
    option_list = list(
        make_option(c("-c", "--cov"), type = "character", help = "Path to .cov file"),
        make_option(c("-n", "--bin.red"), type = "integer", help = "Fold increase in size of coverage bins"),
        make_option(c("-s", "--sub.sample.fold"), type = "integer", help = "Fold by which to sub-sample .cov file"),
        make_option(c("-b", "--bin.num"), type = "integer", default = 50, help = "Number of segments to evaluate per array"),
        make_option(c("-p", "--prob.cutoff"), type = "double", default = 0.01, help = "Cutoff purity ploidy probability value for inclusion in CN probability matrix"),
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

### All of these message() functions just produce text in the .err files and if you run with quiet = FALSE.
### paste0 is short for paste(x, sep=""). Just a time saver.  You can put the space in the quoted text.
### Included message below when trouble shooting.  Had originally put type as integer, fixed it by changing to double

message(paste0("cutoff is ",opt$prob.cutoff))

### Call all of the library here.  Ended up not using d3heatmap b/c of some bugs.  DNAcopy is for the .segment function. gplots is for the histogram.  MASS is for the fitdistr to find the alpha and beta variables of the beta distribution.  pscl is for inverse gamma distribution.

library(skitools)
library(DNAcopy)
library(data.table)
library(pscl)
library(gplots)
library(d3heatmap)
library(MASS)

system(paste('mkdir -p', opt$outdir))

##############################################################
############## Generate input data ###########################
message('Generating input data')

### Stole this code from another script Marcin wrote.  Input is a granges coverage file.  Takes bins and outputs a segment file.  The program was taking too long to run with bin sizes of 200 bp, so we decided to increase the bin size by averaging over a given number of bins, then using this function to call new segments.

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

### This function takes 2 inputs: the fold reduction you want, and a coverage file.  If you set n =2, then the number of bins will decrease by 2 fold.  If n =500, then the number of bins will decrease by 500 fold.
### The input file is a granges object, so I had to do some data wrangling to get it into shape to be an input for .segment.  I take the ratio column from the cov file which is a meta data column. You can access meta data columns from grange objects with $.  Ranges and seqnames (a.k.a. chromosomes) were not metadata so I had to access them differently.
### Mean() doesn't work with Inf and I couldn't find an easy way to get rid of them so I converted them all to NAs which are easy to remove
### In data.table using a period is shorthand for list().  I created a column ratio which was the mean of the old ratio values removing any NAs.  The start was the the first start value in each by segment.  The by is going by every n rows in each chromosome.  So if you have can through 400 rows, but you reach the end of a chromosome, that is still a divider.
### I removed any segments with a ratio that was NA or 0 because those caused errors in the .segment function.
### Had to create a grange object which is the input for .segment
### seq.names is a table with one column of chromosome names and another colum of how many rows there are for each chromosome.  Needed this b/c granges reads in seqnames w/ Rle() which is kind of like a factor
### For ranges I had the start values and then calculated the end values because I know the original bins were all 200 bp and then I know the fold reduction (n).
### The .segment function outputs a table where one column is the mean of the segment, and another column is the number of bins in that segment.  To annotate my data.table new with the segments I created a sequence from 1 to the number of segments and repeated each segment by the number of bins in that segment.

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

### The .segment function has some stochastic element in it so I set the seed for reproducibility

set.seed(1)
y = set.bin.size()

### This function subsamples the final data.table just created.  This parameter is set in the .task file.  The sampling is random and performed with the sample() function.  Increasing the amount of sub-sampling approximates going from WGS to WES to panel

sub.sample <- function(num = opt$sub.sample.fold, z = y){
    x <- z[sample(nrow(z),nrow(z)/num)]
    return(x)
}

### The function is used to present underflow. See link for details https://hips.seas.harvard.edu/blog/2013/01/09/computing-log-sum-exp/

log.sum.exp<- function(x){
    offset <- max(x)
    log(sum(exp(x - offset))) + offset
}

x <- sub.sample()

###############################################################
############ Run this version do determine the SD #############
message('Calculating SD')

### Sloidy is a constant.  See the equation documentation for more detail.  It takes into account bin width, but because the width is the same for all segments, it simplifies to the equation below

sloidy = mean(x$data)

### Hold purity, ploidy and SD constant and look through all values of CN, segment and bin.  For the first run through only search 10 values of purity and ploidy, but search 100 values of SD.  Take the SD of the whole sample as the maximum value to test and then search from that to value/100 by steps of value/100.
### Search through CN from 0 to 10
### Create a matrix with every possible combination of ploidy, purity, SD and CN (use expand.grid to do this)
### num.ppv is the number of unique combinations of purity, ploidy and SD
### ppv is a table with the unique values for purity, ploidy and SD
### opt$bin.num is how many rows of ppv you want to evaluate at a time.  If you have more ram, set this number higher. If less do smaller numer.  If it was 1, then you would loop through ppv one row at a time
### Create a list ahead of time to store output of loop because this is faster in R
### Select only the rows of combo_whole you want by merging with ppv of a certain range
### Create a melted table with all combinations of bin, segment, CN and the selected purity, ploidy and SD.  Create this with careful use of rep() and the condition 'each'.  This is called input.
### NOTE: priors on purity P(alpha)  and ploidy P(tau) are not included in this program anywhere b/c we assume they won't make much difference.
### For the last 3 input steps I leave out a column name in by to perform the function over that variable


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

### Up to this point I have coded in the product, sum product and everything to the right of the equation.  However, it is still in log form


input <- rbindlist(my.list)
setkey(input,SD)
int = input[,.(lse = log.sum.exp(final)), keyby = list(SD)]
out = merge(input, int)[, prob := exp(final - lse)]
MAP = out[,. (MAP = max(lse) + log(densigamma(SD,1,1))),by=SD]
SD.MAP = MAP[MAP == max(MAP)]$SD

saveRDS(out, file = paste0(opt$outdir, "/SD_prob.Rds"))

### int is the integral term over alpha and tau.  This is also still in log form
### out is the final term.  final and lse are both in log form so subtraction is division.
### MAP is the maximum a posteriori
### SD.MAP is the SD with the greatest probability that you will use in the next section of the scrip


##############################################################
########## Run this w/ max SD to get purity and ploidy########
message('Calculating purity and ploidy')

### This code is the same as above, except now there are 100 different purity and ploidy values and only 1 SD.
### setkeyv() allows you to set multiple keys for a data.table
### CN.list is a table of the probabilities for each CN for a given purity ploidy and segment. The probabilities are in log form. There is a column for CN 0 through 10.  This is cast using the dcast function.

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

### Create input file for heatmap using acast.  Save the file as a pdf.  Turn of dendrogram and 

#Make heat map. Turn of Colv and Rowv and dendrogram so that no order is imposed on the heatmap.  It just appears as it is in hm.  The gsub part is just to get the sample name.  It takes the text between the last and second to last forward slash
hm <- acast(out,ploidy~purity,value.var="prob")

pdf(file = paste0(opt$outdir, '/PurityPloidy.pdf'), width = 12, height = 12)
heatmap.2(hm,Colv=NA,Rowv=NA,scale="none",dendrogram="none",revC = TRUE, trace='none',key=FALSE,ylab = "Ploidy", xlab = "Purity",colsep = c(1:100),rowsep=c(1:100),sepwidth = c(0.001,0.001),sepcolor="black", main = paste0("Sample ",gsub(".*/(.*)/.*","\\1",opt$cov), "\nSD = ", SD.MAP))
dev.off()
saveRDS(hm, file = paste0(opt$outdir, "/heatmap.Rds"))
#wij(heatmap(hm,dendrogram=NULL,Rowv=FALSE,Colv=FALSE),filename=paste0(opt$outdir, '/PurityPloidy.html'))

### This filters the CN matrix for purity ploidy combinations that had a probablity above a certain cutoff.  The default is 0.01.
### CN is a table with the probabilities converted back from log to normal values.  In addition, the probability for a given purity ploidy value is normalized so that they sum to 1.  .SDcols is used to defined what columns are called by the term .SD.  Columns 8:18 contain all of the CN probability data. 

#Make CN matrix
out <- out[prob > opt$prob.cutoff]
setkeyv(out,c("purity", "ploidy"))
setkeyv(CN.matrix,c("purity", "ploidy"))
high.prob.CN.matrix <- merge(out, CN.matrix)
CN<-high.prob.CN.matrix[,(exp(.SD)/rowSums(exp(.SD))), .SDcols = 8:18]
CN.output <- cbind(high.prob.CN.matrix[,.(purity, ploidy, segment, prob)], CN)

saveRDS(CN.output, file = paste0(opt$outdir, "/CN_matrix.Rds"))


### The allele fraction is alpha (K* / K).  The probability of that allele fraction is P(alpha) * P(K).  K is the copy state.
### The K* means that for a given copy number k can take on all values from 1 to K.  So if K is 5 then K* can be 1,2,3,4,5.  If K is 1 then K* can only be 1.  This whole process is to calculate the probability that K is > 0.
### Build a histogram where the allele fraction is the coordination on the x axis and the probability is the height on the y axis.


# P(alpha) * P(K)
# alpha (K* / K)
## Create 2 column matrix now: P(alpha)P(k) and alpha(K*/K)


# Create list before running loop
# Not vectorized, definitely a faster way to do this
segments <- unique(CN.output$segment)
hist.list <- vector("list", length(segments))
alpha.list <- vector("list", length(segments))
beta.list <- vector("list", length(segments))

### Loop through segments first

for(z in 1:length(segments)){

    combo_spec <- CN.output[segment==segments[z]]

## Create 2 column matrix now: P(alpha)P(k) and alpha(K*/K)

### Next loop through copy number states
### %o% is the outer product of two matrices.  This multiples every number in one vector by every number in another
### Prob is P(alpha) * P(k).  combo_spec$prob is P(alpha).  copy state 1 starts at row 6, so to get P(K) you have to do i + 5.
### For loc, you have to do do.call(base:::'%o%' because %o% is also a function in skitools so you have to specify which library it's from
### loc is alpha(K*/K).  You have to transpose it t() and then vectorize it c() so that you can cbind it in the next step
 
    num.CN = 10
    my.list<-vector("list",num.CN)
    for(i in 1:10){
        prob <- rep(combo_spec$prob * combo_spec[[i + 5]], each = i)
        loc <- c(t(do.call(base:::'%o%',list(combo_spec$purity, ((1:i) / i)))))
        out <- cbind(prob, loc)
        my.list[[i]] <- out
    }

### prob2 in hist is just prob normalized so it adds up to 1
### hist.matrix is a sample from hist.  You are drawing 10000 values from loc with replacement and the probability that each value is drawn from loc is given by prob2
### You then fit hist.matrix with a beta distribution using the fitdistr() function.  This is stochastic.  For more details see: http://varianceexplained.org/r/empirical_bayes_baseball/
    
    hist <- data.table(do.call(rbind, my.list))
    hist[,prob2 := prob/sum(prob)]
    hist.matrix = hist[, sample(loc, 10000, prob = prob2, replace = TRUE)]
    beta = fitdistr(hist.matrix, dbeta, start = list(shape1 = 1, shape2 = 20))
    hist.list[[z]] <- hist.matrix
    alpha.list[[z]] <- beta$estimate[1]
    beta.list[[z]] <- beta$estimate[2]
    print(z)
}

### Make the final 2 outputs here: a matrix with the hist.matrix values for each segment (1 column per segment) and a granges object with beta distribution alpha and beta parameters for each segment.


hist.out <- data.table(do.call(cbind, hist.list))
colnames(hist.out) <- as.character(segments)
saveRDS(hist.out, file = paste0(opt$outdir, "/hist_all_seg.Rds"))

alpha.out <- data.table(cbind(alpha.list))
beta.out <- data.table(cbind(beta.list))
alpha.beta <- cbind( segments, alpha.out, beta.out)
setkey(alpha.beta, segments)

setkey(x, segment)

# Do this if you want ALL segments annotated with alpha and beta -> this produces NULL for segments that weren't subsampled if subsapmling was performed
x <- alpha.beta[x]

## #Do this if you want only segments with a non NULL alpha and beta
## x <- x[alpha.beta]

seq.names <- as.data.table(table(x$chr))
grange.out <- GRanges(seqnames = Rle(seq.names$V1, seq.names$N), ranges = IRanges(x$start,end=x$start+opt$bin.red*200-1,width=NULL,names=NULL), ratio = x$data, alpha = as.numeric(x$alpha.list), beta = as.numeric(x$beta.list), segment = as.numeric(x$segment))
saveRDS(grange.out, file = paste0(opt$outdir, "/ab_grange.Rds"))

message('done')
