Practice with edgeR
========================================================


```r
rm(list=ls())
source("http://bioconductor.org/biocLite.R")
```

```
## Bioconductor version 2.14 (BiocInstaller 1.14.2), ?biocLite for
##   help
```

```r
source("http://stevenbsmith.net/source/load_R_enviornment_vars.R")
#biocLite("edgeR")
library("edgeR")
```

```
## Loading required package: limma
```

```r
#edgeRUsersGuide()
setwd(paste("/Users/",user,"/bin/Qualifiers/edgeR",sep=""))
datafile =system.file( "extdata/pasilla_gene_counts.tsv", package="pasilla" )
pasillaCountTable =read.table( datafile, header=TRUE, row.names="gene_id" )
head(pasillaCountTable)
```

```
##             untreated1 untreated2 untreated3 untreated4 treated1 treated2
## FBgn0000003          0          0          0          0        0        0
## FBgn0000008         92        161         76         70      140       88
## FBgn0000014          5          1          0          0        4        0
## FBgn0000015          0          2          1          2        1        0
## FBgn0000017       4664       8714       3564       3150     6205     3072
## FBgn0000018        583        761        245        310      722      299
##             treated3
## FBgn0000003        1
## FBgn0000008       70
## FBgn0000014        0
## FBgn0000015        0
## FBgn0000017     3334
## FBgn0000018      308
```

```r
pasillaDesign = data.frame(
  row.names = colnames( pasillaCountTable ),
  condition = c( "untreated", "untreated", "untreated",
                 "untreated", "treated", "treated", "treated" ),
  libType = c( "single-end", "single-end", "paired-end",
               "paired-end", "single-end", "paired-end", "paired-end" ) )
pasillaDesign
```

```
##            condition    libType
## untreated1 untreated single-end
## untreated2 untreated single-end
## untreated3 untreated paired-end
## untreated4 untreated paired-end
## treated1     treated single-end
## treated2     treated paired-end
## treated3     treated paired-end
```

```r
pairedSamples = pasillaDesign$libType == "paired-end"
countTable = pasillaCountTable[ , pairedSamples ]
condition = pasillaDesign$condition[ pairedSamples ]

head(countTable)
```

```
##             untreated3 untreated4 treated2 treated3
## FBgn0000003          0          0        0        1
## FBgn0000008         76         70       88       70
## FBgn0000014          0          0        0        0
## FBgn0000015          1          2        0        0
## FBgn0000017       3564       3150     3072     3334
## FBgn0000018        245        310      299      308
```

```r
group<-factor(c(1,1,2,2))
```
The following is the bare minimum needed to compare between two groups. Each function will be dissected in turn. Also, will need to compare the results to DESeq at some point.


```r
y <- DGEList(counts=countTable,group=group)
y <- calcNormFactors(y)
y <- estimateCommonDisp(y)
y <- estimateTagwiseDisp(y)
et <- exactTest(y)
topTags(et)
```

```
## Comparison of groups:  2-1 
##              logFC logCPM     PValue        FDR
## FBgn0039155 -4.378  5.588 1.989e-183 2.903e-179
## FBgn0003360 -2.961  8.059 2.725e-156 1.989e-152
## FBgn0025111  2.943  7.159 3.007e-154 1.463e-150
## FBgn0026562 -2.447 11.903 1.956e-106 7.138e-103
## FBgn0039827 -4.129  4.281 1.647e-105 4.808e-102
## FBgn0035085 -2.499  5.542  1.618e-96  3.936e-93
## FBgn0029167 -2.226  8.063  4.295e-93  8.958e-90
## FBgn0000071  2.565  5.034  5.474e-79  9.989e-76
## FBgn0029896 -2.546  5.132  1.268e-77  2.057e-74
## FBgn0034897 -2.062  6.097  3.261e-75  4.761e-72
```

#### Calc norm factors:

Let's look under the hood of calcNormFactors:


```r
calcNormFactors
```

```
## function (object, method = c("TMM", "RLE", "upperquartile", "none"), 
##     refColumn = NULL, logratioTrim = 0.3, sumTrim = 0.05, doWeighting = TRUE, 
##     Acutoff = -1e+10, p = 0.75) 
## {
##     if (is(object, "DGEList")) {
##         x <- as.matrix(object$counts)
##         lib.size <- object$samples$lib.size
##     }
##     else {
##         x <- as.matrix(object)
##         lib.size <- colSums(x)
##     }
##     method <- match.arg(method)
##     allzero <- rowSums(x > 0) == 0
##     if (any(allzero)) 
##         x <- x[!allzero, , drop = FALSE]
##     if (nrow(x) == 0 || ncol(x) == 1) 
##         method = "none"
##     f <- switch(method, TMM = {
##         f75 <- .calcFactorQuantile(data = x, lib.size = lib.size, 
##             p = 0.75)
##         if (is.null(refColumn)) refColumn <- which.min(abs(f75 - 
##             mean(f75)))
##         if (length(refColumn) == 0 | refColumn < 1 | refColumn > 
##             ncol(x)) refColumn <- 1
##         f <- rep(NA, ncol(x))
##         for (i in 1:ncol(x)) f[i] <- .calcFactorWeighted(obs = x[, 
##             i], ref = x[, refColumn], libsize.obs = lib.size[i], 
##             libsize.ref = lib.size[refColumn], logratioTrim = logratioTrim, 
##             sumTrim = sumTrim, doWeighting = doWeighting, Acutoff = Acutoff)
##         f
##     }, RLE = .calcFactorRLE(x)/lib.size, upperquartile = .calcFactorQuantile(x, 
##         lib.size, p = p), none = rep(1, ncol(x)))
##     f <- f/exp(mean(log(f)))
##     if (is(object, "DGEList")) {
##         object$samples$norm.factors <- f
##         return(object)
##     }
##     else {
##         return(f)
##     }
## }
## <environment: namespace:edgeR>
```

```r
object<-DGEList(counts=countTable,group=group)
#Defaults:
method = c("TMM", "RLE", "upperquartile", "none")
refColumn = NULL
logratioTrim = 0.3
sumTrim = 0.05
doWeighting = TRUE
Acutoff = -1e+10
p = 0.75
```
The first few lines process/clean and QC the input variables. 
Checks to make sure that the data is  DGEList object, then converts the counts into matrix form. Stores library size, as calculated during he DEGList conversion. The lib size appears to be just the colSums of counts. Check this:
colSums: 8.3584 &times; 10<sup>6</sup>, 9.8413 &times; 10<sup>6</sup>, 9.5718 &times; 10<sup>6</sup>, 1.0344 &times; 10<sup>7</sup>
lib.size: 8.3584 &times; 10<sup>6</sup>, 9.8413 &times; 10<sup>6</sup>, 9.5718 &times; 10<sup>6</sup>, 1.0344 &times; 10<sup>7</sup>
Also note that if the lib size wasnt calculated using the DGEList method, it does it here

```r
if (is(object, "DGEList")) {
        x <- as.matrix(object$counts)
        lib.size <- object$samples$lib.size
    }else {
        x <- as.matrix(object)
        lib.size <- colSums(x)
    }
```
The first line of the following code gets executed implicity during a function call, but had to be modified for the purposes here. The second rgument, choices=method, is not in the function but this is what happens implicilt within the function call. The match.arg simply matches a given method to a list of potential choices. In this case, it takes the first element of method (4 elemtns) matches to the first (TMM) and assigns the signle element TMM as the method variable. The choices=method statement is the default vals of method in the function declaration. Confusing.  

The next lines:
* calcualte the rows (genes) that have zero counts across all samples and keep only the ones with at least one non zero gene count
* then, if there are no genes left or only one sample, set normalization method to none
* execute normaliztion method depending on the value of method variable (switch statement). Store the resulting normlizaton as f. Dig into the TMM and other normalization methods later. 
* finally, dependning on wether the inital input was a DEGobject, return the normalizaion factors either as a vector or within the original object as object$samples$norm.factors (return the whole updated object to the original calcNormFactors call)

The actual calls to each method seem to be hidden or private. will need to look into how to access them.. i.e., .calcFactorQuantile

```r
    method <- match.arg(method,choices=method)
    allzero <- rowSums(x > 0) == 0
    if (any(allzero)) x <- x[!allzero, , drop = FALSE]
    if (nrow(x) == 0 || ncol(x) == 1) method = "none"
#Depending on method...
    f <- switch(method, 
      TMM = {
        f75 <- .calcFactorQuantile(data = x, lib.size = lib.size,p = 0.75)
        if (is.null(refColumn)) refColumn <- which.min(abs(f75 -mean(f75)))
        if (length(refColumn) == 0 | refColumn < 1 | refColumn > ncol(x)) refColumn <- 1
        f <- rep(NA, ncol(x))
        for (i in 1:ncol(x)) f[i] <- .calcFactorWeighted(obs = x[, 
            i], ref = x[, refColumn], libsize.obs = lib.size[i], 
            libsize.ref = lib.size[refColumn], logratioTrim = logratioTrim, 
            sumTrim = sumTrim, doWeighting = doWeighting, Acutoff = Acutoff)
        f
    }, 
      RLE = .calcFactorRLE(x)/lib.size, 
      upperquartile = .calcFactorQuantile(x,lib.size, p = p), 
      none = rep(1, ncol(x))
    ) ## end of switch statment
```

```
## Error: could not find function ".calcFactorQuantile"
```

```r
f <- f/exp(mean(log(f))) #what is this?
```

```
## Error: object 'f' not found
```

```r
## Depending on object type, return within the object or sumply as a vector. 
    if (is(object, "DGEList")) {
        object$samples$norm.factors <- f
        return(object)
    }else {
        return(f)
    }
```

```
## Error: object 'f' not found
```

#### Calc norm factors:

#### estimateCommonDisp:
#### estimateTagwiseDisp:
#### exactTest:

```r
pair = 1:2
dispersion = "auto"
rejection.region = "doubletail"
big.count = 900
prior.count = 0.125
```

```r
    if (!is(object, "DGEList")) 
        stop("Currently only supports DGEList objects as the object argument.")
    if (length(pair) != 2) 
        stop("Pair must be of length 2.")
    rejection.region <- match.arg(rejection.region, c("doubletail", 
        "deviance", "smallp"))
    group <- as.factor(object$samples$group)
    levs.group <- levels(group)
    if (is.numeric(pair)) pair <- levs.group[pair] else pair <- as.character(pair)
    if (!all(pair %in% levs.group)) stop("At least one element of given pair is not a group.\n Groups are: ",paste(levs.group, collapse = " "))
    if (is.null(dispersion)) 
        dispersion <- "auto"
    if (is.character(dispersion)) {
        dispersion <- match.arg(dispersion, c("auto", "common","trended", "tagwise"))
        dispersion <- switch(dispersion, common = object$common.dispersion, 
            trended = object$trended.dispersion, tagwise = object$tagwise.dispersion, 
            auto = getDispersion(object))
        if (is.null(dispersion)) 
            stop("specified dispersion not found in object")
        if (is.na(dispersion[1])) 
            stop("dispersion is NA")
    }
```

```
## Error: specified dispersion not found in object
```

```r
    ldisp <- length(dispersion)
    ntags <- nrow(object$counts)
    if (ldisp != 1 && ldisp != ntags) 
        stop("Dispersion provided by user must have length either 1 or the number of tags in the DGEList object.")
```

```
## Error: Dispersion provided by user must have length either 1 or the number
## of tags in the DGEList object.
```

```r
    if (ldisp == 1) 
        dispersion <- rep(dispersion, ntags)
```
Beggining of the real meat:
* First determine the group to include in the calculation and keep that as y
* Store the library size and normalization factors for that group as respecive vars.

```r
group <- as.character(group)
    j <- group %in% pair
    y <- object$counts[, j, drop = FALSE]
    lib.size <- object$samples$lib.size[j]
    norm.factors <- object$samples$norm.factors[j]
    group <- group[j]
    if (is.null(rownames(y))) 
        rownames(y) <- paste("tag", 1:ntags, sep = ".")
```
* Adjust library size by normalization factors
* Use this adjusted lib size as offset
* Calculate prior count and offset aug... *look into what theyre doing here*.

```r
    lib.size <- lib.size * norm.factors
    offset <- log(lib.size)
    lib.size.average <- exp(mean(offset))
    prior.count <- prior.count * lib.size/mean(lib.size)
    offset.aug <- log(lib.size + 2 * prior.count)
    j1 <- group == pair[1]
    n1 <- sum(j1)
    if (n1 == 0) 
        stop("No libraries for", pair[1])
```
* Calculate abundance? (by some mysterious other function mglmOneGroup, see below)

```r
    y1 <- y[, j1, drop = FALSE]
    abundance1 <- mglmOneGroup(y1 + matrix(prior.count[j1], ntags, 
        n1, byrow = TRUE), offset = offset.aug[j1], dispersion = dispersion)
```

```
## Error: dispersion not floating point number
```

```r
    j2 <- group == pair[2]
    n2 <- sum(j2)
    if (n1 == 0) 
        stop("No libraries for", pair[2])
    y2 <- y[, j2, drop = FALSE]
    abundance2 <- mglmOneGroup(y2 + matrix(prior.count[j2], ntags, 
        n2, byrow = TRUE), offset = offset.aug[j2], dispersion = dispersion)
```

```
## Error: dispersion not floating point number
```

```r
    logFC <- (abundance2 - abundance1)/log(2)
```

```
## Error: object 'abundance2' not found
```

```r
    abundance <- mglmOneGroup(y, dispersion = dispersion, offset = offset)
```

```
## Error: dispersion not floating point number
```

```r
    e <- exp(abundance)
```

```
## Error: object 'abundance' not found
```

```r
    input.mean <- matrix(e, ntags, n1)
```

```
## Error: object 'e' not found
```

```r
    output.mean <- input.mean * lib.size.average
```

```
## Error: object 'input.mean' not found
```

```r
    input.mean <- t(t(input.mean) * lib.size[j1])
```

```
## Error: object 'input.mean' not found
```

```r
    y1 <- q2qnbinom(y1, input.mean = input.mean, output.mean = output.mean, 
        dispersion = dispersion)
```

```
## Error: object 'input.mean' not found
```

```r
    input.mean <- matrix(e, ntags, n2)
```

```
## Error: object 'e' not found
```

```r
    output.mean <- input.mean * lib.size.average
```

```
## Error: object 'input.mean' not found
```

```r
    input.mean <- t(t(input.mean) * lib.size[j2])
```

```
## Error: object 'input.mean' not found
```

```r
    y2 <- q2qnbinom(y2, input.mean = input.mean, output.mean = output.mean, 
        dispersion = dispersion)
```

```
## Error: object 'input.mean' not found
```
P values calcualted using one of several methods. Default is doubletail (inside switch statement) via the exactTestDoubleTail function. Need to look into this. 




