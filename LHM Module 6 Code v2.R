#LHM
#11/7/18
#Advanced Topics in Bioinformatics Module 6

#Packages needed:
library(ggplot2)
library(ggbio)
library(Homo.sapiens)
library(Gviz)
library(GenomicRanges)
library(TxDb.Hsapiens.UCSC.hg19.knownGene)
library(org.Hs.eg.db)

#Regions we're interested in (TADs with 8 CTCF sites):
#GRanges object with 8 ranges and 1 metadata column:
  #      seqnames              ranges strand | orientation
  #         <Rle>           <IRanges>  <Rle> | <character>
  #  [1]     chr1   41960001-41999999      * |        forw
  #  [2]     chr2 220360001-220399999      * |        <NA>
  #  [3]    chr16   54960001-54999999      * |        <NA>
  #  [4]    chr17   61480001-61519999      * |        <NA>
  #  [5]    chr19     3520001-3559999      * |        <NA>
  #  [6]    chr20     3640001-3679999      * |        <NA>
  #  [7]    chr20   43920001-43959999      * |        <NA>
  #  [8]    chr22   39680001-39719999      * |        <NA>
  
#Summary of different options for making these tracks:


#Option 1: ggbio package
#Can get gene names with ease, but doesn't look quite as nice with default settings, or near-default settings, and seems like it will take more effort to make it look nice with additional tracks. The default settings also seem to be based off of gene names, rather than coordinates; if we want to map a coordinate range, we'll have to look more into how to do that. Overall, I'm happier with the Gviz results than with the ggbio ones.

#The chromosome 1 region of interest contains part of a single gene, HIVEP3. Plotting the region just for this gene is shown below. Can expand to other genes on the list, but will need to either (a) look up the genes in each of those regions and use genesymbol [c("first gene", "last gene")] for our regin of interest or (b) see if there is a way to do this using genome coordinates:

data(genesymbol, package = "biovizBase")
wh <- genesymbol[c("HIVEP3")]
wh <- range(wh, ignore.strand = TRUE)
p.txdb <- autoplot(Homo.sapiens, which  = wh, stat="reduce") #Reduce collapses all of the transcripts onto one line
p.txdb


#Option 2: Gviz package
#Can't figure out how to get the gene names to show up in the "dense" version. They only show up (as transcript names) in the non-dense version.

#This is an example of the approach for chromosome 1. The same approach could be used for everything on our list of interest.

#Chromosome ideogram: 
itrack1 <- IdeogramTrack(genome = "hg19", chromosome = 1)

#Genome axis track (shows range of the genome that the next track I'm making is covering):
gaTrack <- GenomeAxisTrack()

#Making the tracks for our genes of interest. As with the approach above, this strategy looks at transcripts, which we can then collapse. Note that names do not seem to display when the transcripts are collapsed. Also, in the non-collapsed version, names only display if the left side of the transcript is visible.
txdb <- TxDb.Hsapiens.UCSC.hg19.knownGene

gTrackDense <- GeneRegionTrack(txdb, chromosome = 1, genome = "hg19", name = "Gene Track", transcriptAnnotation = "symbol", background.title = "darkblue", stacking = "dense")

gTrackOriginal <- GeneRegionTrack(txdb, chromosome = 1, genome = "hg19", name = "Gene Track", transcriptAnnotation = "symbol", background.title = "darkblue")

#Replacing the transcript names with desired gene names (since the txdb object doesn't use gene names):
symbolsDense <- unlist(mapIds(org.Hs.eg.db, gene(gTrackDense), "SYMBOL", "ENTREZID", multiVals = "first"))
symbol(gTrackDense) <- symbolsDense[gene(gTrackDense)]

symbolsOriginal <- unlist(mapIds(org.Hs.eg.db, gene(gTrackOriginal), "SYMBOL", "ENTREZID", multiVals = "first"))
symbol(gTrackOriginal) <- symbolsOriginal[gene(gTrackOriginal)]

#Dense version focused on 40 kb range of interest:
plotTracks(list(itrack1, gaTrack, gTrackDense), from=(41960001), to=(41999999))

#Dense version with 400 kb extra on each side:
plotTracks(list(itrack1, gaTrack, gTrackDense), from=(41960001-400000), to=(41999999+400000))

#Non-dense version focused on 40 kb range of interest:
plotTracks(list(itrack1, gaTrack, gTrackOriginal), from=(41960001), to=(41999999))

#Non-dense version with 400 kb extra on each side:
plotTracks(list(itrack1, gaTrack, gTrackOriginal), from=(41960001-400000), to=(41999999+400000))



#Same approach for each of the other 7 loci of interest:

#  [2]     chr2 220360001-220399999
#Genes within the 40 kb: GMPPA, ASIC4
itrack2 <- IdeogramTrack(genome = "hg19", chromosome = 2)

gTrackDense <- GeneRegionTrack(txdb, chromosome = 2, genome = "hg19", name = "Gene Track", transcriptAnnotation = "symbol", background.title = "darkblue", stacking = "dense")

gTrackOriginal <- GeneRegionTrack(txdb, chromosome = 2, genome = "hg19", name = "Gene Track", transcriptAnnotation = "symbol", background.title = "darkblue")

symbolsDense <- unlist(mapIds(org.Hs.eg.db, gene(gTrackDense), "SYMBOL", "ENTREZID", multiVals = "first"))
symbol(gTrackDense) <- symbolsDense[gene(gTrackDense)]

symbolsOriginal <- unlist(mapIds(org.Hs.eg.db, gene(gTrackOriginal), "SYMBOL", "ENTREZID", multiVals = "first"))
symbol(gTrackOriginal) <- symbolsOriginal[gene(gTrackOriginal)]

#Dense version focused on 40 kb range of interest:
plotTracks(list(itrack2, gaTrack, gTrackDense), from=(220360001), to=(220399999))

#Dense version with 400 kb extra on each side:
plotTracks(list(itrack2, gaTrack, gTrackDense), from=(220360001-400000), to=(220399999+400000))

#Non-dense version focused on 40 kb range of interest:
plotTracks(list(itrack2, gaTrack, gTrackOriginal), from=(220360001), to=(220399999))

#Non-dense version with 400 kb extra on each side:
plotTracks(list(itrack2, gaTrack, gTrackOriginal), from=(220360001-400000), to=(220399999+400000))



#  [3]    chr16   54960001-54999999
#Genes within the 40 kb window: CRNDE (partial, off left side of screen) and IRX5
itrack16 <- IdeogramTrack(genome = "hg19", chromosome = 16)

gTrackDense <- GeneRegionTrack(txdb, chromosome = 16, genome = "hg19", name = "Gene Track", transcriptAnnotation = "symbol", background.title = "darkblue", stacking = "dense")

gTrackOriginal <- GeneRegionTrack(txdb, chromosome = 16, genome = "hg19", name = "Gene Track", transcriptAnnotation = "symbol", background.title = "darkblue")

symbolsDense <- unlist(mapIds(org.Hs.eg.db, gene(gTrackDense), "SYMBOL", "ENTREZID", multiVals = "first"))
symbol(gTrackDense) <- symbolsDense[gene(gTrackDense)]

symbolsOriginal <- unlist(mapIds(org.Hs.eg.db, gene(gTrackOriginal), "SYMBOL", "ENTREZID", multiVals = "first"))
symbol(gTrackOriginal) <- symbolsOriginal[gene(gTrackOriginal)]

#Dense version focused on 40 kb range of interest:
plotTracks(list(itrack16, gaTrack, gTrackDense), from=(54960001), to=(54999999))

#Dense version with 400 kb extra on each side:
plotTracks(list(itrack16, gaTrack, gTrackDense), from=(54960001-400000), to=(54999999+400000))

#Non-dense version focused on 40 kb range of interest:
plotTracks(list(itrack16, gaTrack, gTrackOriginal), from=(54960001), to=(54999999))

#Non-dense version with 400 kb extra on each side:
plotTracks(list(itrack16, gaTrack, gTrackOriginal), from=(54960001-400000), to=(54999999+400000))



#  [4]    chr17   61480001-61519999
#Genes within the 40 kb window: TANC2 (partial, off left side of screen), CYB561
itrack17 <- IdeogramTrack(genome = "hg19", chromosome = 17)

gTrackDense <- GeneRegionTrack(txdb, chromosome = 17, genome = "hg19", name = "Gene Track", transcriptAnnotation = "symbol", background.title = "darkblue", stacking = "dense")

gTrackOriginal <- GeneRegionTrack(txdb, chromosome = 17, genome = "hg19", name = "Gene Track", transcriptAnnotation = "symbol", background.title = "darkblue")

symbolsDense <- unlist(mapIds(org.Hs.eg.db, gene(gTrackDense), "SYMBOL", "ENTREZID", multiVals = "first"))
symbol(gTrackDense) <- symbolsDense[gene(gTrackDense)]

symbolsOriginal <- unlist(mapIds(org.Hs.eg.db, gene(gTrackOriginal), "SYMBOL", "ENTREZID", multiVals = "first"))
symbol(gTrackOriginal) <- symbolsOriginal[gene(gTrackOriginal)]

#Dense version focused on 40 kb range of interest:
plotTracks(list(itrack17, gaTrack, gTrackDense), from=(61480001), to=(61519999))

#Dense version with 400 kb extra on each side:
plotTracks(list(itrack17, gaTrack, gTrackDense), from=(61480001-400000), to=(61519999+400000))

#Non-dense version focused on 40 kb range of interest:
plotTracks(list(itrack17, gaTrack, gTrackOriginal), from=(61480001), to=(61519999))

#Non-dense version with 400 kb extra on each side:
plotTracks(list(itrack17, gaTrack, gTrackOriginal), from=(61480001-400000), to=(61519999+400000))



#  [5]    chr19     3520001-3559999
#Genes within the 40 kb window: FZR1, C19orf71, MFSD12. Something that didn't map to a gene symbol came back as well ("NA")
itrack19 <- IdeogramTrack(genome = "hg19", chromosome = 19)

gTrackDense <- GeneRegionTrack(txdb, chromosome = 19, genome = "hg19", name = "Gene Track", transcriptAnnotation = "symbol", background.title = "darkblue", stacking = "dense")

gTrackOriginal <- GeneRegionTrack(txdb, chromosome = 19, genome = "hg19", name = "Gene Track", transcriptAnnotation = "symbol", background.title = "darkblue")

symbolsDense <- unlist(mapIds(org.Hs.eg.db, gene(gTrackDense), "SYMBOL", "ENTREZID", multiVals = "first"))
symbol(gTrackDense) <- symbolsDense[gene(gTrackDense)]

symbolsOriginal <- unlist(mapIds(org.Hs.eg.db, gene(gTrackOriginal), "SYMBOL", "ENTREZID", multiVals = "first"))
symbol(gTrackOriginal) <- symbolsOriginal[gene(gTrackOriginal)]

#Dense version focused on 40 kb range of interest:
plotTracks(list(itrack19, gaTrack, gTrackDense), from=(3520001), to=(3559999))

#Dense version with 400 kb extra on each side:
plotTracks(list(itrack19, gaTrack, gTrackDense), from=(3520001-400000), to=(3559999+400000))

#Non-dense version focused on 40 kb range of interest:
plotTracks(list(itrack19, gaTrack, gTrackOriginal), from=(3520001), to=(3559999))

#Non-dense version with 400 kb extra on each side:
plotTracks(list(itrack19, gaTrack, gTrackOriginal), from=(3520001-400000), to=(3559999+400000))



#  [6]    chr20     3640001-3679999
#Genes within the 40 kb window: GFRA4 (disappears off left side), ADAM33, SIGLEC1
itrack20 <- IdeogramTrack(genome = "hg19", chromosome = 20)

gTrackDense <- GeneRegionTrack(txdb, chromosome = 20, genome = "hg19", name = "Gene Track", transcriptAnnotation = "symbol", background.title = "darkblue", stacking = "dense")

gTrackOriginal <- GeneRegionTrack(txdb, chromosome = 20, genome = "hg19", name = "Gene Track", transcriptAnnotation = "symbol", background.title = "darkblue")

symbolsDense <- unlist(mapIds(org.Hs.eg.db, gene(gTrackDense), "SYMBOL", "ENTREZID", multiVals = "first"))
symbol(gTrackDense) <- symbolsDense[gene(gTrackDense)]

symbolsOriginal <- unlist(mapIds(org.Hs.eg.db, gene(gTrackOriginal), "SYMBOL", "ENTREZID", multiVals = "first"))
symbol(gTrackOriginal) <- symbolsOriginal[gene(gTrackOriginal)]

#Dense version focused on 40 kb range of interest:
plotTracks(list(itrack20, gaTrack, gTrackDense), from=(3640001), to=(3679999))

#Dense version with 400 kb extra on each side:
plotTracks(list(itrack20, gaTrack, gTrackDense), from=(3640001-400000), to=(3679999+400000))

#Non-dense version focused on 40 kb range of interest:
plotTracks(list(itrack20, gaTrack, gTrackOriginal), from=(3640001), to=(3679999))

#Non-dense version with 400 kb extra on each side:
plotTracks(list(itrack20, gaTrack, gTrackOriginal), from=(3640001-400000), to=(3679999+400000))



#  [7]    chr20   43920001-43959999 #Run #6 first
#Genes within the 40 kb window: MATN4, RBPJL, SDC4. One "NA" transcript shows up as well.

#Dense version focused on 40 kb range of interest:
plotTracks(list(itrack20, gaTrack, gTrackDense), from=(43920001), to=(43959999))

#Dense version with 400 kb extra on each side:
plotTracks(list(itrack20, gaTrack, gTrackDense), from=(43920001-400000), to=(43959999+400000))

#Non-dense version focused on 40 kb range of interest:
plotTracks(list(itrack20, gaTrack, gTrackOriginal), from=(43920001), to=(43959999))

#Non-dense version with 400 kb extra on each side:
plotTracks(list(itrack20, gaTrack, gTrackOriginal), from=(43920001-400000), to=(43959999+400000))



#  [8]    chr22   39680001-39719999
#Genes within the 40 kb window: RPL3, SNORD43
itrack22 <- IdeogramTrack(genome = "hg19", chromosome = 22)

gTrackDense <- GeneRegionTrack(txdb, chromosome = 22, genome = "hg19", name = "Gene Track", transcriptAnnotation = "symbol", background.title = "darkblue", stacking = "dense")

gTrackOriginal <- GeneRegionTrack(txdb, chromosome = 22, genome = "hg19", name = "Gene Track", transcriptAnnotation = "symbol", background.title = "darkblue")

symbolsDense <- unlist(mapIds(org.Hs.eg.db, gene(gTrackDense), "SYMBOL", "ENTREZID", multiVals = "first"))
symbol(gTrackDense) <- symbolsDense[gene(gTrackDense)]

symbolsOriginal <- unlist(mapIds(org.Hs.eg.db, gene(gTrackOriginal), "SYMBOL", "ENTREZID", multiVals = "first"))
symbol(gTrackOriginal) <- symbolsOriginal[gene(gTrackOriginal)]

#Dense version focused on 40 kb range of interest:
plotTracks(list(itrack22, gaTrack, gTrackDense), from=(39680001), to=(39719999))

#Dense version with 400 kb extra on each side:
plotTracks(list(itrack22, gaTrack, gTrackDense), from=(39680001-400000), to=(39719999+400000))

#Non-dense version focused on 40 kb range of interest:
plotTracks(list(itrack22, gaTrack, gTrackOriginal), from=(39680001), to=(39719999))

#Non-dense version with 400 kb extra on each side:
plotTracks(list(itrack22, gaTrack, gTrackOriginal), from=(39680001-400000), to=(39719999+400000))



#Updated coordinates per Cooper's request to match his figure:
#start of tad1 = 35640001
#end of tad1 = 36040000
#start of tad2 = 36080001
#end of tad2 = 36680000
#start of tad3 = 36720001
#end of tad3 = 37360000
itrack5 <- IdeogramTrack(genome = "hg19", chromosome = 5)

gTrackOriginal <- GeneRegionTrack(txdb, chromosome = 5, genome = "hg19", name = "Gene Track", transcriptAnnotation = "symbol", background.title = "darkblue")
symbolsOriginal <- unlist(mapIds(org.Hs.eg.db, gene(gTrackOriginal), "SYMBOL", "ENTREZID", multiVals = "first"))
symbol(gTrackOriginal) <- symbolsOriginal[gene(gTrackOriginal)]

#Want to look at regions covering 3 TADs and their boundaries. Looking at the non-dense version for each of these:

#Full region: #Showing this
plotTracks(list(itrack5, gaTrack, gTrackOriginal), from=(35640001), to=(37360000))
#fullTad1to3 <-plotTracks(list(itrack5, gaTrack, gTrackOriginal), from=(35640001), to=(37360000))

#Boundary 1 (end of TAD1 to start of TAD2): #Showing this
plotTracks(list(itrack5, gaTrack, gTrackOriginal), from=(36040000), to=(36080001)) #UGT3A2

#Boundary 2 (end of TAD2 to start of TAD3): #Showing this
plotTracks(list(itrack5, gaTrack, gTrackOriginal), from=(36680000), to=(36720001)) #Part of SLC1A3

#40 kb region in the middle of TAD2:
mid <- (36080001 + 36680000) / 2
plotTracks(list(itrack5, gaTrack, gTrackOriginal), from=(mid-20000), to=(mid+20000))

#Inside TAD2: #Showing this
plotTracks(list(itrack5, gaTrack, gTrackOriginal), from=(36080001), to=(36680000))

#Another 40 kb region:
plotTracks(list(itrack5, gaTrack, gTrackOriginal), from=(mid-200000), to=(mid-160000))

#Enhancer inside TAD2: 36605831

