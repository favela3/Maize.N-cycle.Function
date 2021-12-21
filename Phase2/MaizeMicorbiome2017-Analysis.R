##Alonso Favela 
#Maize Microbiome 2017 study ====

##Load Programs ==
setwd("/Maize Micorbiome Analysis ")


library(phyloseq)
library(ggplot2)
library(vegan)
library(tidyverse)


######Phyloseq and Plotting:  BELOW===============
biom_file = "Data/otu_table_rarefied_100000.biom"

map_file = "Data/2017MappingFile.txt"
biomOTU = import_biom(biom_file,treefilename = 'Data/rep_seq_aligned_pfiltered.tre', 
                      parseFunction = parse_taxonomy_greengenes)

Map = import_qiime_sample_data("Data/2017MappingFile-R.txt")

OTU16S = merge_phyloseq(biomOTU, Map)


#Removing Outlier files 
 OTU16S <- subset_samples(OTU16S, X.SampleID != "MM17.T2.041")
 OTU16S <- subset_samples(OTU16S, X.SampleID != "MM17.T3.037")
 OTU16S <- subset_samples(OTU16S, X.SampleID != "MM17.T2.027")
 OTU16S <- subset_samples(OTU16S, X.SampleID != "MM17.T2.040")
 OTU16S <- subset_samples(OTU16S, X.SampleID != "MM17.T3.063")
 

##Alpha Diversity of 16S ====
#measurment of alpha diversity
# c("Observed", "Chao1", "ACE", "Shannon", "Simpson", "InvSimpson", "Fisher").
Richness<-plot_richness(OTU16S, x = "Time", measures = "Chao1") +
  theme(text=element_text(family = "Helvetica",size=14))+geom_boxplot()+theme_bw()
Richness

Richness<-plot_richness(OTU16S, x = "Time", measures = "Chao1") +
  theme(text=element_text(family = "Helvetica",size=14))+geom_boxplot()+theme_classic()
Richness

pdf("Output/MC-DiversityinTime.pdf",width=6,height=4,paper='special') 
print(Richness)
dev.off()


Richness<-plot_richness(OTU16S, x = "Genotype", measures = "Chao1") +
  theme(text=element_text(family = "Helvetica",size=14))+geom_boxplot()

pdf("Output/MC-DiversityinGenotypes.pdf",width=6,height=4,paper='special') 
print(Richness)
dev.off()


Richness<-plot_richness(OTU16S, x = "Type", measures = "Chao1") +
  theme(text=element_text(family = "Helvetica",size=14))+geom_boxplot()

pdf("Output/MC-DiversityinTYPE.pdf",width=6,height=4,paper='special') 
print(Richness)
dev.off()

Richness<-plot_richness(OTU16S, x = "GxTFactor", measures = "Chao1") +
  theme(text=element_text(family = "Helvetica",size=14))+geom_boxplot()

pdf("Output/MC-DiversityinGxTfactor.pdf",width=6,height=4,paper='special') 
print(Richness)
dev.off()


Richness<-plot_richness(OTU16S, x = "Year", measures = "Chao1") +
  theme(text=element_text(family = "Helvetica",size=14))+geom_boxplot()
print(Richness)

pdf("Output/MC-Diversityinyearfactor.pdf",width=6,height=4,paper='special') 
print(Richness)
dev.off()

#Beta-diversity of 16S community ====
GP.ord <- ordinate(OTU16S, "NMDS", "bray",)
p1 = plot_ordination(OTU16S, GP.ord, type="Samples", color="Time",  title="Bacterial 16S Ordination")+geom_point(size=3)+geom_text(mapping = aes(label = X.SampleID), size = 2, vjust = 1.5)+theme_bw()
print(p1)
pdf("Output/16STime.pdf",width=6,height=4,paper='special') 
print(p1)
dev.off()

p1 = plot_ordination(OTU16S, GP.ord, type="Samples", color="Genotype",  title="Bacterial 16S Ordination")+geom_point(size=3)+geom_text(mapping = aes(label = X.SampleID), size = 2, vjust = 1.5)+
  theme(legend.position = "none")+theme_bw()
print(p1)
pdf("Output/GenotypeCloud.pdf",width=6,height=4,paper='special') 
print(p1)
dev.off()

##Inbred Time analysis ====
OTU16S <- subset_samples(OTU16S, Type == "Inbred")

GP.ord <- ordinate(OTU16S, "NMDS", "bray",)
p1 = plot_ordination(OTU16S, GP.ord, type="Samples", color="Time",  title="Bacterial 16S Ordination")+geom_point(size=3)+
  geom_text(mapping = aes(label = GxTFactor), size = 2, vjust = 1.5)+
  stat_ellipse(aes(color=Year, group=Year),type = "t")
print(p1)

p1 = plot_ordination(OTU16S, GP.ord, type="Samples", color="Year",  title="Bacterial 16S Ordination")+geom_point(size=3)+
  geom_text(mapping = aes(label = GxTFactor), size = 2, vjust = 1.5)+
  stat_ellipse(aes(color=Time, group=Time),type = "t")
print(p1)

#as time progresses the diversity in the soil microbial community decreses 
p1 = plot_ordination(OTU16S, GP.ord, type="Samples", color="Heterotic.Group",  title="Bacterial 16S Ordination")+geom_point(size=3)+geom_text(mapping = aes(label = X.SampleID), size = 2, vjust = 1.5)
print(p1)


p1 = plot_ordination(OTU16S, GP.ord, type="Samples", color="GxTFactor",  title="Bacterial 16S Ordination")+geom_point(size=3)+
geom_text(mapping = aes(label = GxTFactor), size = 2, vjust = 1.5)+theme(legend.position = "none")
print(p1)

pdf("Output/16SGxT.pdf",width=6,height=4,paper='special') 
print(p1)
dev.off()

p1 = plot_ordination(OTU16S, GP.ord, type="Samples", color="GxTFactor",  title="Bacterial 16S Ordination")+geom_point(size=3)+
  geom_text(mapping = aes(label = GxTFactor), size = 2, vjust = 1.5)+theme(legend.position = "none")
print(p1)

###PLOTTING DIFFERENT FUNCTIONAL GROUP via nitri, meth taxa----
GenotypeAve = merge_samples(OTU16S, "Type", fun = mean) #averages each OTU in all samples belonging to each habitat class
sample_data(GenotypeAve)$Type <- levels(sample_data(OTU16S)$Type)
#relative byt total bumber of genes
GenotypeAve_trans = transform_sample_counts(GenotypeAve, function(x) 100 * x/sum(x) )

# plot_bar(GenotypeAve_trans, "Type", fill = "Phylum")+labs(title = " 16S Abundance Genotypic Ave")+ylab("Relative %")+xlab("Functional Group")+ theme(text=element_text(family = "Helvetica",size=14))+theme_bw()+geom_bar(stat="identity")

####Type Composition Plot Dec.19th====
plot_bar(GenotypeAve_trans, "Type", fill = "Phylum")+
  labs(title = " 16S Abundance Genotypic Ave")+ylab("Relative Percentage")+
  xlab("Functional Group")+ theme(text=element_text(family = "Helvetica",size=14))+
  theme_bw()+
  geom_bar(stat="identity")+##This bit of code should override the black border lines and remove the idenity
  stat_summary(fun=sum,geom="bar")+ 
  theme(legend.position="none")+
  scale_y_continuous(labels = scales::percent_format(scale = 1))
  #scale_y_continuous(labels = scales::percent) ##This Puts scale in precentage



#subset
#NitrifersOTU=subset_taxa(Inbred, Phylum == "Nitrospirae"|Order=="Nitrosomonadales"|Order=="Nitrososphaerales")

NitrifersOTUGenotype=subset_taxa(GenotypeAve_trans, Phylum == "Nitrospirae"|Order=="Nitrosomonadales"|Order=="Nitrososphaerales")

NitrifersOTUGenotype@sam_data
plot_bar(NitrifersOTUGenotype, "Type", fill = "Order", facet_grid = ~Class)+labs(title = "Nitrifers 16S Abundance Type Ave")+ylab("Relative Percentage")+xlab("Functional Group")+ theme(text=element_text(family = "Helvetica",size=14))+theme_bw()#+geom_bar(stat="identity")
##1- inbred
##2- Hybrid
##3- Teosinte
###The only issue is that there is sample unbiasis 
###Nitrosomoas 
plot_bar(NitrifersOTUGenotype, "Genotype", fill = "Order", facet_grid = ~Class)+labs(title = "Nitrifers 16S Abundance Genotypic Ave")+ylab("Relative %")+xlab("Functional Group")+ theme(text=element_text(family = "Helvetica",size=14))+theme_bw()#+geom_bar(stat="identity")

OTU16S
Methanobrevibacter
NitrifersOTUGenotype=subset_taxa(GenotypeAve_trans, Phylum == "Methanobrevibacter"|Order=="Nitrosomonadales"|Order=="Nitrososphaerales")

MethOTUGenotype=subset_taxa(GenotypeAve_trans, Order == "Metholotrpohs"|Order=="Methylophilales"|Family=="Methylobacteriaceae"|Class=="[Methylacidiphilae]")
plot_bar(MethOTUGenotype, "Type", fill = "Order", facet_grid = ~Class)+labs(title = "Methanogen 16S Abundance Genotypic Ave")+ylab("Relative %")+xlab("Functional Group")+ theme(text=element_text(family = "Helvetica",size=14))+theme_bw()#+geom_bar(stat="identity")

#Metholotrpohs
#Methylophilales
#Methylobacteriaceae
#[Methylacidiphilae]


#Average by the GenotypxTime Factor
##Need to apply the correction and divide by the number of samples of the replicates. ``
OTU_merge = merge_samples(OTU16S, "GxTFactor", fun = mean) ##This mean function is only summing the values
GxTAve = transform_sample_counts(OTU_merge, function(x)  x/4 )  ##Need to include this in the data|| This is divdiing the data

sample_data(GxTAve)$GxTFactor <- levels(sample_data(OTU16S)$GxTFactor)##I have to look up this sample levels data 
##Im not sure if these levels are being combined correctly
sample_data(GxTAve)$Genotype <- levels(sample_data(OTU16S)$Genotype)##I have to look up this sample levels data 
sample_data(GxTAve)$Time <- levels(sample_data(OTU16S)$Time)##I have to look up this sample levels data 
sample_data(GxTAve)$Type <- levels(sample_data(OTU16S)$Type)##I have to look up this sample levels data 

GP.ord <- ordinate(GxTAve, "NMDS", "bray",)
p1 = plot_ordination(GxTAve, GP.ord, type="Samples", color="GxTFactor",  title="Bacterial 16S Ordination")+geom_point(size=3)+
  geom_text(mapping = aes(label = GxTFactor), size = 2, vjust = 1.5)+theme_bw()

print(p1)

##Manscript Figure 1B====
##This one may illistrate the point better
#I may need to remove the point labels
p1 = plot_ordination(GxTAve, GP.ord, type="Samples", color="Time",  title="16S rRNA Microbiome: E and P color same")+
  geom_point(size=3)+
  theme_bw()+
  stat_ellipse(aes(color=Time, group=Time),type = "norm")
  #geom_text(mapping = aes(label = GxTFactor), size = 2, vjust = 1.5)
print(p1)

pdf("Output/16SGxTAveSameNoLab.pdf",width=6,height=4,paper='special') 
print(p1)
dev.off()

##This changes the shapes fo the facros
p1 = plot_ordination(GxTAve, GP.ord, type="Samples", shape="Type",  title="16S rRNA Microbiome: Point color type:E color time")+
  geom_point(size=3,shape=1)+
  geom_text(mapping = aes(label = GxTFactor), size = 2, vjust = 1.5)+theme_bw()+
  stat_ellipse(aes(color=Time, group=Time),type = "norm")+scale_shape_identity()

#This has the eliples colored for time point and the Line type coloring the points
#Ploting this only works when one type of variable is mapped for some reason.
p1 = plot_ordination(GxTAve, GP.ord, type="Samples", color="Time",  title="16S rRNA Microbiome: Point color type:E color time")+
  geom_point(size=3,)+
  theme_bw()+
  stat_ellipse(aes(color=Time, group=Time),type = "norm")+scale_shape_identity()
  #geom_text(mapping = aes(label = GxTFactor), size = 2, vjust = 1.5) ##Coode to label points

plot(p1)
##In this figure 
##Need to remove outlier values in LH123

print(p1)
pdf("Output/16SGxTAve2.pdf",width=6,height=4,paper='special') 
print(p1)
dev.off()
#Average by Genotype Factor ====
##This is where I need to pick up =====
OTU_merge = merge_samples(OTU16S, "Genotype", fun = mean) #averages each OTU in all samples belonging to each habitat class
OTU_merge = transform_sample_counts(OTU_merge, function(x)  x/4 )  ##Need to include this in the data|| This is divdiing the data
sample_data(OTU_merge)$Genotype <- levels(sample_data(OTU16S)$Genotype)
#sample_data(OTU_merge)$Type <- levels(sample_data(OTU16S)$Type)

GP.ord <- ordinate(OTU_merge, "NMDS", "bray",)
p1 = plot_ordination(OTU_merge, GP.ord, type="Samples", color="Type",  title="Bacterial 16S Ordination")+
  geom_point(size=3)+geom_text(mapping = aes(label = Genotype), size = 2, vjust = 1.5)+theme_bw()+theme(legend.position = "none")
print(p1)

pdf("Output/16SGenotypeAve2.pdf",width=6,height=4,paper='special') 
print(p1)
dev.off()
###Manusript Figure 1A====
p1 = plot_ordination(OTU_merge, GP.ord, type="Samples", color="Type",  title="16S rRNA")+
  geom_point(size=3)+theme_bw()+
  stat_ellipse(aes(color=Type, group=Type),type = "norm")+geom_text(mapping = aes(label = Type), size = 2, vjust = 1.5)+
  scale_shape_identity()
print(p1)


p1 = plot_ordination(OTU_merge, GP.ord, type="Samples", color="Type",  title="16S+Type+GenoMean")+
  geom_point(size=3)+theme_bw()+
  theme(legend.position = "none")+stat_ellipse(aes(color=Type, group=Type),type = "norm")#+
  #geom_text(mapping = aes(label = Type), size = 2, vjust = 1.5)
print(p1)
pdf("Output/TypeGenotypeMean2.pdf",width=6,height=4,paper='special') 
print(p1)
dev.off()


GP.ord <- ordinate(OTU_merge, "NMDS", "bray",)
p1 = plot_ordination(OTU_merge, GP.ord, type="Samples", color="Year",  title="Bacterial 16S Ordination")+geom_point(size=3)+geom_text(mapping = aes(label = Genotype), size = 2, vjust = 1.5)
print(p1)


#I would like ot plot thje inbred and hybrid genotypic means with out teosinte ====
OTU16S = merge_phyloseq(biomOTU, Map)
OTU16S <- subset_samples(OTU16S, Time != "R")
OTU16S <- subset_samples(OTU16S, Type != "Teosinte")
OTUT1 <- subset_samples(OTU16S, Time == "1")
OTUT2 <- subset_samples(OTU16S, Time == "2")
OTUT3 <- subset_samples(OTU16S, Time == "3")


GP.ord <- ordinate(OTU16S, "NMDS", "bray",)
p1 = plot_ordination(OTU16S, GP.ord, type="Samples", color="Genotype",  title="Bacterial 16S Ordination")+geom_point(size=3)+geom_text(mapping = aes(label = Genotype), size = 2, vjust = 1.5)+
  stat_ellipse(aes(color=Genotype, group=Genotype),type = "t")
print(p1)



GP.ord <- ordinate(OTUT1, "NMDS", "bray",)
p1 = plot_ordination(OTUT1, GP.ord, type="Samples", color="Genotype",  title="Bacterial 16S Ordination")+geom_point(size=3)+geom_text(mapping = aes(label = Genotype), size = 2, vjust = 1.5)+
  geom_polygon(aes(alpha=0.9,fill="Genotype"))
print(p1)

GP.ord <- ordinate(OTUT2, "NMDS", "bray",)
p1 = plot_ordination(OTUT2, GP.ord, type="Samples", color="Genotype",  title="Bacterial 16S Ordination")+geom_point(size=3)+geom_text(mapping = aes(label = Genotype), size = 2, vjust = 1.5)+
  stat_ellipse(aes(color=Genotype, group=Genotype),type = "t")
print(p1)

GP.ord <- ordinate(OTUT3, "NMDS", "bray",)
p1 = plot_ordination(OTUT3, GP.ord, type="Samples", color="Genotype",  title="Bacterial 16S Ordination")+geom_point(size=3)+geom_text(mapping = aes(label = Genotype), size = 2, vjust = 1.5)+
  stat_ellipse(aes(color=Genotype, group=Genotype),type = "t")
print(p1)

#Finding the mean
OTU_merge = merge_samples(OTU16S, "Genotype", fun = mean) #averages each OTU in all samples belonging to each habitat class
sample_data(OTU_merge)$Genotype <- levels(sample_data(OTU16S)$Genotype)

GP.ord <- ordinate(OTU_merge, "NMDS", "bray",)
p1 = plot_ordination(OTU_merge, GP.ord, type="Samples", color="Genotype",  title="Bacterial 16S Ordination")+geom_point(size=3)+geom_text(mapping = aes(label = Genotype), size = 2, vjust = 1.5)
print(p1)

##Lets just look at the taxa that are dominating the micorbial community: This code isnt working and taking too ln
get_taxa_unique(OTU16S, "Phylum")

plot_tree(OTU16S, color = "Treatment", shape = "Season", label.tips = "Genus", 
          size = "abundance", plot.margin = 0.5, ladderize = TRUE)

p2=plot_bar(OTU16S, "Genotype", fill = "Family", facet_grid = ~Class)+
  theme(text=element_text(family = "Helvetica",size=14))
print(p2)


###Here we are sorting and filtering the OTUS====
GenotypeAve = merge_samples(OTU16S, "Type", fun = mean) #averages each OTU in all samples belonging to each habitat class
sample_data(GenotypeAve)$Type <- levels(sample_data(OTU16S)$Type)
#relative byt total bumber of genes
GenotypeAve_trans = transform_sample_counts(GenotypeAve, function(x) x/sum(x) )
##Above we relative the OTU table 

taxa_sums(GenotypeAve_trans) [1:500]
# taxa sums for top 10 taxa (as used below)
sort(taxa_sums(GenotypeAve_trans), TRUE) [1:500]

TopNOTUs = names(sort(taxa_sums(GenotypeAve_trans), TRUE) [1:500])
OTU10 = prune_taxa(TopNOTUs, GenotypeAve_trans) ##This code may need to be altered
#quartz()

##Top 100 OTUs in the Microbiome 
plot_bar(OTU10, "Type", fill = "Order", facet_grid = ~Phylum)+
  theme(legend.position="none")+
  geom_bar(stat="identity")+scale_y_continuous(labels = scales::percent)


OTU_rel_merge = transform_sample_counts(OTU_merge, function(x) 100 * x/sum(x) )
NitrifersOTU=subset_taxa(OTU_rel_merge, Phylum == "Nitrospirae"|Order=="Nitrosomonadales"|Order=="Nitrososphaerales")
plot_bar(NitrifersOTU, "GxTFactor", fill = "Order", facet_grid = ~Class)+labs(title = "nitrfier taxa T2 Ave-Trans")+xlab("Functional Group")+ theme(text=element_text(family = "Helvetica",size=14))+theme_bw()#+geom_bar(stat="identity")
#I'm going to try transforming, averagign, subsetting

p2=plot_bar(OTU_rel_merge, "Genotype", fill = "Class", facet_grid = ~Phylum)+
  theme(text=element_text(family = "Helvetica",size=14))
quartz()
print(p2)
###I should Run DSEQ2

####Nitrogen Cycling Genes======
### Nitrogen Fixation ====
nifH_biom_file = import_biom("Data/otu_table_nifH_UR.biom")
Map = import_qiime_sample_data("Data/2017MappingFile-R.txt")
OTUnifHAll = merge_phyloseq(nifH_biom_file, Map)

###Rarefaction of the functional genes and transformation=====
rarecurve(t(otu_table(OTUnifHAll)), step=5, cex=0.5)

#set.seed(4)
OTUnifH.rarefied= rarefy_even_depth(OTUnifHAll, rngseed=4, sample.size=500, replace=F)
#ps.rarefied = rarefy_even_depth(ps, rngseed=1, sample.size=0.9*min(sample_sums(ps)), replace=F)
rarecurve(t(otu_table(OTUnifH.rarefied)), step=5, cex=0.5)

#Transforming samples
#OTUnifH <- transform_sample_counts(OTUnifH, log)
OTUnifH<-OTUnifH.rarefied
#Outliers samples
#Set1 Outliers
OTUnifH <- subset_samples(OTUnifH, X.SampleID != "MM17.T2.078")
OTUnifH <- subset_samples(OTUnifH, X.SampleID != "MM17.Teo1.08")
OTUnifH <- subset_samples(OTUnifH, X.SampleID != "MM17.T3.041")
OTUnifH <- subset_samples(OTUnifH, X.SampleID != "MM17.Teo3.03")
#Set2 outliers
OTUnifH <- subset_samples(OTUnifH, X.SampleID != "MM17.Teo3.21")
OTUnifH <- subset_samples(OTUnifH, X.SampleID != "MM17.Teo2.04")
OTUnifH <- subset_samples(OTUnifH, X.SampleID != "MM17.T2.047")
#Set 3 Outliers
OTUnifH <- subset_samples(OTUnifH, X.SampleID != "MM17.T2.058")
OTUnifH <- subset_samples(OTUnifH, X.SampleID != "MM17.T1.003")
OTUnifH <- subset_samples(OTUnifH, X.SampleID != "MM17.T2.054")
OTUnifH <- subset_samples(OTUnifH, X.SampleID != "MM17.T2.060")

#Plotting the diversiry 
plot_richness(OTUnifH, x = "Genotype", measures = "Chao1") +
  theme(text=element_text(family = "Helvetica",size=14))+geom_boxplot()

plot_richness(OTUnifH, x = "Type", measures = "Chao1") +
  theme(text=element_text(family = "Helvetica",size=14))+geom_boxplot()

p1<-plot_richness(OTUnifH, x = "Time", measures = "Chao1") +
  theme(text=element_text(family = "Helvetica",size=14))+geom_boxplot()

pdf("Output/nifH-Time-Richness.pdf",width=6,height=4,paper='special') 
print(p1)
dev.off()

p1<-plot_richness(OTUnifH, x = "Genotype", measures = "Chao1") +
  theme(text=element_text(family = "Helvetica",size=14))+geom_boxplot()

pdf("Output/nifH-Richness.pdf",width=6,height=4,paper='special') 
print(p1)
dev.off()


Richness<-plot_richness(OTUnifH, x = "GxTFactor", measures = "Chao1") +
  theme(text=element_text(family = "Helvetica",size=14))+geom_boxplot()

print(Richness)
pdf("Output/nifH-Richness-GxT.pdf",width=6,height=5,paper='special') 
print(Richness)
dev.off()
##Seems like there is modest increases in diversity in Time | 
##Most genotypes seem to be showing the same pattern in the diversity of NifH 

#Lets look at the compositon of these genes
GP.ord <- ordinate(OTUnifH, "NMDS", "bray",)
p1 = plot_ordination(OTUnifH, GP.ord, type="Type", color="Time",  title="NifH Ordination")+geom_point(size=3)+geom_text(mapping = aes(label = X.SampleID), size = 2, vjust = 1.5)
print(p1)
#pdf("Output/nifHTypeRDA.pdf",width=6,height=4,paper='special') 
print(p1)
dev.off()


##
p1 = plot_ordination(OTUnifH, GP.ord, type="Samples", color="Time",  title="NifH Ordination")+geom_point(size=2)
print(p1)


inbrednifH <- subset_samples(OTUnifH, Type == "Inbred")

p1 = plot_ordination(inbrednifH, GP.ord, type="Samples", color="Time",  title="NifH Ordination")+geom_point(size=2)
print(p1)

##Stack Plot
OTUnifAve = merge_samples(OTUnifH, "Genotype", fun = mean) #averages each OTU in all samples belonging to each habitat class
sample_data(OTUnifAve)$Genotype <- levels(sample_data(OTUnifH)$Genotype)
OTUnifAve = transform_sample_counts(OTUnifAve, function(x) 100 * x/sum(x) )

title = " nifH taxa"
p = plot_bar(OTUnifAve, "Genotype", fill = "Rank2", "Abundance", 
             title=title)+ geom_bar(aes(fill=OTU), stat="identity", position="stack")+theme(legend.position = "none")
print(p)

OTUnifAve = merge_samples(OTUnifH, "Type", fun = mean) #averages each OTU in all samples belonging to each habitat class
sample_data(OTUnifAve)$Type <- levels(sample_data(OTUnifH)$Type)
OTUnifAve = transform_sample_counts(OTUnifAve, function(x) 100 * x/sum(x) )

title = " nifH taxa"
p = plot_bar(OTUnifAve, "Type", fill = "Rank2", "Abundance", 
             title=title)+ geom_bar(aes(fill=Rank3), stat="identity", position="stack")+theme(legend.position = "none")
print(p)

##Here I want to plot abudances 
p = plot_bar(OTUnifH, "Genotype", fill = "Rank2", "Abundance", 
             title=title)+ geom_bar(aes(fill=Rank3), stat="identity", position="stack")+scale_fill_manual(values=(rainbow(75)))
print(p)

####### Arch amoA ======
biom_file = import_biom("Data/otu_table_AamoA_UR.biom")
Map = import_qiime_sample_data("Data/2017MappingFile-R.txt")
OTUAamoAAll = merge_phyloseq(biom_file, Map)

#Rarefy the data
rarecurve(t(otu_table(OTUAamoAAll)), step=5, cex=0.5)
#set.seed(4)
##Rareifed 2500
OTUAamoAAll= rarefy_even_depth(OTUAamoAAll, rngseed=4, sample.size=2500, replace=F)

##Alpha Diversity of Archeal amoA ====
#measurment of alpha diversity
Richness<-plot_richness(OTUAamoAAll, x = "Time", measures = "Chao1") +
  theme(text=element_text(family = "Helvetica",size=14))+geom_boxplot()
print(Richness)
pdf("Output/ARCH-DiversityinTime.pdf",width=6,height=4,paper='special') 
print(Richness)
dev.off()

Richness<-plot_richness(OTUAamoAAll, x = "Genotype", measures = "Chao1") +
  theme(text=element_text(family = "Helvetica",size=14))+geom_boxplot()
print(Richness)
pdf("Output/ARCH-DiversityinGenotype.pdf",width=6,height=4,paper='special') 
print(Richness)
dev.off()

Richness<-plot_richness(OTUAamoAAll, x = "GxTFactor", measures = "Chao1") +
  theme(text=element_text(family = "Helvetica",size=14))+geom_boxplot()
print(Richness)
pdf("Output/ARCH-DiversityinGxT.pdf",width=6,height=4,paper='special') 
print(Richness)
dev.off()

Richness<-plot_richness(OTUAamoAAll, x = "Type", measures = "Chao1") +
  theme(text=element_text(family = "Helvetica",size=14))+geom_boxplot()
print(Richness)
pdf("Output/ARCH-DiversityinType.pdf",width=6,height=4,paper='special') 
print(Richness)
dev.off()
##Genotype Average
OTU_merge = merge_samples(OTUAamoAAll, "Genotype", fun = mean) #averages each OTU in all samples belonging to each habitat class
OTU_merge = transform_sample_counts(OTU_merge, function(x)  x/4 )  ##Need to include this in the data
sample_data(OTU_merge)$Genotype <- levels(sample_data(OTUAamoAAll)$Genotype)
#sample_data(OTU_merge)$Type <- levels(sample_data(OTU16S)$Type)

GP.ord <- ordinate(OTU_merge, "NMDS", "bray",)
p1 = plot_ordination(OTU_merge, GP.ord, type="Samples", color="Type",  title="Arch amoA Ordination")+
  geom_point(size=3)+geom_text(mapping = aes(label = Genotype), size = 2, vjust = 1.5)+theme_bw()+theme(legend.position = "none")
print(p1)


#Lets look at the compositon of these genes
#This seems like it will be fruitfull to do some stats on.
GP.ord <- ordinate(OTUAamoAAll, "NMDS", "bray",)
p1 = plot_ordination(OTUAamoAAll, GP.ord, type="Type", color="Time",  title="amoA Ordination")+geom_point(size=3)+geom_text(mapping = aes(label = X.SampleID), size = 2, vjust = 1.5)
print(p1)

#Rarefication removes the outliers 
pdf("Output/amoATimeNMDS.pdf",width=6,height=4,paper='special') 
print(p1)
dev.off()

GP.ord <- ordinate(OTUAamoAAll, "NMDS", "bray",)
p1 = plot_ordination(OTUAamoAAll, GP.ord, type="Type", color="Type",  title="amoA Ordination")+geom_point(size=3)+geom_text(mapping = aes(label = X.SampleID), size = 2, vjust = 1.5)
print(p1)

#Rarefication removes the outliers 
pdf("Output/amoATypeNMDS.pdf",width=6,height=4,paper='special') 
print(p1)
dev.off()


GP.ord <- ordinate(OTUAamoAAll, "NMDS", "bray",)
p1 = plot_ordination(OTUAamoAAll, GP.ord, type="Type", color="Genotype",  title="amoA Ordination")+geom_point(size=3)+geom_text(mapping = aes(label = Genotype), size = 2, vjust = 1.5)
print(p1)
#Rarefication removes the outliers 
pdf("Output/amoAGenotypeNMDS.pdf",width=6,height=4,paper='special') 
print(p1)
dev.off()
###Try correlating NMDS
##Corralted the functional genes to process. 
NMDS<-GP.ord[["points"]]

#Need to average the data points to create a stack plot that makes sense
amoAve = merge_samples(OTUAamoAAll, "Genotype", fun = mean) #averages each OTU in all samples belonging to each habitat class
sample_data(amoAve)$Genotype <- levels(sample_data(OTUAamoAAll)$Genotype)
amoAve = transform_sample_counts(amoAve, function(x) 100 * x/sum(x) )

title = " amoA taxa"
p = plot_bar(amoAve, "Genotype", fill = "Rank2", "Abundance", 
             title=title)+ geom_bar(aes(fill=OTU), stat="identity", position="stack")+theme(legend.position = "none")
print(p)

##Try by type
amoAve = merge_samples(OTUAamoAAll, "Type", fun = mean) #averages each OTU in all samples belonging to each habitat class
sample_data(amoAve)$Type <- levels(sample_data(OTUAamoAAll)$Type)
amoAve = transform_sample_counts(amoAve, function(x) 100 * x/sum(x) )

title = " amoA taxa"
p = plot_bar(amoAve, "Type", fill = "Rank7", "Abundance", 
             title=title)+ geom_bar(aes(fill=OTU), stat="identity", position="stack")+theme(legend.position="none")
print(p)

amoAve = merge_samples(OTUAamoAAll, "Genotype", fun = mean) #averages each OTU in all samples belonging to each habitat class
sample_data(amoAve)$Genotype <- levels(sample_data(OTUAamoAAll)$Genotype)
amoAve = transform_sample_counts(amoAve, function(x) 100 * x/sum(x) )

title = " amoA taxa"
p = plot_bar(amoAve, "Genotype", fill = "Rank7", "Abundance", 
             title=title)+ geom_bar(aes(fill=OTU), stat="identity", position="stack")+theme(legend.position="none")
print(p)
##IF we look within out we see lots pf ot

amoA.TvI <- subset_samples(OTUAamoAAll, Type != "Hybrid")

diagdds = phyloseq_to_deseq2(amoA.TvI, ~ Type)
diagdds = DESeq(diagdds, test="Wald", fitType="parametric")
##Wald Signifance test using in DESEQ2 Analysis
#We are 

res = results(diagdds, cooksCutoff = FALSE)
alpha = 0.5
sigtab = res[which(res$padj < alpha), ]
sigtab = cbind(as(sigtab, "data.frame"), as(tax_table(OTUAamoAAll)[rownames(sigtab), ], "matrix"))
head(sigtab)

dim(sigtab)

theme_set(theme_bw())
scale_fill_discrete <- function(palname = "Set1", ...) {
  scale_fill_brewer(palette = palname, ...)
}
# Phylum order
x = tapply(sigtab$log2FoldChange, sigtab$Phylum, function(x) max(x))
x = sort(x, TRUE)
sigtab$Phylum = factor(as.character(sigtab$Phylum), levels=names(x))
# Genus order
x = tapply(sigtab$log2FoldChange, sigtab$Genus, function(x) max(x))
x = sort(x, TRUE)
sigtab$Genus = factor(as.character(sigtab$Genus), levels=names(x))
#Class order
x = tapply(sigtab$log2FoldChange, sigtab$Class, function(x) max(x))
x = sort(x, TRUE)
sigtab$Class = factor(as.character(sigtab$Class), levels=names(x))
#Order order
x = tapply(sigtab$log2FoldChange, sigtab$Order, function(x) max(x))
x = sort(x, TRUE)
sigtab$Order = factor(as.character(sigtab$Order), levels=names(x))

x = tapply(sigtab$log2FoldChange, sigtab$Family, function(x) max(x))
x = sort(x, TRUE)
sigtab$Family = factor(as.character(sigtab$Family), levels=names(x))

##Plotting the 
ggplot(sigtab, aes(x=Phylum, y=log2FoldChange, color=Phylum)) + geom_point(size=6) + 
  theme(axis.text.x = element_text(angle = -90, hjust = 0, vjust=0.8))

##Here Im going to order the row names
sorted_sigtab <- sigtab2[ (order(row.names(sigtab))), ]

library(gtools)
sorted_sigtab <- sigtab[ (mixedorder(row.names(sigtab))), ]

rownames(sorted_sigtab)
##Array of signifcant rownames
tax_table(OTU16S)

###This ordering is important I need to consider it. 

##Transform it to control for sample number
Type.Sig.Ave = merge_samples(OTU16S, "Type", fun = mean) #averages each OTU in all samples belonging to each habitat class
sample_data(Type.Sig.Ave)$Type <- levels(sample_data(OTU16S)$Type)
##Above we relative the OTU table 
#relative by total number of genes
Type.Sig.Ave.Trans = transform_sample_counts(Type.Sig.Ave, function(x) x/sum(x) )


##This takes the subsetted list from the Inbred vs Teosinte
Type.Sig.Ave <- subset_taxa(Type.Sig.Ave.Trans, rownames(tax_table(Type.Sig.Ave.Trans)) %in% c(rownames(sorted_sigtab)))

scale_fill_discrete <- function(palname = "Set1", ...) {
  scale_fill_brewer(palette = palname, ...)
}

plot_bar(Type.Sig.Ave, "Type", fill = "Phylum", facet_grid = ~Kingdom)+
  theme(legend.position="none")+geom_bar(stat="identity")+scale_fill_manual(values = rainbow(50))

##It seems like I am not ploting the total 
##Here Im removeing all of the NA and Archea
Type.Ave_trans.2=subset_taxa(Type.Sig.Ave, Kingdom == "Bacteria")

###Figure 2 Stack plot =====
##A
plot_bar(Type.Ave_trans.2, "Type", fill = "Phylum", facet_grid = ~Kingdom)+geom_bar(stat="identity")+
  scale_fill_manual(values = rainbow(25))+
  scale_y_continuous(labels = scales::percent_format(scale=100, accuracy = 1))+
  labs(y= "Percent relative abundance", x = "Type")
#Normal Size PDF

##Lots of blast hits, that dont match to anything in the DB
##Bacterial amoA=====
x2 = import_biom("Data/otu_table_BamoA_UR.biom")
Map = import_qiime_sample_data("Data/2017MappingFile-R.txt")
OTUBamoAAll = merge_phyloseq(biom_file, Map)

#Rarefy the data
rarecurve(t(otu_table(OTUBamoAAll)), step=5, cex=0.5)
set.seed(4)
OTUBamoAAll= rarefy_even_depth(OTUAamoAAll, rngseed=4, sample.size=337, replace=F)


Richness<-plot_richness(OTUBamoAAll, x = "Genotype", measures = "Chao1") +
  theme(text=element_text(family = "Helvetica",size=14))+geom_boxplot()
print(Richness)
pdf("Output/BACT-DiversityinTime.pdf",width=6,height=4,paper='special') 
print(Richness)

amoBAve = merge_samples(OTUBamoAAll, "Genotype", fun = mean) #averages each OTU in all samples belonging to each habitat class
sample_data(amoBAve)$Genotype <- levels(sample_data(OTUBamoAAll)$Genotype)
amoBAve = transform_sample_counts(amoBAve, function(x) 100 * x/sum(x) )

title = " BamoA taxa"
p = plot_bar(amoBAve, "Genotype", fill = "Rank2", "Abundance", 
             title=title)+ geom_bar(aes(fill=Rank3), stat="identity", position="stack")
print(p)


amoBAve = merge_samples(OTUBamoAAll, "Type", fun = mean) #averages each OTU in all samples belonging to each habitat class
sample_data(amoBAve)$Type <- levels(sample_data(OTUBamoAAll)$Type)
amoBAve = transform_sample_counts(amoBAve, function(x) 100 * x/sum(x) )

title = " BamoA taxa"
p = plot_bar(amoBAve, "Type", fill = "Rank2", "Abundance", 
             title=title)+ geom_bar(aes(fill=OTU), stat="identity", position="stack")+theme(legend.position = "none")
print(p)

###NirS=====
nirs = import_biom("Data/otu_table_nirS_UR.biom")
Map = import_qiime_sample_data("Data/2017MappingFile-R.txt")
OTUBnirSAll = merge_phyloseq(nirs, Map)

#Rarefy the data
rarecurve(t(otu_table(OTUBnirSAll)), step=5, cex=0.5)
set.seed(4)
OTUBnirSAll= rarefy_even_depth(OTUBnirSAll, rngseed=4, sample.size=1167, replace=F)


Richness<-plot_richness(OTUBnirSAll, x = "Genotype", measures = "Chao1") +
  theme(text=element_text(family = "Helvetica",size=14))+geom_boxplot()
print(Richness)
#pdf("Output/BACT-DiversityinTime.pdf",width=6,height=4,paper='special') 
print(Richness)

nirSAve = merge_samples(OTUBnirSAll, "Genotype", fun = mean) #averages each OTU in all samples belonging to each habitat class
sample_data(nirSAve)$Genotype <- levels(sample_data(OTUBnirSAll)$Genotype)
nirSAve = transform_sample_counts(nirSAve, function(x) 100 * x/sum(x) )

title = " nirS taxa"
p = plot_bar(nirSAve, "Genotype", fill = "Rank2", "Relative %", 
             title=title)+ geom_bar(aes(fill=Rank3), stat="identity", position="stack")
print(p)

p = plot_bar(nirSAve, "Genotype", fill = "Rank2", "Abundance", 
             title=title)+ geom_bar(aes(fill=OTU), stat="identity", position="stack")+theme(legend.position = ("none"))
print(p)
nirSAve = merge_samples(OTUBnirSAll, "Type", fun = mean) #averages each OTU in all samples belonging to each habitat class
sample_data(nirSAve)$Type <- levels(sample_data(OTUBnirSAll)$Type)
nirSAve = transform_sample_counts(nirSAve, function(x) 100 * x/sum(x) )
p = plot_bar(nirSAve, "Type", fill = "Rank2", "Abundance", 
             title=title)+ geom_bar(aes(fill=OTU), stat="identity", position="stack")+theme(legend.position = ("none"))
print(p)

##NirK====
nirk = import_biom("Data/otu_table_nirK_UR.biom")
Map = import_qiime_sample_data("Data/2017MappingFile-R.txt")
OTUBnirKAll = merge_phyloseq(nirk, Map)

#Rarefy the data
rarecurve(t(otu_table(OTUBnirSAll)), step=5, cex=0.5)
set.seed(4)
OTUBnirSAll= rarefy_even_depth(OTUBnirKAll, rngseed=4, sample.size=1948, replace=F)


Richness<-plot_richness(OTUBnirKAll, x = "Genotype", measures = "Chao1") +
  theme(text=element_text(family = "Helvetica",size=14))+geom_boxplot()
print(Richness)
#pdf("Output/BACT-DiversityinTime.pdf",width=6,height=4,paper='special') 
print(Richness)

nirKAve = merge_samples(OTUBnirKAll, "Genotype", fun = mean) #averages each OTU in all samples belonging to each habitat class
sample_data(nirKAve)$Genotype <- levels(sample_data(OTUBnirKAll)$Genotype)
nirKAve = transform_sample_counts(nirKAve, function(x) 100 * x/sum(x) )

title = " nirS taxa"
p = plot_bar(nirKAve, "Genotype", fill = "Rank2", "Relative %", 
             title=title)+ geom_bar(aes(fill=Rank3), stat="identity", position="stack")
print(p)

p = plot_bar(nirSAve, "Genotype", fill = "Rank2", "Abundance", 
             title=title)+ geom_bar(aes(fill=OTU), stat="identity", position="stack")+theme(legend.position = ("none"))
print(p)

                                  
nirKAve = merge_samples(OTUBnirKAll, "Type", fun = mean) #averages each OTU in all samples belonging to each habitat class
sample_data(nirKAve)$Type <- levels(sample_data(OTUBnirKAll)$Type)
nirKAve = transform_sample_counts(nirKAve, function(x) 100 * x/sum(x) )
p = plot_bar(nirSAve, "Type", fill = "Rank2", "Abundance", 
             title=title)+ geom_bar(aes(fill=OTU), stat="identity", position="stack")+theme(legend.position = ("none"))
print(p)

##NosZ
nosZ = import_biom("Data/otu_table_nosZ_UR.biom")
Map = import_qiime_sample_data("Data/2017MappingFile-R.txt")
OTUBnosZKAll = merge_phyloseq(nosZ, Map)

#Rarefy the data
rarecurve(t(otu_table(OTUBnirSAll)), step=5, cex=0.5)
set.seed(4)
OTUBnirSAll= rarefy_even_depth(OTUBnirKAll, rngseed=4, sample.size=1948, replace=F)


Richness<-plot_richness(OTUBnirKAll, x = "Genotype", measures = "Chao1") +
  theme(text=element_text(family = "Helvetica",size=14))+geom_boxplot()
print(Richness)

####### ITS fungi  ======
biom_file = import_biom("Data/otu_table_ITS_UR.biom")
biom_file = import_biom("Data/ITS_table_rarefied_10000.biom")

Map = import_qiime_sample_data("Data/2017MappingFile-R.txt")
OTUITSAll = merge_phyloseq(biom_file, Map)
OTUAITSAll <- subset_samples(OTUITSAll, Time != "R")

#Rarefy the data
# rarecurve(t(otu_table(OTUITSAll)), step=5, cex=0.5)
# #set.seed(4)
# OTUAITSAll= rarefy_even_depth(OTUITSAll, rngseed=4, sample.size=9000, replace=F)
##Alpha Diversity of ITS ====
#measurment of alpha diversity
Richness<-plot_richness(OTUAITSAll, x = "Time", measures = "Chao1") +
  theme(text=element_text(family = "Helvetica",size=14))+geom_boxplot()+theme_classic()
print(Richness)

pdf("Output/ITS-DiversityinTime_rarefied_10000.pdf",width=6,height=4,paper='special') 
print(Richness)
dev.off()

Richness<-plot_richness(OTUAITSAll, x = "Genotype", measures = "Chao1") +
  theme(text=element_text(family = "Helvetica",size=14))+geom_boxplot()
print(Richness)
pdf("Output/ITS-DiversityinGenotype.pdf",width=6,height=4,paper='special') 
print(Richness)
dev.off()

Richness<-plot_richness(OTUAITSAll, x = "GxTFactor", measures = "Chao1") +
  theme(text=element_text(family = "Helvetica",size=14))+geom_boxplot()
print(Richness)
pdf("Output/ITS-DiversityinGxT.pdf",width=6,height=4,paper='special') 
print(Richness)
dev.off()


##Fungal Beta-Diversity ====
GP.ord <- ordinate(OTUAITSAll, "NMDS", "bray",)
p1 = plot_ordination(OTUAITSAll, GP.ord, type="Type", color="Type",  title="Fungi Ordination")+geom_point(size=3)+geom_text(mapping = aes(label = Genotype), size = 2, vjust = 1.5)
print(p1)

pdf("Output/ITSTypeNMDS.pdf",width=6,height=4,paper='special') 
print(p1)
dev.off()

GP.ord <- ordinate(OTUAITSAll, "NMDS", "bray",)
p1 = plot_ordination(OTUAITSAll, GP.ord, type="Type", color="Time",  title="Fungi Ordination")+geom_point(size=3)+geom_text(mapping = aes(label = Genotype), size = 2, vjust = 1.5)
print(p1)

pdf("Output/ITSTimeNMDS.pdf",width=6,height=4,paper='special') 
print(p1)
dev.off()

##Manuscript Figure 1C=====
OTU_merge = merge_samples(OTUAITSAll, "Genotype", fun = mean) #averages each OTU in all samples belonging to each habitat class
OTU_merge = transform_sample_counts(OTU_merge, function(x)  x/4 )  ##Need to include this in the data|| This is divdiing the data
sample_data(OTU_merge)$Genotype <- levels(sample_data(OTUAITSAll)$Genotype)

GP.ord <- ordinate(OTU_merge, "NMDS", "bray",)

p1 = plot_ordination(OTU_merge, GP.ord, type="Samples", color="Type",  title="ITS-Type")+
  geom_point(size=3)+theme_bw()+
  theme(legend.position = "none")+stat_ellipse(aes(color=Type, group=Type),type = "norm")+
  geom_text(mapping = aes(label = Genotype), size = 2, vjust = 1.5)
print(p1)

pdf("Output/ITSTypeGenotype2.pdf",width=6,height=4,paper='special') 
print(p1)
dev.off()
##Manuscript Figure 1D====
OTU_merge = merge_samples(OTUAITSAll, "GxTFactor", fun = mean) ##This mean function is only summing the values.
GxTAve = transform_sample_counts(OTU_merge, function(x)  x/4 )  ##Need to include this in the data|| This is divding the data


sample_data(GxTAve)$GxTFactor <- levels(sample_data(OTUAITSAll)$GxTFactor)
# sample_data(GxTAve)$Genotype <- levels(sample_data(OTU16S)$Genotype)##I have to look up this sample levels data 
#sample_data(GxTAve)$Time <- levels(sample_data(OTUAITSAll)$Time)##I have to look up this sample levels data 
# sample_data(GxTAve)$Type <- levels(sample_data(OTU16S)$Type)##I have to look up this sample levels data 

#Making ordintion
GP.ord <- ordinate(GxTAve, "NMDS", "bray",)

p1 = plot_ordination(GxTAve, GP.ord, type="Samples", color="GxTFactor",  title="ITS Microbiome: E and P color same")+
  geom_point(size=3)+
  theme_bw()#+
  stat_ellipse(aes(color=Time, group=Time),type = "norm")
#geom_text(mapping = aes(label = GxTFactor), size = 2, vjust = 1.5)
print(p1)

p1 = plot_ordination(GxTAve, GP.ord, type="Samples", color="Time",  title="ITS rRNA Microbiome: Point color type:E color time")+
  geom_point(size=3,)+
  theme_bw()+
  stat_ellipse(aes(color=Time, group=Time),type = "norm")+scale_shape_identity()
#geom_text(mapping = aes(label = GxTFactor), size = 2, vjust = 1.5) ##Coode to label points

plot(p1)

pdf("Output/ITS-GxTAveSameNoLab2.pdf",width=6,height=4,paper='special') 
print(p1)
dev.off()


##Here We can use the phylotseq plots
ITSAve = merge_samples(OTUAITSAll, "Type", fun = mean) #averages each OTU in all samples belonging to each habitat class
sample_data(ITSAve)$Type <- levels(sample_data(OTUAITSAll)$Type)
ITSAve = transform_sample_counts(ITSAve, function(x) 100 * x/sum(x) )

title = " Fungal ITS taxa"
p = plot_bar(ITSAve, "Type", fill = "Rank1", "Abundance", 
             title=title)+ geom_bar(aes(fill=Rank3), stat="identity", position="stack")+scale_fill_manual(values=(rainbow(75)))
print(p)
##I Need to look up


-
### Statisical Analysis: 2017 Maize Micorbiome BELOW ======
#These Outliers are MM17.T3.037 and MM17.T2.041 Outlier figure is present in Output folder

library(vegan)
library(reshape2)
library(ggplot2)
library(tidyverse)

#Load data
OTU16SAll <- read.csv("Data/16S-OTU-Meta.csv")

#Outlier Removal 
MC16S<-summary(OTU16SAll)

OTU16SAll= subset(OTU16SAll, X.SampleID != "MM17.T3.037")
OTU16SAll= subset(OTU16SAll, X.SampleID != "MM17.T2.041")
OTU16SAll= subset(OTU16SAll, X.SampleID != "MM17.T2.027")
OTU16SAll= subset(OTU16SAll, X.SampleID != "MM17.T2.040")
OTU16SAll= subset(OTU16SAll, X.SampleID != "MM17.T3.063")




##Hear we are finding the genotypic mean || Here Im calcluatig the mean of a genotype- across time
Run1.results= aggregate(OTU16S[, 9:37604], list(OTU16S$Genotype), mean)

Clust= aggregate(OTU16S[, 9:37604], list(OTU16S$Genotype), mean)

rownames(Clust) <- as.character(unlist(Clust[,1]))

#All.e.std = decostand(OTU16S, method = "total")

#This is a flitering model that removes parts of the microbaial commmunity
filter = function(x){
  temp = which(apply(x,2, function(x) sum(x>0) >= ceiling(0.99 * length(x))))
  x.filtered = x[, temp]
  return(x.filtered)
}

Clust2 = filter(Clust[, 9:37597])

#Here is the average microbiome community with a filter of 99% precentage of samples

#Here we are creating our treatments  
 MM.16S.spec <- OTU16S[,9:37604]
 MM.16S.treatment <- OTU16S[,1:8]
 
#Making a new combined factor thats called row ranf
MM.16S.t.Model<-cbind(MM.16S.treatment,RangeRow=paste(MM.16S.treatment$Range,MM.16S.treatment$Row))

MM.16S.t.Model<-as.data.frame(MM.16S.t.Model)

##CLUSTER SOFT
library(cluster)
MCDistance<-daisy(Clust, metric = c("gower"), stand = FALSE, type = list())

MCDistance<-daisy(Clust2, metric = c("gower"), stand = TRUE, type = list())
library(ape)

Bob<-nj(as.dist(MCDistance)) 
plot(hclust(as.dist(MCDistance)))

adonis(MM.16S.spec ~ Genotype*Time+Row+Range+RangeRow, MM.16S.t.Model)
Model16S<-adonis(MM.16S.spec ~ Genotype*Time+Row+Range+RangeRow, MM.16S.t.Model)

#overall 30% of the community is explained by plant + 7% by Time (seasonality, temp, etc) + 25-30% Space (current location and spatial gradient) of the variance 
#There does not seem to be a signifcant Space*Time interaction here, but for some reason it does explain lots of varaince in the model. 
# Df SumsOfSqs MeanSqs F.Model      R2 Pr(>F)    
#   Genotype       26    4.2653 0.16405  2.3051 0.14695  0.001 ***  
#   Time            2    2.0737 1.03683 14.5690 0.07144  0.001 ***  
#   RangeRow       74    7.8639 0.10627  1.4932 0.27093  0.001 ***  
#   Genotype:Time  52    4.8705 0.09366  1.3161 0.16780  0.001 *** 
###Time:RangeRow 114    8.5293 0.07482  1.0513 0.29385  0.222    
#   Residuals      20    1.4233 0.07117         0.04904           
# Total         288   29.0260                 1.00000           
# ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

# Df SumsOfSqs MeanSqs F.Model      R2 Pr(>F)    
# Genotype       26    4.2653 0.16405  2.2087 0.14695  0.001 ***
#   Time            2    2.0737 1.03683 13.9596 0.07144  0.001 ***
#   Row             1    0.8968 0.89683 12.0748 0.03090  0.001 ***
#   Range           1    0.3204 0.32043  4.3142 0.01104  0.001 ***
#   RangeRow       72    6.6467 0.09231  1.2429 0.22899  0.001 ***
#   Genotype:Time  52    4.8705 0.09366  1.2611 0.16780  0.001 ***
#   Residuals     134    9.9526 0.07427         0.34289           
# Total         288   29.0260                 1.00000           
# ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

adonis(MM.16S.spec ~ Type*Time+Row+Range+RangeRow, MM.16S.t.Model)
##Here 
adonis(MM.16S.spec ~ Genotype*Time+Row+Range+RangeRow, MM.16S.t.Model, strata=MM.16S.t.Model$RangeRow)

adonis(MM.16S.spec ~ Genotype*RangeRow*Time+Row+Range, MM.16S.t.Model)
# Df SumsOfSqs MeanSqs F.Model      R2 Pr(>F)    
# Genotype                26    4.2653 0.16405  2.1773 0.14695  0.001 ***
#   RangeRow                74    7.9016 0.10678  1.4172 0.27223  0.005 ** 
#   Time                     2    2.0360 1.01798 13.5112 0.07014  0.001 ***
#   Genotype:RangeRow        6    0.4879 0.08132  1.0793 0.01681  0.353    
# Genotype:Time           52    4.8655 0.09357  1.2419 0.16763  0.102    
# RangeRow:Time          114    8.5235 0.07477  0.9924 0.29365  0.552    
# Genotype:RangeRow:Time  11    0.7202 0.06547  0.8690 0.02481  0.817    
# Residuals                3    0.2260 0.07534         0.00779           
# Total                  288   29.0260                 1.00000           
# ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1


#Calculate the distance 
adonis(MM.16S.spec ~ Genotype+Time+Rep, MM.16S.treatment)

#Genotype is a more significant factor then time. 
#adonis(formula = MM.16S.spec ~ Genotype + Time + Rep + Type,      data = MM.16S.treatment)
# Df SumsOfSqs MeanSqs F.Model      R2 Pr(>F)    
# Genotype   26    4.4065 0.16948  1.8351 0.14150  0.001 ***
#   Time        2    2.0979 1.04896 11.3579 0.06737  0.001 ***
#   Rep         1    0.2551 0.25514  2.7626 0.00819  0.001 ***
#   Residuals 264   24.3817 0.09235         0.78294           
# Total     293   31.1412                 1.00000           
# ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

#I may want to alter the data via a tranfromation to see what it looks like 
adonis(MM.16S.spec ~ Type+Time+Rep, MM.16S.treatment)
#Type also expains a siginficant proportion of the MC
# Df SumsOfSqs MeanSqs F.Model      R2 Pr(>F)    
#   Type        2    1.8986 0.94931 10.1709 0.06097  0.001 ***
#   Time        2    2.1028 1.05141 11.2647 0.06753  0.001 ***
#   Rep         1    0.2590 0.25904  2.7754 0.00832  0.001 ***
#   Residuals 288   26.8808 0.09334         0.86319           
# Total     293   31.1412                 1.00000           
# ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

###Filtering function for MC =====
#Siginfcants is modestly improved because of the removal of singleton bacteria
All.e.std = decostand(MM.16S.spec, method = "total")
filter = function(x){
  temp = which(apply(x,2, function(x) sum(x>0) >= ceiling(0.99 * length(x))))
  x.filtered = x[, temp]
  return(x.filtered)
}
All.e.filtered = filter(All.e.std)

adonis(All.e.filtered ~ Genotype+Time+Rep, MM.16S.treatment)
#This is the ANOVA when only doing statistics on micorbes that are present in 99% of the samples. 
#Variablity of Taxa do not increase resolution. they infact decrease it. 
#This core of 99% taxa is about 742 OTUs - These are the persistant Core OTUs.

# Df SumsOfSqs MeanSqs F.Model      R2 Pr(>F)    
# Genotype   26    2.3050 0.08865  2.6406 0.18067  0.001 ***
#   Time        2    1.5238 0.76190 22.6939 0.11944  0.001 ***
#   Rep         1    0.1329 0.13285  3.9572 0.01041  0.003 ** 
#   Residuals 262    8.7960 0.03357         0.68947           
# Total     291   12.7577                 1.00000           
# ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

adonis(All.e.filtered ~ Genotype*Time+Row+Range+RangeRow, MM.16S.t.Model)

#So when we remove singleton bacteria from our MC we can see that there is an improvement in our genotype effects AND decreases in the expliantory effects of space
#   Df SumsOfSqs MeanSqs F.Model      R2 Pr(>F)    
#   Genotype       26    2.2793 0.08767  3.4091 0.18112  0.001 ***
#   Time            2    1.5160 0.75799 29.4766 0.12047  0.001 ***
#   Row             1    0.5130 0.51299 19.9490 0.04076  0.001 ***
#   Range           1    0.1694 0.16944  6.5891 0.01346  0.001 ***
#   RangeRow       72    2.6177 0.03636  1.4138 0.20801  0.001 ***
#   Genotype:Time  52    2.0431 0.03929  1.5279 0.16235  0.001 ***
#   Residuals     134    3.4458 0.02572         0.27382           
# Total         288   12.5843                 1.00000           
# ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

#When we strata the block varaibla eand add it ot the interaction it no lonfer becomes signficant.
adonis(All.e.filtered ~ Genotype*Time*RangeRow+Row+Range, MM.16S.t.Model, strata = MM.16S.t.Model$RangeRow)
# Df SumsOfSqs MeanSqs F.Model      R2 Pr(>F)    
# Genotype                26    2.2793 0.08767  3.2368 0.18112  0.001 ***
#   Time                     2    1.5160 0.75799 27.9870 0.12047  0.001 ***
#   RangeRow                74    3.3001 0.04460  1.6466 0.26224  0.030 *  
#   Genotype:Time           52    2.0431 0.03929  1.4507 0.16235  0.108    
# Genotype:RangeRow        6    0.1842 0.03070  1.1335 0.01464  0.373    
# Time:RangeRow          114    2.9441 0.02583  0.9535 0.23395  0.631    
# Genotype:Time:RangeRow  11    0.2363 0.02148  0.7931 0.01877  0.795    
# Residuals                3    0.0813 0.02708         0.00646           
# Total                  288   12.5843                 1.00000           
# ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

adonis(All.e.filtered ~ Genotype/Type*Time+RangeRow+Row+Range, MM.16S.t.Model)

adonis(All.e.filtered ~ Type/Genotype*Time+Type+RangeRow+Row+Range, MM.16S.t.Model)
# Df SumsOfSqs MeanSqs F.Model      R2 Pr(>F)    
#   Type                 2    1.2921 0.64603 25.1227 0.10267  0.001 ***
#   Time                 2    1.5295 0.76475 29.7394 0.12154  0.001 ***
#   RangeRow            93    4.0901 0.04398  1.7103 0.32502  0.001 ***
#   Type:Time            4    0.7477 0.18693  7.2694 0.05942  0.001 ***
#   Type:Genotype:Time  48    1.2953 0.02699  1.0494 0.10293  0.238    
# Residuals          134    3.4458 0.02572         0.27382           
# Total              288   12.5843                 1.00000           
# ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

adonis(All.e.filtered ~ Genotype*Time+Row+Range+RangeRow, MM.16S.t.Model, strata = MM.16S.t.Model$RangeRow)



###Looking at the MC within the INBRED lines ======
Inbred= subset(OTU16SAll, Type == "Inbred")
Inbred= subset(Inbred, Time != "R")


Inbred.treatment <- Inbred[,1:8]
Inbred.spec <- Inbred[,9:37604]

#Making blocking term
Inbred.treatment<-cbind(Inbred.treatment,RangeRow=paste(Inbred.treatment$Range,Inbred.treatment$Row))
Inbred.treatment<-as.data.frame(Inbred.treatment)

adonis(Inbred.spec ~ Genotype+Time+Rep, Inbred.treatment)
##This does not include the rhizosphere
##Rhizosphere increases the statsitcal signfi of a sample
# ##Root core zone doesnt seem to be the best 
# Df SumsOfSqs MeanSqs F.Model      R2 Pr(>F)    
# Genotype   11    0.9965 0.09059  1.0675 0.08693  0.128    
# Time        2    0.9851 0.49257  5.8046 0.08594  0.001 ***
#   Rep         1    0.2318 0.23179  2.7315 0.02022  0.001 ***
#   Residuals 109    9.2496 0.08486         0.80691           
# Total     123   11.4630                 1.00000           
# ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

adonis(Inbred.spec ~ Genotype*Time+Row+Range+RangeRow, Inbred.treatment)

# Df SumsOfSqs MeanSqs F.Model      R2 Pr(>F)    
# Genotype       11    0.9965 0.09059  1.1953 0.08693  0.002 ** 
#   Time            2    0.9851 0.49257  6.4994 0.08594  0.001 ***
#   Row             1    0.5360 0.53601  7.0727 0.04676  0.001 ***
#   Range           1    0.2239 0.22393  2.9548 0.01954  0.001 ***
#   RangeRow       34    3.1161 0.09165  1.2093 0.27184  0.001 ***
#   Genotype:Time  22    1.6645 0.07566  0.9983 0.14520  0.511    
# Residuals      52    3.9409 0.07579         0.34379           
# Total         123   11.4630                 1.00000           
# ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

#This is a pretty interesting model too
adonis(Inbred.spec ~ Genotype*Time+Row+Range+RangeRow, Inbred.treatment, strata=Inbred.treatment$RangeRow)


##Testing time as a removed factor
NoteosinteOTU16S= subset(OTU16S, Type != "Teosinte")
#Teosinte deff drives alot of these patterns

T1= subset(NoteosinteOTU16S, Time == "1")
T2= subset(NoteosinteOTU16S, Time == "2")
T3= subset(NoteosinteOTU16S, Time == "3")

T1.treatment <- T1[,1:8]
T1.spec <- T1[,9:37604]
adonis(T1.spec ~ Genotype+Range+Row, T1.treatment)
# Within timepoint 1 function: 
#This is a huge amount of variance explained
# Df SumsOfSqs  MeanSqs F.Model      R2 Pr(>F)    
#   Genotype  20    1.6153 0.080763  1.1119 0.35056  0.008 ** 
#   Range      1    0.1706 0.170640  2.3493 0.03703  0.001 ***
#   Row        1    0.2068 0.206844  2.8477 0.04489  0.001 ***
#   Residuals 36    2.6148 0.072634         0.56751           
# Total     58    4.6076                  1.00000           
# ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

T2.treatment <- T2[,1:8]
T2.spec <- T2[,9:37604]
adonis(T2.spec ~ Genotype+Range+Row, T2.treatment)
#
#
T3.treatment <- T3[,1:8]
T3.spec <- T3[,9:37604]
adonis(T3.spec ~ Genotype+Range+Row, T3.treatment)

All.e.std = decostand(T2.spec, method = "total")
filter = function(x){
  temp = which(apply(x,2, function(x) sum(x>0) >= ceiling(0.95 * length(x))))
  x.filtered = x[, temp]
  return(x.filtered)
}
All.T2.filtered = filter(All.e.std)
adonis(All.T2.filtered ~ Genotype+Range+Row, T2.treatment)

All.e.std = decostand(T3.spec, method = "total")
filter = function(x){
  temp = which(apply(x,2, function(x) sum(x>0) >= ceiling(0.95 * length(x))))
  x.filtered = x[, temp]
  return(x.filtered)
}
All.T3.filtered = filter(All.e.std)
adonis(All.T3.filtered ~ Genotype+Range+Row, T3.treatment)


##Here I will look at the genotypic means
Run1.results= aggregate(OTU16S[, 9:37604], list(OTU16S$Genotype), mean)

OTU16S.GenotypicMeans<-aggregate(OTU16SAll[, (9:37604)], list(Genotype=OTU16SAll$Genotype,
                                                                          Type=OTU16SAll$Type), mean)
#This is to pull out the data 
Genotype.16S.treatment <- OTU16S.GenotypicMeans[,1:2]
Genotype.16S.spec <- OTU16S.GenotypicMeans[,3:37598]

##Here we are taking the average to look at the effects of the microbiome
adonis(Genotype.16S.spec ~ Type, Genotype.16S.treatment)

# Df SumsOfSqs  MeanSqs F.Model      R2 Pr(>F)    
# Type       2   0.23898 0.119489  7.3614 0.38021  0.001 ***
#   Residuals 24   0.38956 0.016232         0.61979           
# Total     26   0.62854                  1.00000           
# ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1


###ITS Fungal Statsitical Analysis===== 
#Note: within the maize inbred there does not seem to be a difference in microbiome phenotype- but there seems to be a difference acros the hybrid and inbred 
library(vegan)

#Lets start the analysis without any rhizospheres 
OTUITSAll <- read.csv("Data/ITS-OTU-Meta.csv")
OTUITSAll= subset(OTUITSAll, Time != "R")

MM.ITS.spec <- OTUITSAll[,9:2382]
MM.ITS.treatment <- OTUITSAll[,1:8]


#Need to create the block factor
MM.ITS.t.Model<-cbind(MM.ITS.treatment,RangeRow=paste(MM.ITS.treatment$Range,MM.ITS.treatment$Row))
MM.ITS.t.Model<-as.data.frame(MM.ITS.t.Model)


#These are the old model. I want to use an updated more know which has a block location effect and the root and ranges 
adonis(MM.ITS.spec ~ Genotype+Range+Row, MM.ITS.treatment)
adonis(MM.ITS.spec ~ GxTFactor+Range+Row, MM.ITS.treatment)
# Df SumsOfSqs MeanSqs F.Model      R2 Pr(>F)    
# Genotype   26    12.550 0.48268  1.4402 0.11015  0.001 ***
#   Range       1     0.544 0.54443  1.6244 0.00478  0.004 ** 
#   Row         1     0.966 0.96643  2.8836 0.00848  0.001 ***
#   Residuals 298    99.873 0.33514         0.87659           
# Total     326   113.934                 1.00000           
# ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1


#New model| All of my model have an empanhasis on spatial effects within and indivdual plot.
##Whats interesing about these effects are that genotype has a very strong temporal effects with the micorbial community. 
adonis(MM.ITS.spec ~ Genotype*Time+Row+Range+RangeRow, MM.ITS.t.Model)
# Df SumsOfSqs MeanSqs F.Model      R2 Pr(>F)    
#   Genotype       26    12.550 0.48268  1.6315 0.11015  0.001 ***
#   Time            3     5.239 1.74647  5.9031 0.04599  0.001 ***
#   Row             1     0.974 0.97441  3.2935 0.00855  0.001 ***
#   Range           1     0.533 0.53323  1.8023 0.00468  0.003 ** 
#   RangeRow       72    24.074 0.33437  1.1302 0.21130  0.001 ***
#   Genotype:Time  61    22.634 0.37105  1.2542 0.19866  0.001 ***
#   Residuals     162    47.929 0.29586         0.42067           
# Total         326   113.934                 1.00000           
# ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
adonis(MM.ITS.spec ~ Genotype*Time+Row+Range+RangeRow, MM.ITS.t.Model, strata=MM.ITS.t.Model$RangeRow)

#Lets see if there is a Space*Time*Genetic interaction at this scale. 
#This line of thinking kinda makes me wonder if at which scale we will begin to see these kind of spatial genotypeic tradeoffs. 
adonis(MM.ITS.spec ~ Genotype*Time*RangeRow+Row+Range, MM.ITS.t.Model, strata = MM.ITS.t.Model$RangeRow)

#This model has a stronger signfincant space time interaction. This is not the case in our bacterial community. 
#Fungi seem to be doing interesting things in space time AND this is having some interactinve effects
#This model seems to explain large amounts of varance. 
# Df SumsOfSqs MeanSqs F.Model      R2 Pr(>F)    
#   Genotype                26    12.550 0.48268  2.2146 0.11015  0.001 ***
#   Time                     3     5.239 1.74647  8.0129 0.04599  0.001 ***
#   RangeRow                74    25.582 0.34570  1.5861 0.22453  0.001 ***
#   Genotype:Time           61    22.634 0.37105  1.7024 0.19866  0.001 ***
#   Genotype:RangeRow        6     1.442 0.24036  1.1028 0.01266  0.243    
#   Time:RangeRow          140    42.453 0.30323  1.3912 0.37261  0.003 ** 
#   Genotype:Time:RangeRow  12     3.162 0.26350  1.2090 0.02775  0.053 .  
# Residuals                4     0.872 0.21796         0.00765           
# Total                  326   113.934                 1.00000           
# ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1


#With Strata
# Df SumsOfSqs MeanSqs F.Model      R2 Pr(>F)    
# Genotype                26    12.550 0.48268  2.2146 0.11015  0.001 ***
#   Time                     3     5.239 1.74647  8.0129 0.04599  0.001 ***
#   RangeRow                74    25.582 0.34570  1.5861 0.22453  0.012 *  
#   Genotype:Time           61    22.634 0.37105  1.7024 0.19866  0.001 ***
#   Genotype:RangeRow        6     1.442 0.24036  1.1028 0.01266  0.208    
# Time:RangeRow          140    42.453 0.30323  1.3912 0.37261  0.013 *  
#   Genotype:Time:RangeRow  12     3.162 0.26350  1.2090 0.02775  0.060 .  
# Residuals                4     0.872 0.21796         0.00765           
# Total                  326   113.934                 1.00000           
# ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

#Type is inbred, teosinte and hybrids
adonis(MM.ITS.spec ~ Type+Time+Rep, MM.ITS.treatment)
# Df SumsOfSqs MeanSqs F.Model      R2 Pr(>F)    
# Type        2     4.089 2.04465  6.2889 0.03589  0.001 ***
#   Time        3     5.314 1.77149  5.4487 0.04665  0.001 ***
#   Rep         1     0.491 0.49129  1.5111 0.00431  0.020 *  
#   Residuals 320   104.039 0.32512         0.91315           
# Total     326   113.934                 1.00000           
# ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

adonis(MM.ITS.spec ~ Type*Time+Row+Range+RangeRow, MM.ITS.t.Model)

adonis(MM.ITS.spec ~Type/Genotype*Time+Row+Range+RangeRow, MM.ITS.t.Model)


#Is there an interaction effect be
OTUITSAll= subset(OTUITSAll, Time != "R")
MM.ITS.spec <- OTUITSAll[,9:2382]
MM.ITS.treatment <- OTUITSAll[,1:8]

adonis(MM.ITS.spec ~ Genotype*Time+Type+Rep, MM.ITS.treatment)

# Df SumsOfSqs MeanSqs F.Model      R2 Pr(>F)    
# Genotype       26    12.028 0.46260  1.4445 0.11699  0.001 ***
#   Time            2     2.246 1.12316  3.5071 0.02185  0.001 ***
#   Rep             1     0.487 0.48676  1.5199 0.00473  0.014 *  
#   Genotype:Time  52    19.510 0.37520  1.1716 0.18978  0.001 ***
#   Residuals     214    68.534 0.32025         0.66664           
# Total         295   102.805                 1.00000           
# ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

##This comparison only looks at the inred genotypes
Inbred= subset(OTUITSAll, Type == "Inbred")
Inbred= subset(Inbred, Time != "R")

Inbred.treatment <- Inbred[,1:8]
Inbred.spec <- Inbred[,9:2382]

adonis(Inbred.spec ~ Genotype+Time+Rep, Inbred.treatment)

# Df SumsOfSqs MeanSqs F.Model      R2 Pr(>F)    
# Genotype   11     3.850 0.35004  1.0275 0.08637  0.293    
# Time        2     1.150 0.57484  1.6875 0.02579  0.001 ***
#   Rep         1     0.403 0.40316  1.1835 0.00904  0.129    
# Residuals 115    39.175 0.34065         0.87879           
# Total     129    44.578                 1.00000           
# ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

##Testing the space time model in the inbred exclusvisly 
Maize= subset(OTUITSAll, Type != "Teosinte")
Maize= subset(Maize, Time != "R")

Maize.treatment <- Maize[,1:8]
Maize.spec <- Maize[,9:2382]

#Making the new model factor
MM.ITS.Inbred.Model<-cbind(Maize.treatment,RangeRow=paste(Maize.treatment$Range,Maize.treatment$Row))
MM.ITS.Inbred.Model<-as.data.frame(MM.ITS.Inbred.Model)

#New model with time and space! Wowo
adonis(Maize.spec ~ Genotype*Time+Row+Range+RangeRow, MM.ITS.Inbred.Model)
#Without the strata model We see no signficnatnce in the genotype across g
# Df SumsOfSqs MeanSqs F.Model      R2 Pr(>F)    
# Genotype       20     6.959 0.34796 1.05025 0.09124  0.120    
# Time            2     1.488 0.74393 2.24541 0.01951  0.001 ***
#   Row             1     0.989 0.98852 2.98364 0.01296  0.001 ***
#   Range           1     0.524 0.52355 1.58024 0.00686  0.003 ** 
#   RangeRow       61    21.235 0.34812 1.05074 0.27841  0.041 *  
#   Genotype:Time  40    12.612 0.31529 0.95165 0.16535  0.944    
# Residuals      98    32.469 0.33131         0.42568           
# Total         223    76.275                 1.00000           

#Here we use a stata model
adonis(Maize.spec ~ Genotype*Time+Row+Range+RangeRow, MM.ITS.Inbred.Model, strata = MM.ITS.Inbred.Model$RangeRow)

# Df SumsOfSqs MeanSqs F.Model      R2 Pr(>F)    
# Genotype       20     6.959 0.34796 1.05025 0.09124  0.414    
# Time            2     1.488 0.74393 2.24541 0.01951  0.001 ***
#   Row             1     0.989 0.98852 2.98364 0.01296  0.336    
# Range           1     0.524 0.52355 1.58024 0.00686  0.374    
# RangeRow       61    21.235 0.34812 1.05074 0.27841  0.547    
# Genotype:Time  40    12.612 0.31529 0.95165 0.16535  0.958    
# Residuals      98    32.469 0.33131         0.42568           
# Total         223    76.275                 1.00000           

#Here we will test these effects on the interactiv effects 
adonis(Maize.spec ~ Genotype*Time*RangeRow+Row+Range, MM.ITS.Inbred.Model, strata = MM.ITS.Inbred.Model$RangeRow)
adonis(Maize.spec ~ Type*Time+Range+RangeRow, MM.ITS.Inbred.Model, strata = MM.ITS.Inbred.Model$RangeRow)

# Df SumsOfSqs MeanSqs F.Model      R2 Pr(>F)    
# Type        1     0.346 0.34553 1.05595 0.00453  0.001 ***
#   Time        2     1.494 0.74719 2.28341 0.01959  0.001 ***
#   Range       1     0.520 0.52015 1.58957 0.00682  0.001 ***
#   RangeRow   81    28.835 0.35598 1.08789 0.37803  0.001 ***
#   Type:Time   2     0.578 0.28899 0.88316 0.00758  0.837    
# Residuals 136    44.502 0.32722         0.58345           
# Total     223    76.275                 1.00000           
# ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

adonis(Maize.spec ~ Genotype+Time+Rep, Maize.treatment)
# Df SumsOfSqs MeanSqs F.Model      R2 Pr(>F)    
# Genotype   20     7.124 0.35619  1.0759 0.08603  0.040 *  
#   Time        3     3.680 1.22658  3.7049 0.04444  0.001 ***
#   Rep         1     0.492 0.49243  1.4874 0.00595  0.007 ** 
#   Residuals 216    71.511 0.33107         0.86359           
# Total     240    82.807                 1.00000           
# ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

adonis(Maize.spec ~ Type+Time+Rep, Maize.treatment)


# Df SumsOfSqs MeanSqs F.Model      R2 Pr(>F)    
# Type        1     0.501 0.50130  1.5097 0.00605  0.012 *  
#   Time        3     3.783 1.26116  3.7981 0.04569  0.001 ***
#   Rep         1     0.491 0.49073  1.4779 0.00593  0.008 ** 
#   Residuals 235    78.031 0.33205         0.94233           
# Total     240    82.807                 1.00000           
# ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

adonis(Maize.spec ~ Type*Time+Rep, Maize.treatment)

Hybrid= subset(OTUITSAll, Type == "Hybrid")

H.treatment <- Hybrid[,1:8]
H.spec <- Hybrid[,9:2382]

adonis(H.spec ~ Genotype*Time+Rep, H.treatment)



 
 ####
##########Functional Genes=========
 #Here is where the analysis for functional genes begin, apprently I have processed out OTU tables
 library(vegan)
 
#nifH ======
 #I wonder if I can make all of this part of a function to automate this code. 
 #Basically the code would create the files and run the model
 
 nifHALLTreatment <- read.csv("Data/nifH-OTU-Meta.csv")
 #nifHALLTreatment <- read.csv("Data/otu_table_nifH.csv")
 
 
 
 #Ive coded the last row as the taxa, I need to remove this before running the model
 #Removal of factors not needed for the model 
 
 #Remove Taxa line||| I remade the dta frame
 # nifHTaxa<-nifHALLTreatment[324,]
 # nifHTaxa<-t(nifHTaxa)
 # 
 # #nifHALL<-na.omit(nifHALLTreatment)
 # nifHALL<-nifHALLTreatment[-324, ]

 #Removing Rhizosphere and outliers. 
 nifHAll <- subset(nifHALLTreatment, Time != "R")
 #nifHAll <- subset(nifHALLTreatment, Type != "Teosinte")
 
 #droplevels.data.frame(nifHAll)
 #Here we are creating our treatments  
 MM.nifH.spec <- nifHAll[,9:890]
 MM.nifH.treatment <- nifHAll[,1:8]
 
 #This converts the info inot a data frame
 # MM.nifH.spec[] <- lapply(MM.nifH.spec, as.numeric)
 
 #NifH seems to not be very normally distbuted, to correct for this I log transfromed the data
 #MM.nifH.spec<-log(MM.nifH.spec)

 #Making a new combined factor thats called row range
 MM.nifH.t.Model<-cbind(MM.nifH.treatment,RangeRow=paste(MM.nifH.treatment$Range,MM.nifH.treatment$Row))
 MM.nifH.t.Model<-as.data.frame(MM.nifH.t.Model)


 #From these models it seems if if nifH is primaraly sensative to spatial and temporal effects.
 #But not even that sensative, just kindaof
 #The interactive model makes the most sense to me
 #The text is for models with teosinte in them. 
 adonis(MM.nifH.spec ~ Genotype*Time+Row+Range+RangeRow, MM.nifH.t.Model)
# TIme and space seem to be weakly important caluebs
#  Df SumsOfSqs MeanSqs F.Model      R2 Pr(>F)  
#  Genotype       26    11.455 0.44060 0.99163 0.08727  0.569  
#  Time            2     1.174 0.58681 1.32071 0.00894  0.030 *
#    Row             1     0.686 0.68574 1.54337 0.00522  0.041 *
#    Range           1     0.676 0.67585 1.52112 0.00515  0.022 *
#    RangeRow       72    33.330 0.46291 1.04186 0.25391  0.110  
#  Genotype:Time  52    22.631 0.43521 0.97951 0.17240  0.729  
#  Residuals     138    61.315 0.44431         0.46711         
#  Total         292   131.267                 1.00000         
#  ---
#    Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
 adonis(MM.nifH.spec ~ Genotype*Time+Row+Range+RangeRow, MM.nifH.t.Model, strata = MM.nifH.t.Model$RangeRow)
 # Df SumsOfSqs MeanSqs F.Model      R2 Pr(>F)  
 # Genotype       26    11.455 0.44060 0.99163 0.08727  0.594  
 # Time            2     1.174 0.58681 1.32071 0.00894  0.035 *
 #   Row             1     0.686 0.68574 1.54337 0.00522  0.770  
 # Range           1     0.676 0.67585 1.52112 0.00515  0.045 *
 #   RangeRow       72    33.330 0.46291 1.04186 0.25391  0.448  
 # Genotype:Time  52    22.631 0.43521 0.97951 0.17240  0.612  
 # Residuals     138    61.315 0.44431         0.46711         
 # Total         292   131.267                 1.00000         
 # ---
 #   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
 adonis(MM.nifH.spec ~ Genotype*Time*RangeRow+Row+Range, MM.nifH.t.Model,strata=MM.nifH.t.Model$RangeRow)
 
 # Df SumsOfSqs MeanSqs F.Model      R2 Pr(>F)
 # Genotype                26    11.455 0.44060 0.88119 0.08727  0.858
 # Time                     2     1.174 0.58681 1.17362 0.00894  0.198
 # RangeRow                74    34.691 0.46880 0.93760 0.26428  0.578
 # Genotype:Time           52    22.631 0.43521 0.87042 0.17240  0.945
 # Genotype:RangeRow        6     2.795 0.46589 0.93178 0.02130  0.633
 # Time:RangeRow          117    51.699 0.44187 0.88375 0.39385  0.950
 # Genotype:Time:RangeRow  12     5.321 0.44341 0.88683 0.04054  0.811
 # Residuals                3     1.500 0.50000         0.01143       
 # Total                  292   131.267                 1.00000       
 
 ##Here I'm focusing on Type 
 adonis(MM.nifH.spec ~ Type*Time+Row+Range+RangeRow, MM.nifH.t.Model)
 # Df SumsOfSqs MeanSqs F.Model      R2 Pr(>F)   
 # Type        2     1.049 0.52467 1.07541 0.00731  0.130   
 # Time        2     1.176 0.58813 1.20546 0.00819  0.007 **
 #   Row         1     0.617 0.61742 1.26550 0.00430  0.012 * 
 #   Range       1     0.635 0.63503 1.30160 0.00442  0.010 **
 #   RangeRow   91    44.958 0.49405 1.01263 0.31318  0.150   
 # Type:Time   4     1.932 0.48298 0.98995 0.01346  0.545   
 # Residuals 191    93.186 0.48788         0.64913          
 # Total     292   143.554                 1.00000          
 
 adonis(MM.nifH.spec ~ Type/Genotype*Time+Row+Range+RangeRow, MM.nifH.t.Model)
 
 
 ####amoA Bacteria======
 amoBHALLTreatment <- read.csv("Data/BacamoA-OTU-Meta.csv")
 #Bacterial nitrficaion population seem to have signficant relationship with genotype and time
 #But no indivdiaul effects of time. 
 #That Taxa on the last colum seem to be an annomially in the first data set.
 #Nitrfication seems to be functionally most repsonvsive 
 #Remove Factors
 amoBALL <- subset(amoBHALLTreatment, Time != "R")
 #amoBALL <- subset(amoBHALLTreatment, Type != "Teosinte")
 
 MM.amoB.spec <- amoBALL[,9:106]
 MM.amoB.treatment <- amoBALL[,1:8]
 
 MM.amoB.t.Model<-cbind(MM.amoB.treatment,RangeRow=paste(MM.amoB.treatment$Range,MM.amoB.treatment$Row))
 MM.amoB.t.Model<-as.data.frame(MM.amoB.t.Model)
 
 #Remoing teosinte from model increases the model stats-- but weaken the genotype time interaction
 adonis(MM.amoB.spec ~ Genotype*Time+Row+Range+RangeRow, MM.amoB.t.Model)
 # Genotype       26     9.403 0.36164  1.0811 0.08857  0.165    
 # Time            2     1.730 0.86514  2.5862 0.01630  0.001 ***
 #   Row             1     0.533 0.53267  1.5923 0.00502  0.084 .  
 # Range           1     0.396 0.39591  1.1835 0.00373  0.275    
 # RangeRow       72    28.340 0.39362  1.1767 0.26694  0.002 ** 
 #   Genotype:Time  52    19.266 0.37050  1.1076 0.18147  0.050 *  
 #   Residuals     139    46.499 0.33452         0.43798           
 # Total         293   106.166                 1.00000           
 # ---
 #   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
 adonis(MM.amoB.spec ~ Genotype*Time+Row+Range+RangeRow, MM.amoB.t.Model, strata = MM.amoB.t.Model$RangeRow)
 # Df SumsOfSqs MeanSqs F.Model      R2 Pr(>F)    
 # Genotype       26     9.403 0.36164  1.0811 0.08857  0.187    
 # Time            2     1.730 0.86514  2.5862 0.01630  0.001 ***
 #   Row             1     0.533 0.53267  1.5923 0.00502  0.077 .  
 # Range           1     0.396 0.39591  1.1835 0.00373  0.280    
 # RangeRow       72    28.340 0.39362  1.1767 0.26694  0.002 ** 
 #   Genotype:Time  52    19.266 0.37050  1.1076 0.18147  0.069 .  
 # Residuals     139    46.499 0.33452         0.43798           
 # Total         293   106.166                 1.00000           
 adonis(MM.amoB.spec ~ Genotype*Time*RangeRow+Row+Range, MM.amoB.t.Model,strata = MM.amoB.t.Model$RangeRow)

 adonis(MM.amoB.spec ~ Type*Time+Row+Range+RangeRow, MM.amoB.t.Model)
 adonis(MM.amoB.spec ~ Type/Genotype*Time+Row+Range+RangeRow, MM.amoB.t.Model)
 
 #Set working directoy#
 setwd("~/Documents/Projects/Projects/Maize Microbiome/MM2017/Data/RDATA")
 
 ###Loading Important Packages###
 library(ggplot2)
 library(vegan)
 library(tidyverse)
 library(jtools)
 
 
 ###Loading Important Data###
 #Function.2017=read.csv(file="MM2017_Function.csv",head=T)
 Function.2017=read.csv(file="MM2017_Function_Modifed.6.20.20.csv",head=T)
 
 #### amoA Archeal =====
 #amoA group seems to be far more sensative to plant treatment then the bacterial group
 amoAALLTreatment <- read.csv("Data/ArchamoA-OTU-Meta.csv")
 
 amoAALL <- subset(amoAALLTreatment, Time != "R")
 #amoAALL <- subset(amoAALLTreatment, Type != "Teosinte")
 
 MM.amoA.spec <- amoAALL[,9:218]
 MM.amoA.treatment <- amoAALL[,1:8]
 
 MM.amoA.t.Model<-cbind(MM.amoA.treatment,RangeRow=paste(MM.amoA.treatment$Range,MM.amoA.treatment$Row))
 MM.amoA.t.Model<-as.data.frame(MM.amoA.t.Model)
 
 
 adonis(MM.amoA.spec ~ Genotype*Time+Row+Range+RangeRow, MM.amoA.t.Model)
 # Df SumsOfSqs MeanSqs F.Model      R2 Pr(>F)    
 # Genotype       26     3.346  0.1287   1.374 0.08401  0.003 ** 
 #   Time            2     6.798  3.3989  36.285 0.17068  0.001 ***
 #   Row             1     2.708  2.7076  28.905 0.06799  0.001 ***
 #   Range           1     0.406  0.4062   4.337 0.01020  0.002 ** 
 #   RangeRow       72     7.812  0.1085   1.158 0.19614  0.039 *  
 #   Genotype:Time  52     5.550  0.1067   1.139 0.13935  0.072 .  
 # Residuals     141    13.208  0.0937         0.33163           
 # Total         295    39.827                 1.00000           
 # ---
 #   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
 adonis(MM.amoA.spec ~ Genotype*Time+Row+Range+RangeRow, MM.amoA.t.Model, strata = MM.amoA.t.Model$RangeRow)
 # Df SumsOfSqs MeanSqs F.Model      R2 Pr(>F)    
 # Genotype       26     3.346  0.1287   1.374 0.08401  0.001 ***
 #   Time            2     6.798  3.3989  36.285 0.17068  0.001 ***
 #   Row             1     2.708  2.7076  28.905 0.06799  0.001 ***
 #   Range           1     0.406  0.4062   4.337 0.01020  0.001 ***
 #   RangeRow       72     7.812  0.1085   1.158 0.19614  0.001 ***
 #   Genotype:Time  52     5.550  0.1067   1.139 0.13935  0.198    
 # Residuals     141    13.208  0.0937         0.33163           
 # Total         295    39.827                 1.00000           

 adonis(MM.amoA.spec ~ Genotype*Time*RangeRow+Row+Range, MM.amoA.t.Model, strata = MM.amoA.t.Model$RangeRow)
 #In this model only time matters, everyhign else comes back insigndicant. 
 
 adonis(MM.amoA.spec ~ Type*Time+Row+Range+RangeRow, MM.amoA.t.Model)
 
 adonis(MM.amoA.spec ~ Type/Genotype*Time+Row+Range+RangeRow, MM.amoA.t.Model)
 
 
 #### Denitrification Statisitcs ======
 ##nirK gene=====
 nirKALLTreatment <- read.csv("Data/nirK-OTU-Meta.csv")

 nirKALL <- subset(nirKALLTreatment, Time != "R")
 #nirKALL <- subset(nirKALLTreatment, Type != "Teosinte")
 
 MM.nirK.spec <- nirKALL[,9:21030]
 MM.nirK.treatment <- nirKALL[,1:8]
 
 MM.nirK.t.Model<-cbind(MM.nirK.treatment,RangeRow=paste(MM.nirK.treatment$Range,MM.nirK.treatment$Row))
 MM.nirK.t.Model<-as.data.frame(MM.nirK.t.Model)
 
 #Maybe this will imporve from standarization
 #MM.nirK.spec<-log1p(MM.nirK.spec)
 #Statisitcal model| Genotype effects, but no genotype time interactions| Logging changes the calues a bit
 adonis(MM.nirK.spec ~ Genotype*Time+Row+Range+RangeRow, MM.nirK.t.Model)
 # Df SumsOfSqs MeanSqs F.Model      R2 Pr(>F)    
 # Genotype       26    10.984 0.42248  1.0986 0.09233  0.005 ** 
 #   Time            2     3.142 1.57106  4.0854 0.02641  0.001 ***
 #   Row             1     0.556 0.55561  1.4448 0.00467  0.010 ** 
 #   Range           1     0.483 0.48293  1.2558 0.00406  0.070 .  
 # RangeRow       72    29.081 0.40390  1.0503 0.24444  0.016 *  
 #   Genotype:Time  52    20.500 0.39424  1.0252 0.17232  0.180    
 # Residuals     141    54.222 0.38455         0.45577           
 # Total         295   118.968                 1.00000         
 
 adonis(MM.nirK.spec ~ Genotype*Time+Row+Range+RangeRow, MM.nirK.t.Model, strata = MM.nirK.t.Model$RangeRow)
 # Df SumsOfSqs MeanSqs F.Model      R2 Pr(>F)    
 # Genotype       26    10.984 0.42248  1.0986 0.09233  0.002 ** 
 #   Time            2     3.142 1.57106  4.0854 0.02641  0.001 ***
 #   Row             1     0.556 0.55561  1.4448 0.00467  0.752    
 # Range           1     0.483 0.48293  1.2558 0.00406  0.048 *  
 #   RangeRow       72    29.081 0.40390  1.0503 0.24444  0.001 ***
 #   Genotype:Time  52    20.500 0.39424  1.0252 0.17232  0.126    
 # Residuals     141    54.222 0.38455         0.45577           
 # Total         295   118.968                 1.00000           
 # ---
 #   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
 adonis(MM.nirK.spec ~ Genotype*Time*RangeRow+Row+Range, MM.nirK.t.Model, strata = MM.nirK.t.Model$RangeRow)
 
 adonis(MM.nirK.spec ~ Type*Time+Row+Range+RangeRow, MM.nirK.t.Model)
 
 adonis(MM.nirK.spec ~ Type/Genotype*Time+Row+Range+RangeRow, MM.nirK.t.Model)
 
 
 #nirS Genes=====
 
 nirSALLTreatment <- read.csv("Data/nirS-OTU-Meta.csv")
 
 nirSALL <- subset(nirSALLTreatment, Time != "R")
 #nirSALL <- subset(nirSALLTreatment, Type != "Teosinte")
 
 MM.nirS.spec <- nirSALL[,9:2615]
 MM.nirS.treatment <- nirSALL[,1:8]
 
 MM.nirS.t.Model<-cbind(MM.nirS.treatment,RangeRow=paste(MM.nirS.treatment$Range,MM.nirS.treatment$Row))
 MM.nirS.t.Model<-as.data.frame(MM.nirS.t.Model)
 
 #Running models and Results 
 adonis(MM.nirS.spec ~ Genotype*Time+Row+Range+RangeRow, MM.nirS.t.Model)
 # Df SumsOfSqs MeanSqs F.Model      R2 Pr(>F)    
 # Genotype       26    10.581 0.40697 1.08047 0.09333  0.031 *  
 #   Time            2     1.409 0.70462 1.87070 0.01243  0.001 ***
 #   Row             1     1.015 1.01503 2.69481 0.00895  0.001 ***
 #   Range           1     0.798 0.79755 2.11741 0.00703  0.001 ***
 #   RangeRow       72    28.471 0.39543 1.04984 0.25112  0.059 .  
 # Genotype:Time  52    17.994 0.34603 0.91868 0.15871  0.996    
 # Residuals     141    53.109 0.37666         0.46843           
 # Total         295   113.377                 1.00000           
 # ---
 #   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
 
 adonis(MM.nirS.spec ~ Genotype*Time+Row+Range+RangeRow, MM.nirS.t.Model, strata = MM.nirS.t.Model$RangeRow)
 # Df SumsOfSqs MeanSqs F.Model      R2 Pr(>F)    
 # Genotype       26    10.581 0.40697 1.08047 0.09333  0.322    
 # Time            2     1.409 0.70462 1.87070 0.01243  0.001 ***
 #   Row             1     1.015 1.01503 2.69481 0.00895  0.942    
 # Range           1     0.798 0.79755 2.11741 0.00703  0.249    
 # RangeRow       72    28.471 0.39543 1.04984 0.25112  0.852    
 # Genotype:Time  52    17.994 0.34603 0.91868 0.15871  0.996    
 # Residuals     141    53.109 0.37666         0.46843           
 # Total         295   113.377                 1.00000           
 # ---
 #   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
 
 adonis(MM.nirS.spec ~ Genotype*Time*RangeRow+Row+Range, MM.nirS.t.Model, strata = MM.nirS.t.Model$RangeRow)
#Fully interacive time is the only predictive model 
 adonis(MM.nirS.spec ~ Type*Time+Row+Range+RangeRow, MM.nirS.t.Model, strata = MM.nirS.t.Model$RangeRow)
 
 adonis(MM.nirS.spec ~ Type/Genotype*Time+Row+Range+RangeRow, MM.nirS.t.Model, strata = MM.nirS.t.Model$RangeRow)
 
 
 ##nosZ Gene====
 nosZALLTreatment <- read.csv("Data/nosZ-OTU-Meta.csv")
 
 nosZALL <- subset(nosZALLTreatment, Time != "R")
 #nosZALL <- subset(nosZALLTreatment, Type != "Teosinte")
 
 MM.nosZ.spec <- nosZALL[,9:7303]
 MM.nosZ.treatment <- nosZALL[,1:8]
 
 MM.nosZ.t.Model<-cbind(MM.nosZ.treatment,RangeRow=paste(MM.nosZ.treatment$Range,MM.nosZ.treatment$Row))
 MM.nosZ.t.Model<-as.data.frame(MM.nosZ.t.Model)
 
 #Model Testing
 adonis(MM.nosZ.spec ~ Genotype*Time+Row+Range+RangeRow, MM.nosZ.t.Model)
 # Df SumsOfSqs MeanSqs F.Model      R2 Pr(>F)    
 # Genotype       26    12.225 0.47020  1.0274 0.08870  0.133    
 # Time            2     1.856 0.92801  2.0277 0.01347  0.001 ***
 #   Row             1     0.786 0.78573  1.7168 0.00570  0.001 ***
 #   Range           1     0.546 0.54594  1.1929 0.00396  0.064 .  
 # RangeRow       72    33.931 0.47126  1.0297 0.24619  0.048 *  
 #   Genotype:Time  52    23.951 0.46060  1.0064 0.17378  0.351    
 # Residuals     141    64.530 0.45766         0.46820           
 # Total         295   137.825                 1.00000           
 # ---
 #   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
 adonis(MM.nosZ.spec ~ Genotype*Time+Row+Range+RangeRow, MM.nosZ.t.Model, strata = MM.nosZ.t.Model$RangeRow)
 # Df SumsOfSqs MeanSqs F.Model      R2 Pr(>F)    
 # Genotype       26    12.225 0.47020  1.0274 0.08870  0.015 *  
 #   Time            2     1.856 0.92801  2.0277 0.01347  0.001 ***
 #   Row             1     0.786 0.78573  1.7168 0.00570  0.136    
 # Range           1     0.546 0.54594  1.1929 0.00396  0.044 *  
 #   RangeRow       72    33.931 0.47126  1.0297 0.24619  0.036 *  
 #   Genotype:Time  52    23.951 0.46060  1.0064 0.17378  0.302    
 # Residuals     141    64.530 0.45766         0.46820           
 # Total         295   137.825                 1.00000           
 # ---
 #   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
 adonis(MM.nosZ.spec ~ Genotype*Time*RangeRow+Row+Range, MM.nosZ.t.Model, strata = MM.nosZ.t.Model$RangeRow)
 ##Time is only slightly explaintory!
 
 adonis(MM.nosZ.spec ~ Type*Time+Row+Range+RangeRow, MM.nosZ.t.Model, strata = MM.nosZ.t.Model$RangeRow)
 
 adonis(MM.nosZ.spec ~ Type/Genotype*Time+Row+Range+RangeRow, MM.nosZ.t.Model, strata = MM.nosZ.t.Model$RangeRow)
 
 ##Most of the change sin function are occuring as a result of teosinte. In the functional community it does seem that 
 
 
 ##Plotting the r2 of the functional ANOVA results =====
 library(ggplot2)
 library(tidyverse)
 
 FunctionalGene.ANOVA <- read.csv("Data/Functional.Gene.ANOVA.csv")
 
 ggplot(FunctionalGene.ANOVA) +
   aes(x = FunctionalGenes, fill = factor(R2)) +
   geom_bar(position = "fill")
 
 ggplot(FunctionalGene.ANOVA, aes(x=as.numeric(FunctionalGenes), y=(R2), fill=Treatment)) +geom_area()
 

ggplot(data=FunctionalGene.ANOVA, aes(x=FunctionalGenes, y=R2, fill=Treatment))  + geom_bar(stat = "identity")

ggplot(data=FunctionalGene.ANOVA, aes(x=FunctionalGenes, y=R2, fill=Treatment, order=Treatment))  +
  geom_bar(stat = "identity")+geom_col(position = position_stack(reverse = TRUE))+theme_bw()



#####DESEQ2 Analysis. Determing Differeing OTUS=====

library("DESeq2")
##DESEQ Does a 2 way two comparions
OTU16S.TvI <- subset_samples(OTU16S, Type != "Hybrid")

OTU16S.TvI@sam_data

diagdds = phyloseq_to_deseq2(OTU16S.TvI, ~ Type)
diagdds = DESeq(diagdds, test="Wald", fitType="parametric")
##Wald Signifance test using in DESEQ2 Analysis
#We are 

res = results(diagdds, cooksCutoff = FALSE)
alpha = 0.05
sigtab = res[which(res$padj < alpha), ]
sigtab = cbind(as(sigtab, "data.frame"), as(tax_table(OTU16S)[rownames(sigtab), ], "matrix"))
head(sigtab)

dim(sigtab)

theme_set(theme_bw())
scale_fill_discrete <- function(palname = "Set1", ...) {
  scale_fill_brewer(palette = palname, ...)
}
# Phylum order
x = tapply(sigtab$log2FoldChange, sigtab$Phylum, function(x) max(x))
x = sort(x, TRUE)
sigtab$Phylum = factor(as.character(sigtab$Phylum), levels=names(x))
# Genus order
x = tapply(sigtab$log2FoldChange, sigtab$Genus, function(x) max(x))
x = sort(x, TRUE)
sigtab$Genus = factor(as.character(sigtab$Genus), levels=names(x))
#Class order
x = tapply(sigtab$log2FoldChange, sigtab$Class, function(x) max(x))
x = sort(x, TRUE)
sigtab$Class = factor(as.character(sigtab$Class), levels=names(x))
#Order order
x = tapply(sigtab$log2FoldChange, sigtab$Order, function(x) max(x))
x = sort(x, TRUE)
sigtab$Order = factor(as.character(sigtab$Order), levels=names(x))

x = tapply(sigtab$log2FoldChange, sigtab$Family, function(x) max(x))
x = sort(x, TRUE)
sigtab$Family = factor(as.character(sigtab$Family), levels=names(x))

##Plotting the 
ggplot(sigtab, aes(x=Phylum, y=log2FoldChange, color=Phylum)) + geom_point(size=6) + 
  theme(axis.text.x = element_text(angle = -90, hjust = 0, vjust=0.8))

##Here Im going to order the row names
sorted_sigtab <- sigtab2[ (order(row.names(sigtab))), ]

library(gtools)
sorted_sigtab <- sigtab[ (mixedorder(row.names(sigtab))), ]

rownames(sorted_sigtab)
##Array of signifcant rownames
tax_table(OTU16S)

###This ordering is important I need to consider it. 

##Transform it to control for sample number
Type.Sig.Ave = merge_samples(OTU16S, "Type", fun = mean) #averages each OTU in all samples belonging to each habitat class
sample_data(Type.Sig.Ave)$Type <- levels(sample_data(OTU16S)$Type)
##Above we relative the OTU table 
#relative by total number of genes
Type.Sig.Ave.Trans = transform_sample_counts(Type.Sig.Ave, function(x) x/sum(x) )


##This takes the subsetted list from the Inbred vs Teosinte
Type.Sig.Ave <- subset_taxa(Type.Sig.Ave.Trans, rownames(tax_table(Type.Sig.Ave.Trans)) %in% c(rownames(sorted_sigtab)))

scale_fill_discrete <- function(palname = "Set1", ...) {
  scale_fill_brewer(palette = palname, ...)
}

plot_bar(Type.Sig.Ave, "Type", fill = "Phylum", facet_grid = ~Kingdom)+
  theme(legend.position="none")+geom_bar(stat="identity")+scale_fill_manual(values = rainbow(50))

##It seems like I am not ploting the total 
##Here Im removeing all of the NA and Archea
Type.Ave_trans.2=subset_taxa(Type.Sig.Ave, Kingdom == "Bacteria")

#Significant OTUs from types|| Relative abudnace needs to be corrected
###Figure 2 Stack plot =====
##A
plot_bar(Type.Ave_trans.2, "Type", fill = "Phylum", facet_grid = ~Kingdom)+geom_bar(stat="identity")+
  scale_fill_manual(values = rainbow(25))+
  scale_y_continuous(labels = scales::percent_format(scale=100, accuracy = 1))+
  labs(y= "Percent relative abundance", x = "Type")
#Normal Size PDF
##B
plot_bar(Type.Ave_trans.2, "Type", fill = "Order", facet_grid = ~Phylum)+
  theme(legend.position="none")+geom_bar(stat="identity")+scale_fill_manual(values = rainbow(120))+
  theme(strip.text.x = element_text(size = 8, colour = "black", angle = 90))+
  scale_y_continuous(labels = scales::percent_format(scale=100, accuracy = 1))+
  labs(y= "Percent relative abundance", x = "Type")
#5x10 Size

##ITs a similar Differences between the acido vs actino bacteria differences. 

#Subset grop and removing NA from the Data
Acidobacteria=subset_taxa(Type.Ave_trans.2, Phylum == "Acidobacteria")
Acidobacteria.NONA=subset_taxa(Acidobacteria, Order != "NA")
##I can remove the other facet if I wanted to 
##Here we will plot with in group
plot_bar(Acidobacteria, "Type", fill = "Order", facet_grid = ~Class, title="Acidobacteria")+
  #theme(legend.position="none")+
  geom_bar(stat="identity")+scale_fill_manual(values = rainbow(15))+
  theme(strip.text.x = element_text(size = 8, colour = "black", angle = 90))+
  scale_y_continuous(labels = scales::percent_format(scale=100, accuracy = 0.1))+labs(y= "Percent Abundance", x = "Type")

#Removed NA
plot_bar(Acidobacteria.NONA, "Type", fill = "Order", facet_grid = ~Class, title="Acidobacteria")+
  #theme(legend.position="none")+
  geom_bar(stat="identity")+scale_fill_manual(values = rainbow(15))+
  theme(strip.text.x = element_text(size = 8, colour = "black", angle = 90))+
  scale_y_continuous(labels = scales::percent_format(scale=100, accuracy = 0.1))+labs(y= "Percent Abundance", x = "Type")


#Acidobacteria-6

Acidobacteria6=subset_taxa(Acidobacteria, Class == "Acidobacteria-6")

plot_bar(Acidobacteria6, "Type", fill = "Order", facet_grid = ~Order, title="Acidobacteria")+
  #theme(legend.position="none")+
  geom_bar(stat="identity")+scale_fill_manual(values = rainbow(15))+
  theme(strip.text.x = element_text(size = 8, colour = "black", angle = 90))+
  scale_y_continuous(labels = scales::percent)+labs(y= "Percent Abundance", x = "Type")

#Solibacterales
Solibacterales=subset_taxa(Acidobacteria, Class == "Solibacteres")

plot_bar(Solibacterales, "Type", fill = "Family", facet_grid = ~Family, title="Solibacterales")+
  #theme(legend.position="none")+
  geom_bar(stat="identity")+scale_fill_manual(values = rainbow(15))+
  theme(strip.text.x = element_text(size = 8, colour = "black", angle = 90))+
  scale_y_continuous(labels = scales::percent)+labs(y= "Percent Abundance", x = "Type")

##Here some blanks in taxaonmy that may be useful with an updated datebase. 


###Actinobacteria Figures 
Actinobacteria=subset_taxa(Type.Ave_trans.2, Phylum == "Actinobacteria")
Actinobacteria.NONA=subset_taxa(Actinobacteria, Order != "NA")

plot_bar(Actinobacteria, "Type", fill = "Order", facet_grid = ~Class, title="Actinobacteria")+
  #theme(legend.position="none")+
  geom_bar(stat="identity")+scale_fill_manual(values = rainbow(10))+
  theme(strip.text.x = element_text(size = 8, colour = "black", angle = 90))+
  scale_y_continuous(labels = scales::percent)+labs(y= "Percent Abundance", x = "Type")

#Removed NA
plot_bar(Actinobacteria.NONA, "Type", fill = "Order", facet_grid = ~Class, title="Actinobacteria")+
  #theme(legend.position="none")+
  geom_bar(stat="identity")+scale_fill_manual(values = rainbow(10))+
  theme(strip.text.x = element_text(size = 8, colour = "black", angle = 90))+
  scale_y_continuous(labels = scales::percent)+labs(y= "Percent Abundance", x = "Type")

##Here we are trying to see differences in actinomyces known common soil bacteria that causes network like strcutres
Actinomycetales=subset_taxa(Actinobacteria, Order == "Actinomycetales")
Actinomycetales.NONA=subset_taxa(Actinomycetales, Family != "NA")

plot_bar(Actinomycetales.NONA, "Type", fill = "Family", facet_grid = ~Family, title="Actinomycetales")+
  theme(legend.position="none")+
  geom_bar(stat="identity")+scale_fill_manual(values = rainbow(33))+
  theme(strip.text.x = element_text(size = 8, colour = "black", angle = 90))+
  scale_y_continuous(labels = scales::percent)+labs(y= "Percent Abundance", x = "Type")

###Proteobacteria 
Proteobacteria=subset_taxa(Type.Ave_trans.2, Phylum == "Proteobacteria")

##These Guys are really large portions of the microbial community. 
plot_bar(Proteobacteria.NA, "Type", fill = "Order", facet_grid = ~Class, title="Proteobacteria")+
  #theme(legend.position="none")+
  geom_bar(stat="identity")+scale_fill_manual(values = rainbow(30))+
  theme(strip.text.x = element_text(size = 8, colour = "black", angle = 90))+
  scale_y_continuous(labels = scales::percent)+labs(y= "Percent Abundance", x = "Type")

##Family
Proteobacteria.NA=subset_taxa(Proteobacteria, Order != "NA")

plot_bar(Proteobacteria.NA, "Type", fill = "Order", facet_grid = ~Order, title="Proteobacteria")+
  #theme(legend.position="none")+
  geom_bar(stat="identity")+scale_fill_manual(values = rainbow(45))+
  theme(strip.text.x = element_text(size = 8, colour = "black", angle = 90))+
  scale_y_continuous(labels = scales::percent)+labs(y= "Percent Abundance", x = "Type")

#K>P>C>O>F>G>S
##This is very cool! Im kinda excited by this. 
Burkholderiales=subset_taxa(Proteobacteria, Order == "Burkholderiales")
Burkholderiales=subset_taxa(Proteobacteria, Order == "Burkholderiales")

plot_bar(Burkholderiales, "Type", fill = "Family", facet_grid = ~Family, title="Burkholderiales")+
  #theme(legend.position="none")+
  geom_bar(stat="identity")+scale_fill_manual(values = rainbow(10))+
  theme(strip.text.x = element_text(size = 8, colour = "black", angle = 90))+
  scale_y_continuous(labels = scales::percent)+labs(y= "Percent Abundance", x = "Type")

#Bdellovibrionales

Bdellovibrionales<-subset_taxa(Proteobacteria, Order == "Bdellovibrionales")
plot_bar(Bdellovibrionales, "Type", fill = "Family", facet_grid = ~Family, title="Bdellovibrionales")+
  #theme(legend.position="none")+
  geom_bar(stat="identity")+scale_fill_manual(values = rainbow(3))+
  theme(strip.text.x = element_text(size = 8, colour = "black", angle = 90))+
  scale_y_continuous(labels = scales::percent)+labs(y= "Percent Abundance", x = "Type")

###I'll just check for nitrosomoas 
Acidobacteria=subset_taxa(Type.Ave_trans.2, Phylum == "Acidobacteria")

NitrifersDESEQ2=subset_taxa(Type.Ave_trans.2, Phylum == "Nitrospirae"|Order=="Nitrosomonadales"|Order=="Nitrososphaerales")

plot_bar(NitrifersDESEQ2, "Type", fill = "Family", facet_grid = ~Family, title="Nitroso")+
  #theme(legend.position="none")+
  geom_bar(stat="identity")+scale_fill_manual(values = rainbow(3))+
  theme(strip.text.x = element_text(size = 8, colour = "black", angle = 90))+
  scale_y_continuous(labels = scales::percent_format(scale=100, accuracy = 0.01))+labs(y= "Percent Abundance", x = "Type")

##Methanotrph
MethOTU.DESEQ=subset_taxa(Type.Ave_trans.2, Order == "Metholotrpohs"|Order=="Methylophilales"|Family=="Methylobacteriaceae"|Class=="[Methylacidiphilae]")

plot_bar(MethOTU.DESEQ, "Type", fill = "Family", facet_grid = ~Family, title="Meth")+
  #theme(legend.position="none")+
  geom_bar(stat="identity")+scale_fill_manual(values = rainbow(3))+
  theme(strip.text.x = element_text(size = 8, colour = "black", angle = 90))+
  scale_y_continuous(labels = scales::percent_format(scale=100, accuracy = 0.01))+labs(y= "Percent Abundance", x = "Type")
