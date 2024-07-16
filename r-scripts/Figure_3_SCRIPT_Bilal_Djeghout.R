#' ---
#' title: Campylobacter detection analysis
#' self-contained: true
#' ---
#+ setup, warning=F

#' ## Load data and packages
knitr::opts_chunk$set(warning = FALSE, out.width = '100%', dpi=300, fig.cap = "")
#' Load libraries that we will need

library(data.table) # For cleaning and coding data
library(readxl)     # For reading data from excel files
library(ggplot2)    # For making graphs
library(ggbeeswarm) # For making beeswarm graphs.
#' Load data
dat <- read_excel(path = "MetaGen_Dataset_BD.xlsx", skip = 1, .name_repair = "universal")[1:66,]

#' Change dataset into a data.table for easier data wrangling
setDT(dat)

#' This is how I make a simple cross-tab on a subset of the data. Here I'm comparing the
#' Bristol Stool Scale with a variable defined by `.__Campylobacter_isolates`>0

dat[ EPA_results=="Ca+" & Type=="Raw", table(Bristol__scale , `.__Campylobacter_isolates`>0 )]

#' We could try to visualise this with a stacked bar but I don't think it's helpful.
dat[ EPA_results=="Ca+" ] |>
  ggplot() +
  aes(x=Bristol__scale , fill=`.__Campylobacter_isolates`>0) +
  geom_bar(position = "stack", col="black") +
  facet_wrap(~Type)

#' Anyway, from here there's no evidence that the rate of culture detection
#' changes with Bristol stool scale. (There's a typo in the table (BS=4,
#' Type=Raw, 2/2 != 50%))

#' ## Encode the data

#' I'm using data.table syntax to create the variables within the data that I will use for plotting.


#' I will leave the original data as it is, and make new variables when I want to change something.

dat[, detectByCulture := `.__Campylobacter_isolates`>0]
dat[, detectByqPCR    := qPCR_ct != "Not_detected"]
dat[, detectBySequencing := as.numeric(Number_of_Contigs)>0][is.na(detectBySequencing), detectBySequencing:=0]
dat[, detectBySequencing := as.numeric(Completeness_of_genome_...)>0][is.na(detectBySequencing), detectBySequencing:=0]

dat[, goldStandard := EPA_results == "Ca+"]

dat[, includeSeq := Platform != "Library_Fail"]

# For these, where the value is not numeric it will be replaced with NA which is OK.
dat[, completeness := as.numeric(`Completeness_of_genome_...`)]
dat[, coverage := as.numeric(`Coverage.._X`)]
dat[, totalReads := as.numeric(Total_._reads_.millions.)]
dat[, campyReads := as.numeric(`._of_Campylobacter_Genus_reads`)]
dat[, propCampyReads := as.numeric(`Percentage.of.Campylobacter.reads.out.of.total.reads.sequenced`)]
dat[, propCampyReads := campyReads / totalReads]

# Species
dat[, campySpecies := `Campy._MAG_species`]
dat[campySpecies=="C. Jejuni", campySpecies := "C. jejuni"]
# This variable will hold the species.
dat[, campySpeciesDetail := campySpecies]
# This species variable is now just for whether sample is resolved to species level.
dat[campySpecies=="C. coli", campySpecies := "C. jejuni"]
dat[campySpecies=="C. jejuni", campySpecies := "Detected"]

# Sequence type
dat[, mlst := !is.na(as.numeric(Campy._MAG_ST))]

dat[, MAGst := as.numeric(Campy._MAG_ST)]
dat[, cultureST := `pubMLST_.ST.`]

# Now a variable for the level to which the MAGs are resolved.
dat[, resolution := "Not detected or resolved"]
dat[campySpecies=="Detected", resolution := "Species level"]
dat[mlst==TRUE, resolution := "MLST level"]

dat[ , N50_old := N50]
dat[ , N50 := as.numeric(N50)]
dat[ , CT := as.numeric(qPCR_ct)]

#' ## Table 1

#' Now Table 1 as a graph.  I think this is cleaner, allows deeper understanding.

dat[goldStandard==TRUE, .(Patient_No, Type, goldStandard, detectBySequencing, detectByCulture, detectByqPCR, Bristol__scale)] |>
  melt(id.vars = c("Patient_No", "Type", "goldStandard", "Bristol__scale")) |>
  ggplot() +
    aes(y=factor(Patient_No), x=variable, fill=value) +
    geom_tile(col="grey") +
    facet_grid(Bristol__scale~Type, space="free", scale="free", switch = "y") +
    scale_fill_manual(values=c("blue", "red")) +
  theme_bw() + labs(x="Method", y="Bristol Stool Scale", fill="Detected")+
  theme(axis.text.x = element_text(angle=30, hjust=1)) +
  theme(axis.text.y = element_blank()) +
  theme(axis.ticks.y=element_blank())


#' You can make a cross-tab of the detection by qPCR vs detection by culture (here restricted to raw samples that are known to the Ca+)
dat[Type=="Raw" & goldStandard==TRUE , table(detectByCulture, detectByqPCR)]

#' And if you want a statistical test of the difference between methods use a McNemar test:
dat[Type=="Raw" & goldStandard==TRUE , table(detectByCulture, detectByqPCR) |> mcnemar.test()]


#' ## Table 2
#'
#' Table 2 as a set of graphs.
#'
#' First I just wanted to understand the relationships between the variables.
#' What is the relationship between coverage and the number of campylobacter reads?
#' For most it seems to be linear but for a few it doesn't correspond?

plot(dat$coverage, dat$campyReads, log="xy")

#' Now work out which  This doesn't exactly correspond to table 2 so will need to refine definition.
dat[, table(includeSeq, Type, goldStandard)]

#' We want to show the number of reads, coverage, proportion of reads, completeness stratified by characterised.

#' First do a single outcome:
dat[(includeSeq)] |>
  ggplot() +
    aes(x=campySpecies, y=campyReads,shape=Type) +
    #geom_beeswarm() +
    geom_point(position=position_dodge(width=0.3)) +
    scale_y_log10() +
    theme_bw() +
    labs(y="Campylobacter Genus Reads")

#' Now look at all the quality outcomes together.  This works but we don't
#' really want 'completeness' on a log scale so need to do something about that.
#' To get all the measures into one graph we can 'melt' the data..
#' We can use the shape and the colour of the points to compare different groups.

dat[(includeSeq) , .(campySpecies,totalReads, campyReads, propCampyReads, N50,completeness, coverage, goldStandard, Type, Patient_No)] |>
  melt(id.vars=c("goldStandard", "Type", "Patient_No","campySpecies")) |>
  ggplot() +
  aes(x=campySpecies, y=value,shape=goldStandard, color=Type) +
  geom_point(position=position_dodge(width=0.3))+
  #geom_beeswarm() +
  scale_y_log10() +
  scale_color_manual(values=c("black", "red")) +
  theme_bw() +
  facet_wrap(~variable,scale="free_y", nrow=2)

#' ## Compare the filtered with the unfiltered results.

#' Here I make a similar graph but connect the filtered with the raw samples in
#' the cases where both are included.  Where filtering worked, the coverage and
#' completeness was substantially higher

dat[(includeSeq & goldStandard) , .(campySpecies,totalReads, campyReads, propCampyReads, N50,completeness=(completeness), coverage, goldStandard, Type, Patient_No)] |>
  melt(id.vars=c("goldStandard", "Type", "Patient_No","campySpecies")) |>
  ggplot() +
  aes(x=Type, y=value, group=Patient_No) + geom_point(aes(col=campySpecies)) + geom_line() +
  facet_wrap(~variable, scale="free_y") +
  scale_y_log10() +
  scale_color_manual(values=c("black", "red")) +
  theme_bw()


####### Start of code for new plot

library(patchwork)
library(ggtext)

blanklabelplot<-ggplot()+labs(x="Type of stool processing method")+theme_classic()+
  guides(x = "none", y = "none")
labels = c("totalReads"="Number of total reads",
           "campyReads"="Number of *Campylobacter* reads",
           "propCampyReads"="*Campylobacter* reads proportion",
           "N50"="Mean contig length",
           "completeness"="*Campylobacter* genome<br>completeness (%)",
           "coverage"="*Camylobacter* coverage (x)")
dat2 <- dat[(includeSeq & goldStandard) , .(campySpecies,totalReads, campyReads, propCampyReads, N50,completeness=(completeness), coverage, goldStandard, Type, Patient_No)]
plots <- lapply(c("totalReads", "campyReads", "propCampyReads", "N50","completeness","coverage") , \(v) {
   ggplot(dat2) +
     aes(y=get(v) , x=Type, col=campySpecies) +
     scale_color_manual(values=c("black", "red"),
                        labels=c("*Campylobacter* species identified",
                                 "*Campylobacter* species not identified"),
                        name="MAG *Campylobacter* species detection") +
     geom_point() +
     scale_x_discrete(labels=c("Filter"="Filtered", "Raw"="Unfiltered"))+
     labs(x=NULL,y=labels[v])+
     scale_y_log10() +
     geom_line(aes(group=Patient_No)) + theme_bw() +
     theme(legend.position = "bottom", axis.title.y=element_markdown(),
           legend.text = element_markdown(),
           legend.title = element_markdown())
 } )
plots[[3]] <- plots[[3]] + scale_y_continuous()
(wrap_plots(plots) + plot_layout(tag_level = "new")) / blanklabelplot +
   plot_layout(guides="collect", heights=c(100,1), tag_level = "new") +
   plot_annotation(tag_levels=list(c("",""),LETTERS[1:6])) &
   theme(legend.position="right")

ggsave("graph.png", width=10,height=6,dpi="retina")
ggsave("graph.tiff", width=10,height=6,dpi="retina")


##### End of code for new plot.


#' Look at the correlation between some of the measures:
#'

#' I like the next picture. MAG completeness and identification success strongly linked to read count.
dat[(includeSeq )] |> ggplot() +
  aes(x=campyReads, y=completeness, col=resolution, shape=goldStandard) +
  geom_point() +
  scale_color_manual(values=c("black", "red","blue")) +
  scale_x_log10() +
  theme_bw()

## What is different about the points that don't conform to the line?
dat[(includeSeq )] |> ggplot() +
  aes(x=campyReads, y=coverage, col=resolution, shape=goldStandard) +
  geom_point() +
  scale_color_manual(values=c("black", "red","blue")) +
  scale_x_log10() + scale_y_log10()+
  theme_bw()

#' Same as above but with detection of AMR as the output represented by the shape.
dat[(includeSeq )] |> ggplot() +
  aes(x=campyReads, y=completeness, col=resolution, shape=(Campy._MAG_AMR_genes!="ND")) +
  geom_point() +
  scale_color_manual(values=c("black", "red","blue")) +
  scale_x_log10() +
  theme_bw()

#' You could facet by type if you wanted.
dat[(includeSeq )] |> ggplot() +
  aes(x=campyReads, y=completeness, col=resolution, shape=(Campy._MAG_AMR_genes=="ND")) +
  geom_point() +
  scale_color_manual(values=c("black", "red","blue")) +
  scale_x_log10() +
  theme_bw()+
  facet_wrap(~Type)


#' ## Comparison between culture and sequencing

#' Here I make a short table just showing the samples we are interested in.

dat[(includeSeq & goldStandard & Type=="Raw" & resolution!="Not detected or resolved"),
     .(MAGst, cultureST, Campy._MAG_AMR_genes, starAMR.._Campy._genome_AMR_genes)] |> knitr::kable()

dat[(includeSeq & goldStandard & Type=="Raw" & resolution!="Not detected or resolved"),
    .(Campy._MAG_AMR_genes, starAMR.._Campy._genome_AMR_genes)] |> knitr::kable()


dat[(includeSeq & goldStandard & Type=="Filter" & resolution!="Not detected or resolved"),
    .(Campy._MAG_AMR_genes, starAMR.._Campy._genome_AMR_genes)] |> knitr::kable()


#' Other quality measures

#' N50 and CT (need to add the 'not detected' CTs back in.)

dat[(goldStandard)] |>
  ggplot() +
  aes(x=N50, y=CT, color=resolution) + geom_point() + scale_x_log10() + facet_wrap(~Type) +
  theme_bw()  +
  scale_color_manual(values=c("black", "red", "blue"))

dat[(goldStandard)] |>
  ggplot() +
  aes(x=N50, y=campyReads, color=resolution) + geom_point() +
  scale_x_log10() +
  scale_y_log10() +
  facet_wrap(~Type) +
  theme_bw()  +
  scale_color_manual(values=c("black", "red", "blue"))

