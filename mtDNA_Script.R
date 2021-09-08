avc <- read.csv('mtDNA_Variants/AvC - Sheet1.csv')
## For each position, select row indexes in which the exon isn't an empty
## value. 
uniques <- c()
for(x in levels(factor(avc[,'BP']))){
  y <- rownames(avc[avc[,'BP']==x & avc[,'EXON']!="-",])
  uniques <- append(uniques,y)
}
## Subset the original dataframe with the chosen indexes
avc <- avc[as.numeric(uniques),]

## To display the variant name, we want 
## [Position]:[Reference]>[Substitution]. The mutation column will be used
## later
avc$Mutation <- paste(avc$BP,avc$A2,sep = ':')
avc$Mutation <- paste(avc$Mutation,avc$A1, sep = '>')

## Order the dataframe by the order of the positions (BP)
avc <- avc[order(avc$BP, decreasing = F),]

## For the purposes of the figure, change "SYMBOL" to "MT_Locus"
colnames(avc)[17] <- 'MT_Locus'

## For the figure legend, we want the loci to be ordered in the order 
## found in the MT chromosome. 
avc$MT_Locus <- factor(avc$MT_Locus, levels = unique(avc$MT_Locus))

## We'll display the variants' differential representation with log of 
## odds ratio, using Haldane correction to make sure we're not dividing by
## 0
avc$logOR <- log((avc$Frequency_Active+0.5)/(avc$Frequency_Control+0.5))

library(ggplot2)
library(ggrepel)

## Manually assign the colors of the mito loci
genecols <- c('black','purple','orange','red','green','yellow','grey',
              'pink','skyblue','tomato','sienna','maroon')
## Graph!
ggplot(avc, aes(x=BP, y=logOR, label = Mutation, fill = MT_Locus,
                color = Consequence, shape = Consequence)) + 
  geom_point(size = 3) + theme_classic() + 
  geom_hline(yintercept = 0) + 
  geom_text_repel(size = 2.5) +
  theme(legend.text = element_text(size = 7), 
        legend.title = element_text(size = 9),
        legend.spacing.y = unit(0, "cm")) +
  scale_shape_manual(values = c(23,22,21)) +
  scale_fill_manual(values = genecols) +
  guides(fill=guide_legend(override.aes=list(color = genecols)),
         col = guide_legend(ncol = 3))

## Do the same for inactive vs control
ivc <- read.csv('mtDNA_Variants/mito_IvC_spreadsheet.csv')
uniques <- c()
for(x in levels(factor(ivc[,'BP']))){
  y <- rownames(ivc[ivc[,'BP']==x & ivc[,'EXON']!="-",])
  uniques <- append(uniques,y)
}
ivc <- ivc[as.numeric(uniques),]

ivc$Mutation <- paste(ivc$BP,ivc$A2,sep = ':')
ivc$Mutation <- paste(ivc$Mutation,ivc$A1, sep = '>')
ivc <- ivc[order(ivc$BP, decreasing = F),]
colnames(ivc)[17] <- 'MT_Locus'
ivc$MT_Locus <- factor(ivc$MT_Locus, levels = unique(ivc$MT_Locus))
ivc$logOR <- log((ivc$Frequency_Inactive+0.5)/(ivc$Frequency_Control+0.5))

genecols <- c('black','purple','orange','red')
ggplot(ivc, aes(x=BP, y=logOR, label = Mutation, fill = MT_Locus,
                color = Consequence, shape = Consequence)) + 
  geom_point(size = 3) + theme_classic() + 
  geom_hline(yintercept = 0) + 
  geom_text_repel(size = 2.5) +
  theme(legend.text = element_text(size = 7), 
        legend.title = element_text(size = 9),
        legend.spacing.y = unit(0, "cm")) +
  scale_shape_manual(values = c(23,22,21)) +
  scale_fill_manual(values = genecols) +
  guides(fill=guide_legend(override.aes=list(color = genecols)),
         col = guide_legend(ncol = 3))

## Let's say now you want to feature both active and inactive variants 
## in the same figure.. How would you go about doing that?