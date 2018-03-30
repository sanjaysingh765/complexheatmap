#load required libraries
 library(extrafont)
 #font_import() # only one time required when first time use the      library       extrafont
 #y
 fonts() 
 loadfonts()
 library(ComplexHeatmap)
 library(circlize)
library(RColorBrewer)


#upload expression data
rawdata<- read.csv("input.fpkm",header=TRUE,row.names=1,check.names = F)
head(rawdata)
dim(rawdata)

#choose the clumn you want to use for heatmap
data <- as.matrix(rawdata[, c(3:14)]) #choose the correct sample
head(data)

#convert value in log scale
heatdata <- log2(data+1)

#scale row
rm <- rowMeans(heatdata)
sx <- apply(heatdata, 1, sd)
zz <- sweep(heatdata,1,rm)
zz <- sweep(zz, 1, sx, "/")
zzz <- t(scale(t(heatdata)))

#bottom annotation
long_cn = do.call("paste0", rep(list(colnames(data)), 1))  # just to construct long text
ha_rot_cn = HeatmapAnnotation(text = anno_text(long_cn, rot = 45,  just = "right",offset = unit(18, "mm"))) #change value to change the angle and position of text


#Top annotation
df2 = data.frame(Genotype = c(rep("Sensitive", 6), rep("Resistant", 6)),
    Treatment = c(rep("Untreated", 3), rep("Treated", 3),rep("Untreated", 3),rep("Treated", 3)))
ha = HeatmapAnnotation(df = df2, col = list(Genotype = c("Sensitive" =  "#c53333", "Resistant" = "#2d72bf"),Treatment = c("Untreated" = "#d87e29","Treated" =  "#ad3ab4")),annotation_legend_param = list(title_gp = gpar(fontsize = 6), labels_gp = gpar(fontsize = 6)))

#change color for row label
row_label_col = rawdata$color


#global legend annotation
ht_global_opt(heatmap_legend_title_gp = gpar(fontsize = 6), heatmap_legend_title_position="topcenter",heatmap_legend_labels_gp = gpar(fontface = "plain",fontsize = 6))



# plot the pie chart
png("blog.png", units="in", family="Times",  width=6, height=4, res=300, pointsize = 9) #pointsize is font size| increase image size to see the key

Heatmap(zzz, 
#km = 2, for kmean clustering
name = "Color key", # legend title
rect_gp = gpar(col = "black"), # cell border
#col = colorRamp2(c(-2, 0,2),c("#259AFF","#FFFFBF", "#E81D00")),
#col = colorRamp2(c(-2, 0, 2), c("springgreen", "black", "violet")),
col = colorRamp2(c(-2, 0, 2), c("green", "black", "red")),
#heatmap_legend_param=list(color_bar="continuous", legend_direction="vertical", legend_width=unit(5,"cm"), title_position="topcenter", title_gp=gpar(fontsize=8, fontface="bold"),labels_gp = gpar(fontface = "plain",fontsize = 6)),

#Row annotation configurations
cluster_rows=TRUE,
row_title="Gene",
row_title_side="left",
row_title_gp=gpar(fontsize=10, fontface="bold"),
row_title_rot=90,
show_row_names=TRUE,
row_names_side="left",
show_row_dend=TRUE,
row_dend_side = "left",  #should the row cluster be put on the left or right of the heatmap?
row_names_gp = gpar(col = row_label_col,fontsize = 8,fontface = "italic"),  #change row label color

#Column annotation configuratiions
cluster_columns=FALSE,
column_title="Samples",
column_title_side="top",
column_title_gp=gpar(fontsize=10, fontface="bold"),
#column_title_rot=45,
show_column_names=F,
column_names_gp=gpar(fontsize=14),
show_column_dend=TRUE,


#Dendrogram configurations: columns
clustering_distance_columns=function(x) as.dist(1-cor(t(x))),
clustering_method_columns="ward.D2",  #"euclidean", "maximum", "manhattan", "canberra", "binary", "minkowski", #"pearson",     #"spearman", "kendall"  
column_dend_height=unit(30,"mm"),

#Dendrogram configurations: rows
clustering_distance_rows=function(x) as.dist(1-cor(t(x))),
clustering_method_rows="ward.D2",  #"euclidean", "maximum", "manhattan", "canberra", "binary", "minkowski", #"pearson",     #"spearman", "kendall"  
row_dend_width=unit(30,"mm"),


 #Annotations
top_annotation_height=unit(0.5,"cm"),
top_annotation=ha,
bottom_annotation_height=unit(1, "cm"),
bottom_annotation=ha_rot_cn)+
#another heatmap
Heatmap(rawdata$Group, 
name = "Group",
column_names_gp=gpar(fontsize=8),
show_column_names= TRUE,
width = unit(4, "mm"),
col =brewer.pal(3, "Set3")
#col = c("red", "green"),
)


#Add annotation names
for(an in colnames(df2)) {
    decorate_annotation(an, {
        # annotation names on the left
        grid.text(an, gp=gpar(fontsize=6,fontface = "bold"), unit(0, "npc") - unit(2, "mm"), 0.5, default.units = "npc", just = "right")
    })
}

dev.off()  



