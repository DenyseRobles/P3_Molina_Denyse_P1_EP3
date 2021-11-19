
#CARGAR LIBRERIAS
library("Biostrings")
library("msa")
library("ape")
library("seqinr")
library("ggtree")
library("treeio")
#SUBIR SECUENCIAS 
insulinas<-readAAStringSet("01_Secuencias/Insulinas.fasta")
#ALINEAMIENTO CON ALGORITMOS
insulinasMuscle<-msa(insulinas, method = "Muscle")
insulinasMuscle
insulinasClustal<-msa(insulinas, method = "ClustalW")
insulinasClustal
#INFERIR ÁRBOLES
#CONVERTIR SECUENCIAS
insulinasClustalc<-msaConvert(insulinasClustal, type = "seqinr::alignment") 
insulinasMusclec<-msaConvert(insulinasMuscle, type = "seqinr::alignment") 
#DISTANCIAS
insulinasClustald<-dist.alignment(insulinasClustalc, "identity")
insulinasMuscled<-dist.alignment(insulinasMusclec, "identity")
#MATRIZ
insulinasClustalm<-as.matrix(insulinasClustald)
insulinasMusclem<-as.matrix(insulinasMuscled)
#ÁRBOL NJ (Neighbor joining)
ClustalTree<- nj (insulinasClustalm)
plot(ClustalTree, main= "ÁRBOL FILOGENÉTICO EN CLUSTALW INSULINAS")
MuscleTree<- nj (insulinasMusclem)
plot(MuscleTree, main= "ÁRBOL FILOGENÉTICO EN MUSCLE INSULINAS")

ggtree(ClustalTree) + geom_tiplab(as_ylab=TRUE, color='black')
ggtree(MuscleTree) + geom_tiplab(as_ylab=TRUE, color='black')


