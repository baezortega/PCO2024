# ANALYSIS OF SOMATIC SUBSTITUTIONS IN CHOLANGIOCYTE ORGANOIDS (RE-NanoSeq)
# Adrian Baez-Ortega, Wellcome Sanger Institute, 2023-24

# This analysis is part of the publication:
# Petrus-Reurer, Lei, Tysoe, et al. (2025).
# Immunogenicity of autologous and allogeneic human primary cholangiocyte organoid cellular therapies


# Diploid total and coding genome size (from dNdScv package)
GENOME.SIZE = 6.2e9
CODING.SIZE = 69553910


# NanoSeq sample ID equivalences
SAMPLES = c("Primary"   = "PD61632b_ds0001",
            "P0 (0 d)"  = "PD61632c_ds0001",
            "P1 (28 d)" = "PD61632e_ds0001",
            "P5 (49 d)" = "PD61632g_ds0001")


# (1) Load NanoSeq mutation burdens
burdens = as.data.frame(t(sapply(SAMPLES, function(id) {
    as.matrix(read.table(list.files("data", paste0(id, "_vs_PD.*.mut_burden.tsv"), full.names=T),
                         sep="\t", header=T, as.is=T))["corrected", ]
})))


# (2) Calculate mutation burdens per genome
burdens$burden_coding = burdens$burden * CODING.SIZE
burdens$burden_genome = burdens$burden * GENOME.SIZE
burdens$burden_genome_lci = burdens$burden_lci * GENOME.SIZE
burdens$burden_genome_uci = burdens$burden_uci * GENOME.SIZE
burdens$burden_genome_post_primary = burdens$burden_genome - burdens$burden_genome[1]
burdens$burden_coding_post_primary = burdens$burden_coding - burdens$burden_coding[1]


# (2) Plot and output mutation burdens

# (a) Barplot of mutation burdens
cairo_pdf(paste0(Sys.Date(), "_NanoSeq_Burdens.pdf"), 7, 5)
par(mar=c(3.5, 6.5, 5.5, 2))
b = barplot(burdens$burden_genome,
            names.arg=rownames(burdens), ylab="Mutations per cell\n",
            col="steelblue4", border=NA, ylim=c(0, 3000), las=1, cex.lab=1.2, cex.names=1.15)
barplot(rep(burdens$burden_genome[1], nrow(burdens)), add=T, col="skyblue3", border=NA, axes=F)
abline(h=c(0, burdens$burden_genome[1]), col=c("black", "grey"), lwd=2)
segments(x0=b, y0=burdens$burden_genome_uci, y1=burdens$burden_genome_lci, lwd=2, xpd=NA)
text(b, burdens$burden_genome_uci + 330, xpd=NA, font=2,
     labels=paste0("Total: ", round(burdens$burden_genome), "\n\n"))
text(b, burdens$burden_genome_uci + 330, xpd=NA,
     labels=paste0("\n(", round(burdens$burden_genome_lci), "–", round(burdens$burden_genome_uci),
                   ")\nCoding: ", round(burdens$burden_coding)))
dev.off()


# (b) Table of mutations acquired during culture (i.e. since Primary)
burden.culture = data.frame("Sample" = rownames(burdens)[-1],
                            "Total mutations per cell acquired in culture" = 
                                round(burdens$burden_genome_post_primary[-1]),
                            "Coding mutations per cell acquired in culture" = 
                                round(burdens$burden_coding_post_primary[-1]),
                            check.names=F)

write.table(burden.culture,
            file=paste0(Sys.Date(), "_NanoSeq_Mutations_Culture.tsv"),
            sep="\t", quote=F, row.names=F)

