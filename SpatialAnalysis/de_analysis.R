library(edgeR)
library(limma)
library(stringr)
library(rlist)
library(plyr)
library(EnhancedVolcano)
library(fgsea)
library(qusage)

# Get a list of significant gene names to be shown in the volcano plot.
# Parameter n specifies number of genes in both up- and down-regulation to be shown.
# Parameter extra_genes gives a list of gene names in addition to be shown.
top_genes <- function(my_top_table, n=10, q_thresh=0.05, logFC_thresh=1, extra_genes=c()) {
    df_fc_up <- subset(my_top_table, logFC>=logFC_thresh & adj.P.Val<q_thresh)
    up_names_all <- rownames(df_fc_up[order(-df_fc_up$logFC),])
    n_select <- n
    if (n_select > length(up_names_all)) {
        n_select <- length(up_names_all)
    }
    up_names <- up_names_all[0:n_select]

    df_fc_down <- subset(my_top_table, logFC<=-logFC_thresh & adj.P.Val<q_thresh)
    down_names_all <- rownames(df_fc_down[order(df_fc_down$logFC),])
    n_select <- n
    if (n_select > length(down_names_all)) {
        n_select <- length(down_names_all)
    }
    down_names <- down_names_all[0:n_select]

    key_genes <- append(append(up_names, down_names), extra_genes)
    return(key_genes)
}

# Find the p-value cutoff related to the q-value specified in a given DE result.
determine_p_cutoff <- function(my_top_table, q_thresh=0.05) {
    q_diff <- abs(my_top_table$adj.P.Val - q_thresh)
    idx <- which(q_diff == min(q_diff))[1]
    print(paste0("Set p-value cutoff at ", my_top_table[idx, 'P.Value'], "."))
    return(my_top_table[idx, 'P.Value'])
}

# Generate .csv files of up- and down-regulated genes of a given DE result.
gen_de_csv <- function(my_top_table, prefix, all_genes=F) {
    df_table <- as.data.frame(my_top_table)
    if (!all_genes) {
        df_table <- df_table[df_table$adj.P.Val <= 0.05,]
    }

    if (nrow(df_table) > 0) {
        df_table_up <- df_table[df_table$logFC > 0,]
        if (nrow(df_table_up) > 0) {
            df_table_up <- df_table_up[order(df_table_up$adj.P.Val),]
            write.csv(df_table_up, file=paste0(prefix, ".up.csv"), quote=FALSE)
        } else {
            message(paste(prefix, "has no up-regulated gene."))
        }

        df_table_down <- df_table[df_table$logFC < 0,]
        if (nrow(df_table_down) > 0) {
            df_table_down <- df_table_down[order(df_table_down$adj.P.Val),]
            write.csv(df_table_down, file=paste0(prefix, ".down.csv"), quote=FALSE)
        } else {
            message(paste(prefix, "has no down-regulated gene."))
        }

    } else {
        message(paste(prefix, "has no significant gene after FDR control at 5%."))
    }
}

# Perform GSEA analysis based on a given DE result, and generate GSEA results into .csv files.
gen_pathway_csv <- function(my_top_table, prefix, seed) {
    pathway_record <- data.frame(cbind(c("H", "C2-CP", "C5-BP"), c("./utils/h.all.v7.1.symbols.gmt", "./utils/c2.cp.v7.1.symbols.gmt", "./utils/c5.bp.v7.1.symbols.gmt")))
    colnames(pathway_record) <- c("gene_set", "path")

    for (i in 1:nrow(pathway_record)) {
        print(paste0(i, ". Processing gene set ", pathway_record$gene_set[i], "..."))

        pathway <- read.gmt(pathway_record$path[i])
        my_ranks <- my_top_table$logFC
        names(my_ranks) <- row.names(my_top_table)

        set.seed(seed)
        fgseaRes <- fgsea(pathway, my_ranks, minSize=15, maxSize=500, nproc=1)

        if (nrow(fgseaRes) > 0) {
            df.out <- data.frame(fgseaRes[order(-NES),])
            res <- sapply(df.out$leadingEdge, function(x) {
                        paste(unlist(x), collapse=',')
                   })
            df.out$leadingEdge <- res
            write.csv(df.out, file=paste0(prefix, ".", pathway_record$gene_set[i], ".csv"))
        } else {
            message(paste(prefix, "has no gene enrichment result."))
        }

    }
}

# DE analysis using limma.
# df_seg gives AOI information; df_cnt gives gene-count matrix over AOIs.
# conditions specifies selection condition on AOIs to be considered in the analysis.
# field_name specifies which categorical variable to be considered in the analysis.
# field_levels gives all the levels of field_name variable to be used for building the linear model,
# and the contrast to be estimated is the last level minus the first level in field_levels.
analyze_by_contrasts <- function(df_seg, df_cnt, conditions, field_name, field_levels) {
    df_seg_cond <- subset(df_seg, segment %in% conditions)
    print(table(df_seg_cond[,field_name]))

    snames <- rownames(df_seg_cond)
    cnt_mat <- df_cnt[, snames]

    my.dge <- DGEList(cnt_mat)
    my.dge <- calcNormFactors(my.dge)
    field_vec <- factor(df_seg_cond[snames, field_name], levels=field_levels)

    my.model <- model.matrix(~0 + field_vec)
    colnames(my.model) <- paste(rep(field_name, length(levels(field_vec))), levels(field_vec), sep='')
    y <- voom(my.dge, my.model, plot=F)
    my.fit <- lmFit(y, my.model)

    contr.str <- paste0(tail(colnames(coef(my.fit)), n=1), "-", colnames(coef(my.fit))[1])
    print(contr.str)
    my.contr <- makeContrasts(contrasts=contr.str, levels=colnames(coef(my.fit)))
    tmp <- contrasts.fit(my.fit, my.contr)
    tmp <- eBayes(tmp)
    top.table <- topTable(tmp, sort.by='P', n=Inf)

    return(top.table)
}
