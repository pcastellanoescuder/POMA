
PomaVolcano <- function(data,
                        pval = c("raw", "adjusted"),
                        Pval_cutoff = 0.05,
                        FC_cutoff = 1.5,
                        xlim = 2){

  if (missing(pval)) {
    pval <- "raw"
    warning("pval argument is empty! Raw p-value will be used")
  }
  if (!(pval %in% c("raw", "adjusted"))) {
    stop(crayon::red(clisymbols::symbol$cross, "Incorrect value for pval argument!"))
  }

  names <- rownames(data)

  if(pval == "raw"){
    df <- data.frame(P.Value = data$P.Value, FC = data$Fold_Change_Ratio, names = names)
  }
  else{
    df <- data.frame(P.Value = data$adj.P.Val, FC = data$Fold_Change_Ratio, names = names)
  }


  df <- mutate(df, threshold = as.factor(ifelse(df$P.Value >= Pval_cutoff,
                                                yes = "none",
                                                no = ifelse(df$FC < FC_cutoff,
                                                            yes = ifelse(log2(df$FC) < -log2(FC_cutoff),
                                                                         yes = "Down-regulated",
                                                                         no = "none"),
                                                            no = "Up-regulated"))))

  ggplot(data = df, aes(x = log2(FC), y = -log10(P.Value), colour = threshold)) +
    geom_point(size=1.75) +
    xlim(c(-(xlim), xlim)) +
    xlab("log2 fold change") +
    ylab("-log10 p-value") +
    scale_y_continuous(trans = "log1p")+
    ggtitle("Comparisson: Group2/Group1") +
    ggrepel::geom_label_repel(data = df[df$P.Value < Pval_cutoff & (df$FC > FC_cutoff | log2(df$FC) < -log2(FC_cutoff)),],
              aes(x = log2(FC), y = -log10(P.Value), label = names), show.legend = FALSE) +
    geom_vline(xintercept = -log2(FC_cutoff), colour = "black") +
    geom_vline(xintercept = log2(FC_cutoff), colour = "black") +
    geom_hline(yintercept = -log10(Pval_cutoff), colour = "black") +
    theme(legend.position = "none") +
    theme_minimal() +
    scale_color_manual(values = c("Down-regulated" = "#E64B35", "Up-regulated" = "#3182bd", "none" = "#636363"))

  }

