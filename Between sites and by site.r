### Distance to centroid between sites/by site

pacman::p_load("tidyverse", "phyloseq", "microbiome", "ggpubr")
library(vegan)
dist.3col <- getFromNamespace("dist.3col", "iCAMP")

load("data/kegg.rdata")

meta <- meta(ps_rare)
grp1 <- ifelse(sample_data(physeq_rel)$g == "C", "AS", "NS")
grp2 <- sample_data(physeq_rel) %>%
data.frame %>%
mutate(g = paste(g, site, sep="")) %>%
pull(g)


disp.env.grp1 <- select(meta, pH, NO3, NH4, TC, TN, TP, AP, TK, AK) %>% scale %>% dist %>% betadisper(grp1)
disp.env.grp2 <- select(meta, pH, NO3, NH4, TC, TN, TP, AP, TK, AK) %>% scale %>% dist %>% betadisper(grp2)

disp.ph.grp1 <- meta["pH"] %>% dist %>% betadisper(grp1)
disp.no.grp1 <- meta["NO3"] %>% dist %>% betadisper(grp1)
disp.nh.grp1 <- meta["NH4"] %>% dist %>% betadisper(grp1)
disp.tc.grp1 <- meta["TC"] %>% dist %>% betadisper(grp1)
disp.tn.grp1 <- meta["TN"] %>% dist %>% betadisper(grp1)
disp.tp.grp1 <- meta["TP"] %>% dist %>% betadisper(grp1)
disp.ap.grp1 <- meta["AP"] %>% dist %>% betadisper(grp1)
disp.tk.grp1 <- meta["TK"] %>% dist %>% betadisper(grp1)
disp.ak.grp1 <- meta["AK"] %>% dist %>% betadisper(grp1)

disp.ph.grp2 <- meta["pH"] %>% dist %>% betadisper(grp2)
disp.no.grp2 <- meta["NO3"] %>% dist %>% betadisper(grp2)
disp.nh.grp2 <- meta["NH4"] %>% dist %>% betadisper(grp2)
disp.tc.grp2 <- meta["TC"] %>% dist %>% betadisper(grp2)
disp.tn.grp2 <- meta["TN"] %>% dist %>% betadisper(grp2)
disp.tp.grp2 <- meta["TP"] %>% dist %>% betadisper(grp2)
disp.ap.grp2 <- meta["AP"] %>% dist %>% betadisper(grp2)
disp.tk.grp2 <- meta["TK"] %>% dist %>% betadisper(grp2)
disp.ak.grp2 <- meta["AK"] %>% dist %>% betadisper(grp2)

ec <- read.delim("data/EC.txt", row.names=1) %>%
mutate(EC = scale(EC))
g1 <- ec$grp
g2 <- rep(c(1:54), each=5)
dis <- ec["EC"] %>% dist
disp.ec.grp1 <- betadisper(dis, g1)
disp.ec.grp2 <- betadisper(dis, g2)

p_disp.ec.grp1 <- disp.ec.grp1$distances %>%
as.data.frame %>%
rownames_to_column("samp") %>%
mutate(group=g1, name="EC")
colnames(p_disp.ec.grp1) <- c("samp", "x", "group", "name")

p_disp.ec.grp2 <- disp.ec.grp2$distances %>%
as.data.frame %>%
rownames_to_column("samp") %>%
mutate(group=g1, name="EC")
colnames(p_disp.ec.grp2) <- c("samp", "x", "group", "name")

get_p <- function(x, name) {
    a <- x$distances %>%
    as.data.frame %>%
    rownames_to_column("samp") %>%
    mutate(group=grp1, name=name)
    colnames(a) <- c("samp", "x", "group", "name")
    a
}

p_disp.grp1 <- rbind(
    get_p(disp.ph.grp1, "pH"),
    get_p(disp.tc.grp1, "TC"),
    get_p(disp.tn.grp1, "TN"),
    get_p(disp.no.grp1, "NO"),
    get_p(disp.nh.grp1, "NH"),
    get_p(disp.tp.grp1, "TP"),
    get_p(disp.ap.grp1, "AP"),
    get_p(disp.tk.grp1, "TK"),
    get_p(disp.ak.grp1, "AK"),
    p_disp.ec.grp1
) %>%
mutate(
    # name = fct_relevel(name, c("pH", "TC", "TN", "TP", "TK", "NO", "AP", "AK", "NH")),
    name = fct_relevel(name, c("pH", "TC", "TN", "NO", "NH", "EC", "TP", "AP", "TK", "AK")),
    name = fct_recode(name,
                      EC="EC", pH = "pH", `Total carbon` = "TC", `Total nitrogen` = "TN",
                     `Total phosphorus` = "TP", `Total potassium` = "TK", `Nitrate nitrogen` = "NO",
                     `Available phosphorus` = "AP", `Available potassium` = "AK", `Ammonia nitrogen` = "NH"))

p_disp.grp2 <- rbind(
    get_p(disp.ph.grp2, "pH"),
    get_p(disp.tc.grp2, "TC"),
    get_p(disp.tn.grp2, "TN"),
    get_p(disp.no.grp2, "NO"),
    get_p(disp.nh.grp2, "NH"),
    get_p(disp.tp.grp2, "TP"),
    get_p(disp.ap.grp2, "AP"),
    get_p(disp.tk.grp2, "TK"),
    get_p(disp.ak.grp2, "AK"),
    p_disp.ec.grp2
) %>%
mutate(
    
    # name = fct_relevel(name, c("pH", "TC", "TN", "TP", "TK", "NO", "AP", "AK", "NH")),
    name = fct_relevel(name, c("pH", "TC", "TN", "NO", "NH", "EC", "TP", "AP", "TK", "AK")),
    name = fct_recode(name,
                      EC="EC", pH = "pH", `Total carbon` = "TC", `Total nitrogen` = "TN",
                     `Total phosphorus` = "TP", `Total potassium` = "TK", `Nitrate nitrogen` = "NO",
                     `Available phosphorus` = "AP", `Available potassium` = "AK", `Ammonia nitrogen` = "NH"))


ggplot(p_disp.grp1, aes(group, x, color=group)) +
geom_boxplot() +
ggpubr::stat_compare_means(label = "p.signif", comparisons = list(c("AS", "NS")), show.legend = FALSE) +
labs(x="", y="Distance to centroid between sites\nSoil chemistry") +
facet_wrap(~name, scale="free", nrow=2) +
theme_classic(base_size=20) +
theme(legend.position="none")



disp.env.grp1$distances %>%
as.data.frame %>%
rownames_to_column("samp") %>%
mutate(group=grp1) -> p_disp.env.grp1
colnames(p_disp.env.grp1) <- c("samp", "x", "group")
p1 <- ggplot(p_disp.env.grp1, aes(group, x, color=group)) +
geom_boxplot() +
ggpubr::stat_compare_means(label = "p.signif", comparisons = list(c("AS", "NS")), show.legend = FALSE) +
labs(x="", y="Distance to centroid between sites\nSoil chemistry") +
theme_classic(base_size=20) +
theme(legend.position="none")

disp.env.grp2$distances %>%
as.data.frame %>%
rownames_to_column("samp") %>%
mutate(group=grp1) -> p_disp.env.grp2
colnames(p_disp.env.grp2) <- c("samp", "x", "group")
p2 <- ggplot(p_disp.env.grp2, aes(group, x, color=group)) +
geom_boxplot() +
ggpubr::stat_compare_means(label = "p.signif", comparisons = list(c("AS", "NS")), show.legend = FALSE) +
labs(x="", y="Distance to centroid by site\nSoil chemistry") +
theme_classic(base_size=20) +
theme(legend.position="none")

ggpubr::ggarrange(p1, p2)