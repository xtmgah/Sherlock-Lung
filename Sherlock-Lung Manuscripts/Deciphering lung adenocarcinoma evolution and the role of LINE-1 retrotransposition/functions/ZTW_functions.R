set_wd()
libztw()

load('BBsolution_final3_short.RData')
load('Clinical/clinical_data.RData')
unique(wgs_groups_info$SP_Group)

barcode2spg <- wgs_groups_info$SP_Group
names(barcode2spg) <- wgs_groups_info$Tumor_Barcode

sp_group_color_new <- pal_npg()(4)
names(sp_group_color_new) <- c('AS_N','EU_N','EU_S','Others')
show_col(sp_group_color_new['AS_N'])
show_col(sp_group_color_new['EU_N'])
show_col(sp_group_color_new['EUS'])


sp_group_color <- pal_npg()(4)
names(sp_group_color) <- c('N_A','N_U','S_U','Others')

sp_group_lift <- names(sp_group_color_new)
names(sp_group_lift) <- names(sp_group_color)

sp_group_data <- tibble(SP_Group=names(sp_group_lift),SP_Group_New=sp_group_lift)

sp_group_major_new <- c('AS_N','EU_N','EU_S')
sp_group_major <- c('N_A','N_U','S_U')


save(barcode2spg,sp_group_color,sp_group_color_new,sp_group_lift,sp_group_data,sp_group_major_new,sp_group_major,file='ZTW_functions.RData')



## interested gene list
bar_genes <- c('UNG','TDG','TDG1','SMUG1','MBD4','NEIL1','NEIL2','NEIL3','NTHL1','OGG1','MUTYH','MPG','APE1','SETX','ADAR1','ADAR2','MLH1','MSH2','MHS6','PMS2', 'POLE', 'POLD1','PNK','XRCC1','NTH1')
hrd_genes <- c('BRCA1','BRCA2','RAD51B','RAD51C','RAD51D','XRCC2','XRCC3','BARD1','BRIP1','CHEK1','CHEK2','FAM175A','NBN','PALB2','MSH2','MSH6','GEN1','PMS2','ATR','ATM','MRE11A')
pol_gene <- c('POLB')

dnagenelist <- unique(c(bar_genes,hrd_genes))







