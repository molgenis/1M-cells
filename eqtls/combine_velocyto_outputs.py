import loompy as loompy

looms_v2 = ['/groups/umcg-bios/tmp01/projects/1M_cells_scRNAseq/ongoing/splicing/lanes/HG19/180920_lane1/velocyto/180920_lane1.loom', '/groups/umcg-bios/tmp01/projects/1M_cells_scRNAseq/ongoing/splicing/lanes/HG19/180925_lane1/velocyto/180925_lane1.loom', '/groups/umcg-bios/tmp01/projects/1M_cells_scRNAseq/ongoing/splicing/lanes/HG19/180925_lane2/velocyto/180925_lane2.loom', '/groups/umcg-bios/tmp01/projects/1M_cells_scRNAseq/ongoing/splicing/lanes/HG19/180926_lane1/velocyto/180926_lane1.loom', '/groups/umcg-bios/tmp01/projects/1M_cells_scRNAseq/ongoing/splicing/lanes/HG19/180926_lane2/velocyto/180926_lane2.loom', '/groups/umcg-bios/tmp01/projects/1M_cells_scRNAseq/ongoing/splicing/lanes/HG19/180927_lane1/velocyto/180927_lane1.loom', '/groups/umcg-bios/tmp01/projects/1M_cells_scRNAseq/ongoing/splicing/lanes/HG19/180927_lane2/velocyto/180927_lane2.loom', '/groups/umcg-bios/tmp01/projects/1M_cells_scRNAseq/ongoing/splicing/lanes/HG19/180928_lane1/velocyto/180928_lane1.loom', '/groups/umcg-bios/tmp01/projects/1M_cells_scRNAseq/ongoing/splicing/lanes/HG19/180928_lane2/velocyto/180928_lane2.loom', '/groups/umcg-bios/tmp01/projects/1M_cells_scRNAseq/ongoing/splicing/lanes/HG19/181003_lane1/velocyto/181003_lane1.loom', '/groups/umcg-bios/tmp01/projects/1M_cells_scRNAseq/ongoing/splicing/lanes/HG19/181003_lane2/velocyto/181003_lane2.loom', '/groups/umcg-bios/tmp01/projects/1M_cells_scRNAseq/ongoing/splicing/lanes/HG19/181003_lane3/velocyto/181003_lane3.loom', '/groups/umcg-bios/tmp01/projects/1M_cells_scRNAseq/ongoing/splicing/lanes/HG19/181011_lane1/velocyto/181011_lane1.loom', '/groups/umcg-bios/tmp01/projects/1M_cells_scRNAseq/ongoing/splicing/lanes/HG19/181011_lane2/velocyto/181011_lane2.loom', '/groups/umcg-bios/tmp01/projects/1M_cells_scRNAseq/ongoing/splicing/lanes/HG19/181022_lane1/velocyto/181022_lane1.loom', '/groups/umcg-bios/tmp01/projects/1M_cells_scRNAseq/ongoing/splicing/lanes/HG19/181022_lane2/velocyto/181022_lane2.loom', '/groups/umcg-bios/tmp01/projects/1M_cells_scRNAseq/ongoing/splicing/lanes/HG19/181023_lane1/velocyto/181023_lane1.loom', '/groups/umcg-bios/tmp01/projects/1M_cells_scRNAseq/ongoing/splicing/lanes/HG19/181023_lane2/velocyto/181023_lane2.loom', '/groups/umcg-bios/tmp01/projects/1M_cells_scRNAseq/ongoing/splicing/lanes/HG19/181024_lane1/velocyto/181024_lane1.loom', '/groups/umcg-bios/tmp01/projects/1M_cells_scRNAseq/ongoing/splicing/lanes/HG19/181024_lane2/velocyto/181024_lane2.loom', '/groups/umcg-bios/tmp01/projects/1M_cells_scRNAseq/ongoing/splicing/lanes/HG19/181024_lane3/velocyto/181024_lane3.loom', '/groups/umcg-bios/tmp01/projects/1M_cells_scRNAseq/ongoing/splicing/lanes/HG19/181030_lane1/velocyto/181030_lane1.loom', '/groups/umcg-bios/tmp01/projects/1M_cells_scRNAseq/ongoing/splicing/lanes/HG19/181030_lane2/velocyto/181030_lane2.loom', '/groups/umcg-bios/tmp01/projects/1M_cells_scRNAseq/ongoing/splicing/lanes/HG19/181106_lane2/velocyto/181106_lane2.loom', '/groups/umcg-bios/tmp01/projects/1M_cells_scRNAseq/ongoing/splicing/lanes/HG19/181107_lane1/velocyto/181107_lane1.loom', '/groups/umcg-bios/tmp01/projects/1M_cells_scRNAseq/ongoing/splicing/lanes/HG19/181107_lane2/velocyto/181107_lane2.loom', '/groups/umcg-bios/tmp01/projects/1M_cells_scRNAseq/ongoing/splicing/lanes/HG19/181108_lane1/velocyto/181108_lane1.loom', '/groups/umcg-bios/tmp01/projects/1M_cells_scRNAseq/ongoing/splicing/lanes/HG19/181108_lane2/velocyto/181108_lane2.loom', '/groups/umcg-bios/tmp01/projects/1M_cells_scRNAseq/ongoing/splicing/lanes/HG19/181122_lane1/velocyto/181122_lane1.loom', '/groups/umcg-bios/tmp01/projects/1M_cells_scRNAseq/ongoing/splicing/lanes/HG19/181122_lane2/velocyto/181122_lane2.loom', '/groups/umcg-bios/tmp01/projects/1M_cells_scRNAseq/ongoing/splicing/lanes/HG19/181213_lane1/velocyto/181213_lane1.loom', '/groups/umcg-bios/tmp01/projects/1M_cells_scRNAseq/ongoing/splicing/lanes/HG19/181213_lane2/velocyto/181213_lane2.loom', '/groups/umcg-bios/tmp01/projects/1M_cells_scRNAseq/ongoing/splicing/lanes/HG19/181213_lane3/velocyto/181213_lane3.loom', '/groups/umcg-bios/tmp01/projects/1M_cells_scRNAseq/ongoing/splicing/lanes/HG19/181214_lane1/velocyto/181214_lane1.loom', '/groups/umcg-bios/tmp01/projects/1M_cells_scRNAseq/ongoing/splicing/lanes/HG19/181214_lane2/velocyto/181214_lane2.loom', '/groups/umcg-bios/tmp01/projects/1M_cells_scRNAseq/ongoing/splicing/lanes/HG19/181218_lane1/velocyto/181218_lane1.loom', '/groups/umcg-bios/tmp01/projects/1M_cells_scRNAseq/ongoing/splicing/lanes/HG19/181218_lane2/velocyto/181218_lane2.loom', '/groups/umcg-bios/tmp01/projects/1M_cells_scRNAseq/ongoing/splicing/lanes/HG19/181219_lane1/velocyto/181219_lane1.loom', '/groups/umcg-bios/tmp01/projects/1M_cells_scRNAseq/ongoing/splicing/lanes/HG19/181219_lane2/velocyto/181219_lane2.loom', '/groups/umcg-bios/tmp01/projects/1M_cells_scRNAseq/ongoing/splicing/lanes/HG19/181219_lane3/velocyto/181219_lane3.loom']

loom_v2_merged = '/groups/umcg-bios/tmp01/projects/1M_cells_scRNAseq/ongoing/splicing/objects/loom_v2_UT_CA_and_others.loom'

loompy.combine(looms_v2, loom_v2_merged, key="Accession")

looms_v3 = ['/groups/umcg-bios/tmp01/projects/1M_cells_scRNAseq/ongoing/splicing/lanes/HG19/190109_lane1/velocyto/190109_lane1.loom', '/groups/umcg-bios/tmp01/projects/1M_cells_scRNAseq/ongoing/splicing/lanes/HG19/190109_lane2/velocyto/190109_lane2.loom', '/groups/umcg-bios/tmp01/projects/1M_cells_scRNAseq/ongoing/splicing/lanes/HG19/190110_lane1/velocyto/190110_lane1.loom', '/groups/umcg-bios/tmp01/projects/1M_cells_scRNAseq/ongoing/splicing/lanes/HG19/190115_lane1/velocyto/190115_lane1.loom', '/groups/umcg-bios/tmp01/projects/1M_cells_scRNAseq/ongoing/splicing/lanes/HG19/190115_lane2/velocyto/190115_lane2.loom', '/groups/umcg-bios/tmp01/projects/1M_cells_scRNAseq/ongoing/splicing/lanes/HG19/190122_lane1/velocyto/190122_lane1.loom', '/groups/umcg-bios/tmp01/projects/1M_cells_scRNAseq/ongoing/splicing/lanes/HG19/190123_lane2/velocyto/190123_lane2.loom', '/groups/umcg-bios/tmp01/projects/1M_cells_scRNAseq/ongoing/splicing/lanes/HG19/190124_lane1/velocyto/190124_lane1.loom', '/groups/umcg-bios/tmp01/projects/1M_cells_scRNAseq/ongoing/splicing/lanes/HG19/190124_lane2/velocyto/190124_lane2.loom', '/groups/umcg-bios/tmp01/projects/1M_cells_scRNAseq/ongoing/splicing/lanes/HG19/190130_lane1/velocyto/190130_lane1.loom', '/groups/umcg-bios/tmp01/projects/1M_cells_scRNAseq/ongoing/splicing/lanes/HG19/190130_lane2/velocyto/190130_lane2.loom', '/groups/umcg-bios/tmp01/projects/1M_cells_scRNAseq/ongoing/splicing/lanes/HG19/190204_lane1/velocyto/190204_lane1.loom', '/groups/umcg-bios/tmp01/projects/1M_cells_scRNAseq/ongoing/splicing/lanes/HG19/190204_lane2/velocyto/190204_lane2.loom']

loom_v3_merged = '/groups/umcg-bios/tmp01/projects/1M_cells_scRNAseq/ongoing/splicing/objects/loom_v3_UT_CA_and_others.loom'

loompy.combine(looms_v3, loom_v3_merged, key="Accession")
