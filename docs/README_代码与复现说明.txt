代码与复现说明

1. 原始 TCGA/GEO/BeatAML/GSE116256/STRING 数据需要从公共数据库获取。
2. 本投稿包不包含公共数据库原始大文件。
3. 本包包含派生结果表、关键统计结果表和脚本，可用于复核主要结果。
4. AS 分析未执行，因为当前项目没有 PSI/event-level 输入。
5. de novo single-cell reanalysis 未执行，因为没有 Seurat/h5ad/10x/count object 或 UMAP 坐标。
6. formal immune deconvolution 不作为主结论；bulk 结果按 signature-based inference 表述。
7. BeatAML 结果为 ex vivo pharmacogenomic association，不是 clinical response 或 treatment guidance。
8. 固定六基因模型 CLCN5, ITGB2, ARHGEF5, TRIM32, SAT1, ACOX2 不可随意替换。
