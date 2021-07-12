# mqpipe

Pipeline for internship at Gabriela Chiosis Lab. Extracts Protein IDs, gene names, peptide intensities, and peptide IDs from raw Maxquant data in TSV format. Then joins the clean protein data with the sample manifest and filters sample entries by percentage or by raw NA count. The pipeline performs an imputation on NA values. Finally, it adds gene names and a description to their respective proteins.
Work in progress.
