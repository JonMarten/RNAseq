# Updating .bgen files.
The Limix pipeline requires bgen v1.2 as input. This embeds the sample ids in the file itself. Data is currently stored in an older format, requiring the .bgen files to be remade.

Additionally, the current files have no proper alternate id, using a duplicate of the RSID column instead (this may be a consequnce of the old bgen format, I'm unsure). Variants with no RSID assigned use '.' rather than a unique identifier. Additionally, there are several duplicated variants in the data which cause the pipeline to crash. 

The new .bgen files are updated to version 1.2 with sample ids embedded, duplicates removed and every variant assigned a unique alternate identifier. Genomic positions are updated to b38 of the human genome to match the mapping of the RNA seq data. **Alleles are NOT updated.** This is because aligning strands is time consuming.

The new alternate identifier is: `[chromosome]_[position_b38]_[A1]_[A2]:[rsid]`
Where:
* A1 is alphabeticaly before A2 for SNPs or
* A1 is the shorter allele for indels

Where a variant has no rsid, the alternate identifier is `[chromosome]_[position_b38]_[A1]_[A2]` and the RSID field is also set to this value.

