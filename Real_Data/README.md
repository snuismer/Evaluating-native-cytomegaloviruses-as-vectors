# Data key for raw data file

## Data description by column:
- Animal ID. A unique identifier for each animal.
- Town. The town upon which sampling was centered. All entries are Bafodia.
- Date. The date of sampling and rodent capture.
- Season. Either "Wet" or "Dry" with Wet indicating sampling between May-November and Dry indicating sampling between December-April
- MitoGenomeSpecies. Species ID from mtDNA. All entries are *M. natalensis*.
- Sex. Either "Male" or "Female". Determined at capture.
- Weight_g. The weight of the animal in grams. Not used but included for completeness.
- Tail_mm. The length of the tail in mm.
- Body_W_Tail_mm. The length of the body from snout tip to tail tip in mm.
- Body_WO_Tail_mm. The length of the body from snout tip to tail base in mm.
- Pred_ELW. Eye lens weight predicted from body length as described in the manuscript.
- Pred_Age. Animal age in days predicted using Eye lens weight as described in the manuscript.
- TotalReads. The total number of sequencing reads.
- CMV1M44_Reads. The total number of reads aligning to the M44 segment of MnatCMV1.
- CMV1M69_Reads. The total number of reads aligning to the M69 segment of MnatCMV1.
- CMV1M128_Reads. The total number of reads aligning to the M128 segment of MnatCMV1.
- CMV2M44_Reads. The total number of reads aligning to the M44 segment of MnatCMV2.
- CMV2M69_Reads. The total number of reads aligning to the M69 segment of MnatCMV2.
- CMV2M128_Reads. The total number of reads aligning to the M128 segment of MnatCMV2.
- CMV3M44_Reads. The total number of reads aligning to the M44 segment of MnatCMV3.
- CMV3M69_Reads. The total number of reads aligning to the M69 segment of MnatCMV3.
- CMV3M128_Reads. The total number of reads aligning to the M128 segment of MnatCMV3.
- CMV1_Binary. Infection status of the animal by MnatCMV1. "1" indicates the animal had more than 1x10^{-7} standardized reads aligning to all three segments of MnatCMV1. "0" indicates the animal had fewer than 1x10^{-7} standardized reads aligning to at least one of the three segments of MnatCMV1.
- CMV2_Binary. Infection status of the animal by MnatCMV2. "1" indicates the animal had more than 1x10^{-7} standardized reads aligning to all three segments of MnatCMV2. "0" indicates the animal had fewer than 1x10^{-7} standardized reads aligning to at least one of the three segments of MnatCMV2.
- CMV3_Binary. Infection status of the animal by MnatCMV3. "1" indicates the animal had more than 1x10^{-7} standardized reads aligning to all three segments of MnatCMV3. "0" indicates the animal had fewer than 1x10^{-7} standardized reads aligning to at least one of the three segments of MnatCMV3.
- CoI_Binary. Indicates a case of co-infection by MnatCMV2 and MnatCMV3. "1" indicates co-infection where both CMV2_Binary and CMV3_Binary were "1" whereas "0" indicates an absence of co-infection where either or both CMV2_Binary and CMV3_Binary were "0".


