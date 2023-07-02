#  Authors 
#  Alexander Freudenberg, alexander.freudenberg@stads.de

#  Copyright (C) 2023 Alexander Freudenberg, Torsten Pook

#  Licensed under the Apache License, Version 2.0 (the "License");
#  you may not use this file except in compliance with the License.
#  You may obtain a copy of the License at

#     http://www.apache.org/licenses/LICENSE-2.0

#  Unless required by applicable law or agreed to in writing, software
#  distributed under the License is distributed on an "AS IS" BASIS,
#  WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
#  See the License for the specific language governing permissions and
#  limitations under the License.


# =====================
# Load packages
# =====================

set.seed(123)
data_dir <- "./data/"
if(!dir.exists(data_dir)){
    stop("Modify variable data_dir to point to an existing directory, which will be used for storing the simulated data.")
}

if(!require(MoBPS)){
    devtools::install_github("tpook92/MoBPS", subdir="pkg")
}
if(!require(MoBPSmaps)){
    devtools::install_github("tpook92/MoBPS", subdir="pkg-maps")
}
if(!require(stringr)){
    install.packages("stringr")
}

# =====================
# Auxiliary function
# =====================

# Hash table for looking up PLINK byte values in packed 2bit format
conversion_table_2bit_plink <- as.raw(c(0, 2, 3, 85, 8, 10, 11, 85, 12, 14, 15, 85, 85, 85, 85, 85, 32, 34, 35, 85, 40, 42, 43, 85, 44, 46, 47, 85, 85, 85, 85, 85, 48, 50, 51, 85, 56, 58, 59, 85, 60, 62, 63, 85, 85, 85, 85, 85, 85, 85, 85, 85, 85, 85, 85, 85, 85, 85, 85, 85, 85, 85, 85, 85, 128, 130, 131, 85, 136, 138, 139, 85, 140, 142, 143, 85, 85, 85, 85, 85, 160, 162, 163, 85, 168, 170, 171, 85, 172, 174, 175, 85, 85, 85, 85, 85, 176, 178, 179, 85, 184, 186, 187, 85, 188, 190, 191, 85, 85, 85, 85, 85, 85, 85, 85, 85, 85, 85, 85, 85, 85, 85, 85, 85, 85, 85, 85, 85, 192, 194, 195, 85, 200, 202, 203, 85, 204, 206, 207, 85, 85, 85, 85, 85, 224, 226, 227, 85, 232, 234, 235, 85, 236, 238, 239, 85, 85, 85, 85, 85, 240, 242, 243, 85, 248, 250, 251, 85, 252, 254, 255))

#' Export to PLINK binary format
#'
#' This function exports a genotype matrix created by MoBPS to PLINK binary format (.bed).
#'
#' @param geno Genotype dataset returned by the get.geno function in MoBPS.
#' @param pheno Phenotype vector of the population. Defaults all individuals to be controls.
#' @param filename Root filename of the output files.
#'
#' @return Nothing, except if filename.bed already exists, in which case it returns the character vector.
mobps_to_bed <- function(geno, pheno = rep("1", nrow(geno)), filename = "simulation"){

    nindiv <- ncol(geno)
    nsnps <- nrow(geno)

    # Convert MoBPS genotype data to byte format
    geno_raw <- matrix(as.raw(geno),ncol = nindiv)

    # Set-up byte-sized matrix in PLINK dimensions
    twobit_packed <- matrix(as.raw(0), ncol = nindiv/4, nrow = nsnps)

    # Shift each column by the number of bits that is required to store four individuals in one byte
    # Afterwards, do a bitwise or operation to pack the individuals and store the result in the respective column of the new matrix  
    cat("Pack values.\n")
    for(col in seq(1,nindiv, by = 4)){
        # This will holds four columns packed into one
        col_packed <- raw(nsnps)
        for(i in 0:3){
            col_packed <- col_packed | rawShift(geno_raw[, col + i], i * 2)
        }
        new_col_index <- ceiling(col/4)
        twobit_packed[,new_col_index] <- col_packed   
    }
    # Convert twobit to PLINK format
    # FIXME: the as.integer instructions casts all values into 32bit integers - yet, subscripts with raws are not implemented in R... 
    cat("Convert to PLINK format.\n")
    plink_format <- conversion_table_2bit_plink[as.integer(twobit_packed)+1] 
    plink_format <- matrix(plink_format, nrow = nsnps)
    # Sanity checks
    # Are any column of plink_format fully zero?
    if(any(apply(plink_format == as.raw(0), 2, all))){
        error("At least one column in packed format is all zero.")
    }
    # Were some values mistakenly coded as missings
    # In this case, the hash table above is false as MoBPS doesn't return NAs
    if(any(plink_format == as.raw(85))){
        stop("Some genotypes were coded as missing.")
    }

    # Check if file in binary format already exists 
    filename_bed <- paste0(filename, ".bed")
    if(file.exists(filename_bed)){
        warning(paste("File", filename, "already exists. Returning packed character matrix."))
        return(plink_format)
    }

    # Create file and open connection
    cat("Write files.\n")
    file.create(filename_bed)
    con <- file(filename_bed, "wb")

    header_bytes <- as.raw(c(0x6c, 0x1b, 0x01))
    writeBin(header_bytes, con)

    # Transpose plink_char as it is column-major and PLINK is SNP-major
    plink_write_format <- as.vector(t(plink_format))
    writeBin(plink_write_format, con, useBytes=TRUE)

    close(con)

    # Write bim and fam files
    filename_bim <- paste0(filename, ".bim")
    filename_fam <- paste0(filename, ".fam")
    
    # Replace all underscores by tildes (see PLINK docs)
    ids <- str_replace_all(colnames(geno),"_", "~")

    # Extract and replace sex of all individuals - coded as F and M in indiv IDs
    sex <- str_sub_all(colnames(geno), 1, 1)
    sex <- str_replace_all(unlist(sex), c("M" = "1", "F" = "2"))

    # Create tables for writing to files
    bim_table <- data.frame(chr = 1, id = rownames(geno), pos = 0, bp = 1:nsnps, al1 = "D", al2 = "d")
    fam_table <- data.frame(fam_id = ids, own_id = ids, fat_id = 0, mot_id = 0, sex = sex, pheno = pheno)

    # Write to files
    write.table(bim_table, filename_bim, quote = F, sep = "\t", row.names = F, col.names = F)
    write.table(fam_table, filename_fam, quote = F, sep = " ", row.names = F, col.names = F)
    
    return(NULL)
}

# =====================
# Main
# =====================

RandomFieldsUtils::RFoptions(warn_parallel=FALSE, install="no")

# Number of individuals to be simulated for phenotyping 
# Note that setting this value high requires a lot of memory
indiv <- 5e4

# Simulate a population based on the map_cattle2 map
dataset <- founder.simulation(nindi = 1e2, sex.quota = 0.5, map = map_cattle2, display.progress = FALSE, verbose = FALSE)
population <- creating.diploid(dataset = dataset, map = map_cattle2, n.additive = c(1000,1000), verbose=FALSE)
population <- breeding.diploid(population, breeding.size = indiv, heritability = c(0.5,0.5),
                                phenotyping = "all",
                                display.progress = FALSE, verbose=FALSE)
population <- breeding.diploid(population, breeding.size = 1e2, heritability = c(0.5,0.5),
                                phenotyping = "all",
                                display.progress = FALSE, verbose=FALSE)

# Extract genotypes and phenotypes
geno <- get.geno(population, gen=2) 
pheno <- get.pheno(population, gen=2)
bv <- get.bv(population, gen=2) 

# Write to PLINK binary format
mobps_to_bed(geno, pheno[1,], paste0(data_dir, "mobps_simulation"))
# Write true breeding values for reference
write.table(data.frame(bv=bv), paste0(data_dir,"mobps_simulation.bv"), quote = F, sep = " ", row.names = F, col.names = F)
