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

if(!require(MoBPS)){
    devtools::install_github("tpook92/MoBPS", subdir="pkg")
}
if(!require(MoBPSmaps)){
    devtools::install_github("tpook92/MoBPS", subdir="pkg-maps")
}
indiv <- 1e3

dataset <- founder.simulation(nindi=indiv, sex.quota = 0.5, map = map_cattle2, display.progress = !FALSE, verbose= !FALSE)
population <- creating.diploid(dataset = dataset, map = map_cattle2, n.additive = c(1000,1000), verbose=FALSE)
population <- breeding.diploid(population, heritability = c(0.5,0.5),
                                phenotyping.database = cbind(1,1),
                                n.observation = c(1,0), display.progress = FALSE, verbose=FALSE)
population <- breeding.diploid(population, heritability = c(0.5,0.5),
                                phenotyping.database = cbind(1,2),
                                n.observation = c(1,1), display.progress = FALSE, verbose=FALSE)

geno <- get.geno(population, gen=1) 
pheno <- get.pheno(population, gen=1)
bv <- get.bv(population, gen=1) 