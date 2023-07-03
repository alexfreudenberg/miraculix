#  Authors 
#  Alexander Freudenberg, alexander.freudenberg@stads.de

#  Copyright (C) 2022-2023 Alexander Freudenberg

#  Licensed under the Apache License, Version 2.0 (the "License");
#  you may not use this file except in compliance with the License.
#  You may obtain a copy of the License at

#     http://www.apache.org/licenses/LICENSE-2.0

#  Unless required by applicable law or agreed to in writing, software
#  distributed under the License is distributed on an "AS IS" BASIS,
#  WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
#  See the License for the specific language governing permissions and
#  limitations under the License.

import Dates;
include("benchmark_suite.jl")

# =====================
# Get arguments for benchmarking
# =====================

VALID_MODES = Set(["miraculix","PLINK","cuBLAS","GCTA"])

if length(ARGS) > 0
    mode = ARGS[1]
    if mode ∉ VALID_MODES
        error("First command-line argument needs to be `GPU` or `PLINK`, specifying the Benchmark Suite to run.")
    end
    if length(ARGS) > 1
        arg = ARGS[2]
        tag = @tagged mode && eval(arg)
    else
        tag = @tagged mode
    end
else
    error("No command-line arguments for execution mode provided.")
end

if !isdir(LOG_DIR)
    mkdir(LOG_DIR)
end
date = Dates.today()

# =====================
# Start benchmarks
# =====================
results_file = "results_$(join(ARGS,' '))-$date.json"
println(results_file)

results = run(suite[tag], verbose = true, samples = 5, evals = 1)

BenchmarkTools.save("$LOG_DIR/" * results_file,results)