setClass("CESAnalysis", representation(maf = "data.table", trinucleotide_mutation_weights = "list",
          progressions = "character", mutrates_list = "list", dndscv_out_list = "list",
          excluded = "data.table", selection_results = "data.table", gene_epistasis_results = "data.table", coverage = "list",
          genome = "BSgenome", advanced = "list", genome_data_dir = "character", status = "list", samples = "data.table", 
          mutations = "list"))

setMethod("show", "CESAnalysis", 
  function(object) {
    steps = names(object@status)
    for (step in steps) {
      cat(paste0(step,": ", object@status[[step]], "\n"))
      if (step == "MAF data") {
        num_excluded_records = object@excluded[, .N]
        if(num_excluded_records > 0) {
          cat(paste0("(View ", num_excluded_records, " excluded MAF records with excluded_maf_records())\n"))
        }
      }
    }
    cat(paste0("[Created in cancereffectsizeR, version ", object@advanced$version, ".]"))
  }
)


setValidity("CESAnalysis",
    function(object) {
     # add validation later
     TRUE
    }
)
