
#' Get MAF data from TCGA cohort
#' 
#' This convenience function queries the Genomic Data Commons API to get MAF data
#' generated with the Aliquot Ensemble Somatic Variant Merging and Masking workflow for
#' the specified project, and writes an MAF file. The API always provides data from the
#' latest data release. This function might work with non-TCGA MAF data hosted on GDC
#' (e.g., TARGET and GENIE-MSK), but it hasn't been tested and users should proceed with
#' caution.
#' 
#' TCGA cohort MAFs will be structured as downloaded, with a Unique_Patient_Identifier
#' column generated from the first 12 characters of Tumor_Sample_Barcode. When passed to
#' preload_maf() or load_maf(), this column will supersede Tumor_Sample_Barcode. In the
#' handful of patients with multiple Tumor_Sample_Barcodes (essentially replicated
#' sequencing, with very high variant overlap), these functions will effectively take the
#' union of these samples for each patient. Relatedly, the small number of TCGA
#' non-primary tumor samples should not be handled this way (and such samples are by
#' default removed by this function).
#' 
#' Temporary aliquot MAF files downloaded by this function are deleted after they are read.
#' 
#' @param project TCGA project name (e.g., "TCGA-BRCA").
#' @param filename Output filename where MAF data should be saved. Must end in '.maf'
#'   (plaintext) or '.maf.gz' (gzip compressed).
#' @param exclude_TCGA_nonprimary Default TRUE. For TCGA projects, exclude samples not
#'   associated with a patient's initial primary tumor. (In many TCGA projects, a small
#'   handful of patients have metastatic, recurrent, or additional primary samples.)
#' @param test_run Default FALSE. When TRUE, gets MAF data for a few samples instead of the whole cohort.
#' @param allow_failure Default FALSE. When TRUE, if data download fails for some samples,
#'   still return data from whatever samples did have complete, successful data download.
#'   The data returned will be valid, but the cohort may be incomplete.
#' @export
get_TCGA_project_MAF = function(project = NULL, filename = NULL, test_run = FALSE,
                               exclude_TCGA_nonprimary = TRUE, allow_failure = FALSE) {
  if (! is.character(project) || length(project) != 1) {
    stop("project should be 1-length character (e.g., TCGA-BRCA).")
  }
  
  if (! is.character(filename) || length(filename) != 1) {
    stop("filename should be a pathname for the output MAF file.")
  }
  if (file.exists(filename)) {
      stop('Specified output filename already exists.')
  }
  dir = dirname(filename)
  if(! dir.exists(dir)) {
    stop("The directory specified in the output file path does not exist.")
  }
  
  if(! file_test('-w', dir)) {
    if (dir == '.') {
      msg = paste0("You don't have write permissions in your working directory.",
                   " Change directories or specify a different path for your output filename.")
      stop(pretty_message(msg, emit = F))
    } else {
      stop("The directory specified in the output file path is not writeable.")
    }
  }
  
  if(endsWith(filename, '.maf.gz')) {
    gzipped = TRUE
  } else if(endsWith(filename, '.maf')) {
    gzipped = FALSE
  } else {
    stop("filename must end in .maf or .maf.gz")
  }
  
  if (file.create(filename, showWarnings = FALSE)) {
    file.remove(filename)
  } else {
    stop("The specified filename couldn't be used.")
  }
  
  if(! is.logical(exclude_TCGA_nonprimary) || length(exclude_TCGA_nonprimary) != 1) {
    stop("remove_nonprimary should be T/F.")
  }
                   
  
  if(! is.logical(allow_failure) || length(allow_failure) != 1) {
    stop("allow_failure should be T/F.")
  }
  
  if(! is.logical(test_run) || length(test_run) != 1) {
    stop("test_run should be T/F.")
  }
  
  projects_endpt = 'https://api.gdc.cancer.gov/projects'
  files_endpt = 'https://api.gdc.cancer.gov/files'
  data_endpt = 'https://api.gdc.cancer.gov/data'
  versions_endpt = 'https://api.gdc.cancer.gov/files/versions'
  
  response = httr::GET(projects_endpt, 
                       query = list(fields = "project_id", size = '10000', format = 'JSON'))
  if (response$status_code != 200) {
    msg = paste0("GDC API query (for names of valid TCGA projects) failed ", "(got response ", 
                  response$status_code, "). Perhaps there's a problem with your network;",
                 " alternatively, check the GDC site to ensure that the data portal is online.")
    stop(pretty_message(msg, emit = F))
  }
  tcga_projects = sapply(rjson::fromJSON(rawToChar(response$content))$data$hits, '[', 'project_id')
  if(! project %in% tcga_projects) {
    fixed_project = paste0('TCGA-', project)
    if (fixed_project %in% tcga_projects) {
      project = fixed_project
    } else {
      pretty_message(paste0("TCGA projects: ", paste(tcga_projects, collapse = ", "), "."))
      stop("project is not a valid TCGA project name.")
    }
  }
  
  # We'll only remove the non-primary tumors for the user
  is_tcga_project = startsWith(project, 'TCGA')
  
  if (! exclude_TCGA_nonprimary & is_tcga_project) {
    stop("Setting exclude_TCGA_nonprimary = FALSE doesn't do anything when the project is not a TCGA project.")
  }
  
  exclude_TCGA_nonprimary = exclude_TCGA_nonprimary && is_tcga_project
  
  user_call = paste0("# Generated with cancereffectsizeR v", packageVersion('cancereffectsizeR'),
                   ': get_TCGA_project_MAF(project = ', deparse(project), ', test_run = ', deparse(test_run),
                   ', allow_failure = ', deparse(allow_failure))
  if(is_tcga_project) {
    user_call = paste0(user_call, ", exclude_TCGA_nonprimary = ", deparse(exclude_TCGA_nonprimary))
  }
  user_call = paste0(user_call, ')')
  
  filters = sprintf(
  '{
    "op": "and",
    "content": [
      {
        "op": "in",
        "content":{
          "field": "cases.project.project_id",
          "value": ["%s"]
        }
      },
      {
        "op": "in",
        "content":{
          "field": "analysis.workflow_type",
          "value": ["Aliquot Ensemble Somatic Variant Merging and Masking"]
        }
      },
      {
        "op": "in",
        "content":{
          "field": "access",
          "value": ["open"]
        }
      }
    ]
  }', project)
  
  num_files = ifelse(test_run, '5', '100000')
  response = httr::GET(files_endpt, 
                       query = list(filters = filters, 
                                    fields = "file_name,md5sum,release.version", size = num_files, format = 'JSON'))
  
  if (response$status_code != 200) {
    msg = paste0("Could not get list of ", project, " cases (GDC API query failed with status code ", 
         response$status_code, ")." )
    stop(pretty_message(msg, emit = F))
  }

  content = rjson::fromJSON(rawToChar(response$content))
  files = rbindlist(lapply(content$data$hits, '[', c("id", "file_name", "md5sum")))
  if(files[, .N] == 0) {
    stop("No matching TCGA participants (patients).")
  }
  
  message("Downloading ", files[, .N] , " temporary MAF files...")
  tmp_dir = tempdir()
  files[, path := paste0(tmp_dir, '/', file_name)]
  files[, url := paste0(data_endpt, '/', id)]
  
  pbapply::pbmapply(download.file, files$url, files$path, MoreArgs = list(quiet = TRUE, mode = 'wb'))
  
  message("Verifying files...")
  actual_sums = tools::md5sum(files$path)
  files[names(actual_sums), obs_sum := actual_sums, on = 'path']
  files[, failed := TRUE]
  files[obs_sum == md5sum, failed := FALSE]
  
  num_failed = files[failed == TRUE, .N]
  if(num_failed > 0) {
    failing_urls = files[failed == TRUE, url]
    if(length(failing_urls) > 20) {
      failing_urls = c(failing_urls[1:15], paste0("and ", length(failing_urls) - 15, " more."))
    }
    failing_urls = paste(failing_urls, collapse = ",\n")
    if(allow_failure) {
      unlink(files[failed == TRUE, path])
      files = files[failed == FALSE]
      warning("Some files did not download correctly:\n", failing_urls, immediate. = TRUE)
    } else {
      unlink(files$path)
      stop("Some files did not download correctly:\n", failing_urls)
    }
  }
  
  to_read = files$path
  names(to_read) = files$id
  cohort_maf = rbindlist(lapply(to_read, fread, skip = 'Hugo'), idcol = "source_file_id") # column headers start with Hugo_Symbol (comment lines precede)
  
  if(is_tcga_project) {
    cohort_maf[, Unique_Patient_Identifier := substr(Tumor_Sample_Barcode, 1, 12)]
    cohort_maf[, c("V1", "V2", "V3", "type_vial", "portion_analyte", "plate", "center") := tstrsplit(Tumor_Sample_Barcode, split = "-")]
    cohort_maf[, c("V1", "V2", "V3", "portion_analyte", "plate", "center") := NULL] # already extracted participant ID
    cohort_maf[, tissue_type := substr(type_vial, 1, 2)]
    
    # 01 = primary solid tumor; 03 = primary blood-derived cancer (sample of peripheral blood)
    # (from https://gdc.cancer.gov/resources-tcga-users/tcga-code-tables/sample-type-codes)
    cohort_maf = cohort_maf[, primary_sample := tissue_type %in% c("01", "03")]
    
    num_not_primary = cohort_maf[primary_sample == FALSE, uniqueN(Tumor_Sample_Barcode)]
    if(num_not_primary > 0) {
      if (exclude_TCGA_nonprimary) {
        pretty_message(paste0("Removing ", num_not_primary, " samples that are not of tissue types 01 (primary solid tumor)",
                              " or 03 (primary blood-derived cancer). Usually, there are other primary samples from the same patient(s),",
                              " and these will remain in the MAF data."))
        cohort_maf = cohort_maf[primary_sample == TRUE]
      } else {
        pretty_message("Some samples are not of tissue types 01 (primary solid tumor) or 03 (primary blood-derived cancer).")
      }
    }
    cohort_maf[, c("type_vial", "tissue_type", "primary_sample") := NULL]
  } else {
    msg = paste0("Since this isn't a TCGA project, sample ID parsing is not implemented. Multiple samples (as indicated by Tumor_Sample_Barcode)", 
                 " may come from the same patient, and no Unique_Patient_Identifier column has been created.")
    warning(msg)
  }

  message("Deleting temporary files...")
  unlink(files$path)
  
  if(is_tcga_project) {
    num_samples = uniqueN(cohort_maf$Tumor_Sample_Barcode)
    num_participants = uniqueN(cohort_maf$Unique_Patient_Identifier)
    pretty_message(paste0("Writing MAF file covering ", num_samples , " samples from ",
                          num_participants, " patients."))
  } else {
    message("Writing MAF file covering ", num_samples, " samples.")
  }
  
  # Get current GDC release (one would think there would be a simpler way)
  response = httr::GET(paste0(versions_endpt, '/', cohort_maf$source_file_id[1], "?format=JSON"))
  if (response$status_code != 200) {
    warning("Failed to add GDC release to MAF file header, due to a failed API call (weird!).")
    latest_release = "Unknown (at least 33.1)"
  } else {
    latest_release = tail(rjson::fromJSON(rawToChar(response$content))[[1]]$latest_release, n = 1)
  }
  
  if (gzipped) {
    out = gzfile(filename)
  } else {
    out = file(filename)
  }
  
  headers = c(user_call, 
              paste0("# MAF data from Genomic Data Commons data release ", latest_release, 
                     " (Aliquot Ensemble Somatic Variant Merging and Masking workflow)."),
              paste(names(cohort_maf), collapse = "\t"))
  writeLines(headers, out)
  close(out)
  fwrite(cohort_maf, filename, append = TRUE, sep = "\t")
  pretty_message(paste0("MAF file saved to ", filename, "."))
}



