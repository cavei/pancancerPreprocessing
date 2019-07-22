#' From maf format prepare the mutation table
#'
#' @param maf a maf format table
#' @param impact select by impact
#' @param mutType select by mutation type
#' @param filterByThisEntrez filter to keep only entred listed
#' @param patients keep only the patients listed
#' @param select_type logical to perform selection on the type
#' @param tumor_site_type keep only tumor samples of the type. 01 primary solid tumor.
#' @param select_columns named numeric to specyfy the columns to select
#'
#' @rdname mutation-handler
#'
#' @return a matrix entrez x samples
#'
#' @export
prepareMutations <- function(maf, impact=NULL, mutType=NULL, filterByThisEntrez=NULL, patients=NULL,
                             select_type=TRUE, tumor_site_type="01",
                             select_columns=c(hugo=1,entrez=2, impact=94, type=9, patient=16)) {

  if (is.null(names(select_columns)))
    stop("you need names associated to the columns to select")

  if (!is.numeric(select_columns))
    stop("select_columns must be culumn index")

  if (NCOL(maf) < max(select_columns))
    stop("your imput does not appear to be a valid maf file")

  if (select_type) {
    bcode <- maf$Tumor_Sample_Barcode
    select <- substr(bcode, 14,15) == tumor_site_type
    maf <- maf[select, , drop=F]
    if (nrow(maf)==0){
      warning("No lines for maf. Maybe no sample for primary tumor")
      stop(paste0("No samples with primary solid tumor. Choose among the following.\n",
        paste(unique(substr(bcode, 14,15)), collapse=", "))
      )
    }
  }

  allPatientsMeasured <- unique(extractTCGAPatientsName(maf[[16]]))
  smallMaf <- maf[, select_columns]
  colnames(smallMaf) <- names(select_columns)
  smallMaf$patient <- extractTCGAPatientsName(smallMaf$patient)

  if (!is.null(impact))
    smallMaf <- selectImpact(smallMaf, impact)

  if (!(is.null(mutType)))
    smallMaf <- selectMutation(maf, mutType)

  mut <- summarizeMutation(smallMaf)
  mut$entrez <- as.character(mut$entrez)

  if (!is.null(filterByThisEntrez))
    mut <- mut[mut$entrez %in% filterByThisEntrez, , drop=F]

  if (!is.null(patients)) {
    patients <- intersect(unique(mut$patient), patients)
    mut <- mut[mut$patient %in% patients, , drop=F]
  }
  entrez <- unique(mut$entrez)
  patients <- unique(mut$patient)
  mutations <- matrix(0, nrow=length(entrez), ncol=length(patients), dimnames = list(entrez, patients))
  for (i in seq_along(mut$entrez)) {
    x = as.character(mut$entrez[i])
    y = mut$patient[i]
    mutations[x,y] <- 1
  }
  if (length(setdiff(patients, allPatientsMeasured)) > 0)
    warning("some patients were excluded beacuse they do not have any of your favourite mutations.
            Consider adding an 0 column to prevent patients loss.")
  list(data=mutations, allPatients=allPatientsMeasured, mutationsTypes=mut)
}


smallMaf_min_columns <- c("hugo", "entrez", "patient", "type", "impact")

checkSmallMaf_format <- function(smallMaf) {
  missing <- setdiff(smallMaf_min_columns, colnames(smallMaf))
  if (length(missing)) {
    stop(paste0("Invalid smallMaf format: columns \"", paste(missing, collapse = ", ", "\" are missing")))
  }
}

#' Collapse mutations
#'   Given the same gene and patients mutation impact and type are collapsed.
#'
#' @inheritParams selectMutation
#'
#' @rdname mutation-handler
#'
#' @export
summarizeMutation <- function(smallMaf) {
  if (!is.data.frame(smallMaf))
    stop("smallMaf need to be a data.frame")

  checkSmallMaf_format(smallMaf)

  patient2gene2mutation <- tapply(1:NROW(smallMaf), paste(smallMaf$entrez, smallMaf$patient, sep="_"),
                                  function(idx) {
                                    mat <- smallMaf[idx, , drop=F]
                                    type <- paste(smallMaf$type[idx], collapse=";")
                                    impact <- paste(smallMaf$impact[idx], collapse=";")
                                    hugo <- mat[1, "hugo"]
                                    entrez <- mat[1, "entrez"]
                                    patient <- mat[1, "patient"]
                                    data.frame(hugo, entrez, impact, type, patient, stringsAsFactors = F)
                                  })
  do.call(rbind,patient2gene2mutation)
}

#' Select mutation by type
#'
#' @param smallMaf reduced version of a maf file
#' @param type mutation type
#'
#' @return smallMaf subset (matrix or data.frame)
#' @rdname mutation-handler
#'
#' @export
selectMutation <- function(smallMaf, type){
  checkSmallMaf_format(smallMaf)
  mutationSelection <- smallMaf$type %in% type
  # hugo <- smallMaf[[1]][mutationSelection]
  # entrez <- smallMaf[[2]][mutationSelection]
  # samples <- smallMaf[[16]][mutationSelection]
  # patients <- substr(samples, 1,12)
  # type <- smallMaf[[9]][mutationSelection]
  # data.frame(hugo=hugo, entrez=entrez, patient=patients, type=type, stringsAsFactors = FALSE)
  smallMaf[mutationSelection, ,drop=F]
}

#' Select mutation by impact
#'
#' @inheritParams selectMutation
#'
#' @rdname mutation-handler
#'
#' @export
selectImpact <- function(smallMaf, type){
  checkSmallMaf_format(smallMaf)
  impactSelection <- smallMaf$impact %in% type
  # hugo <- smallMaf$hugo[impactSelection]
  # entrez <- smallMaf$entrez[[2]][impactSelection]
  # samples <- smallMaf$1[[16]][impactSelection]
  # patients <- substr(samples, 1,12)
  # type <- smallMaf[[9]][impactSelection]
  # impact <- smallMaf[[94]][impactSelection]
  # data.frame(hugo=hugo, entrez=entrez, patient=patients, type=type, impact=impact, stringsAsFactors = FALSE)
  smallMaf[impactSelection, , drop=F]
}
