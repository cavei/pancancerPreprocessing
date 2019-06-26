assignCGToGenePromoter <- function(data, promoterUP=-2000, promoterDW=500) {
  discretePos <- rep("geneBody", length(data$Position_to_TSS))
  discretePos[data$Position_to_TSS >= promoterUP & data$Position_to_TSS <= promoterDW] <- "Promoter"
  discretePos
}

#' @importFrom stats na.omit
#' @importFrom AnnotationDbi mapIds
#' @importFrom SummarizedExperiment rowData
#' @importFrom org.Hs.eg.db org.Hs.eg.db
#'
createEntrezMethylationMatrix <- function(methylationAssay) {
  requireNamespace("org.Hs.eg.db")

  meta = rowData(methylationAssay)

  expandedMeta <- set_expand(meta, cols=c("Gene_Symbol", "Position_to_TSS"))
  expandedMeta$Position_to_TSS <- as.numeric(expandedMeta$Position_to_TSS)

  expandedMeta <- stats::na.omit(expandedMeta)
  expandedMeta$discretePos <- assignCGToGenePromoter(expandedMeta)

  promoterMeta <- unique(expandedMeta[expandedMeta$discretePos=="Promoter",c("Composite.Element.REF", "Gene_Symbol")])

  # symbol2entrez <- mapIds(org.Hs.eg.db, keys=promoterMeta$Gene_Symbol, column="ENTREZID", keytype="SYMBOL", multiVals="list")
  symbol2entrez <- AnnotationDbi::mapIds(org.Hs.eg.db::org.Hs.eg.db, keys=promoterMeta$Gene_Symbol, column="ENTREZID", keytype="SYMBOL")

  promoterMeta$entrez <- symbol2entrez
  promoterMeta <- stats::na.omit(promoterMeta)

  methylationPromoter <- methylationAssay[promoterMeta$Composite.Element.REF] ## filter with valid CG
  assayPromoter <- methylationPromoter@assays[[1]]

  summarized <- tapply(row.names(assayPromoter), promoterMeta$entrez, function(cg){
    colMeans(assayPromoter[cg, ,drop=F], na.rm = T)
  })

  summarized <- do.call(rbind, summarized)
  colnames(summarized) <- substr(colnames(summarized), 1, 12)
  keep <- apply(summarized, 1, function(x) { !all(is.na(x) )})
  list(assay=summarized[keep, , drop=F], annotation=promoterMeta)
}

createDiscreteVersionOfMethylation <- function(data, breaks= c(0, 0.2, 0.8, 1), labels=c(0,0.5,1)) {
  if (!(length(breaks)-1 == length(labels)))
    stop("Labels length must be the same of length minus 1.")

  data[is.na(data)] <- NA

  binary <- t(apply(data, 1, function(x) {
    as.character(cut(x, breaks=breaks, labels=labels))
  }))

  colnames(binary) <- colnames(data)
  binary
}

# Clinical Function
resolveDuplicatedAndMakeDataFrame<- function(m) {
  if(NCOL(m)!=3)
    stop("Incorrect format. 'm' must have 3 columns: barcode, status, days")
  lu <- tapply(seq_len(NROW(m)), m[,1], function(idx) {
    r <- unique(m[idx,, drop=F])
    if (NROW(r)==1)
      return(r)

    if (all(r[,2]=="0")){
      j <- which.max(as.numeric(r[,3]))
      return(r[j, ,drop=F])
    }

    r <- r[which(r[,2]=="1"),,drop=F]
    j <- which.min(as.numeric(r[,3]))
    r[j, , drop=F]
  })
  lu <- do.call(rbind, lu)
  data.frame(status=as.numeric(lu[,2]),
             days = as.numeric(lu[,3]),
             row.names=lu[,1], stringsAsFactors=F)
}

createFollowUp <- function(followup, newtumor) {
  dropped <- setdiff(followup$bcr_patient_barcode,
                     newtumor$bcr_patient_barcode)
  recDaysNA <- rep(NA, length(dropped))

  names(recDaysNA) <- dropped

  recDays <- newtumor$days_to_new_tumor_event_after_initial_treatment
  names(recDays) <- newtumor$bcr_patient_barcode
  recDays <- c(recDays, recDaysNA)[followup$bcr_patient_barcode]

  barcode = followup$bcr_patient_barcode
  status  = followup$vital_status
  fup     = followup$days_to_last_followup
  death   = followup$days_to_death
  rec     = followup$new_tumor_events
  recDays = recDays
  remission = followup$followup_treatment_success
  remission1 = followup$primary_therapy_outcome_success
  cancerStatus = followup$person_neoplasm_cancer_status

  os <- t(sapply(seq_len(length(barcode)), function(idx) {
    if (status[idx]=="Dead") {
      os_binary = 1
      days = death[idx]
    } else if (status[idx]=="Alive"){
      os_binary = 0
      days = fup[idx]
    } else {
      stop("unrecognized stauts")
    }
    c(barcode[idx],os_binary, days)
  }))

  pfs <- t(sapply(seq_len(length(barcode)), function(idx) {
    if (rec[idx]!="" & rec[idx]!="NA" & rec[idx]!="NO"){
      pfs_binary = 1
      days = recDays[idx]
    } else if (status[idx]=="Dead") {
      # old condition for former data was only remission[idx] == "Progressive Disease"
      if (remission[idx] == "Progressive Disease" |
          remission1[idx] == "Progressive Disease" |
          cancerStatus[idx] == "WITH TUMOR") {
        pfs_binary=1
      } else {
        pfs_binary = 0
      }
      days = death[idx]

    } else if (status[idx]=="Alive"){
      pfs_binary = 0
      days = fup[idx]
    } else {
      stop("unrecognized stauts")
    }
    c(barcode[idx], pfs_binary, days)
  }))

  survAnnot.os  <- resolveDuplicatedAndMakeDataFrame(os)
  survAnnot.pfs <- resolveDuplicatedAndMakeDataFrame(pfs)
  return(list(os=survAnnot.os, pfs=survAnnot.pfs))
}

extractTCGAPatientsName <- function(nms){
  substr(nms, 1, 12)
}
