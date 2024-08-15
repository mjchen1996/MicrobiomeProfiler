## Extend readLines() for parsing a remote .gz format file.
#' @importFrom utils download.file
#' @keywords internal
readLines_gzfile <- function(url, ...) {
    gzf <- tempfile(fileext = ".gz")
    download.file(url, gzf)
    con <- gzfile(gzf)
    on.exit(close(con))
    content <- readLines(con, ...)
    return(content)
}

## PubChem Resource Description Framework Api
#' @importFrom yulab.utils yread
pubchem_rest <- function(rest_url) {
    message("Reading PubChem RDF online: \"", rest_url, "\"...")
    content <- yread(rest_url, reader = readLines_gzfile)
    content <- content[grep("^[^@]", content)]
    content <- content |>
        strsplit("\t") |>
        do.call(rbind, args = _)
    res <- data.frame(from = content[, 1], linktype = content[, 2],
                      to = content[, 3])
    return(res)
}

## Download PubChem annotation of the latest version of PubChem pathway and
## stored in a 'GSON' object.
# @param rm.protein remove protein nodes in PubChem pathway, defualt TRUE.
gson_PubChemPathway <- function(rm.protein = TRUE) {
    url <- "https://ftp.ncbi.nlm.nih.gov/pubchem/RDF/pathway/pc_pathway.ttl.gz"
    pc_pw <- pubchem_rest(url)

    gsid2name <- pc_pw[pc_pw$linktype == "dcterms:title", c("from", "to")]
    gsid2name <- setNames(gsid2name, c("gsid", "name"))
    gsid2name$gsid <- sub("pathway:", "", gsid2name$gsid)
    gsid2name$name <- gsub("\"| \\.$", "", gsid2name$name)

    if (rm.protein) {
        patt <- "compound:"
    } else {
       patt <- "compound:|protein:"
    }
    gsid2gene <- pc_pw[grep(patt, pc_pw$to), c("from", "to")]
    gsid2gene <- setNames(gsid2gene, c("gsid", "gene"))
    i <- which(gsid2gene$gsid != "")
    times <- diff(c(i, nrow(gsid2gene) + 1))
    gsid2gene$gsid <- rep(gsid2gene$gsid[i], times)
    gsid2gene$gsid <- sub("pathway:", "", gsid2gene$gsid)
    gsid2gene$gene <- gsub("[compound:|protein:| \\.| ,]", "", gsid2gene$gene)

    raw_html <- readLines("https://ftp.ncbi.nlm.nih.gov/pubchem/RDF/pathway/")
    raw_html <- raw_html[grep("pc_pathway.ttl.gz", raw_html)]
    version <- sub("^.+(\\d{4}-\\d{2}-\\d{2} \\d{2}:\\d{2}).+$", "\\1", raw_html)
    version <- paste0("Release ", version)
    keytype <- "PubChem Compound ID (CID)"
    if (!rm.protein) {
        keytype <- paste0(keytype, " or Protein Accession Number")
    }
    gson_pc_pw <- gson(gsid2gene = gsid2gene,
                       gsid2name = gsid2name,
                       species = "Homo sapiens",
                       gsname = "PubChem Pathway",
                       version = version,
                       keytype = keytype,
                       accessed_date = as.character(Sys.Date()))
    return(gson_pc_pw)
}
