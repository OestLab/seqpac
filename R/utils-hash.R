#' @keywords internal
#' @noRd
.md5_hash <- function(x) {
  tmp <- tempfile()
  on.exit(unlink(tmp), add = TRUE)
  saveRDS(x, tmp, compress = FALSE)
  unname(tools::md5sum(tmp))
}