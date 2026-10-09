# Saves the STRING responses that test-backend-string.R replays: version,
# get_string_ids, and network (physical, functional, regulatory) for TP53,
# MDM2, CDKN1A, and MAPK1, from STRING 12.5. The get_string_ids request also
# sends P0DTC2 (SARS-CoV-2 spike), which has no human STRING protein, so
# the saved response has no row for it, as STRING leaves out unmapped IDs.
#
# Run from the package root:
#   Rscript tests/testthat/_fixtures/make_string_fixtures.R
stable_address <- "https://version-12-5.string-db.org"
fixture_dir <- "tests/testthat/_fixtures/string"
dir.create(fixture_dir, showWarnings = FALSE)

save_response <- function(file, url, body) {
    Sys.sleep(1)
    response <- httr::POST(url, body = body, encode = "form")
    httr::stop_for_status(response)
    writeLines(httr::content(response, as = "text", encoding = "UTF-8"),
               file.path(fixture_dir, file))
}

save_response("version.json", "https://string-db.org/api/json/version",
              list(caller_identity = "MSstatsBioNet"))
save_response("get_string_ids.json",
              paste0(stable_address, "/api/json/get_string_ids"),
              list(identifiers = paste("P04637", "Q00987", "P38936",
                                       "P28482", "P0DTC2", sep = "\r"),
                   species = "9606", echo_query = 1, limit = 1,
                   caller_identity = "MSstatsBioNet"))
string_ids <- c("9606.ENSP00000269305", "9606.ENSP00000258149",
                "9606.ENSP00000384849", "9606.ENSP00000215832")
for (network_type in c("physical", "functional", "regulatory")) {
    save_response(paste0("network_", network_type, ".json"),
                  paste0(stable_address, "/api/json/network"),
                  list(identifiers = paste(string_ids, collapse = "\r"),
                       species = "9606", network_type = network_type,
                       required_score = 0,
                       caller_identity = "MSstatsBioNet"))
}
