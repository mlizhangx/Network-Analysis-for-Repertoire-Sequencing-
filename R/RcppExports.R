

.hamAdjacencyMatSparse <- function(strings, maxdist, drop_deg_zero) {
    .Call(`_NAIR_hamAdjacencyMatSparse`, strings, maxdist, drop_deg_zero)
}

hamDistBounded <- function(a, b, k) {
    .Call(`_NAIR_hamDistBounded`, a, b, k)
}

.levAdjacencyMatSparse <- function(strings, maxdist, drop_deg_zero) {
    .Call(`_NAIR_levAdjacencyMatSparse`, strings, maxdist, drop_deg_zero)
}

levDistBounded <- function(a, b, k) {
    .Call(`_NAIR_levDistBounded`, a, b, k)
}

