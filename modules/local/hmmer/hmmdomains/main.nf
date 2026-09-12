// HMMER_HMMRANK picks the best profile per sequence; this picks the best profile per stretch
// of sequence, so a multi-domain sequence keeps one hit per domain instead of collapsing to one.
process HMMER_HMMDOMAINS {
    tag "$meta.id"
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/mulled-v2-b2ec1fea5791d428eebb8c8ea7409c350d31dada:a447f6b7a6afde38352b24c30ae9cd6e39df95c4-1' :
        'quay.io/biocontainers/mulled-v2-b2ec1fea5791d428eebb8c8ea7409c350d31dada:a447f6b7a6afde38352b24c30ae9cd6e39df95c4-1' }"

    input:
    tuple val(meta), path(domtblouts)   // HMMER_HMMSEARCH.out.domain_summary
    val max_overlap                     // fraction of the shorter envelope two hits may share and still both be kept

    output:
    tuple val(meta), path("*.hmmdomains.tsv.gz"),       emit: domains
    tuple val(meta), path("*.hmmarchitectures.tsv.gz"), emit: architectures
    path "versions.yml", emit: versions, topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def prefix = task.ext.prefix ?: "${meta.id}"

    """
    #!/usr/bin/env Rscript
    library(readr)
    library(dplyr)
    library(tidyr)
    library(stringr)

    # Columns follow HMMER's --domtblout layout; the unused ones still need a name each to
    # keep the rest aligned, hence d0..d8. Splitting a table produces 23 character columns
    # for every row of it, so each file is cut down to the dozen typed columns that survive
    # before the next is read -- holding all of them at full width at once is what makes this
    # run out of memory on a search with many profiles.

    read_domtbl <- function(fname) {
        read_fwf(
            fname, fwf_cols(content = c(1, NA)),
            col_types = cols(content = col_character()), comment = '#'
        ) %>%
            filter(! str_detect(content, '^ *#')) %>%
            separate(
                content,
                c(
                    'accno', 'd0', 'tlen', 'model', 'd1', 'qlen', 'd2', 'd3', 'd4', 'd5', 'd6',
                    'cevalue', 'ievalue', 'score', 'd7',
                    'hmm_from', 'hmm_to', 'ali_from', 'ali_to', 'env_from', 'env_to', 'd8', 'rest'
                ),
                '\\\\s+', extra = 'merge', convert = FALSE
            ) %>%
            transmute(
                profile = basename(fname) %>% str_remove('^${prefix}\\\\.') %>% str_remove('\\\\.domtbl\\\\.gz\$'),
                accno, model,
                across(c(tlen, qlen, hmm_from, hmm_to, ali_from, ali_to, env_from, env_to), as.integer),
                across(c(cevalue, ievalue, score), as.double)
            )
    }

    domains <- bind_rows(lapply(c('${domtblouts.join("','")}'), read_domtbl)) %>%
        # One hmm file may hold several models, and then the file name alone doesn't say what
        # matched; only in that case is the model name worth carrying into the label.
        group_by(profile) %>%
        mutate(label = if (n_distinct(model) > 1) paste(profile, model, sep = ':') else profile) %>%
        ungroup()

    # Best-scoring first, keeping a hit only where it stays clear of every hit already kept.
    # A little overlap is allowed, since neighbouring domains commonly share a few residues,
    # measured against the shorter of the two envelopes so the tolerance means the same for a
    # short profile as a long one. Envelope, not ali, coordinates: ali bounds stop at the
    # aligned core and would under-count how much of the sequence a domain really occupies.
    #
    # Sorting first leaves each sequence's hits in one contiguous block, already in the order
    # they need to be considered, so the scan can walk index ranges over plain vectors. Handing
    # each sequence to a function instead costs milliseconds per sequence, which turns into
    # hours once a search covers a metagenome's worth of them.

    domains <- domains[
        order(
            domains[['accno']], -domains[['score']], domains[['ievalue']],
            domains[['profile']], domains[['model']]
        ),
    ]

    from   <- domains[['env_from']]
    to     <- domains[['env_to']]
    len    <- to - from + 1L
    keep   <- logical(nrow(domains))
    blocks <- rle(domains[['accno']])[['lengths']]
    ends   <- cumsum(blocks)
    starts <- ends - blocks + 1L

    for (b in seq_along(starts)) {
        if (blocks[b] == 1L) {
            keep[starts[b]] <- TRUE
            next
        }
        kept <- integer(0)
        for (i in starts[b]:ends[b]) {
            if (!length(kept) || all(
                pmax(0L, pmin(to[i], to[kept]) - pmax(from[i], from[kept]) + 1L) / pmin(len[i], len[kept]) <= ${max_overlap}
            )) {
                kept <- c(kept, i)
            }
        }
        keep[kept] <- TRUE
    }

    resolved <- domains[keep, ] %>%
        arrange(accno, env_from) %>%
        group_by(accno) %>%
        mutate(i = row_number(), n = n()) %>%
        ungroup() %>%
        transmute(
            accno, i, n, profile, model, label, score, cevalue, ievalue, tlen, qlen,
            hmm_from, hmm_to, hmm_len = hmm_to - hmm_from + 1L,
            ali_from, ali_to, ali_len = ali_to - ali_from + 1L,
            env_from, env_to, env_len = env_to - env_from + 1L
        )

    write_tsv(resolved, '${prefix}.hmmdomains.tsv.gz')

    # A schematic of where the domains sit, for reading rather than parsing -- anything
    # computed on belongs in the table above. Each dash stands for a twentieth of the
    # sequence left uncovered, so the picture is comparable between sequences whatever their
    # length, and any real gap keeps at least one dash rather than rounding away to nothing.
    # Domains kept despite overlapping leave no gap to draw and simply sit next to each other.

    dashes <- function(gap, tlen) {
        ifelse(gap > 0L, pmax(1L, as.integer(round(gap / tlen * 20))), 0L)
    }

    resolved %>%
        group_by(accno) %>%
        mutate(
            piece = paste0(
                strrep('-', dashes(env_from - if_else(i == 1L, 0L, lag(env_to)) - 1L, tlen)),
                '<', label, '>',
                if_else(i == n, strrep('-', dashes(tlen - env_to, tlen)), '')
            )
        ) %>%
        summarise(
            tlen = tlen[1], n_domains = n(), covered = sum(env_len),
            architecture = paste(label, collapse = '|'),
            sketch = paste(piece, collapse = ''),
            .groups = 'drop'
        ) %>%
        write_tsv('${prefix}.hmmarchitectures.tsv.gz')

    writeLines(
        c(
            "\\"${task.process}\\":",
            paste0("    r-base: ", paste0(R.Version()[c("major","minor")], collapse = ".")),
            paste0("    r-tidyverse: ", packageVersion('tidyverse'))
        ),
        "versions.yml"
    )
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"

    """
    echo 'accno\ti\tn\tprofile\tmodel\tlabel\tscore\tcevalue\tievalue\ttlen\tqlen\thmm_from\thmm_to\thmm_len\tali_from\tali_to\tali_len\tenv_from\tenv_to\tenv_len' > ${prefix}.hmmdomains.tsv
    gzip ${prefix}.hmmdomains.tsv

    echo 'accno\ttlen\tn_domains\tcovered\tarchitecture' > ${prefix}.hmmarchitectures.tsv
    gzip ${prefix}.hmmarchitectures.tsv

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        r-base: \$(Rscript -e "cat(strsplit(R.version[['version.string']], ' ')[[1]][3])")
        r-tidyverse: \$(Rscript -e "cat(as.character(packageVersion('tidyverse')))")
    END_VERSIONS
    """
}
