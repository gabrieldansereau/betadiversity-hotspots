## List files to version control manually

# Raw files
# Normally these should not be changed
_raw_files = [
    "data/raw/ebd_sample.txt",
    "data/raw/ebd_warblers_cut.csv",
    "data/raw/ebd_warblers_head.csv",
    "data/raw/ebd_warblers.csv",
]

# Processed files
# The corresponding placeholder should also be created anywhere these files are written.
_jld_files = [
    # "data/jld2/bart-distributions_xtras.jld2",
    # "data/jld2/bart-distributions.jld2",
    # "data/jld2/bio-distributions.jld2",
    "data/jld2/comparison-results.jld2",
    # "data/jld2/raw-distributions.jld2",
    # "data/jld2/rf-distributions.jld2",
]
_proc_files = [
    # "data/proc/bart_predictions_lower.csv",
    # "data/proc/bart_predictions_pres.csv",
    # "data/proc/bart_predictions_prob.csv",
    # "data/proc/bart_predictions_upper.csv",
    # "data/proc/bart_summaries.csv",
    # "data/proc/bart_varimps.csv",
    "data/proc/ebd_warblers_prep.csv"
]
_raster_files = [
    "data/raster/bart_xtras_lower-distrib.tif",
    "data/raster/bart_xtras_prob-distrib.tif",
    "data/raster/bart_xtras_upper-distrib.tif",
]

## Create placeholders
_placeholder_files = vcat(_jld_files, _proc_files, _raster_files, _raw_files)
_placeholder_files = replace.(_placeholder_files, "." => "_placeholder.")

function _create_placeholder_files!(files)
    for f in files
        touch(f)
        open(f, "w") do io
            write(io, string(Dates.now()))
        end;
    end
end
# _create_placeholder_files!(_placeholder_files)

## Create function to manually version control files (separately for raw & proc)
"""
    verify_raw_files(; touch_placeholders=false)
    verify_proc_files(; touch_placeholders=false)

Verifies if the data files are more recent than their corresponding placeholder
files. Placeholder files were created for important data files which are too big
to be version controlled.

This function will throw a warning when the data files are more recent than
their placeholder files, indicating a possibly unwanted and dangerous change for
which verifications should be performed.

Using `touch_placeholders=true` will update the placeholders' timestampsso they
are more recent than the data files. The placeholders should then be committed
with a description of the changes made to the data files (or discarded if it's
confirmed not important change was made).

"""
function verify_raw_files(; touch_placeholders=false)
    files = _raw_files
    placeholder_files = filter(startswith("data/raw/"), _placeholder_files)
    more_recent = mtime.(files) .> mtime.(placeholder_files)
    if any(more_recent)
        @warn """\n
        The following raw data files are more recent than their corresponding placeholders.
        Raw data files should not normally change, so make sure these changes were made on
        purpose. These files are not version controlled.

            $(files[more_recent])

        Run the following command after verification to update the placeholder files
        timestamps and remove this warning.

            verify_raw_files(; touch_placeholders=true)
        """
        if touch_placeholders
            @info "Re-writing placeholder files to update timestaps and disable the warning"
            for p in placeholder_files[more_recent]
                open(p, "w") do io
                    write(io, string(Dates.now()))
                end;
            end
        end
    end
end
# verify_raw_files()

function verify_proc_files(; touch_placeholders=false)
    files = vcat(_proc_files, _jld_files, _raster_files)
    placeholder_files = filter(!startswith("data/raw/"), _placeholder_files)
    more_recent = mtime.(files) .> mtime.(placeholder_files)
    if any(more_recent)
        @warn """\n
        The following processed data files are more recent than their corresponding placeholders.
        Processed data files are expected to changed in the analyses, but make sure the
        changes were made on purpose. These files are not version controlled.

            $(files[more_recent])

        Run the following command after verification to update the placeholder files
        timestamps and remove this warning.

            verify_proc_files(; touch_placeholders=true)
        """
        if touch_placeholders
            @info "Re-writing placeholder files to update timestaps and disable the warning"
            for p in placeholder_files[more_recent]
                open(p, "w") do io
                    write(io, string(Dates.now()))
                end;
            end
        end
    end
end
# verify_proc_files()
