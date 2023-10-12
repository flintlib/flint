using GMP_jll, MPFR_jll, Artifacts
import Pkg

function gmp_artifact_dir()
    artifacts_toml = joinpath(dirname(dirname(Base.pathof(GMP_jll))), "StdlibArtifacts.toml")

    # If this file exists, it's a stdlib JLL and we must download the artifact ourselves
    if isfile(artifacts_toml)
        meta = artifact_meta("GMP", artifacts_toml)
        hash = Base.SHA1(meta["git-tree-sha1"])
        if !artifact_exists(hash)
            dl_info = first(meta["download"])
            Pkg.Artifacts.download_artifact(hash, dl_info["url"], dl_info["sha256"])
        end
        return artifact_path(hash)
    end

    # Otherwise, we can just use the artifact directory given to us by GMP_jll
    return GMP_jll.find_artifact_dir()
end

function mpfr_artifact_dir()
    artifacts_toml = joinpath(dirname(dirname(Base.pathof(MPFR_jll))), "StdlibArtifacts.toml")

    # If this file exists, it's a stdlib JLL and we must download the artifact ourselves
    if isfile(artifacts_toml)
        meta = artifact_meta("MPFR", artifacts_toml)
        hash = Base.SHA1(meta["git-tree-sha1"])
        if !artifact_exists(hash)
            dl_info = first(meta["download"])
            Pkg.Artifacts.download_artifact(hash, dl_info["url"], dl_info["sha256"])
        end
        return artifact_path(hash)
    end

    # Otherwise, we can just use the artifact directory given to us by MPFR_jll
    return MPFR_jll.find_artifact_dir()
end
